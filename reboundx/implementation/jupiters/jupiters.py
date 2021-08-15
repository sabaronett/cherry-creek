import numpy as np
import os
import psutil
import rebound
import reboundx
import sys
import time

# initialize constants
init_a = float(sys.argv[1]) # au
T0 = 1.2324996566303585e10  # Sun's age ~50 Myr pre-TRGB
M0 = 0.98984348238912501    # initial mass of Sun

def makesim(init_a):
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.units = ('yr', 'AU', 'Msun')
    sim.add(m=M0, hash='Sun')
    sim.add(m=9.547919e-4, a=init_a, r=5.11347118e-9, hash='Jupiter')
    sim.move_to_com()
    sim.dt = 0.05
    sim.collision = "direct"
    rebx = reboundx.Extras(sim)
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)
    sim.status()
    return sim, rebx, tides

def makesubdir(name):
    if not os.path.exists(name):
        os.makedirs(name)

def memory_usage_psutil():
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20)
    return mem # return the memory usage in MB

def writetxt(times, values, path='data.txt'):
    with open(path, 'w') as f:           # will overwrite existing file
        for i in range(times.size):
            f.write('%.16E\t%.16E\n' % (times[i], values[i]))

# load MESA data
data = np.loadtxt('input/eta_0.5/m.txt') # return (N, 2) array
mtimes = data[:, 0]                      # return only 1st col
masses = data[:, 1]                      # return only 2nd col
data = np.loadtxt('input/eta_0.5/r.txt')
rtimes = data[:, 0]
Rsuns = data[:, 1]                       # data in Rsun units
data = np.loadtxt('input/eta_0.5/l.txt')
ltimes = data[:, 0]
Lsuns = data[:, 1]                       # data in Lsun units

# conversions and precalculations
radii = np.zeros(Rsuns.size)        # convert Rsun to AU
for i,r in enumerate(Rsuns):
    radii[i] = r * 0.00465047       # 215 Rsun ~ 1 AU
watts = np.zeros(Lsuns.size)        # convert Lsun to W (MKS units)
for i,l in enumerate(Lsuns):
    watts[i] = l * 3.828e26         # IAU Resolution B3 conversion
lumins = np.zeros(watts.size)       # convert W to sim units
for i,w in enumerate(watts):
    lumins[i] = (w *((6.7e-12)**2)*(5e-31))/((3.2e-8)**3)
t_fs = np.zeros(lumins.size)        # precalculate t_f (Eq. 1)
for i,l in enumerate(lumins):
    t_fs[i] = np.cbrt(masses[i]*radii[i]**2/l)
taus = np.zeros(t_fs.size)          # precalc tau (Eq. 2)
G = 4*np.pi**2                      # units of AU, yr, and Msun
for i,t_f in enumerate(t_fs):
    taus[i] = 2.*radii[i]**3/G/masses[i]/t_f

# initialize sim and create Interpolator objects
timer_start = time.perf_counter()
sim, rebx, tides = makesim(init_a)
sim.ri_whfast.safe_mode = 0  # boost WHFast performance (advanced)
sim.ri_whfast.corrector = 11 # increase WHFast accuracy (advanced)
starmass = reboundx.Interpolator(rebx, mtimes, masses, 'spline')
starradius = reboundx.Interpolator(rebx, rtimes, radii, 'spline')
startau = reboundx.Interpolator(rebx, ltimes, taus, 'spline')

# update Sun's mass and radius accordingly
ps = sim.particles
ps[0].m = starmass.interpolate(rebx, t=T0)
ps[0].r = starradius.interpolate(rebx, t=T0)
ps[0].params["tctl_k2"] = 0.038 # ~ lambda_2, Schroder & Smith (2008)
ps[0].params["tctl_tau"] = startau.interpolate(rebx, t=T0)
ps[0].params["Omega"] = 0 # explicitly set to 0 (would be 0 by default)

# initialize main sim
tmax = 60e6                # max sim integration time
Nup = 10000                # 6-kyr param update interval
ts = np.linspace(0., tmax, Nup)
a = np.zeros(Nup)          # record semiaxis
mem_psutil = np.zeros(Nup) # mem usage tracking

try:
    for j,t in enumerate(ts):
        mem_psutil[i] = memory_usage_psutil()
        sim.move_to_com()   
        a[j] = ps[1].a                                   # record
        ps[0].m = starmass.interpolate(rebx, t=T0+sim.t) # update
        ps[0].r = starradius.interpolate(rebx, t=T0+sim.t)
        ps[0].params["tctl_tau"] = startau.interpolate(rebx,t=T0+sim.t)
        sim.ri_whfast.recalculate_coordinates_this_timestep = 1
        sim.integrator_synchronize()                     # corrector
        sim.integrate(t)                                 # til next Nup
except rebound.Collision as error:
    print(error)

# write semiaxis vs sim.t
makesubdir('output')  # create file output directory
fname = 'output/{:.2f}au.txt'.format(init_a)
writetxt(ts, a, path=fname)
# performance
max_mem = np.amax(mem_psutil)
timer_stop = time.perf_counter() 
runtime = timer_stop - timer_start
# cout
hh = runtime // 3600
remainder = runtime - hh*3600
mm = remainder // 60
ss = remainder - mm*60
print('________________________________')
print('Job Resource Usage Summary\n')
print('    Real Memory Used : %dmb'%(max_mem))
print('    Walltime Used    : %02d:%02d:%02d'%(hh, mm, ss))
print('________________________________')
