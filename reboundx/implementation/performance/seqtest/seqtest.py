import time
import psutil
import os
import numpy as np
import rebound
import reboundx

# initialize constants
M0 = 0.8646552426064663 # initial mass of star
timer_start = time.perf_counter()
# for sequential trials
Nloops = 3
ptimes = np.zeros(Nloops)

def memory_usage_psutil():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20)
    return mem

def makesim():
    """
    Main REBOUND sim setup.
    """
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    sim.add(m=M0, hash='Sun')
    sim.add(m=3.e-6, a=1., hash='Earth')
    sim.collision = 'direct' # check if RGB Sun engulfs planet
    sim.integrator = 'whfast'
    sim.dt = 0.1*sim.particles[1].P
    sim.move_to_com()
    rebx = reboundx.Extras(sim)
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)
    return sim, rebx, tides

def engulf(sim):
    """
    Custom engulf (close encounter) handler.
    """
    ps = sim.particles
    ps[1].m = 0
    ps[1].x, ps[1].y, ps[1].z = 0, 0, 0

def writetxt(times, values, path='data.txt'):
    """
    Function to output time-value series data to a two-column text file.
        
    Parameters
    ----------
    times : numpy.ndarray
        Monotonic array of times to be written out to the first
        column of the data file.
    values : numpy.ndarray
        Array of values (particular stellar property) to be written out to the
        second column of the data file.
    path : str
        Path and filename of the data file to be outputted.
        Default path set to working directory and filename "data.txt"
    """
    with open(path, 'w') as f: # will overwrite existing file
        for i in range(times.size):
            f.write('%.16E\t%.16E\n' % (times[i], values[i]))

for j in range(Nloops):
    # load MESA data
    data = np.loadtxt('input/m.txt')    # return (N, 2) array
    mtimes = data[:, 0]                 # return only 1st col
    masses = data[:, 1]                 # return only 2nd col
    data = np.loadtxt('input/r.txt')
    rtimes = data[:, 0]
    Rsuns = data[:, 1]                  # data in Rsun units
    data = np.loadtxt('input/l.txt')
    ltimes = data[:, 0]
    Lsuns = data[:, 1]                  # data in Lsun units

    # conversions and precalculations
    radii = np.zeros(Rsuns.size)        # convert Rsun to AU
    for i, r in enumerate(Rsuns):
        radii[i] = r * 0.00465047       # 215 Rsun ~ 1 AU
        radii[i] += 0.15                # upshift to match Schroder
    watts = np.zeros(Lsuns.size)        # convert Lsun to W (MKS units)
    for i, l in enumerate(Lsuns):
        watts[i] = l * 3.828e26         # IAU Resolution B3 conversion
    lumins = np.zeros(watts.size)       # convert W to sim units
    for i, w in enumerate(watts):
        lumins[i] = (w * ((6.7e-12)**2) * (5e-31)) / ((3.2e-8)**3)
    t_fs = np.zeros(lumins.size)        # precalculate t_f (Eq. 1)
    for i, l in enumerate(lumins):
        t_fs[i] = np.cbrt(masses[i]*radii[i]**2/l) 
    pretaus = np.zeros(t_fs.size)       # precalc pretau (Eq. 2)
    G = 4*np.pi**2                      # units of AU, yrs and solar masses
    for i, t_f in enumerate(t_fs):
        pretaus[i] = 2.*radii[i]**3/G/masses[i]/t_f

    # initialize sim and create Interpolator objects
    sim, rebx, tides = makesim()
    starmass = reboundx.Interpolator(rebx, mtimes, masses, 'spline')
    starradius = reboundx.Interpolator(rebx, rtimes, radii, 'spline')
    starptau = reboundx.Interpolator(rebx, ltimes, pretaus, 'spline')

    # update Sun's mass and radius accordingly
    ps = sim.particles
    T0 = 1.23895e10 # Sun's age ~ 4 Myr pre-TRGB (sim start)
    ps[0].m = starmass.interpolate(rebx, t=T0)
    ps[0].r = starradius.interpolate(rebx, t=T0)
    ps[0].params["tctl_k1"] = 0.038 # ~ lambda_2, Schroder & Smith (2008)
    ps[0].params["tctl_tau"] = starptau.interpolate(rebx, t=T0)
    ps[0].params["Omega"] = 0 # explicitly set to 0 (would be 0 by default)
    sim.move_to_com()

    # main sim
    Nout = 1000
    mass = np.zeros(Nout)
    radius = np.zeros(Nout)
    a = np.zeros(Nout)
    ts = np.linspace(0., 4.e6, Nout)
    proc_time, mem_psutil = np.zeros(Nout), np.zeros(Nout) # performance tracking

    for i, t in enumerate(ts):
        try:
            sim.integrate(t)
        except:
            engulf(sim)
        ps[0].m = starmass.interpolate(rebx, t=T0+sim.t)
        ps[0].r = starradius.interpolate(rebx, t=T0+sim.t)
        sim.move_to_com() # lost mass had momentum, thus recenter to COM frame
        # record for post-sim plots
        mass[i] = sim.particles[0].m
        radius[i] = sim.particles[0].r
        a[i] = sim.particles[1].a
        # update tidal parameter
        d = ps[0] - ps[1] # componentwise difference between particles
        r = np.sqrt(d.x**2 + d.y**2 + d.z**2)
        ps[0].params["tctl_tau"] = starptau.interpolate(rebx, t=T0+sim.t)*r # Eq. 2
        # record current memory usage (MB)
        proc_time[i] = time.perf_counter() - timer_start
        mem_psutil[i] = memory_usage_psutil()
        
    writetxt(ts, mass, 'output/m.txt')
    writetxt(ts, radius, 'output/r.txt')
    writetxt(ts, a, 'output/a.txt')
    writetxt(proc_time, mem_psutil, 'output/mem.txt')

    # performance timer
    timer_stop = time.perf_counter() 
    runtime = timer_stop - timer_start
    h = runtime // 3600
    remainder = runtime - h*3600
    m = remainder // 60
    s = remainder - m*60
    if h != 0:
        print('Wall time: %dh %dmin %ds'%(h, m, s))
    elif m != 0:
        print('Wall time: %dmin %ds'%(m, s))
    else:
        print('Wall time: %ds'%(s))
    
    # for sequential trials only
    timer_start = time.perf_counter()
    ptimes[j] = runtime
writetxt(np.arange(1, Nloops+1), ptimes, 'output/seqtimes.txt')