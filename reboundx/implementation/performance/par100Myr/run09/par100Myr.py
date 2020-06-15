import time
import psutil
import os
import numpy as np
import rebound
import reboundx

# initialize constants
T0 = 12293.5e6        # Sun's age ~ 100 Myr pre-TRGB (sim start)
M0 = 0.98506175042481 # initial mass of star
names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars']
symbols = ['☉', '☿', '♀︎', '⊕', '♂︎']
timer_start = time.perf_counter()
# for sequential trials
Nloops = 1
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
    sim.add(m=M0, hash=names[0])
    sim.add(m=0.166e-6, a=0.39, hash=names[1])
    sim.add(m=2.45e-6, a=0.723, hash=names[2])
    sim.add(m=3.e-6, a=1., hash=names[3])
    sim.add(m=0.323e-6, a=1.524, hash=names[4])
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

for k in range(Nloops):
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
    a = np.zeros([Nout, sim.N])
    ts = np.linspace(0., 100.e6, Nout) # 100 Myr sim
    cp = 1                             # index of closest survivng planet
    emass = 0.                         # mass of engulfed planets
    proc_time, mem_psutil = np.zeros(Nout), np.zeros(Nout) # performance tracking

    for i, t in enumerate(ts):
        sim.integrate(t)
        d = ps[0] - ps[cp]             # componentwise difference to nearest planet
        r = np.sqrt(d.x**2 + d.y**2 + d.z**2)
        
        if r <= ps[0].r:               # nearest planet engulfed
            emass += ps[cp].m          # add engulfed planet mass
            ps[cp].m = 0               # zero planet mass and move to COM
            ps[cp].x, ps[cp].y, ps[cp].z = 0, 0, 0
            cp += 1                    # next closest surviving planet
            sim.dt = 0.1*ps[cp].P      # adjust timestep accordingly
            d = ps[0] - ps[cp]         # update distance to nearest survivor
            r = np.sqrt(d.x**2 + d.y**2 + d.z**2)
        
        # evolve Sun and recenter to COM
        ps[0].m = starmass.interpolate(rebx, t=T0+sim.t) + emass
        ps[0].r = starradius.interpolate(rebx, t=T0+sim.t)
        sim.move_to_com()
        
        # record values for post-sim plots
        mass[i] = sim.particles[0].m
        radius[i] = sim.particles[0].r
        for j in range(1, sim.N):
            a[i, j] = ps[j].a
        
        # update tidal parameter relative to nearest surviving planet
        ps[0].params["tctl_tau"] = starptau.interpolate(rebx, t=T0+sim.t)*r # Eq. 2
        # record current memory usage (MB)
        proc_time[i] = time.perf_counter() - timer_start
        mem_psutil[i] = memory_usage_psutil()
        
    writetxt(ts, mass, 'output/m_5Myr.txt')
    writetxt(ts, radius, 'output/r_5Myr.txt')
    for j in range(1, sim.N):
        fname = 'output/p' + str(j) + '_tides_5Myr.txt'
        writetxt(ts, a[:, j], path=fname)
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
    ptimes[k] = runtime
writetxt(np.arange(1, Nloops+1), ptimes, 'output/seqtimes.txt')
