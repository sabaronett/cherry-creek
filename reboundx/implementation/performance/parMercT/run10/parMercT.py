import time
import psutil
import os
import numpy as np
import rebound
import reboundx

# initialize constants
M0 = 0.8868357536545315  # initial mass of star
R0 = 0.33436215847158252 # initial radius of star
names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars']
symbols = ['☉', '☿', '♀︎', '⊕', '♂︎']
timer_start = time.perf_counter()
Nloops = 1 # for sequential trials
ptimes = np.zeros(Nloops)

def memory_usage_psutil():
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20)
    return mem # return the memory usage in MB

def makesim():
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    sim.add(m=M0, r=R0, hash=names[0])
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

def writetxt(times, values, path='data.txt'):
    with open(path, 'w') as f: # will overwrite existing file
        for i in range(times.size):
            f.write('%.16E\t%.16E\n' % (times[i], values[i]))

for k in range(Nloops):
    # initialize sim and create Interpolator objects
    sim, rebx, tides = makesim()

    # update Sun's mass and radius accordingly, and set tidal parameters
    ps = sim.particles
    sim.move_to_com()
    ps[0].params["tctl_k1"] = 0.038 # ~ lambda_2, Schroder & Smith (2008)
    ps[0].params["tctl_tau"] = 0.002643238893883989
    ps[0].params["Omega"] = 0 # explicitly set to 0 (would be 0 by default)

    # initialize main sim
    Nout = 1000
    mass = np.zeros(Nout)
    radius = np.zeros(Nout)
    a = np.zeros([Nout, sim.N])
    ts = np.linspace(0., 9.2e5, Nout)  # 920 Kyr sim
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
        
        sim.move_to_com()
        
        # record values for post-sim plots
        mass[i] = sim.particles[0].m
        radius[i] = sim.particles[0].r
        for j in range(1, sim.N):
            a[i, j] = ps[j].a
        
        # record current memory usage (MB)
        proc_time[i] = time.perf_counter() - timer_start
        mem_psutil[i] = memory_usage_psutil()
        
    writetxt(ts, mass, 'output/m.txt')
    writetxt(ts, radius, 'output/r.txt')
    for j in range(1, sim.N):
        fname = 'output/p' + str(j) + '_T.txt'
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
