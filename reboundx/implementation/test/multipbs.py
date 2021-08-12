import subprocess
import numpy as np

pwd = '/home/barons2/github/sabaronett/cherry-creek/reboundx/implementation/test'
aus = np.arange(1, 1.2, 0.1) # in AU

with open('processes.txt', 'w') as f:
    for au in aus:
        rc = ['qsub',
              '-j', 'oe',
              '-l', 'select=1:ncpus=1:mem=1gb',
              '-l', 'walltime=0:01:00',
              '-m', 'abe',
              '-M', 'barons2@unlv.nevada.edu',
              '-N', 'test_{:.1f}au'.format(au),
              '-q', 'small',
              '-S', '/home/barons2/.conda/envs/rebx-3.4.1/bin/python',
              pwd+'/args.py', '{:.1f}'.format(au)]
        complete = subprocess.run(rc)
        f.write('args = %s\n'%(complete.args))
        f.write('  returncode = %s\n'%(complete.returncode))
        f.write('  stdout = %s\n'%(complete.stdout))
        f.write('  stderr = %s\n\n'%(complete.stderr))
