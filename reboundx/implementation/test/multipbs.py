import subprocess
import numpy as np

pwd = '/home/barons2/github/sabaronett/cherry-creek/reboundx/implementation/test'
aus = np.arange(1, 1.2, 0.1) # in AU

for au in aus:
    rc = ['qsub',
            '-j', 'oe',
            '-l', 'select=1:ncpus=1:mem=1gb',
            '-l', 'walltime=0:01:00',
            '-m', 'abe',
            '-M', 'barons2@unlv.nevada.edu',
            '-N', 'test_{:.1f}au'.format(au),
            '-q', 'small',
            '--', pwd+'./test.sh', '{:.1f}'.format(au)]
    subprocess.run(rc)
