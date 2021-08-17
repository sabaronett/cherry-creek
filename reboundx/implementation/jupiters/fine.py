import subprocess
import numpy as np

pwd = '/home/barons2/github/sabaronett/cherry-creek/reboundx/implementation/jupiters/'
init_as = np.arange(1.60, 2.02, 0.02) # in au

for a in init_as:
    rc = ['qsub',
          '-j', 'oe',
          '-l', 'select=1:ncpus=1:mem=1gb',
          '-l', 'walltime=00:16:00',
          '-m', 'abe',
          '-M', 'barons2@unlv.nevada.edu',
          '-N', '{:.2f}au'.format(a),
          '-q', 'small',
          '--', pwd+'run.sh', '{:.2f}'.format(a)]
    subprocess.run(rc)
