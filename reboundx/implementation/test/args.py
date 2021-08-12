import sys

arg1 = sys.argv[1]

print('argv[1] =', arg1)
with open('{:s}au.txt'.format(arg1), 'w') as f:
    f.write('args = %s\n'%(arg1))
