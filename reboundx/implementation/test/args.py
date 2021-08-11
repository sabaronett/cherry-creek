import sys

print('argv[1] =', sys.argv[1])
with open('{:.1f}au.txt', 'w') as f:
    f.write('args = %s\n'%(sys.argv[1]))
