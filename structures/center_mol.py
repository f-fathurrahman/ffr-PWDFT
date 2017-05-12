from __future__ import print_function
import ase.io
import sys
from ase.units import Bohr

Narg = len(sys.argv)
if Narg != 3:
    print('Error: exactly two arguments are needed')
    print('This Narg = ', Narg)
    exit()

filename = sys.argv[1]
L = float(sys.argv[2])*Bohr

atoms = ase.io.read(filename)
atoms.set_pbc([True,True,True])
atoms.set_cell([L,L,L])

atoms.center()
atoms.write(filename)

prefix = filename.replace('.xyz','')
atoms.write(prefix+'.xsf')
