#Need mayavi mayavi2 and ipython
#Exemple script : http://code.enthought.com/projects/mayavi/docs/development/html/mayavi/mlab.html#simple-scripting-with-mlab
#1. Lancer ipython  : ipython -pylab -wthread
#2. Executer script : %run my_script

import numpy
from enthought.mayavi.mlab import *
from pylab import *

DATAFILE  = loadtxt('config_3D/config_python_DIMER_10000020.dat')
x         = DATAFILE[:,0]
y         = DATAFILE[:,1]
z         = DATAFILE[:,2]
u         = DATAFILE[:,3]
v         = DATAFILE[:,4]
w         = DATAFILE[:,5]

from enthought.mayavi import mlab
s = mlab.quiver3d(x, y, z, u, v, w, line_width=2, scale_factor=1)
mlab.outline()
mlab.show()
