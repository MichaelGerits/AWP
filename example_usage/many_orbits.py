'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Many orbits script
'''

# Python standard libraries
from sys import path
path.append( '../src/python_tools' )

# AWP libraries
from Spacecraft import Spacecraft as SC
from planetary_data import earth
import plotting_tools as pt

# 3rd party libraries
import numpy as np

aops   = np.arange( 0, 360, 30 )
inc   = 60
coes   = [ earth[ 'radius' ] + 10000, 0.05, 0.0, 0.0, 0.0, 0.0 ]
scs    = []
config = {
	'tspan': '1',

}

print( len( aops ))
max_steps = 99999999999999999999999999999
#TODO simulate simultaniosly to get rid of the problem
if __name__ == '__main__':
	for aop in aops:
		coes[ 2 ] = inc
		coes[ 5 ] = aop
		config= {
			'coes': coes,
			'tspan': '1',
			'propagate': True

		}
		sc = SC( config )
		if len(sc.ets) < max_steps:
			max_steps = len(sc.ets)
		scs.append( sc )

	rs = [ sc.states[ :, :3 ] for sc in scs ]
	vs = [ sc.states[ :, 3:6 ] for sc in scs ]
	quats = [ sc.states[ :, 6:10 ] for sc in scs ]
	#shape = (sc_amount, amount of simulated points, amount oflogged values)
	pt.animate_orbits( max_steps, rs, vs, quats, args = { 'show': True, 'ani_name': 'mult_orbit.gif', 'lb_axes': False })