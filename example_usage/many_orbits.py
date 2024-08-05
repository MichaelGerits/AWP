'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Many orbits script
'''

# Python standard libraries
from sys import path
path.append("C:\\Users\\MichaelGerits\\SABInternship\\PythonProjectFiles\\AWP\\src\\python_tools")
# AWP libraries
from Spacecraft import Spacecraft as SC
from planetary_data import earth
import plotting_tools as pt

# 3rd party libraries
import numpy as np

aops   = np.arange( 0, 360, 10 )
inc   = 60
coes   = [ earth[ 'radius' ] + 11000, 0.0, 0.0, 0.0, 0.0, 0.0 ]
scs    = []
config = {
	'tspan': '1',

}

print( len( aops ))
max_steps = 1e10
stp_index = 0
#TODO simulate simultaniosly to get rid of the problem
if __name__ == '__main__':
	for aop in aops:
		coes[ 2 ] = inc
		coes[ 5 ] = aop
		config= {
			'coes': coes,
			'actuators'	 : [0., 0., 1e-3, 0.1e-6, 0.1e-6, -0.5e-6], #these arethe forces and torques that act with respect to the body axis
			'mass0'		 : 100.,
			'inertia0'	 : np.array([[100., 0., 0.],
						   			 [0., 100., 0.], 
									 [0., 0., 100.],]),
			#Tspan is either the amount or seconds. 
			#If it is a string,it is the amount of orbits
			'tspan'      : '2', 
			#this decides at which points the integrator stores points to be plotted
			'dt' : 100,
			'orbit_perts': { 'J2': True },

		}
		sc = SC( config )
		if sc.n_steps < max_steps:
			max_steps = sc.n_steps
		scs.append( sc )

	rs = [ sc.states[ :, :3 ] for sc in scs ]
	vs = [ sc.states[ :, 3:6 ] for sc in scs ]
	quats = [ sc.states[ :, 6:10 ] for sc in scs ]
	times = scs[0].ets - len(scs[0].ets) * [scs[0].ets[0]] #TODO: try to break it
	#shape = (sc_amount, amount of simulated points, amount oflogged values)
	pt.animate_orbits( max_steps, rs, vs, quats, times, args = { 'show': False, 'ani_name': 'mult_orbit.gif', 'lb_axes': False, 'or_axes': True })