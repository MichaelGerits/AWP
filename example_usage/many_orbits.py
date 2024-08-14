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

max_steps = np.inf
stp_index = 0

coes = []
scs = []
coes.append( [ 26600, 0.74, 63.4, 0.0, 270, 0.0 ] )#Molniya orbit
coes.append( [ 42164, 0.24, 63.4, 0.0, 270, 100 ] )#Tundra orbit
coes.append( [ 42164.0, 0,  0.,  0, 0, 100 ] )#GEO orbit

if __name__ == '__main__':
	for i in coes: #itterates over all the spacecraft
		config= {
			'date0'		 : '2020-01-01',
			'coes'       : i,
			'actuators'	 : [0., 0., 0., 0., 0., 0.], #these are the forces and torques that act with respect to the body axis
			'mass0'		 : 100.,
			'inertia0'	 : np.array([[0.1, 0., 0.],
						   			 [0., 1., 0.], 
									 [0., 0., 1.]]),
			'drag_Cp'		 : np.array([-1., 0., 0.]), #position of the Cp's in the attitude body fixed frame
			'solarPress_Cp'	 : np.array([-1., 0., 0.]),
			'tspan'      : '1', #Tspan is either the amount or seconds. If it is a string,it is the amount of orbits
			'dt' : 500, #this decides at which points the integrator STORES points to be plotted
			'orbit_perts': {'J2': True, 
				   			#'n_bodies': [pd.moon, pd.jupiter, pd.saturn, pd.sun],
							'grav_grad': True,
							#'atmos_drag': {'CD': 2.2, 'A':10},
							#'solar_press': {'ref': 1, 'A': 10}
							}
			}
		sc = SC( config ) #initializes the spacecraft
		sc.calc_latlons() #calculates the longitude and lattitude positions
		#checks for the orbit with the smalles amount of steps, so during animation, no errors are given
		if sc.n_steps < max_steps:
			max_steps = sc.n_steps
			times = sc.ets - (len(sc.ets) * [sc.et0])
		scs.append( sc )

	#shape = (sc_amount, amount of simulated points, logged values)
	rs = [ sc.states[ :, :3 ] for sc in scs ]
	vs = [ sc.states[ :, 3:6 ] for sc in scs ]
	quats = [ sc.states[ :, 6:10 ] for sc in scs ]
	latlons = [ sc.latlons for sc in scs ]
	
	
	#pt.plot_groundtracks (latlons, args = {'show': True, 'grid': False, 'labels': [ 'Molniya', 'Tundra', 'Geosynchronous' ],})
	pt.plot_orbits(rs, {'show': True, 'labels': [ 'Molniya', 'Tundra', 'Geosynchronous' ]})
	pt.animate_orbits( max_steps, rs, vs, quats, times, args = { 'ani_name': 'mult_orbit.gif', 'lb_axes': True, 'or_axes': False, 'labels': [ 'Molniya', 'Tundra', 'Geosynchronous' ], 'frames': None })