'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

'''

# Python standard libraries
from sys import path
import numpy as np
import spiceypy as spice


#adding in a new path for the custom files to work too
path.append("C:\\Users\\MichaelGerits\\SABInternship\\PythonProjectFiles\\AWP\\src\\python_tools")

# AWP libraries
from Spacecraft import Spacecraft as SC
from planetary_data import earth
import spice_data as sd
import plotting_tools as pt
import numerical_tools as nt
import planetary_data as pd




if __name__ == '__main__':
	
	#coes = [semi-major axis(km) ,eccentricity ,inclination (deg) , ture anomaly, aop(deg), raan]
	#state = state values are in unit km and km/s, rad and rad/s

	coes = [ earth[ 'radius' ] + 5000, 0., 30, 0., 10., 20. ]
	#coes = [ earth[ 'radius' ] + 1260, 0., 100.7, 0., 0., 90. ] #sun synchronous orbit
	#coes = [ 42164, 0.24, 63.4, 0.0, 270, 100 ] #Tundra orbit
	#coes = [ 26600, 0.74, 63.4, 0.0, 270, 0.0 ] #Molniya orbit
	#state = [earth[ 'radius' ] + 5000,  0., 0., 0.,  5.91874728, 0., np.cos(0), 0., 0., np.sin(0), 0., 0., 0.00]
	sc   = SC(
			{
			'cb' 		 : earth,
			'date0'		 : '2020-01-01',
			'coes'       : coes,
			#'orbit_state': state,
			'actuators'	 : [0., 0., 0., 0., 0., 0.], #these are the forces and torques that act with respect to the body axis
			'mass0'		 : 100.,
			'inertia0'	 : np.array([[0.1, 0., 0.],
						   			 [0., 1., 0.], 
									 [0., 0., 1.]]),
			'drag_Cp'		 : np.array([-1., 0., 0.]), #position of the Cp's in the attitude body fixed frame
			'solarPress_Cp'	 : np.array([-1., 0., 0.]),
			'tspan'      : '1', #Tspan is either the amount or seconds. If it is a string,it is the amount of orbits
			'dt' : 50, #this decides at which points the integrator STORES points to be plotted
			'orbit_perts': {#'J2': True, 
				   			#'n_bodies': [pd.moon, pd.sun, pd.jupiter, pd.saturn ],
							#'grav_grad': True,
							#'atmos_drag': {'CD': 2.2, 'A':10},
							#'solar_press': {'ref': 1, 'A': 10},
							'mag_torque':{'di_moment': np.array([1., 0., 0.])}
							}
			} )
	sc.plot_states(args = {'show': True, 'time_unit': 'days'})
	#sc.plot_coes(args = {'show': True, 'time_unit': 'days'})
	sc.plot_sun_dirs(args = {'show': True, 'time_unit': 'days'})
	#sc.plot_eclipse_array(args = {'show': True, 'time_unit': 'days'})
	sc.plot_altitudes()
	#sc.plot_groundtracks( args = {'show': True, 'surface_body': 'earth'})
	sc.plot_3d(ani = True,  args = { 'show': True, 'ani_name': 'orbit.gif', 'frames': None, 'showTime': True, 'fps': 5}) 	#frames decides how many frames that are stored are shown, none shows them all
	



