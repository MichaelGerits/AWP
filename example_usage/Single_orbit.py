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
#make sure to change this when running your own files
path.append("C:\\Users\\MichaelGerits\\OneDrive - Alchemich Farma BVBA\\Documenten\\coding\\python\\earospace_programs\\AWP_Attitude_addition\\src\\python_tools")

# AWP libraries
from Spacecraft import Spacecraft as SC
import spice_data as sd
import plotting_tools as pt
import numerical_tools as nt
import planetary_data as pd




if __name__ == '__main__':
	
	#coes = [semi-major axis(km) ,eccentricity ,inclination (deg) , ture anomaly, aop(deg), raan]
	#state = state values are in unit km and km/s, rad and rad/s

	#coes = [ pd.mercury[ 'radius' ] + 750, 0., 90, 0., 0., 0. ]
	#coes = [ pd.earth[ 'radius' ] + 400, 0., 51.6, 0., 0., 51.6 ] #sun synchronous orbit
	#coes = [ 42164, 0.24, 63.4, 0.0, 270, 100 ] #Tundra orbit
	#coes = [ 26600, 0.74, 63.4, 0.0, 270, 0.0 ] #Molniya orbit
	state = [-3037.482918000000,  5692.380440000000, -2137.245501000000, -5.26287783523500,  -0.74973955380600, 5.51613905818000, np.cos(0), 0., 0., np.sin(0), 0., 0., 0.]
	sc   = SC(
			{
			'cb' 		 : pd.earth,
			'date0'		 : '2024-10-04',
			#'coes'       : coes,
			'orbit_state': state,
			'actuators'	 : [0., 0., 0., 0., 0., 0.], #these are the forces and torques that act with respect to the body axis
			'mass0'		 : 472851.3,
			'inertia0'	 : np.array([[10276978., -1084837., 597098.],      
						   			 [-1084837., 31940398., -614081.], 
									 [597098., -614081., 40019058.]]),
			'drag_Cp'		 : np.array([0.00E+00, 3.82E-02, -2.56E+00]), #position of the Cp's in the attitude body fixed frame
			'solarPress_Cp'	 : np.array([0.00E+00, 3.82E-02, -2.56E+00]),
			'tspan'      : '100', #Tspan is either the amount or seconds. If it is a string,it is the amount of orbits
			'dt' : 100, #this decides at which points the integrator STORES points to be plotted
			'orbit_perts': {'J2': True, 
				   			'n_bodies': [pd.sun, pd.jupiter, pd.moon ],
							'grav_grad': True,
							'atmos_drag': {'CD': 5.20, 'A':1702.82},
							'solar_press': {'ref': 0, 'A': 7.5, 'eclipse_bodies': [pd.moon]},
							#'mag_torque':{'di_moment': np.array([1., 0., 0.])}
							}
			} )
	sc.plot_states(args = {'show': True, 'time_unit': 'hours'})
	sc.plot_coes(args = {'show': True, 'time_unit': 'days'})
	#sc.plot_sun_dirs(args = {'show': True, 'time_unit': 'days'})
	#sc.plot_eclipse_array( args = {'show': True, 'time_unit': 'days'})

	sc.plot_altitudes()
	#sc.plot_groundtracks( args = {'show': True, 'surface_body': 'earth'})
	sc.plot_3d(ani = False,  args = {'cb_radius': pd.earth['radius'], 'show': True, 'ani_name': 'orbit.gif', 'frames': None, 'showTime': True, 'fps': 5}) 	#frames decides how many frames that are stored are shown, none shows them all
	#sc.plot_perts()	



