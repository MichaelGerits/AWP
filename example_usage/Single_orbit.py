'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Hello world of Spacecraft class
Two-body propagation with J2 perturbation for 100 periods
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




if __name__ == '__main__':
	#import the spice kernel to calc latlons in the earth_IAU frame and the atmosphere data
	spice.furnsh( sd.pck00010 )

	#coes = [semi-major axis(km) ,eccentricity ,inclination (deg) , ture anomaly, aop(deg), raan]
	#state = state values are in unit km and km/s, rad and rad/s

	coes = [ earth[ 'radius' ] + 5000, 0., 0, 0., 0., 0. ]
	state = [earth[ 'radius' ] + 2000,  0., 0., 0.,  5.91874728, 0., np.cos(np.pi/2), 0., 0., np.sin(np.pi/2), 0., 0., 0., 0.]
	#coes = [ 26600, 0.64, 63.4, 0.0, 0.0, 0.0 ] #Molniya
	sc   = SC(
			{
			#'coes'       : coes,
			'orbit_state': state,
			'actuators'	 : [0., 0., 0., 0., 0., 0.], #these are the forces and torques that act with respect to the body axis
			'mass0'		 : 100.,
			'inertia0'	 : np.array([[1, 0., 0.],
						   			 [0., 1, 0.], 
									 [0., 0., 1.],]),
			#position of the Cp's in the attitude body fixed frame
			'drag_Cp'		 : np.array([0., 5., 0.]),
			'solarPress_Cp'	 : np.array([1., 0., 0.]),
			#Tspan is either the amount or seconds. 
			#If it is a string,it is the amount of orbits
			'tspan'      : '5', 
			#this decides at which points the integrator STORES points to be plotted
			'dt' : 200,
			'orbit_perts': {'atmos_drag': {'CD': 2.2, 'A': 2} }
			} )
	#frames decides how many frames that are stored are shown, none shows them all
	print(sc.state0)
	sc.plot_states()
	sc.plot_3d(ani = True,  args = { 'show': True, 'ani_name': 'orbit.gif', 'frames': None, 'showTime': True, 'fps': 5})
	#sc.plot_coes()
	#sc.calc_latlons()
	#pt.plot_groundtracks([sc.latlons], {'show': True})
	



