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

#adding in a new path for the custom files to work too
path.append("C:\\Users\\MichaelGerits\\SABInternship\\PythonProjectFiles\\AWP\\src\\python_tools")

# AWP libraries
from Spacecraft import Spacecraft as SC
from planetary_data import earth

if __name__ == '__main__':
	#coes = [semi-major axis(km) ,eccentricity ,inclination (deg) , raan(deg), aop(deg), ta]
	
	#coes = [ earth[ 'radius' ] + 25000, 0.1, 30.0, 0.0, 0.0, 0.0 ]
	state = [2.82402000e+04,  0., 0., 0., 3.41236793, 1.97013154, np.cos(np.pi/2), 0., 0., np.sin(np.pi/2), 0.]
	coes = [ 26600, 0.64, 63.4, 0.0, 0.0, 0.0 ]
	sc   = SC(
			{
			'coes'       : coes,
			#'orbit_state': state,
			#Tspan is either the amount or seconds. 
			#If it is a string,it is the amount of orbits
			'tspan'      : '1', 
			'orbit_perts': { 'J2': True }
			} )
	sc.plot_3d(ani = True,  args = { 'show': False, 'ani_name': 'orbit.gif', 'frames': None})
