'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Hello world of Spacecraft class
Two-body propagation with J2 perturbation for 100 periods
'''

# Python standard libraries
from sys import path
#adding in a new path for the custom files to work too
path.append("C:\\Users\\MichaelGerits\\SABInternship\\PythonProjectFiles\\AWP\\src\\python_tools")

# AWP libraries
from Spacecraft import Spacecraft as SC
from planetary_data import earth

if __name__ == '__main__':
	# coes = [semi-major axis(km) ,eccentricity ,inclination (deg) , raan(deg), aop(deg), ta]
	coes = [ earth[ 'radius' ] + 1000, 0.05, 30.0, 0.0, 0.0, 0.0 ]
	#coes = [ 26600, 0.64, 63.4, 0.0, 0.0, 0.0 ]
	sc   = SC(
			{
			'coes'       : coes,
			#Tspan is either the amount or seconds. 
			#If it is a string,it is the amount of orbits
			'tspan'      : '1', 
			'dt'         : 100.0,
			'orbit_perts': { 'J2': True }
			} )
	sc.plot_3d(ani = True)
