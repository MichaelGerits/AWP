'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

SPICE data filepaths
'''

import os

base_dir = os.path.join(
	os.path.dirname( os.path.realpath( __file__ ) ),
	os.path.join( '..', '..', 'data', 'spice' ) 
	)

leapseconds_kernel = os.path.join( base_dir, 'lsk/naif0012.tls' ) #kernel for the leapseconds for ephemeris time
de432              = os.path.join( base_dir, 'spk/de432s.bsp'   ) #kernel for the ephemeris positions of the barycentre of diff bodies
de421              = os.path.join( base_dir, 'spk/de421.bsp'   ) #kernel for the ephemeris positions of bodies that interacted with Juno
pck00010           = os.path.join( base_dir, 'pck/pck00010.tpc' ) #kernel for planetary constants

