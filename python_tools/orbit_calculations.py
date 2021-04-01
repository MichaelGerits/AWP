'''
AWP
Orbit calculation tools
'''

# 3rd party libraries
import numpy as np
import math
import spiceypy as spice

# personal libraries
import numerical_tools as nt
import planetary_data  as pd

def esc_v( r, mu = pd.earth[ 'mu' ] ):
	'''
	Calculate escape velocity at given radial distance from body
	'''
	return math.sqrt( 2 * mu / r )
 
def state2coes( state, et = 0, mu = pd.earth[ 'mu' ], deg = True, print_results = False ):
	rp,e,i,raan,aop,ma,t0,mu,ta,a,T = spice.oscltx( state, et, mu )

	if deg:
		i    *= nt.r2d
		ta   *= nt.r2d
		aop  *= nt.r2d
		raan *= nt.r2d

	if print_results:
		print( 'a'   , a    )
		print( 'e'   , e    )
		print( 'i'   , i    )
		print( 'RAAN', raan )
		print( 'AOP' , aop  )
		print( 'TA'  , ta   )
		print()

	return [ a, e, i, ta, aop, raan ]

def state2period( state, mu = pd.earth['mu'] ):

	# specific mechanical energy
	epsilon = nt.norm( state[ 3:6 ] ) ** 2 / 2.0 - mu / nt.norm( state[ :3 ] )

	# semi major axis
	a = -mu / ( 2.0 * epsilon )

	# period
	return 2 * math.pi * math.sqrt( a ** 3 / mu )

def coes2state( coes, mu = pd.earth[ 'mu' ], deg = True ):
	a, e, i, ta, aop, raan = coes
	if deg:
		i    *= nt.d2r
		ta   *= nt.d2r
		aop  *= nt.d2r
		raan *= nt.d2r

	rp = a * ( 1 - e )

	return spice.conics( [ rp, e, i, raan, aop, ta, 0, mu], 0 )

def state2ap( state, mu = pd.earth[ 'mu' ] ):
	h       = nt.norm( np.cross( state[ :3 ], state[ 3: ] ) )
	epsilon = nt.norm( state[ 3: ] ) ** 2 / 2.0 - mu / nt.norm( state[ :3 ] )
	e       = math.sqrt( 2 * epsilon * h ** 2 / mu ** 2 + 1 )
	a       = h ** 2 / mu / ( 1 - e ** 2 )
	ra      = a * ( 1 + e )
	rp      = a * ( 1 - e )
	return  ra, rp