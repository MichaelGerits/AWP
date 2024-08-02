'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Spacecraft class definition
'''

# Python standard libraries
import os
import math as m

# 3rd party libraries
from scipy.integrate import solve_ivp
import spiceypy          as spice
import numpy             as np
import matplotlib.pyplot as plt
plt.style.use( 'dark_background' )

# AWP libraries
import orbit_calculations as oc
import numerical_tools    as nt
import plotting_tools     as pt
import planetary_data     as pd
import spice_data         as sd

def null_config():
	return {
		'cb'             : pd.earth,
		'date0'          : '2021-04-01',
		'et0'            : None,
		'frame'          : 'J2000',
		'dt'			 : 100,
		'orbit_state'    : [], 	#orbit defined by position and velocity vector quaternion and angular velocity
		'coes'           : [], #[semi-major axis(km) ,eccentricity ,inclination (deg) , , ,]
		'orbit_perts'    : {}, #defines a list of what pertubations are to be included
		'propagator'     : 'LSODA', #defines which ODE solver is used
		'atol'           : 1e-6, #absolute max error
		'rtol'           : 1e-6, #relative max error
		'stop_conditions': {}, #list of condistions to stop propagations
		'print_stop'     : True,
		'dense_output'   : False,
		'mass0'          : 0,
		'output_dir'     : '.',
		'propagate'      : True
	}

class Spacecraft:

	def __init__( self, config ):
		self.config = null_config()
		for key in config.keys():
			self.config[ key ] = config[ key ]

		self.orbit_perts = self.config[ 'orbit_perts' ]
		self.cb          = self.config[ 'cb' ]

		if self.config[ 'coes' ]:
			# 'oc' is for the orbit calculations lib
			self.config[ 'orbit_state' ] = oc.coes2state( self.config[ 'coes' ], mu = self.config[ 'cb' ][ 'mu' ] )
			#adds 0 rotation quaternion to the orbital state
			self.config[ 'orbit_state' ] = np.append(self.config['orbit_state'], np.zeros(7))

		#convert the optional amount of orbits to amount of seconds to simulate
		if type( self.config[ 'tspan' ] ) == str:
			self.config[ 'tspan' ] = float( self.config[ 'tspan'] ) *\
				oc.state2period( self.config[ 'orbit_state' ], self.cb[ 'mu' ] )

		#allocates memory for the initial orbital state
		self.state0 = np.zeros( 14 )
		#adds initial pos (3d) and vel (3d) + a 4d quaternoin + 3d angular rate
		self.state0[ :13 ] = self.config[ 'orbit_state' ][ :13 ]
		
		#adds the initial mass
		self.state0[ 13 ] = self.config[ 'mass0' ]

		self.coes_calculated      = False
		self.latlons_calculated   = False
		self.altitudes_calculated = False
		self.ra_rp_calculated     = False
		self.eclipses_calculated  = False

		self.assign_stop_condition_functions()
		self.assign_orbit_perturbations_functions()
		self.load_spice_kernels()

		if self.config[ 'propagate' ]:
			self.propagate_orbit()

	def assign_stop_condition_functions( self ):
		'''
		The stop conditions methods are passed into scipy's solve_ivp
		function as events. The methods can have 2 attributes added to
		them that the solve_ivp function will use to stop propagation
		if the method returns a 0 crossing. In order for the propagation to
		stop, the "terminal" attribute must be set to True.
		For 0 crossings, the "direction" attribute will dictate which direction
		will trigger the stop condition.
		Positive    --> negative to positive 0 crossing
		Negative    --> positive to negative 0 crossing
		0 (default) --> both
		'''

		self.check_min_alt.__func__.direction   = -1
		self.check_max_alt.__func__.direction   =  1
		self.check_enter_SOI.__func__.direction = -1

		self.check_min_alt.__func__.terminal = True
		self.stop_condition_functions        = [ self.check_min_alt ]

		if 'min_alt' not in self.config[ 'stop_conditions' ].keys():
			self.config[ 'stop_conditions' ][ 'min_alt' ] =\
				self.cb[ 'deorbit_altitude' ]

		self.stop_conditions_map = {
			'min_alt'  : self.check_min_alt,
			'max_alt'  : self.check_max_alt,
			'enter_SOI': self.check_enter_SOI
			}

		for key in self.config[ 'stop_conditions' ].keys():
			method                   = self.stop_conditions_map[ key ]
			method.__func__.terminal = True
			self.stop_condition_functions.append( method )

	def assign_orbit_perturbations_functions( self ):
	
		self.orbit_perts_funcs_map = {
			'J2'      : self.calc_J2,
			'n_bodies': self.calc_n_bodies
		}
		self.orbit_perts_funcs = []

		for key in self.config[ 'orbit_perts' ]:
			self.orbit_perts_funcs.append( 
				self.orbit_perts_funcs_map[ key ] )

	def load_spice_kernels( self ):
		spice.furnsh( sd.leapseconds_kernel )
		self.spice_kernels_loaded = [ sd.leapseconds_kernel ]

		if self.config[ 'et0' ] is not None:
			self.et0 = self.config[ 'et0' ]
		else:
			self.et0 = spice.str2et( self.config[ 'date0' ] )

	def check_min_alt( self, et, state ):
		return nt.norm( state[ :3 ] ) -\
			      self.cb[ 'radius' ] -\
			      self.config[ 'stop_conditions' ][ 'min_alt' ]

	def check_max_alt( self, et, state ):
		return nt.norm( state[ :3 ] ) -\
			      self.cb[ 'radius' ] -\
			      self.config[ 'stop_conditions' ][ 'max_alt' ]

	def check_enter_SOI( self, et, state ):
		body      = self.config[ 'stop_conditions' ][ 'enter_SOI' ]
		r_cb2body = spice.spkgps( body[ 'SPICE_ID' ], et,
						self.config[ 'frame' ], self.cb[ 'SPICE_ID' ] )[ 0 ]
		r_sc2body = r_cb2body - state[ :3 ]

		return nt.norm( r_sc2body ) - body[ 'SOI' ]

	def print_stop_condition( self, parameter ):
		print( f'Spacecraft has reached {parameter}.' )

	def calc_n_bodies( self, et, state ):
		a = np.zeros( 3 )
		for body in self.config[ 'orbit_perts' ][ 'n_bodies' ]:
			r_cb2body  = spice.spkgps( body[ 'SPICE_ID' ], et,
				self.config[ 'frame' ], self.cb[ 'SPICE_ID' ] )[ 0 ]
			r_sc2body = r_cb2body - state[ :3 ]

			a += body[ 'mu' ] * (\
				 r_sc2body / nt.norm( r_sc2body ) ** 3 -\
				 r_cb2body / nt.norm( r_cb2body ) ** 3 )
		return a

	def calc_J2( self, et, state ):
		z2     = state[ 2 ] ** 2
		norm_r = nt.norm( state[ :3 ] )
		r2     = norm_r ** 2
		tx     = state[ 0 ] / norm_r * ( 5 * z2 / r2 - 1 )
		ty     = state[ 1 ] / norm_r * ( 5 * z2 / r2 - 1 )
		tz     = state[ 2 ] / norm_r * ( 5 * z2 / r2 - 3 )
		return 1.5 * self.cb[ 'J2' ] * self.cb[ 'mu' ] *\
			   self.cb[ 'radius' ] ** 2 \
			 / r2 ** 2 * np.array( [ tx, ty, tz ] )

	def diffy_q( self, et, state ):
		'''
		initialises the original states and "derives" it for the ODE
		'''
		rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, w1, w2, w3, mass = state
		r         = np.array( [ rx, ry, rz ] )
		v         = np.array( [ vx, vy, vz ] )
		q 		  = np.array( [q0, q1,q2, q3] )
		w 		  = np.array( [w1, w2, w3] ) #BODY AXIS rotational rate
		mass_dot  = 0.0 #time derivative of the mass
		state_dot = np.zeros( 14 )
		et       += self.et0

		#TODO: add angular accelerations + petrubutions
		a = -r * self.cb[ 'mu' ] / nt.norm( r ) ** 3

		for pert in self.orbit_perts_funcs:
			a += pert( et, state )

		#the angular rate matrix to get the dot of the position quaternion
		w_matrix = np.array([[0, w1, w2, w3],
					   		 [-w1, 0, -w3, w2],
							 [-w2, w3, 0, -w1],
							 [-w3, -w2, w1, 0]])

		state_dot[ :3  ] = v #r dot
		state_dot[ 3:6 ] = a #v dot
		state_dot[ 6:10 ] = np.transpose(0.5 * np.dot(np.transpose(q), w_matrix)) #q dot
		state_dot[ 10:13 ] = np.zeros(3) #omega dot #TODO change this to be the angular acceleration
		state_dot[ 13 ] = mass_dot
		return state_dot

	def propagate_orbit( self ):
		print( 'Propagating orbit..' )

		self.ode_sol = solve_ivp(
			t_eval		 = np.arange(self.et0, self.et0 + self.config['tspan'], self.config[ 'dt']), #desides at which timesteps the values should be stored
			fun          = self.diffy_q, #array with vel and acc (time derivatives of state)
			t_span       = ( self.et0, self.et0 + self.config[ 'tspan' ] ),#time span (from et0 to et0+time span)
			y0           = self.state0, #initial state
			method       = self.config[ 'propagator' ], #what ODE solver is used
			events       = self.stop_condition_functions, #stopping conditions
			rtol         = self.config[ 'rtol' ], #relative accuracy lim
			atol         = self.config[ 'atol' ], #absolute accuracy lim
			dense_output = self.config[ 'dense_output' ] )

		self.states  = self.ode_sol.y.T
		self.ets     = self.ode_sol.t
		self.n_steps = self.states.shape[ 0 ]

	def calc_altitudes( self ):
		self.altitudes = np.linalg.norm( self.states[ :, :3 ], axis = 1 ) -\
						self.cb[ 'radius' ]
		self.altitudes_calculated = True

	def calc_coes( self ):
		print( 'Calculating COEs..' )
		self.coes = np.zeros( ( self.n_steps, 6 ) )

		for n in range( self.n_steps ):
			self.coes[ n, : ] = oc.state2coes( 
				self.states[ n, :6 ], { 'mu': self.cb[ 'mu' ] } )
			
		self.coes_rel        = self.coes[ : ] - self.coes[ 0, : ]
		self.coes_calculated = True

	def calc_apoapses_periapses( self ):
		if not self.coes_calculated:
			self.calc_coes()

		self.apoapses  = self.coes[ :, 0 ] * ( 1 + self.coes[ :, 1 ] )
		self.periapses = self.coes[ :, 0 ] * ( 1 - self.coes[ :, 1 ] )

		self.ra_rp_calculated = True

	def calc_latlons( self ):
		self.latlons = nt.cart2lat( self.states[ :, :3 ],
			self.config[ 'frame' ], self.cb[ 'body_fixed_frame' ], self.ets )
		self.latlons_calculated = True

	def calc_eclipses( self, method = 'either', v = False, vv = False ):
		self.eclipse_array = oc.calc_eclipse_array(
			self.ets, self.states[ :, :3 ],
			self.cb, self.config[ 'frame'] )
		self.eclipses = oc.find_eclipses( self.ets, self.eclipse_array,
			method, v, vv )
		self.eclipses_calculated = True

	def plot_eclipse_array( self, args = { 'show': True } ):
		if not self.eclipses_calculated:
			self.calc_eclipse_array()

		pt.plot_eclipse_array( self.ets, self.eclipse_array, args )

	def plot_3d( self, args, ani= True):
		if ani == True:
			pt.animate_orbits( self.n_steps,[ self.states[ :, :3 ] ],[self.states[:, 3:6]], [ self.states[:, 6:10] ], args )
		pt.plot_orbits( [ self.states[ :, :3 ] ], args )

	def plot_groundtracks( self, args = { 'show': True } ):
		if not self.latlons_calculated:
			self.calc_latlons()

		pt.plot_groundtracks( [ self.latlons ], args )

	def plot_coes( self, args = { 'show': True }, step = 1 ):
		if not self.coes_calculated:
			self.calc_coes()

		pt.plot_coes( self.ets[ ::step ], [ self.coes[ ::step ] ],
			args )

	def plot_states( self, args = { 'show': True } ):
		pt.plot_states( self.ets, self.states[ :, :6 ], args )

	def plot_altitudes( self, args = { 'show': True } ):
		if not self.altitudes_calculated:
			self.calc_altitudes()

		pt.plot_altitudes( self.ets, [ self.altitudes ], args )
