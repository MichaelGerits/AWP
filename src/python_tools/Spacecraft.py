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
from pyatmos import expo
import progressbar
import ppigrf
import warnings


from datetime import datetime
plt.style.use( 'dark_background' )

# AWP libraries
import orbit_calculations as oc
import numerical_tools    as nt
import plotting_tools     as pt
import planetary_data     as pd
import spice_data         as sd
from tinyQuaternion import Quaternion

def null_config():
	'''
	this function defines a "base" configuration for the spacecraft class 
	such that not all parameters have to be specified every time
	'''
	return {
		'cb'             : pd.earth, #defines the central body
		'date0'          : '2020-04-01', #defines the start date
		'et0'            : None, #define the start ephemeris time (if no date)
		'frame'          : 'J2000', #defines the inertial reference frame of the simulation
		'dt'			 : 200, #defines the dt at which points are LOGGED
		'tspan'			 : '1',
		'orbit_state'    : [], 	#orbit defined by position and velocity vector quaternion and angular velocity
		'actuators'		 : [0., 0., 0., 0., 0., 0.], #implement a body axis centered linear and rotaional force
		'coes'           : [], #[semi-major axis(km) ,eccentricity ,inclination (deg) , ture anomaly, aop(deg), raan(deg)]
		'orbit_perts'    : {}, #defines a list of what pertubations are to be included
		'propagator'     : 'DOP853', #defines which ODE solver is used
		'atol'           : 1e-6, #absolute max error
		'rtol'           : 1e-6, #relative max error
		'stop_conditions': {}, #list of condistions to stop propagations
		'print_stop'     : True,
		'dense_output'   : False,
		'mass0'          : 1, #initial mass
		'inertia0'		 : np.array([[1., 0., 0.], #moment of inertia tensor
						   			 [0., 1., 0.], 
									 [0., 0., 1.],]),
		'drag_Cp'		 : np.zeros(3), #position of centre of pressure relative to body fixed axes
		'solarPress_Cp'	 : np.zeros(3), #position of centre of pressure relative to body fixed axes
		'output_dir'     : '.',
		'propagate'      : True #bool for propagation on initialisation
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
			self.config[ 'orbit_state' ] = np.append(self.config['orbit_state'], [np.cos(0), 0., 0., 0., 0., 0., 0.])

		#convert the optional amount of orbits to amount of seconds to simulate
		if type( self.config[ 'tspan' ] ) == str:
			self.config[ 'tspan' ] = float( self.config[ 'tspan'] ) *\
				oc.state2period( self.config[ 'orbit_state' ], self.cb[ 'mu' ] )

		if np.linalg.norm(self.config['actuators']) != 0:
			self.actuators = config['actuators']
		#allocates memory for the initial orbital state
		self.state0 = np.zeros( 14 )
		#adds initial pos (3d) and vel (3d) + a 4d quaternoin + 3d angular rate
		self.state0[ :13 ] = self.config[ 'orbit_state' ]
		
		#adds the initial mass
		self.state0[ 13 ] = self.config[ 'mass0' ]

		self.coes_calculated      = False
		self.latlons_calculated   = False
		self.altitudes_calculated = False
		self.ra_rp_calculated     = False
		self.eclipses_calculated  = False
		self.sun_dirs_calculated  = False

		#assigns specifics of the orbit integrator
		self.assign_stop_condition_functions()
		self.assign_orbit_perturbations_functions()
		#loads the spice kernels from the spice_data.py file
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
			'n_bodies': self.calc_n_bodies,
			'grav_grad': self.calc_grav_gradient,
			'atmos_drag': self.calc_atmos_drag,
			'mag_torque':self.calc_mag_torque,
			'solar_press': self.calc_solar_press
		}
		self.orbit_perts_funcs = []

		for key in self.config[ 'orbit_perts' ]: #adds whichever perts were initialized to be computed
			self.orbit_perts_funcs.append( 
				self.orbit_perts_funcs_map[ key ] )

	def load_spice_kernels( self ):
		print('\nLOADING SPICE KERNELS...\n')
		spice.furnsh( sd.leapseconds_kernel )
		spice.furnsh( sd.pck00010 )
		spice.furnsh( sd.de432 )
		self.spice_kernels_loaded = [ sd.leapseconds_kernel, sd.de432, sd.pck00010 ]

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
		alpha = np.zeros( 3 )
		q = Quaternion(q=state[6:10])
		_q = q.conjugate
		#get theposition vector of each body relative to the sc  
		for body in self.config[ 'orbit_perts' ][ 'n_bodies' ]:
			r_cb2body  = spice.spkgps( body[ 'SPICE_ID' ], et,
				self.config[ 'frame' ], self.cb[ 'SPICE_ID' ] )[ 0 ]
			r_sc2body = r_cb2body - state[ :3 ]

			#calc acceleration
			a += body[ 'mu' ] * (\
				 r_sc2body / nt.norm( r_sc2body ) ** 3 -\
				 r_cb2body / nt.norm( r_cb2body ) ** 3 )

			#calc attitude effect
			_r = _q.rotatePoint(r_sc2body) #rotate position vector to the body axis frame
			norm_r = nt.norm( _r )

			mult = np.matmul(self.config[ 'inertia0' ], _r)
			cross = np.cross(_r, mult)
			#same as grav gradient torque
			T = 3*self.cb[ 'mu' ]/(norm_r**5) * cross 

			alpha += np.matmul(np.linalg.inv(self.config[ 'inertia0' ]), T)
		return (a, alpha)

	def calc_J2( self, et, state ):
		'''
		calc the J2 effect on acceleration in km/s^2
		'''
		z2     = state[ 2 ] ** 2
		norm_r = nt.norm( state[ :3 ] )
		r2     = norm_r ** 2
		tx     = state[ 0 ] / norm_r * ( 5 * z2 / r2 - 1 )
		ty     = state[ 1 ] / norm_r * ( 5 * z2 / r2 - 1 )
		tz     = state[ 2 ] / norm_r * ( 5 * z2 / r2 - 3 )
		return (1.5 * self.cb[ 'J2' ] * self.cb[ 'mu' ] * self.cb[ 'radius' ] ** 2 / r2 ** 2 * np.array( [ tx, ty, tz ] ), 
		  		np.zeros(3))

	def calc_grav_gradient( self, et, state):
		'''
		calc the gravity gradiens effect.
		'''
		q = Quaternion(q=state[6:10])
		_q = q.conjugate 
		r =  state[ :3 ]  
		_r = _q.rotatePoint(r) #rotate position vector to the body axis frame
		norm_r = nt.norm( _r )

		mult = np.matmul(self.config[ 'inertia0' ], _r)
		cross = np.cross(_r, mult)
		T = 3*self.cb[ 'mu' ]/(norm_r**5) * cross #the units work themselves out

		alpha = np.matmul(np.linalg.inv(self.config[ 'inertia0' ]), T)

		return(np.zeros(3), alpha)
	
	def calc_atmos_drag(self, et, state):
		'''
		calculates the atmospheric drag by making use of an atmosphere model up to 1000km. 
		Anything above is consodered no drag.
		includes the relative airspeed using the rotational velocity of the atmosphere
		'''

		r = state[:3] 
		v = state[3:6] 
		mass = state[13]
		q = Quaternion(q=state[6:10])
		_q = q.conjugate

		alt = nt.norm(r) - self.cb[ 'radius' ]
		

		#if the spacecraft is to high for the atmosphere calc
		if alt > 1000:
			return (np.zeros(3), np.zeros(3))
		
		#load the exponential model at alt
		expo_geom = expo(alt)

		#get the density
		rho = expo_geom.rho
		CD = self.config['orbit_perts']['atmos_drag']['CD']
		A = self.config['orbit_perts']['atmos_drag']['A']
		
		#get the relative velocity with the atmosphere (atmosphere rotates)
		v_rel = v*1000-np.cross(self.cb['atm_rotation_vec'],r*1000) #change to si units

		#calc drag force
		force = -v_rel * nt.norm(v_rel) * 0.5 * rho * CD * A
		_force = _q.rotatePoint(force) #convert to body fixed frame to calc torque
		torque = np.cross(self.config[ 'drag_Cp' ], _force) #calc torque

		alpha = np.matmul(np.linalg.inv(self.config[ 'inertia0' ]), torque)
		a = force/mass/1000 #convert to km/s^2

		return (a, alpha)
	
	def calc_mag_torque(self, et, state): #TODO: speed up
		r = state[:3]
		q = Quaternion(q=state[6:10]) 
		D = self.config['orbit_perts']['mag_torque']['di_moment']
		

		#gets the longitude, lattitude and the radius
		rad, lon, lat = nt.cart2lat( np.array([r]),self.config[ 'frame' ], self.cb[ 'body_fixed_frame' ], np.array([et]))[0]
		alt = rad-self.cb['radius']

		date = spice.et2utc(et, format_str='ISOC', prec=0).replace('T', ' ')[:10].split('-')
		date = datetime(int(date[0]), int(date[1]), int(date[2]))

		with warnings.catch_warnings(action="ignore"): #gets the magnetic field
			#TODO:Newer magnetic field model will be needed due to deprecated function
			Be, Bn, Bu = ppigrf.igrf(lon, lat, alt, date) 
		ENU_vec = np.array([Be[0,0], Bn[0,0], Bu[0,0]])*1e-9

		lon_rad = lon*nt.r2d
		lat_rad = lat*nt.r2d

		R = np.array([
    				 [-np.sin(lon_rad),                	  np.cos(lon_rad),                	  0				 ], 
    				 [-np.sin(lat_rad) * np.cos(lon_rad), -np.sin(lat_rad) * np.sin(lon_rad), np.cos(lat_rad)],
    				 [np.cos(lat_rad) * np.cos(lon_rad),  np.cos(lat_rad) * np.sin(lon_rad),  np.sin(lat_rad)]
					 ])
		
		IAU_vec = np.matmul(R, ENU_vec) #transforms the ENU vec to IAU vector

		J200_vec = nt.frame_transform(np.array([IAU_vec]), self.cb[ 'body_fixed_frame' ], self.config[ 'frame' ], np.array([et]))[0] #tranforms to the J200 inertial frame

		body_vec = q.rotatePoint(J200_vec) #transforms to the body frame

		Torque = np.cross(D, body_vec)

		alpha = np.matmul(np.linalg.inv(self.config[ 'inertia0' ]), Torque)
		
		return (np.zeros(3), alpha)
	
	def calc_solar_press(self, et,  state): 
		r = state[:3] 
		mass = state[13]
		
		if oc.check_eclipse(et, r, self.config['cb'],self.config['orbit_perts']['solar_press']['eclipse_bodies'],  self.config['frame']) != -1: #checks if it goes trough an eclipse
			return(np.zeros(3), np.zeros(3))
		
		q = Quaternion(q=state[6:10])
		_q = q.conjugate #rotation vector to vonvert to body axis

		r_cb2body  = spice.spkgps( pd.sun[ 'SPICE_ID' ], et, self.config[ 'frame' ], self.cb[ 'SPICE_ID' ] )[ 0 ] #get the vector form central body to sun
		r_body2sc = r - r_cb2body
		d_sun = nt.norm(r_body2sc)
		

		a = (1+self.config['orbit_perts']['solar_press']['ref'])*pd.sun['G1']*self.config['orbit_perts']['solar_press']['A']*1e-6/(d_sun**3)/mass*r_body2sc
		torque = np.cross(self.config['solarPress_Cp'], _q.rotatePoint(1000*a*mass)) #rotate to body frame and multiply by 1000 to go to SI units (Nm)

		alpha = np.matmul(np.linalg.inv(self.config['inertia0']), torque) #rad/s^2
		return (a, alpha)
	
	def diffy_q( self, et, state ):
		'''
		initialises the original states and "derives" it for the ODE
		'''
		rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, w1, w2, w3, mass = state
		r         = np.array( [ rx, ry, rz ] )
		v         = np.array( [ vx, vy, vz ] )
		q 		  = nt.normed(np.array( [q0, q1, q2, q3] ))
		w 		  = np.array( [w1, w2, w3] ) #BODY AXIS rotational rate
		Force		  = np.array( self.config['actuators'][:3]) #body axis force (applied)
		Torque		  = np.array( self.config['actuators'][3:]) #body axis torque (applied)

		inertiaTens = self.config['inertia0']

		#adjusting torque and internal force to the body axis
		_q = Quaternion(q=np.array([q0, q1, q2, q3]))
		a_b = Force/ mass / 1000 #convert to km/s^2
		a_g = _q.rotatePoint(a_b) #total acceleration is in inertial frame, so we convert
		alpha = np.matmul(np.linalg.inv(inertiaTens), np.transpose(Torque))

		mass_dot  = 0.0 #time derivative of the mass
		state_dot = np.zeros( 14 )

		a = -r * self.cb[ 'mu' ] / nt.norm( r ) ** 3 + a_g #km/s^2 
		
		#add up all the pertubation effects
		for pert in self.orbit_perts_funcs:
			effect = pert( et, state )
			a += effect[0]
			alpha += effect[1]

			#print(effect)
				

		#the angular rate matrix to get the dot of the position quaternion
		w_matrix = np.array([[0, w1, w2, w3],
					   		 [-w1, 0, -w3, w2],
							 [-w2, w3, 0, -w1],
							 [-w3, -w2, w1, 0]])
		
		H = np.matmul(inertiaTens, w)

		#get the time derivative of the state
		state_dot[ :3  ] = v #r dot
		state_dot[ 3:6 ] = a #v dot
		state_dot[ 6:10 ] = np.transpose(0.5 * np.dot(q, w_matrix)) #q dot
		state_dot[ 10:13 ] = alpha - np.matmul(np.cross(w, H), np.linalg.inv(inertiaTens)) #consider rotational dynamics
		state_dot[ 13 ] = mass_dot

		self.bar.update((et - self.et0)/(self.config['dt']))

		return state_dot

	def propagate_orbit( self ):
		#sets a format for the loading bar
		self.widgets = [' [',
        	progressbar.Timer(format= 'Propagating Orbit: %(elapsed)s'),
        	'] ',
            progressbar.Bar('*'),' (',
            progressbar.ETA(), ') ',
        ]
		self.bar = progressbar.ProgressBar(max_value=(self.config['tspan']+1)/(self.config['dt']), widgets=self.widgets).start()

		self.ode_sol = solve_ivp(
			t_eval		 = np.arange(self.et0, self.et0 + self.config['tspan'], self.config[ 'dt']), #desides at which timesteps the values should be stored
			fun          = self.diffy_q, #array with vel and acc (time derivatives of state)
			t_span       = ( self.et0, self.et0 + self.config[ 'tspan' ] ),#time span (from et0 to et0+time span)
			y0           = self.state0, #initial state
			method       = self.config[ 'propagator' ], #what ODE solver is used
			events       = self.stop_condition_functions, #stopping conditions
			rtol         = self.config[ 'rtol' ], #relative accuracy lim
			atol         = self.config[ 'atol' ], #absolute accuracy lim
			#max_step 	 = 50, #TODO: fix the integrator/diffy_q to get rid of the errors
			dense_output = self.config[ 'dense_output' ] )

		self.states  = self.ode_sol.y.T
		self.ets     = self.ode_sol.t
		self.n_steps = self.states.shape[ 0 ]

	def calc_altitudes( self ):
		print('\nCalculating Altitudes...')
		self.altitudes = np.linalg.norm( self.states[ :, :3 ], axis = 1 ) -\
						self.cb[ 'radius' ]
		self.altitudes_calculated = True
	
	def calc_coes( self ):
		print( '\nCalculating COEs..' )
		self.coes = np.zeros( ( self.n_steps, 6 ) )

		for n in range( self.n_steps ):
			self.coes[ n, : ] = oc.state2coes(self.states[ n, :6 ], { 'mu': self.cb[ 'mu' ], 'et': self.ets[n] } )
			rdif = []
			vdif = []
			angdif = []
			rdif.append(nt.norm(self.states[n, :3]) - nt.norm(self.state0[:3]))
			vdif.append(nt.norm(self.states[n, 3:6]) - nt.norm(self.state0[3:6]))
			angdif.append(nt.vecs2angle(self.states[n, :3], self.states[n, 3:6]) - 90)
		print(f'max position delta: {max(rdif)} (km)')
		print(f'max velocity delta: {max(vdif)} (km/s)')
		print(f'max angular delta (from 90°): {max(angdif)} (°)')
			
		self.coes_rel        = self.coes[ : ] - self.coes[ 0, : ]
		self.coes_calculated = True

	def calc_sun_dirs( self ):
		'''
		calculates the directional angles between the body fixed frame and the direction of the sun in degrees
		'''
		print( '\nCalculating sun directions..' )
		self.sun_dirs = np.zeros( ( self.n_steps, 3 ) ) #allocates the memory
		axes = np.array([[1,0,0],
						 [0,1,0],  #defines the body axes
						 [0,0,1]])

		for n in range( self.n_steps ):
			q = Quaternion(q=self.states[n, 6:10])
			_q = q.conjugate #rotation vector to convert to body axis

			r_cb2sun  = spice.spkgps( pd.sun[ 'SPICE_ID' ], self.ets[n], self.config[ 'frame' ], self.cb[ 'SPICE_ID' ] )[ 0 ] #get the vector form central body to sun
			r_sc2sun = (r_cb2sun - self.states[n, :3])
			sun_dir = _q.rotatePoint(nt.normed(r_sc2sun)) 

			alpha = nt.vecs2angle(axes[0], sun_dir)
			beta = nt.vecs2angle(axes[1], sun_dir)
			gamma = nt.vecs2angle(axes[2], sun_dir)

			self.sun_dirs[n] = [alpha, beta, gamma]
		self.sun_dirs_calculated = True
	
	def calc_apoapses_periapses( self ):
		print("\nCalculating rp and ap....")
		if not self.coes_calculated:
			self.calc_coes()

		self.apoapses  = self.coes[ :, 0 ] * ( 1 + self.coes[ :, 1 ] )
		self.periapses = self.coes[ :, 0 ] * ( 1 - self.coes[ :, 1 ] )

		self.ra_rp_calculated = True

	def calc_latlons( self ):
		print("\nCalculating lattitudes and longitudes...")
		#reason we change the frame here is to account for the earth's rotation
		self.latlons = nt.cart2lat( self.states[ :, :3 ],
			self.config[ 'frame' ], self.cb[ 'body_fixed_frame' ], self.ets )
		self.latlons_calculated = True

	def calc_eclipses( self, bodies=[], method = 'either', v = False, vv = False ):
		self.eclipse_array = oc.calc_eclipse_array(self.ets, self.states[ :, :3 ], self.cb, bodies=bodies ,frame = self.config[ 'frame'] ) #TODO: add in multiple bodies
		self.eclipses = oc.find_eclipses( self.ets, self.eclipse_array,
			method, v, vv )
		self.eclipses_calculated = True

	def plot_eclipse_array( self, bodies=[],args = { 'show': True } ):
		if not self.eclipses_calculated:
			self.calc_eclipses(bodies=bodies)

		pt.plot_eclipse_array( self.ets, self.eclipse_array, args )

	def plot_3d( self, args, ani= True):
		pt.plot_orbits( [ self.states[ :, :3 ] ], args)
		if ani == True:
			pt.animate_orbits( self.n_steps,[ self.states[ :, :3 ] ],[self.states[:, 3:6]], [ self.states[:, 6:10] ], self.ets - len(self.ets) * [self.ets[0]], args )

	def plot_groundtracks( self, args = { 'show': True } ):
		if not self.latlons_calculated:
			self.calc_latlons()

		pt.plot_groundtracks( [ self.latlons ], args )

	def plot_coes( self, args = { 'show': True }, step = 1 ):
		if not self.coes_calculated:
			self.calc_coes()

		pt.plot_coes( self.ets[ ::step ], [ self.coes[ ::step ] ],
			args )

	def plot_sun_dirs(self, args = { 'show': True}, step = 1):
		if not self.sun_dirs_calculated:
			self.calc_sun_dirs()

		pt.plot_sun_dirs( self.ets[ ::step ], self.sun_dirs[ ::step ], args )

	def plot_states( self, args = { 'show': True } ):
		pt.plot_states( self.ets, self.states[ :, :13 ], args )

	def plot_altitudes( self, args = { 'show': True } ):
		if not self.altitudes_calculated:
			self.calc_altitudes()

		pt.plot_altitudes( self.ets, [ self.altitudes ], args )
