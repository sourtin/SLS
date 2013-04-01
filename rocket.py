from interpolation import interpolate
from atmosphere import atmosphere

constants = {
	'G' : 6.67e-11,
	'c' : 299792458,
	'g0' : 9.80665,
}

class body:
	mass = 0
	radius = 0
	atm = None
	
	def __init__(self, mass, radius, atmcsv=None):
		self.mass = mass
		self.radius = radius
		if atmcsv:
			self.atm = atmosphere(atmcsv)

class mission:
	altitude = [0,0]
	velocity = 0
	host = None
	energy = 0
	
	def potential(self,alt):
		return -constants['G'] * self.host.mass / (alt + self.host.radius)
	
	# mission from alt_from to alt_to, attaining a velocity of vel or else orbital
	# launching from and orbiting round host, with extra energy E
	def __init__(self, alt_from, alt_to, vel=False, host=None, E=0):
		self.altitude = [alt_from, alt_to]
		if host: self.host = host
		self.velocity = constants['G'] * self.host.mass / (self.altitude[1] + self.host.radius)
		
		if vel: self.velocity = vel**2
		self.energy = .5*self.velocity
		self.energy += self.potential(self.altitude[1]) - self.potential(self.altitude[0])
		self.energy += E

# Defaults and Predefined Data
Sun = body(1.9891e30, 6.96342e8)
Earth = body(5.972e24, 6373319, 'earth.csv')
Moon = body(7.3477e22, 1.7371e6)

CD=interpolate([[0,0.3],[0.2,0.25],[0.5,0.2],[0.78,0.25],[0.84,0.30],[1.0,0.38],[1.2,0.5],
				[1.5,0.56],[2.0,0.52],[2.5,0.45],[3.0,0.43],[3.5,0.42],[4.0,0.42],[1e26,0.42]])

mission.host = Earth
LEO = mission(0, 1160000)
