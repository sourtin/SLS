from math import sqrt, exp, cos, sin, pi, ceil, atan2
from interpolation import interpolate, interpolate2, diff, digitise, hyp, integrate
from atmosphere import atmosphere
from rocket import body, Earth, CD, constants
import copy

'''
coordinates:
 x: +ve in right direction
 y: +ve in up direction
 theta: +ve in clockwise direction,
 		0deg in right direction

computer:
	pitch
	area
	ae
	mass
	thrust
	data recorder
pitch:
	modes = 0 time
		    1 altitude
		    2 distance
	points = [[t/h/s,pitch],...]
stages:
	id
	mass
	time of pickup (assumed 0)
	time of ejection
	cross sectional area
	side area
burns:
	stage id
	start time
	quantity or time
	mdot or mdot range
	isp or isp range
	exhaust area
	angle (offset from current pitch)

NB if two stages are next to each other, assign the side areas such that the first one to
be ejected has area max(0,a1-a2), and the second has area a2
'''

def sign(x):
	if (x < 0): return -1
	elif (x == 0): return 0
	elif (x > 0): return 1
	else: return float('nan')
def interpolate3(val,t0,dt):
	if type(val) is type([]):
		b = (val[1] - val[0]) / dt
		a = val[0] - b * t0
		return lambda t: a + b * t
	return lambda t: val
def rot_vector(v,t):
	x = v['x']; y = v['y']
	c = cos(t); s = sin(t)
	x2 = y*s + x*c
	y2 = y*c - x*s
	return {'x':x2, 'y':y2}
def cart2polar(v):
	r = hyp(v['x'], v['y'])
	t = atan2(v['y'], v['x'])
	return {'r':r, 't':t}
def polar2cart(v):
	r = v['r']; t = v['t']
	x = r*cos(t); y = r*sin(t)
	return {'x':x, 'y':y}

class sim:
	def __init__(self,tres=.5,host=Earth,cd=CD):
		self.dt = tres		# time resolution
		self.host = host
		self.cd = cd		# drag coeff
		self.ppm = 0		# pitch mode
		self.ppp = [[0.],[0.]]  # pitch points
		self.data = []
		self.stages = {}; self.astages = {}
		self.burns = [];  self.aburns = []
		
		# state:
		self.t = 0
		self.s = 0
		self.pos = {'x':0, 'y':host.radius, 'alt':0}
		self.v = {'x':0, 'y':0}
		self.A = {'top':0, 'side':0, 'exh':0}
		self.m = 0
		self.burning = False
		self.pitch = 0
		self.reps = 0
	def drag(self):
		atm = self.host.atm; A = self.A; alt = self.pos['alt']
		if self.burning: dr_th = atm.pres(alt) * A['exh']
		else: dr_th = 0
		
		vh = self.v['x']; vv = self.v['y']; v = hyp(vh,vv)
		cc = .5 * self.cd(v) * atm.dens(alt)
		cp = cos(self.pitch); sp = sin(self.pitch)
		dt = cc*A['top']*(vv*cp+vh*sp)**2 + dr_th
		ds = cc*A['side']*(vv*sp-vh*cp)**2
		
		return {'x':(dt*sp + ds*cp)*sign(vh), 'y':(dt*cp - ds*sp)*sign(vv)}
	def prog_pitch(self, mode, points):
		self.ppm = mode
		x, y = zip(*points)
		x = [float(x) for x in x]
		y = [float(y) for y in y]
		dx = diff(x)
		if any([k < 0 for k in dx]):
			raise ValueError('invalid pitch program')
		self.ppp = [x, y]
	def update_pitch(self):
		if self.ppm == 0: t = self.t
		elif self.ppm == 1: t = self.pos['alt']
		elif self.ppm == 2: t = self.s
		else: t = 0
		
		x, y = self.ppp
		if t <= x[0]: return y[0] * pi / 180
		elif t >= x[-1]: return y[-1] * pi / 180
		j = digitise(t, x)[0]
		f = interpolate2(x[j],y[j], x[j+1],y[j+1])
		return f(t) * pi / 180 # radians
	def stage(self,id,mass,csa=0,side_area=0,t0=0,teject=float('inf')):
		s = {'m':mass, 't':[t0,teject], 'area':{'top':csa,'side':side_area}}
		self.stages[id] = s
	def burn(self,id,t0,mdot,isp,Ae,quantity,time=False,angle=0):
		if id not in self.stages:
			raise ValueError("No stage with id '{0}' exists".format(id))
		if quantity > self.stages[id]['m']:
			raise ValueError("Cannot burn {0}kg of fuel in stage '{1}'".format(quantity,id))
		if t0 < self.stages[id]['t'][0]:
			raise ValueError("Cannot begin burn of stage '{0}' before it is loaded".format(id))
		
		mda = mdot
		if type(mdot) is type([]):
			mda = .5 * (mdot[0] + mdot[1])
		if not time: time = quantity / mda
		mf = interpolate3(mdot, t0, time)
		sf = interpolate3(isp, t0, time)
		
		if t0+time > self.stages[id]['t'][1]:
			raise ValueError("Cannot continue burning '{0}' after it is ejected".format(id))
		
		b = {'id':id, 't':[t0,t0+time], 'mdot':mf, 'isp':sf, 'area':Ae, 'angle':angle}
		self.burns += [b]
	def update_state(self):
		t = self.t
		dels = []
		for k,v in self.stages.iteritems():
			if v['t'][0] <= t:
				self.astages[k] = v
				dels += [k]
		for k in dels: del self.stages[k]
		dels = []
		for i in xrange(len(self.burns)-1,-1,-1):
			v = self.burns[i]
			if v['t'][0] <= t:
				self.aburns += [v]
				del self.burns[i]
		for k,v in self.astages.iteritems():
			if v['t'][1] <= t:
				dels += [k]
		for k in dels: del self.astages[k]
		self.aburns = [v for v in self.aburns if v['t'][1] >= t]
		m = 0; csa = 0; sidea = 0
		for k,v in self.astages.iteritems():
			m += v['m']
			csa = max(csa, v['area']['top'])
			sidea += v['area']['side']
		A = sum(v['area'] for v in self.aburns)
		self.m = m
		self.A['top'] = csa
		self.A['side'] = sidea
		self.A['exh'] = A
		self.pitch = self.update_pitch()
	def step(self,thrust,dt):
		host = self.host; m = self.m; drag = self.drag()
		g = constants['G'] * host.mass / (host.radius + self.pos['alt'])**2
		vx = self.v['x']; vy = self.v['y']; vv = hyp(vx,vy)
		x = self.pos['x']; y = self.pos['y']; r = hyp(x,y)
		ax = (thrust['x'] - drag['x'])/m - g*x/r
		ay = (thrust['y'] - drag['y'])/m - g*y/r
		dsx = vx*dt + .5*ax*dt**2; dsy = vy*dt + .5*ay*dt**2
		x += dsx; y += dsy; alt = hyp(x,y) - host.radius
		#print m, thrust, drag, ax, ay, alt, ",", g, x, y, r
		if alt < 0:
			y = host.radius
			thrust=drag={'x':0, 'y':0}
			alt=x=av=ah=dsx=dsy=vx=vy=v=0
		self.s += hyp(dsx, dsy); self.pos['alt'] = alt
		self.v['x'] = vx+ax*dt; self.v['y'] = vy+ay*dt
		self.pos['x'] = x; self.pos['y'] = y
		return drag
	def full_report(self,datum,header=False,n=8):
		d = datum; a = d['areas']
		ni = "{:.%de},"%(n-6); si = "{:^%ds} "%n
		def print_cart(v):
			x = v['x']; y = v['y']; r = hyp(x,y)
			return (ni*3).format(x,y,r)
		if header: print (si*22).format("t","distance","x","y","alt","v_x","v_y","v","A_cs","A_side",
							"A_e","m0","dm","pitch","T_x","T_y","T","D_x","D_y","D","r","theta")[:-1]
		s  = (ni*2).format(d['t'],d['distance'])
		s += (ni*3).format(d['pos']['x'],d['pos']['y'],d['pos']['alt'])
		s += print_cart(d['vel'])
		s += (ni*3).format(a['top'],a['side'],a['exh'])
		s += (ni*3).format(d['mass'][0], d['mass'][1], d['pitch']*180/pi)
		s += print_cart(d['thrust'])
		s += print_cart(d['drag'])
		s += (ni*2).format(d['polar']['r'],d['polar']['t'])[:-1]
		print s
	def report(self,datum,header=False,n=14):
		d = datum; a = d['areas']
		ni = "{:.%de},"%(n-6); si = "{:^%ds} "%n
		def print_hyp(v):
			x = v['x']; y = v['y']; r = hyp(x,y)
			return ni.format(r)
			return (ni*3).format(x,y,r)
		def print_cart(v):
			x = v['x']; y = v['y']; r = hyp(x,y)
			return (ni*3).format(x,y,r)
		if header: print (si*13).format("t","distance","x","y","alt","v_x","v_y","v","m","pitch","T","D","F")[:-1]
		s  = (ni*2).format(d['t'],d['distance'])
		s += (ni*3).format(d['pos']['x'],d['pos']['y'],d['pos']['alt'])
		s += print_cart(d['vel'])
		s += (ni*2).format(d['mass'][0], d['pitch']*180/pi)
		T = d['thrust']; D = d['drag'];
		F = {'x':T['x']-D['x'],'y':T['y']-D['y']}
		s += print_hyp(d['thrust'])
		s += print_hyp(d['drag'])
		s += print_hyp(F)[:-1]
		print s
	def log(self):
		head = True; i = 0
		for datum in self.data[self.reps:]:
			if i%50==0:
				#if i%15==0: print; head = True
				self.report(datum, head)
				head = False;
			i+=1
		self.reps = len(self.data)
	def compute(self):
		t = self.t; dt = self.dt
		self.update_state()
		m0 = self.m
		thrust = {'x':0, 'y':0}
		self.burning = False
		for burn in self.aburns:
			id = burn['id']
			bt = t - burn['t'][0]
			if bt > dt: bt = dt
			dm = bt * burn['mdot'](t)
			Ft = dm * burn['isp'](t) * constants['g0'] / dt
			self.astages[id]['m'] -= dm
			self.m -= dm
			Ft = rot_vector({'x':0, 'y':Ft}, burn['angle'])
			thrust['x'] += Ft['x']
			thrust['y'] += Ft['y']
			self.burning = True
		thrust = rot_vector(thrust, self.pitch)
		drag = self.step(thrust, dt)
		datum = {
			't': t,
			'distance': self.s,
			'pos': self.pos,
			'vel': self.v,
			'areas': self.A,
			'mass': [m0,m0-self.m],
			'pitch': self.pitch,
			'thrust':thrust,
			'drag':drag,
			'polar':cart2polar(self.pos),
		}
		self.data += copy.deepcopy([datum])
		self.t += dt
	def run(self,t):
		while self.t <= t:
			self.compute()
		self.log()
	def analyse(self):
		dat=[[x['distance'], hyp(x['drag']['x'],x['drag']['y'])] for x in self.data]
		dat2=[[x['distance'], hyp(x['thrust']['x'],x['thrust']['y'])] for x in self.data]
		ED = integrate(copy.deepcopy(dat))
		EE = integrate(dat2)
		print "Drag Energy = %f GJ" % (ED/1e9)
		print "Thrust Energy = %f GJ" % (EE/1e9)
		return dat