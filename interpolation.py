from math import sqrt
# Turn '1 2; 3 4' into matrix
def matrix(s):
	M = [x.split() for x in s.split(';')]
	M = [[float(x) for x in y] for y in M]
	return M

# Crude matrix pretty-printing
def matprint(M):
	for row in M:
		print row

# GaussJordan Elimination Algorithm:
def gje(M):
	n = len(M)
	for p in range(n):
		for k in range(p+1, n):
			if abs(M[k][p]) > abs(M[p][k]):
				M[p], M[k] = M[k], M[p]
		k = M[p][p]
		if k == 0:
			print 'Errror: Matrix has no reduced form'
			raise
		M[p] = [x/k for x in M[p]]
		for i in range(n):
			if i != p:
				k = M[i][p]
				M[i] = [a - k * b for a, b in zip(M[i], M[p])]
	return M

# Solve for coefficients of a linear system
# Pass equations in the form [a, b, c, ... n, p]
#  where ax + by + cz + ... nw = p
# There must be n equations
# Returns [x, y, z, ..., w]
def solvesystem(*eqns):
	solution = gje(list(eqns))
	return [x[-1] for x in solution]
	

def diff(x):
	return [x[i+1]-x[i] for i in range(len(x)-1)]
	
def digitise(x,l):
	if x<l[0] or x>l[-1]:
		raise ValueError('argument {0} out of interpolation bounds'.format(x))
	for i,v in enumerate(l[:-1]):
		if x>=v and x<=l[i+1]:
			return i,x-v
	return len(l)-2,x-l[-1]

def interpolate2(x1, y1, x2, y2):
	dx = x2 - x1
	m = (y2 - y1) / dx
	c = (y1*x2 - y2*x1) / dx
	return lambda x: m * x + c

# Akima interpolation
def interpolate(data):
	x = [float(i[0]) for i in data]
	y = [float(i[1]) for i in data]
	n = len(data)
	
	if n < 2:
		return {
			0 : lambda x: 0,
			1 : lambda x: y[0],
		}[n]
	elif n == 2: return interpolate2(x[0],y[0], x[1],y[1])
	
	dx = diff(x)
	if any([i < 0 for i in dx]):
		raise ValueError('invalid independent variable values')
	
	m = [j / i for i, j in zip(dx, diff(y))]
	m1 = [3.*m[0] - 2.*m[1], 2.*m[0] - m[1]] + m + [2.*m[n-2] - m[n-3], 3.*m[n-2] - 2.*m[n-3]]
	dm = [abs(i) for i in diff(m1)]
	
	f1 = dm[2:n+2]; f2 = dm[0:n]
	f12 = [i + j for i, j in zip(f1, f2)]
	
	m2 = m1[1:n+1]; m3 = m1[2:n+2]; tv = 1e-9 * max(f12)
	b = [(i3 * i1 + i4 * i2) / i5 if i5 > tv else i1 for i1,i2,i3,i4,i5 in zip(m2,m3,f1,f2,f12)]
	
	m2 = b[0:n-1]; m3 = b[1:n]
	c = [(3.*i1 - 2.*i2 - i3)/i4 for i1,i2,i3,i4 in zip(m,m2,m3,dx)]
	d = [(i2 + i3 - 2.*i1)/i4**2 for i1,i2,i3,i4 in zip(m,m2,m3,dx)]
	
	def s(x_):
		i,dx = digitise(x_,x)
		return ((dx * d[i] + c[i]) * dx + b[i]) * dx + y[i]
	return s

def integrate(data):
	trapeze = lambda a,b: .5 * (b[0] - a[0]) * (a[1] + b[1])
	return sum([trapeze(a[0],a[1]) for a in zip(data[:-1],data[1:])])
	
def drange(start, stop, step):
	r = start
	while r < stop:
		yield r
		r += step

def hyp(*x):
	x=[x**2 for x in list(x)]
	return sqrt(sum(x))
		
def mncd(step=0.1):
     for mn in drange(0,3.5,step):
             print mn,",",CD(mn)