import csv
from math import sqrt
from interpolation import interpolate

class atmosphere:
	temp = 0		# temperature
	dens = 0		# density
	pres = 0		# pressure
	molm = 0		# molar mass
	sphr = 0		# specific heat ratio
	
	def __init__(self, path):
		data = []
		with open(path, 'rb') as f:
			r = csv.reader(f)
			for row in r:
				data += [row]
				#[float(x) for x in row]
		data += [[4.4e26, 2.7, 0, 0, 0, 1.46]] #universe
		alt,temp,dens,pres,molm,sphr = zip(*data)
		
		self.temp = interpolate(zip(alt, temp))
		self.dens = interpolate(zip(alt, dens))
		self.pres = interpolate(zip(alt, pres))
		self.molm = interpolate(zip(alt, molm))
		self.sphr = interpolate(zip(alt, sphr))
	
	def cair(self,alt):
		return sqrt(8314.462 * self.sphr(alt) * self.temp(alt) / self.molm(alt)) 