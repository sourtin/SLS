#Space Launch Simulator
from sim import sim

pp = [[31, 0], [60, 26], [100, 46.85], [160, 68.8], [240, 60.8],
		[390, 72.436], [540, 79.744], [620, 81.3], [700, 91.4]]

satv = sim()
satv.stage('payload',1939+45693,0,425)
satv.stage('LES',4042)
satv.stage(1,2239795,113,0,0,163)
satv.stage(1j,5206,0,0,0,193)
satv.stage(2,479964,79,0,0,550)
satv.stage(2j,3663,0,0,0,561)
satv.stage(3,119119,34,0,0)
satv.prog_pitch(0,pp)
satv.burn(1,.3,[13231,13348],304,49,1793760)
satv.burn(1,135.20,10638,304,39,281151)
satv.burn(1j,163,176.3,237,0,617,angle=10)
satv.burn(2,164,1224,424,3.*5,363053)
satv.burn(2,460.62,979,424,3.*4,38560)
satv.burn(2,500,728,434,3.*4,35111)
satv.burn(2j,549,13.38,231,0,53.524,angle=10)
satv.burn(3,552.2,213.416,431,3.,31400)
satv.burn(3,9856.2,204.91,445,3.,71068)

satv.run(712) # run 712s of simulation
dat=satv.analyse() # print energies


#Use this to generate a data file for gnuplot

s = '# distance/m, drag/N\n'
for x in dat:
	s += '%f %f\n' % (x[0],x[1])
s = s[:-1]

f = open('sim.dat', 'w')
f.write(s)
f.close()