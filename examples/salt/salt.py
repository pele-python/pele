import potentials.salt as salt
import numpy as np
import basinhopping as bh
import take_step.adaptive_step as adaptive_step
import take_step.random_displacement as ts
import storage.savenlowest
        
natoms = 196

# create PatchyParticle potential
pot = salt.salt()    

tmp=np.zeros(4*3+6)
tmp[-6:-3] = 10000.
tmp[3]=0.5
tmp[7]=0.5
tmp[11]=0.5

gx=np.arange(1,20,0.01)
gy=gx.copy()
import pylab as pl
for i in xrange(gx.size):
    tmp[3]=gx[i]/10000.
    gy[i] = pot.getEnergy(tmp)
pl.plot(gx, gy)
tmp[3] = 0.5
for i in xrange(gx.size):
    tmp[7]=gx[i]/10000.
    gy[i] = pot.getEnergy(tmp)
pl.plot(gx, gy)
tmp[7]=tmp[11]=0.
tmp[0]=0.25
for i in xrange(gx.size):
    tmp[9]=gx[i]/10000.
    gy[i] = pot.getEnergy(tmp)
pl.plot(gx, gy)
pl.ylim((-17,17))
pl.show()

# load coords array
x = np.random.random(3*natoms+6)
x[-6:-3] = 22
x[-3:]=0.

#for i in xrange(0,4):
#    for j in xrange(0,4):
#        for k in xrange(0,4):
#            x[(i*4*4 + j*4 +k)*3+0] = 1.0*i/4.
#            x[(i*4*4 + j*4 +k)*3+1] = 1.0*j/4.
#            x[(i*4*4 + j*4 +k)*3+2] = 1.0*k/4.


print "start"
#e,g = pot.getEnergyGradient(x)
#print g[1:12]
#print "numerical"
#gn = pot.NumericalDerivative(x, 1e-8)
#print gn[1:12]
#print "difference"
#print (g-gn)[1:12]
#exit()

# use adaptive step size, 0.3 start, acceptance rate 0.5, adjust every 20
manstep = adaptive_step.manageStepSize (.1, 0.3, 20)    
step = ts.takeStep( getStep=manstep.getStepSize)

# store the lowest 10 minma
minima = storage.savenlowest.SaveN(nsave=10)

# start a basin hopping run
opt = bh.BasinHopping(x, pot,                      
                      temperature=1.,
                      takeStep=step.takeStep,
                      event_after_step=[manstep.insertStepWrapper],
                      storage = minima.insert
                      )

# do 100 mc steps
opt.run(20)

#print minima.data[0][1]
#print pot.getEnergy(minima.data[0][1])
#print pot.NumericalDerivative(minima.data[0][1], 1e-8)

a = minima.data[0][1].copy()
pot.toReal(a)
print a
# save the minima
with open("pylowest", "w") as fout:
  for i in minima.data:
      fout.write( str(natoms) + "\n")
      fout.write( "Energy = " + str(i[0]) + "\n")
      tmp=i[1].copy()
      pot.toReal(tmp)
      for atom in tmp.reshape(x.size/3, 3)[0:natoms/2]:
          fout.write("A " + str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + "\n")
      for atom in tmp.reshape(x.size/3, 3)[natoms/2:natoms]:
          fout.write("B " + str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + "\n")
