import potentials.salt as salt
import numpy as np
import basinhopping as bh
import take_step.adaptive_step as adaptive_step
import take_step.random_displacement as ts
import storage.savenlowest
        

        
natoms = 4

# create PatchyParticle potential
pot = salt.salt()    

# load coords array
x = np.random.random(3*natoms+6)
x[-6:-3] = 3
x[-3:]=0.
x[-2]=3.

#for i in xrange(0,4):
#    for j in xrange(0,4):
#        for k in xrange(0,4):
#            x[(i*4*4 + j*4 +k)*3+0] = 1.0*i/4.
#            x[(i*4*4 + j*4 +k)*3+1] = 1.0*j/4.
#            x[(i*4*4 + j*4 +k)*3+2] = 1.0*k/4.
#exit()
#x,e,t1,t2 = quench.fire(x, pot.getEnergyGradient)
#print "start"
#x[1]=0.0
e,g = pot.getEnergyGradient(x)
print e,pot.getEnergy(x)
print g
#print "numerical"
gn = pot.NumericalDerivative(x, 1e-6)
print gn
#print "difference"
print (g-gn)/g
#exit()
a=x.copy()
pot.toReal(a)
print a
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
                      #,quenchRoutine=quench.fire
                      )

# do 100 mc steps
opt.run(100)

#print minima.data[0][1]
#print pot.getEnergy(minima.data[0][1])
#print pot.NumericalDerivative(minima.data[0][1], 1e-8)

a = minima.data[0][1].copy()
pot.toReal(a)
print a
# save the minima
with open("pylowest.xyz", "w") as fout:
  for i in minima.data:
      fout.write( str(natoms) + "\n")
      fout.write( "Energy = " + str(i[0]) + "\n")
      tmp=i[1].copy()
      pot.toReal(tmp)
      for atom in tmp.reshape(x.size/3, 3)[0:natoms/2]:
          fout.write("A " + str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + "\n")
      for atom in tmp.reshape(x.size/3, 3)[natoms/2:natoms]:
          fout.write("B " + str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + "\n")
