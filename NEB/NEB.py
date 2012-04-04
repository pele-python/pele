import numpy as np

def distance_cart(x1,x2):
    return x2-x1
     
class NEB:
  def __init__(self, initial, final, potential, distance=distance_cart, nimages=10, k=1.0):
    self.distance=distance
    self.potential=potential
    self.k = k
    
    #initialiye coordinate&gradient array
    self.coords=np.zeros([nimages,initial.size])
    self.grad=np.zeros([nimages,initial.size])
    
    #interpolate initial points
    self.interpolate(initial, final, nimages)    
    
    # copy initial and final structure
    self.initial = initial
    self.final = final
    
  
  # do a stupid steepest descent, have to add scipy interface for minimizer first
  def run(self,nsteps=100):
    for i in xrange(nsteps):
      self.UpdateForces()
      self.coords = self.coords + 0.01*self.grad
    
    
  # update forces for all images
  def UpdateForces(self):
    coords = self.coords
    for i in xrange(0,10):
      pl=self.initial
      pr=self.final
      if i!=0:
        pl = coords[i-1,:]
      if i!=9:
        pr = coords[i+1,:]
      E,g = self.NEBForce(coords[i,:],pl,pr)
      self.grad[i,:] = g
    #print self.grad
      
  
  # update force for one image
  def NEBForce(self,p,pl,pr):
      d1 = p - pl
      d2 = pr - p
      t = d1/np.linalg.norm(d1) + d2/np.linalg.norm(d2)
      t = t / np.linalg.norm(t)
      
      E,g = self.potential.getEnergyGradient(p)
      
      gpar = g - np.dot(g,t)*t
      gspring = np.dot(self.k*(d2-d1),t)*t
      
      return E,gpar+gspring
  
  # initial interpolation  
  def interpolate(self,initial, final, nimages):
    delta=(final-initial)/(nimages+1)
    for i in xrange(0,nimages):
      self.coords[i,:] = initial + delta*(i+1)

      
# From here on testing

import math
import numpy
    
# LEPS 2d potential
class LEPS:
    def getEnergy( self, r ):       
        '''
        potential energy as a function of position
        for the LEPS potential on a line
        python version
        '''
        x=r[0]
        y=r[1]
        a = 0.05
        b = 0.3
        c = 0.05
        alpha = 1.942
        r0 = 0.742
        dAB = 4.746
        dBC = 4.746
        dAC = 3.445

        def Q( d, r ):
            return d*( 3*numpy.exp(-2*alpha*(r-r0))/2 - numpy.exp(-alpha*(r-r0)) )/2
               
        def J( d, r ):
            return d*( numpy.exp(-2*alpha*(r-r0)) - 6*numpy.exp(-alpha*(r-r0)) )/4

        
        rAB = x;
        rBC = y;
        rAC = rAB + rBC;
               
        JABred = J(dAB, rAB)/(1+a)
        JBCred = J(dBC, rBC)/(1+b)
        JACred = J(dAC, rAC)/(1+c)
                              
        return Q(dAB, rAB)/(1+a) + \
               Q(dBC, rBC)/(1+b) + \
               Q(dAC, rAC)/(1+c) - \
               numpy.sqrt( JABred*JABred + \
                           JBCred*JBCred + \
                           JACred*JACred - \
                           JABred*JBCred - \
                           JBCred*JACred - \
                           JABred*JACred )
                           
    def getEnergyGradient( self, r ):
        '''
        force as a function of position
        for the LEPS potential on a line
        python version
        '''
        x=r[0]
        y=r[1]
        a = 0.05
        b = 0.3
        c = 0.05
        alpha = 1.942
        r0 = 0.742
        dAB = 4.746
        dBC = 4.746
        dAC = 3.445


        def Q( d, r ):
            return d*( 3*numpy.exp(-2*alpha*(r-r0))/2 - numpy.exp(-alpha*(r-r0)) )/2
               
        def J( d, r ):
            return d*( numpy.exp(-2*alpha*(r-r0)) - 6*numpy.exp(-alpha*(r-r0)) )/4
                 
        def dQ( d, r ):
            return alpha*d*( -3*numpy.exp(-2*alpha*(r-r0)) + numpy.exp(-alpha*(r-r0)) )/2;
               
        def dJ( d, r ):
            return alpha*d*( -2*numpy.exp(-2*alpha*(r-r0)) + 6*numpy.exp(-alpha*(r-r0)) )/4;
        
        rAB = x;
        rBC = y;
        rAC = rAB + rBC;
               
        JABred = J(dAB, rAB)/(1+a);
        JBCred = J(dBC, rBC)/(1+b);
        JACred = J(dAC, rAC)/(1+c);

        dJABred = dJ(dAB, rAB)/(1+a);
        dJBCred = dJ(dBC, rBC)/(1+b);
        dJACred = dJ(dAC, rAC)/(1+c);
                              
        Fx = dQ(dAB, rAB)/(1+a) + \
             dQ(dAC, rAC)/(1+c) - \
             ( 2*JABred*dJABred + \
               2*JACred*dJACred - \
               dJABred*JBCred - \
               JBCred*dJACred - \
               dJABred*JACred - \
               JABred*dJACred ) / \
             ( 2 * numpy.sqrt( JABred*JABred + \
                               JBCred*JBCred + \
                               JACred*JACred - \
                               JABred*JBCred - \
                               JBCred*JACred - \
                               JABred*JACred ))


        Fy = dQ(dBC, rBC)/(1+b) + \
             dQ(dAC, rAC)/(1+c) - \
             ( 2*JBCred*dJBCred + \
               2*JACred*dJACred - \
               JABred*dJBCred - \
               dJBCred*JACred - \
               JBCred*dJACred - \
               JABred*dJACred ) / \
             ( 2 * numpy.sqrt( JABred*JABred + \
                               JBCred*JBCred + \
                               JACred*JACred - \
                               JABred*JBCred - \
                               JBCred*JACred - \
                               JABred*JACred ))

        return self.getEnergy,[ -Fx, -Fy ]
  
if __name__ == "__main__":
  import pylab as pl
  x = np.arange(.5,3.,.1)
  y = np.arange(.5,3.,.1)
  z = np.zeros([len(x),len(y)])
  potential = LEPS()
  for i in range(0, len(x)):
    for j in range(0, len(y)):
        z[j,i] =  potential.getEnergy([x[i],y[j]])
  print "done"
  #pl.imshow(z)
  #pl.show()
  initial = np.array([1.,2.]) #np.random.random(3)
  final = np.array([2.,.7]) #np.random.random(3)
  print "Initial: ", initial
  print "Final: ", final
  #pl.imshow(z)
  
  neb = NEB(initial, final, potential)
  tmp = neb.coords
  
  pl.contourf(x,y,z)
  pl.colorbar()
  pl.plot(tmp[:,0],tmp[:,1], 'o-')
  neb.run()

    
  tmp = neb.coords
  pl.plot(tmp[:,0],tmp[:,1], 'rx-')
  pl.show()
  
    