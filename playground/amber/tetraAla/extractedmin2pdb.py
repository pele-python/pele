#!/usr/bin/env python

import sys
import string

###############################################################
##                                                            #
## Edyta Malolepsza                                           #
## David Wales' group, University of Cambridge                #
## in case of problems please send email: em427@cam.ac.uk     #
##                                                            #
###############################################################
##
## transform path.info into path_all.pdb
## ./path2pdb.py
##
###############################################################

prmtop = open("coords.prmtop").read()
f = string.split(prmtop, "\n")

q0 = f.index("%FLAG POINTERS                                                                  ")
q1 = f.index("%FLAG ATOM_NAME                                                                 ")
q2 = f.index("%FLAG RESIDUE_LABEL                                                             ")
q3 = f.index("%FLAG RESIDUE_POINTER                                                           ")

## names of tables are related to names in prmtop file

atomNumber = int(string.split(f[q0+2])[0])

print atomNumber 

atomName = []
residueLabel = []

an = 0
line = 0
while (an<atomNumber):
  for j in range(20):
    if (an<atomNumber):
      an += 1
      atomName.append(f[q1+2+line][j*4:(j+1)*4].strip())
    else:
      break
  line += 1

for i in range(q3-q2-2):
  for j in range((len(f[q2+2+i])+1)/4):
    residueLabel.append(string.strip(f[q2+2+i][j*4:4*(j+1)]))

info = open("path.info").read()
ff = string.split(info, "\n")

xyz = open('path_all.pdb','w')

for i in range(len(ff)/(atomNumber)):
    m = atomNumber  # number of lines for each stationary points
    l = 0       # number of lines before coordinates
    mm = 1                 # number of residue
    print i 
    for j in range(atomNumber):
#        print j
        x = float(string.split(ff[m*i+l+j])[0])
        y = float(string.split(ff[m*i+l+j])[1])
        z = float(string.split(ff[m*i+l+j])[2])
        # bug fix on 12 Jan 2012 
        if (atomName[j]=='N' and atomName[j+1]=='H'):
            mm += 1

        xyz.write("%4s%7d %-4s %4s%5d%12.3f%8.3f%8.3f\n" % ('ATOM', j+1, atomName[j], 'ALA', mm, x, y, z))
    xyz.write("END\n")

xyz.write("\n")
xyz.close()

