import numpy as np
import os
#
strl = "HZRTA_CaI_PRD_"
strlcv = "HZRTA_CaI_PRD_"
fllam = 0#1 only required for strong fields
#
f1 = open('param.dat','r')
header = f1.readline()
line = f1.readline()
line = line.split()
tgm = int(line[0])
mu = [0.0 for x in range(tgm)]
chi = [0.0 for x in range(tgm)]
ct = int(0)
for line in f1:
 line = line.split()
 if ct%2 == 0:
  mu[ct//2] = float(line[0])
 else:
  chi[ct//2] = float(line[0])
 ct = ct+1
 if ct//2 >= tgm:
  break
f1.close()
f2 = open('lambda.dat','r')
line = f2.readline()
line = line.split()
ivac = int(line[0])
line = f2.readline()
line = line.split()
wl0 = float(line[0])
line = f2.readline()
line = line.split()
nwl = int(line[0])
wavl = [0.0 for x in range(nwl)]
ct = 0
for line in f2:
 line = line.split()
 wavl[ct] = float(line[0])
 ct = ct+1
f2.close()
#
stok = [[[0.0 for z in range(tgm)] for y in range(nwl)] for x in range(0,4)]
ff = open('field_input.dat','r')
fline = ff.readline()
for fline in ff:
 fline = fline.split()
 Bint = float(fline[0])
 thetb = float(fline[1])
 chib = float(fline[2])
 strr = str(fline[3])
 fm = open('magnetic.dat','w')
 fm.write("%9.4f \n %1.16f \n %1.16f \n %2d" % (float(Bint),thetb,chib,fllam))
 fm.close()
 os.system("../hz_rta") 
 for im in range(tgm):
  ststrl = 'stokes_'
  ststrr = '.res'
  stfile = ststrl+str(im+1)+ststrr
  fs = open(stfile,'r')
  line = fs.readline()
  ct = 0
  for line in fs:
   line = line.split()
   for i in range(0,4):
    stok[i][ct][im] = float(line[i+1])
   ct = ct+1
  fs.close()
#
 filout = strl+str(Bint)+strr
 np.savez(filout,mu=mu,chi=chi,wavl=wavl,stok=stok)
 os.system("rm *.res")
