import numpy, sharedmem, time
from shape_corr import *

ptype = 4
nmin = 10
dirno = '085'
mmin = 10.0
mmax = 15.0
subgrpbool = 1
logbool = 1
ndim = 3
gdim = str(ndim) + 'd'
rpno = 25
boxlen = 1000.0
permax = boxlen/2.0
Rmin = 0.1
Rmax = permax
bins = numpy.logspace(numpy.log10(Rmin),numpy.log10(Rmax),rpno,base=10)
rpbins = numpy.zeros(rpno+1,dtype=numpy.float32)
rpbins[0:rpno] = bins[:]
rpbins[rpno] = finalRmax

if (subgrpbool == 1):
   
   rdstr1 = ''
   center_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/centerdata_' + str(ndim) + 'dDM',dtype = (numpy.float,3))
   grplen_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/grplen_' + str(ndim) + 'dDM',dtype = (numpy.int,3))
   pno_DM = grplen_DM[:,2]
   eigvec_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/eigvecdata_' + str(ndim) + 'dDM',dtype = (numpy.float,(ndim,ndim)))
   mass_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/massdata_' + str(ndim) + 'dDM',dtype = numpy.float)
   mass_DM = mass_DM*pow(10,10)
   q_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/qval_' + str(ndim) + 'dDM',dtype = numpy.float)
   idCent_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/idcc_' + str(ndim) + 'dDM',dtype = (numpy.int))

   center_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/centerdata_' + str(ndim) + 'dstar',dtype = (numpy.float,3))
   grplen_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/grplen_' + str(ndim) + 'dstar',dtype = (numpy.int,3))
   pno_star = grplen_star[:,2]
   eigvec_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/eigvecdata_' + str(ndim) + 'dstar',dtype = (numpy.float,(ndim,ndim)))
   mass_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/massdata_' + str(ndim) + 'dstar',dtype = numpy.float)
   mass_star = mass_star*pow(10,10)
   q_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/qval_' + str(ndim) + 'dstar',dtype = numpy.float)
   idCent_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/idcc_' + str(ndim) + 'dstar',dtype = (numpy.int))

   piv = numpy.searchsorted(idCent_DM,idCent_star,side='left')
   chk1 = ((piv < len(idCent_DM)))
   piv3d = piv[chk1]

   center_DM2 = center_DM[piv3d]
   grplen_DM2 = grplen_DM[piv3d]
   pno_DM2 = pno_DM[piv3d]
   eigvec_DM2 = eigvec_DM[piv3d]
   mass_DM2 = mass_DM[piv3d]
   idCent_DM2 = idCent_DM[piv3d]

   center_star2 = center_star[chk1]
   grplen_star2 = grplen_star[chk1]
   pno_star2 = pno_star[chk1]
   eigvec_star2 = eigvec_star[chk1]
   mass_star2 = mass_star[chk1]
   idCent_star2 = idCent_star[chk1]

   chk2 = ((grplen_DM2[:,2] > nmin) & (grplen_star2[:,2] > nmin) & (idCent_DM2 == idCent_star2) & (mass_DM2 >= pow(10,mmin)) & (mass_DM2 < pow(10,mmax)))

   center_DM3 = center_DM2[chk2]
   grplen_DM3 = grplen_DM2[chk2]
   pno_DM3 = pno_DM2[chk2]
   eigvec_DM3 = eigvec_DM2[chk2]
   mass_DM3 = mass_DM2[chk2]
   idCent_DM3 = idCent_DM2[chk2]

   center_star3 = center_star2[chk2]
   grplen_star3 = grplen_star2[chk2]
   pno_star3 = pno_star2[chk2]
   eigvec_star3 = eigvec_star2[chk2]
   mass_star3 = mass_star2[chk2]
   idCent_star3 = idCent_star2[chk2]

   print [len(center_DM3),len(grplen_DM3),len(pno_DM3),len(eigvec_DM3), len(mass_DM3), len(idCent_DM3)]
   print [len(center_star3),len(grplen_star3),len(pno_star3),len(eigvec_star3), len(mass_star3), len(idCent_star3)]

   Dlen = 5000000
   center_pivD = numpy.fromfile('/home/twopointdata/randdmpos' + dirno + '.rand' + str(Dlen) + 'alltracer',dtype = (numpy.float64,3))

gnp = 16
Blen = Dlen/gnp
fDMcen_z = fDMcenter[:,2]

nsubl = 2
lsub = boxlen/nsubl
str_nsubl = str(nsubl)
if (nsubl == 1):
   lsub = 0
if (nsubl == 2):
   str_nsub = ''
sboxno = numpy.zeros(nsubl**3)
sxmin = numpy.zeros(nsubl**3)
sxmax = numpy.zeros(nsubl**3)
symin = numpy.zeros(nsubl**3)
symax = numpy.zeros(nsubl**3)
szmin = numpy.zeros(nsubl**3)
szmax = numpy.zeros(nsubl**3)

for zi in range(0,nsubl):
       for yi in range(0,nsubl):
              for xi in range(0,nsubl):
                     sxmin[(nsubl**2)*zi + (nsubl)*yi + xi] = xi*lsub
                     sxmax[(nsubl**2)*zi + (nsubl)*yi + xi] = (xi+1)*lsub
                     symin[(nsubl**2)*zi + (nsubl)*yi + xi] = yi*lsub
                     symax[(nsubl**2)*zi + (nsubl)*yi + xi] = (yi+1)*lsub
                     szmin[(nsubl**2)*zi + (nsubl)*yi + xi] = zi*lsub
                     szmax[(nsubl**2)*zi + (nsubl)*yi + xi] = (zi+1)*lsub
                     sboxno[(nsubl**2)*zi + (nsubl)*yi + xi] = (nsubl**2)*zi + (nsubl)*yi + xi 

gbn = sboxno[0:len(sboxno)]
etaR = 0.87
gnp = 16
rpoints = numpy.zeros(len(rpbins)-1)
for ari in range(0,len(rpbins)-1):
        rpoints[ari] = 0.5*(rpbins[ari] + rpbins[ari+1])
galPcorrpoints_jk = numpy.zeros((len(gbn),len(rpbins)-1))
galPcorrpoints_mjk = numpy.zeros(len(rpbins)-1)

for jki in range(0,len(gbn)):
        Hsd = numpy.zeros(len(rpbins)-1)
        Hsr = numpy.zeros(len(rpbins)-1)
        Hrdrr = numpy.zeros(len(rpbins)-1)
        rjklen = numpy.fromfile('/home/EEEDdata_errors/' + rdstr2 + '_'+ str(Dlen) + 'jklen5' + '.' + str(jki) + clrstr1 + str_nsubl  + '.' + dirno + '.' + str(mmin) + '.' + typ1 + '.alltracer',dtype=numpy.int)
        rDlen = numpy.double(rjklen[0])
        rDlen1 = numpy.double(rjklen[1])     
        rDlen2 = numpy.double(rjklen[2])
        rDlen3 = numpy.double(rjklen[3])
        rDlen4 = numpy.double(rjklen[4])
        rjklen_12 = numpy.fromfile('/home/EEEDdata_errors/' + rdstr2 + '_'+ str(Dlen) + 'jklen5' + '.' + str(0) + '' + str_nsubl  + '.' + dirno + '.' + str(mmin) + '.' + typ1 + '.alltracer',dtype=numpy.int)
        rDlen_12 = numpy.double(rjklen_12[0])
        rDlen1_12 = numpy.double(rjklen_12[1])
        rDlen2_12 = numpy.double(rjklen_12[2])
        rDlen3_12 = numpy.double(rjklen_12[3])
        rDlen4_12 = numpy.double(rjklen_12[4])
        for sij in range(0,gnp):
                si1 = sij
                sj1 = sij + 1 
                Hsd = Hsd + numpy.fromfile('/home/EEEDdata_errors/' + rdstr2 + '_'+ str(Dlen) + 'ed' + '.' + str(jki) + clrstr1 + str_nsubl  + '.' + dirno + '.' + str(mmin) + '.' + str(si1) + str(sj1) + typ1 + '.alltracer',dtype = numpy.float)
                Hsr = Hsr + numpy.fromfile('/home/EEEDdata_errors/' + rdstr2 + '_'+ str(Dlen) + 'rhisted' + '.' + str(jki) + clrstr1 + str_nsubl  + '.' + dirno + '.' + str(mmin) + '.' + str(si1) + str(sj1) + typ1 + '.alltracer',dtype = numpy.float)
                Hrdrr = Hrdrr + numpy.fromfile('/home/EEEDdata_errors/' + rdstr2 + '_'+ str(Dlen) + 'anged' + '.' + str(jki) + clrstr1 + str_nsubl  + '.' + dirno + '.' + str(mmin) + '.' + str(si1) + str(sj1) + typ1 + '.alltracer',dtype = numpy.float)                         
        for ari in range(0,len(rpbins)-1):    
                for rj in range(0,len(zpbins)-1):
                       galPcorrpoints_jk[jki,ari] = galPcorrpoints_jk[jki,ari] + (Hsd[ari]/Hsr[ari])
        print jki
rpoints = np.array(rpoints)/1000.0
galPcorrpoints_jk = galPcorrpoints_jk - 1.0/3.0
errvals_jk = numpy.zeros(len(rpbins)-1)
for mjk in range(0,len(rpbins)-1):
        galPcorrpoints_mjk[mjk] = numpy.mean(galPcorrpoints_jk[:,mjk])
        errvals_jk[mjk] = pow((nsubl**3 - 1)*(numpy.mean((galPcorrpoints_jk[:,mjk] - galPcorrpoints_mjk[mjk])**2)),1.0/2.0)
 
for jki in range(0,len(gbn)):
        matplotlib.pyplot.plot(rpoints[:],(galPcorrpoints_jk[jki,:][:])/bgf)
matplotlib.pyplot.show()

matplotlib.pyplot.errorbar(rpoints[:],galPcorrpoints[:]/bgf,yerr =errvals_jk[:]/bgf,color = 'k',lw = 2.5) 
matplotlib.pyplot.xscale('log')
matplotlib.pyplot.yscale('log')
matplotlib.pyplot.legend(loc=0,fontsize=10)
matplotlib.pyplot.show()

