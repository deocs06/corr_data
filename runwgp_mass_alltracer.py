import numpy, sharedmem, time
from shape_corr import *

ptype = 4
nmin = 10
dirno = '085'
mmin = 10.0
mmax = 15.0
subgrpbool = 1
logbool = 1
ndim = 2
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

arbins = numpy.zeros(len(rpbins)-1)
for aif in range(0,len(arbins)-1):
        arbins[aif] = (numpy.pi)*(pow(rpbins[aif+1],2.0) - pow(rpbins[aif],2.0))*Dlen*len(fq)/(boxlen*boxlen)
arbins[len(arbins) - 1] = ((boxlen*boxlen) - (numpy.pi)*(pow(rpbins[aif+1],2.0)))*Dlen*len(fq)/(boxlen*boxlen)
slices1 = [slice(i,i+1) for i in range(0,gnp)]

Hfin = numpy.zeros((len(rpbins)-1,len(zpbins)-1))
for sij in range(0,gnp):
       si1 = sij
       sj1 = sij + 1
       Hfin = Hfin + numpy.fromfile(dirname1+'/wgpdata_alltracer/' + rdstr2 + '_'+ str(Dlen) + 'Sgp' + clrstr1 + '.' + dirno + '.' + str(mmin) + '.' + str(si1) + str(sj1) + typ1 + '.alltracer',dtype = (numpy.float,len(rpbins)-1))
      
rpoints = numpy.zeros(len(arbins))
galPcorrpoints = numpy.zeros(len(arbins))
galPcorr2points = numpy.zeros(len(arbins))
galPcorrcpoints = numpy.zeros(len(arbins))
e_cut2 = (((1 - fq**2)/(1 + fq**2))**2).sum()/(2.0*len(fq))
etaR = 1 - e_cut2
bgf = 1.0
typ3 = ''

for ari in range(0,len(arbins)):
      rpoints[ari] = 0.5*(rpbins[ari] + rpbins[ari+1])
      for rj in range(0,len(arbins)):
             galPcorrpoints[ari] = galPcorrpoints[ari] + ((Hfin[ari,rj:rj+1]).sum())*(2*diff(zpbins)[0])*(boxlen/(2.0*(zpbins[rj+1] - zpbins[rj])))
             galPcorr2points[ari] = galPcorr2points[ari] + ((Hfin2[ari,rj:rj+1]).sum())*(2*diff(zpbins)[0])*(boxlen/(2.0*(zpbins[rj+1] - zpbins[rj])))
             galPcorrcpoints[ari] = galPcorrcpoints[ari] + ((Hfinc[ari,rj:rj+1]).sum())
      errPcorrpoints[ari] = boxlen*pow((galPcorr2points[ari]/boxlen - (pow(galPcorrpoints[ari]/boxlen,2.0)/galPcorrcpoints[ari])),1.0/2.0)
      galPcorrpoints[ari] = galPcorrpoints[ari]/arbins[ari]
      galPcorr2points[ari] = galPcorr2points[ari]/arbins[ari]
      errPcorrpoints[ari] = errPcorrpoints[ari]/pow(arbins[ari],1.0)

rpoints = np.array(rpoints)/1000.0
galPcorrpoints = (1.0/(2.0*etaR*bgf))*np.array(galPcorrpoints)/1000.0
errPcorrpoints = (1.0/(2.0*etaR*bgf))*np.array(errPcorrpoints)/1000.0

fig1 = matplotlib.pyplot.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(rpoints[:],galPcorrpoints[:],color = 'b',lw = 2.0,label= 'wg')
matplotlib.pyplot.xscale('log')
matplotlib.pyplot.yscale('log')
matplotlib.pyplot.legend(loc=0)
matplotlib.pyplot.xlim([0.01,100])
matplotlib.pyplot.xlabel(r'$r_{p}$ ($h^{-1}$Mpc)',fontsize=20)
matplotlib.pyplot.ylabel(r'$w_{\delta } (r_{p})$',fontsize=20)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.show()

