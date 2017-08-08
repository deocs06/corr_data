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
gini_ij = 0
gfin_ij = len(gbn)

for jki in range(gini_ij,gfin_ij):

      jk_chk1 = ~((cen_S[:,0]>=sxmin[gbn[jki]])&(cen_S[:,0]<sxmax[gbn[jki]])&(cen_S[:,1]>=symin[gbn[jki]])&(cen_S[:,1]<symax[gbn[jki]])&(cen_S[:,2]>=szmin[gbn[jki]])&(cen_S[:,2]<szmax[gbn[jki]]))  
      fDMcenter_jk = fDMcenter[jk_chk1]
      fgrplen_jk = fgrplen[jk_chk1]
      fpno_jk = fpno[jk_chk1] 
      feigvec_jk = feigvec[jk_chk1]
      fmass_jk = fmass[jk_chk1]
      fq_jk = fq[jk_chk1]

      jk_chk2 = ~((cen_D[:,0]>=sxmin[gbn[jki]])&(cen_D[:,0]<sxmax[gbn[jki]])&(cen_D[:,1]>=symin[gbn[jki]])&(cen_D[:,1]<symax[gbn[jki]])&(cen_D[:,2]>=szmin[gbn[jki]])&(cen_D[:,2]<szmax[gbn[jki]]))
      cen_D_jk = cen_D[jk_chk2]      

      jk_chk3 = ~((cen_randR[:,0]>=sxmin[gbn[jki]])&(cen_randR[:,0]<sxmax[gbn[jki]])&(cen_randR[:,1]>=symin[gbn[jki]])&(cen_randR[:,1]<symax[gbn[jki]])&(cen_randR[:,2]>=szmin[gbn[jki]])&(cen_randR[:,2]<szmax[gbn[jki]]))
      cen_randR_jk = cen_randR[jk_chk3]

      jk_chk4 = ~((cen_randD[:,0]>=sxmin[gbn[jki]])&(cen_randD[:,0]<sxmax[gbn[jki]])&(cen_randD[:,1]>=symin[gbn[jki]])&(cen_randD[:,1]<symax[gbn[jki]])&(cen_randD[:,2]>=szmin[gbn[jki]])&(cen_randD[:,2]<szmax[gbn[jki]]))
      cen_randD_jk = cen_randD[jk_chk4]
      
      slices1 = [slice(i,i+1) for i in range(0,gnp)]
      gnp = 16
      Dlen1 = len(fq_jk)
      Dlen2 = len(cen_D_jk)
      Blen2 = Dlen2/gnp
      Dlen3 = len(cen_randR_jk)
      Blen3 = Dlen3/gnp
      Dlen4 = len(cen_randD_jk)
      Blen4 = Dlen4/(pnD)
   
      jklen = numpy.array([Dlen,Dlen1,Dlen2,Dlen3,Dlen4])
      jklen.tofile('/home/EEEDdata_errors/' + rdstr2 + '_'+ str(Dlen) + 'jklen5' + '.' + str(jki) + clrstr1 + str_nsubl  + '.' + dirno + '.' + str(mmin) + '.alltracer')
    
      with sharedmem.Pool() as pool:
           def work(slice):
               si1 = slice.start
               sj1 = slice.stop
               seCent = fDMcenter_jk
               seeigvec = feigvec_jk
               sq = fq_jk
               if (sj1 < 16):
                  sDcent = cen_D_jk[si1*Blen2 : sj1*Blen2]
                  sRcent = cen_randR_jk[si1*Blen3 : sj1*Blen3]
               if (sj1 == 16):
                  sDcent = cen_D_jk[si1*Blen2 : Dlen2]
                  sRcent = cen_randR_jk[si1*Blen3 : Dlen3]
               hij1,hr1 = ed_corr(seCent,seeigvec,sDcent,rpbins,si1,sj1,boxlen)                        
               hij1.tofile('/home/EEEDdata_errors/' + rdstr2 + '_'+ str(Dlen) + 'ed' + '.' + str(jki) + clrstr1 + str_nsubl  + '.' + dirno + '.' + str(mmin) + '.' + str(si1) + str(sj1) + '.alltracer')
               hr1.tofile('/home/EEEDdata_errors/' + rdstr2 + '_'+ str(Dlen) + 'rhisted' + '.' + str(jki) + clrstr1 + str_nsubl  + '.' + dirno + '.' + str(mmin) + '.' + str(si1) + str(sj1) + '.alltracer')     
           pool.map(work,slices1) 
