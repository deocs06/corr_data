import numpy, sharedmem, time

def wgp_shapcorr(eCent1,eeigvec1,eq1,DCent1,rpbins1,zpbins1,i1,j1,boxlen):

    permax = boxlen/2.0
    eCent2d = eCent1
    DCent2d = DCent1

    Ntote = len(eCent1)
    NtotD = len(DCent1)
   
    pareige = []

    print [i1,j1,Ntote,NtotD]

    for ii1 in range(0,Ntote):
              pareige.append(eeigvec1[ii1][0])             
    pareige = numpy.array(pareige)
    
    Hgp = numpy.zeros((len(rpbins1)-1,len(zpbins1)-1)) 
    for j in range(0,len(eCent1)):   
          Rdiff = DCent2d[0:NtotD] - eCent2d[j]
          Rmask = 1 - (numpy.abs(Rdiff/permax)).astype(numpy.int32)
          Roffset = numpy.ma.masked_array((-1)*boxlen*((Rdiff/permax).astype(numpy.int32)),Rmask)
          Rdiff = Rdiff + Roffset.filled(0)
          Rmeas = pow(Rdiff,2.0)
          Rmeas = Rmeas[:,0] + Rmeas[:,1]
          zmeas = Rdiff[:,2]
          Rmeas = pow(Rmeas,1.0/2.0)
          dij_c = numpy.float32((pareige[j,0]*Rdiff[:,0] + pareige[j,1]*Rdiff[:,1]))
          chkeD = (Rmeas > 0)
          Rdiff = Rdiff[chkeD]
          Rmeas = Rmeas[chkeD] 
          dij_c = dij_c[chkeD]
          zmeas = numpy.abs(zmeas[chkeD])
          dji_c = numpy.float32(dij_c/Rmeas)
          SeD = (2*pow(dji_c,2.0) - 1.0)*((1 - pow(eq1[j],2.0))/(1 + pow(eq1[j],2.0)))
          Hgp_j = numpy.histogram2d(Rmeas,zmeas,bins=[rpbins1,zpbins1],weights = SeD)
          Hgp = Hgp + Hgp_j[0]
    return Hgp

def wgp_countcorr(eCent1,DCent1,rpbins1,zpbins1,i1,j1,boxlen):

    permax = boxlen/2.0

    eCent2d = eCent1
    DCent2d = DCent1

    Ntote = len(eCent1)
    NtotD = len(DCent1)
   
    print [i1,j1,Ntote,NtotD]

    Hgp = numpy.zeros((len(rpbins1)-1,len(zpbins1)-1))
    ntc = len(eCent1)
    for j in range(0,len(eCent1)):  
          Rdiff = DCent2d[0:NtotD] - eCent2d[j]
          Rmask = 1 - (numpy.abs(Rdiff/permax)).astype(numpy.int32)
          Roffset = numpy.ma.masked_array((-1)*boxlen*((Rdiff/permax).astype(numpy.int32)),Rmask)
          Rdiff = numpy.array(Rdiff + Roffset.filled(0)).copy()
          Rmeas = pow(Rdiff,2.0)
          Rmeas = Rmeas[:,0] + Rmeas[:,1]
          zmeas = Rdiff[:,2]
          Rmeas = pow(Rmeas,1.0/2.0)
          chkeD = (Rmeas > 0)
          Rdiff = Rdiff[chkeD]
          Rmeas = Rmeas[chkeD] 
          zmeas = numpy.abs(zmeas[chkeD])
          Hgp_j = numpy.histogram2d(Rmeas,zmeas,bins=[rpbins1,zpbins1])
          Hgp = Hgp + Hgp_j[0]
          print ['cc',j,ntc]
    return Hgp

def ed_corr(eCent1,eeigvec1,DCent1,rpbins1,i1,j1,boxlen):

    permax = boxlen/2.0

    eCent2d = eCent1
    DCent2d = DCent1

    Ntote = len(eCent1)
    NtotD = len(DCent1)
   
    pareige = []

    print [i1,j1,Ntote,NtotD]

    for ii1 in range(0,Ntote):
              pareige.append(eeigvec1[ii1][0])             
    pareige = numpy.array(pareige)

    Hgp = numpy.zeros(len(rpbins1)-1)
    Hr = numpy.zeros(len(rpbins1)-1)
    Hang = numpy.zeros(len(rpbins1)-1)
    for j in range(0,len(eCent1)):
          Rdiff = DCent2d[0:NtotD] - eCent2d[j]
          Rmask = 1 - (numpy.abs(Rdiff/permax)).astype(numpy.int32)
          Roffset = numpy.ma.masked_array((-1)*boxlen*((Rdiff/permax).astype(numpy.int32)),Rmask)
          Rdiff = numpy.array(Rdiff + Roffset.filled(0)).copy()
          Rmeas = pow(Rdiff,2.0)
          Rmeas = Rmeas[:,0] + Rmeas[:,1] + Rmeas[:,2]
          Rmeas = pow(Rmeas,1.0/2.0)
          dij_c = numpy.float32((pareige[j,0]*Rdiff[:,0] + pareige[j,1]*Rdiff[:,1] + pareige[j,2]*Rdiff[:,2]))
          chkeD = (Rmeas > 0)
          Rdiff = Rdiff[chkeD]
          Rmeas = Rmeas[chkeD] 
          dij_c = dij_c[chkeD]
          dji_c = numpy.float32(dij_c/Rmeas)
          dji_c[numpy.abs(dji_c) > 1.0] = 1.0
          SeD = pow(dji_c,2.0)
          Hgp_j = numpy.histogram(Rmeas,bins=rpbins1,weights = SeD)
          Hgp = Hgp + Hgp_j[0]
          Hr_j = numpy.histogram(Rmeas,bins=rpbins1)
          Hr = Hr + Hr_j[0]
    return Hgp, Hr
