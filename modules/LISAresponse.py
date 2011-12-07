from synthlisa import *
import matplotlib
matplotlib.use('Agg')
from pylab import *
from spharm import Spharmt, getspecindx
import numpy


class planeGWresponse(object):
    def __init__(self,LISA,elat,elon,pol):
        self.Lisa = LISA
        self.wave = SimpleMonochromatic(1e-3,0,1,2,elat,elon,pol)
        self.elat = elat
        self.elon = elon
        self.pol = pol
        return
    def getk(self):
        return array(self.wave.putk())
    def getepol(self,polarisation='plus'):
        try:
            if polarisation == 'plus':
                return array(self.wave.putep(self.elat,self.elon,self.pol))
            elif polarisation == 'cross':
                return array(self.wave.putec(self.elat,self.elon,self.pol))
            else:
                raise ValueError
        except ValueError:
            print "polarisation must be either 'plus' or 'cross'"
        return

    def getp(self,n,t=0):
        return array(self.Lisa.putp(n,t)) 
    def getn(self,n,t=0):
        return array(self.Lisa.putn(n,t)) 
    def armlength(self,n,t=0):
        return self.Lisa.armlength(n,t)
    def getfn(self,n,t=0):
        return 1/(2*pi*self.Lisa.armlength(n,t))

    def getlink_arm(self,link):
        link_dictionary = {'123':2,'312':1,'231':3,
                           '1-32':-3,'2-13':-1,'3-21':-2}
        return link_dictionary[link]
    def getlink_receiver(self,link):
        return int(link[-1])
        
    def linkF(self,f,link='123',t=0,polarisation='plus'):
        """
        The response function of a one-way Doppler measurement.
        The event of reception is taken to be the reference.
        """
        l = self.getlink_arm(link)
        r = self.getlink_receiver(link)
        n = self.getn(l,t)
        L = self.armlength(l,t)
        fl = self.getfn(l,t)
        p = self.getp(r,t)
        k = self.getk()
        E = self.getepol(polarisation)
        kn = dot(k,n)
        kp = dot(k,p)
        nEn = dot(n,dot(E,n))
        argu = (1-kn)*f/fl/2
        Famp = -1j*(f/fl/2)*nEn*sinc(argu/pi) #fractional frequency
#        Famp = -1j*(f/fl/2)*nEn*sinc(argu/pi) / (1j*2*pi*f) #phase difference       
        Fphse = -argu
        F = Famp*exp(1j*Fphse)
        return F

    def timedelay(self,tdi,t=0):
        Delay = 0
        for delay in tdi:
            Delay = Delay + delay[0]*self.armlength(delay[1],t)
        return Delay

    """
    def __TDIname = ( ( { 'lnk':'slr' , 'sgn':+/-1 , 'tdi':[ (+/-1,Ll) , (...) , ... ] } ,
                        {...} , ... ) ,
                      ref )
        return
    """

    "////////// Sagnac G1 |L|-closed TDI observables \\\\\\\\\\"
    def __TDIalpha1(self,ref):
        alpha1 = (({'lnk':'123','sgn':1,'tdi':[(-1,3),(-1,1)]},
                   {'lnk':'312','sgn':1,'tdi':[(-1,3)]},
                   {'lnk':'231','sgn':1,'tdi':[]},
                   {'lnk':'3-21','sgn':-1,'tdi':[]},
                   {'lnk':'2-13','sgn':-1,'tdi':[(-1,-2)]},
                   {'lnk':'1-32','sgn':-1,'tdi':[(-1,-2),(-1,-1)]})
                  ,ref)
        return alpha1
    def __TDIalpha2(self,ref):
        alpha2 = (({'lnk':'231','sgn':1,'tdi':[(-1,1),(-1,2)]},
                  {'lnk':'123','sgn':1,'tdi':[(-1,1)]},
                  {'lnk':'312','sgn':1,'tdi':[]},
                  {'lnk':'1-32','sgn':-1,'tdi':[]},
                  {'lnk':'3-21','sgn':-1,'tdi':[(-1,-3)]},
                  {'lnk':'2-13','sgn':-1,'tdi':[(-1,-3),(-1,-2)]})
                 ,ref)
        return alpha2
    def __TDIalpha3(self,ref):
        alpha3 = (({'lnk':'312','sgn':1,'tdi':[(-1,2),(-1,3)]},
                  {'lnk':'231','sgn':1,'tdi':[(-1,2)]},
                  {'lnk':'123','sgn':1,'tdi':[]},
                  {'lnk':'2-13','sgn':-1,'tdi':[]},
                  {'lnk':'1-32','sgn':-1,'tdi':[(-1,-1)]},
                  {'lnk':'3-21','sgn':-1,'tdi':[(-1,-1),(-1,-3)]})
                 ,ref)
        return alpha3

    def alpha1(self,f,t=0,polarisation='plus',ref=1):
        TDIinfo = self.__TDIalpha1(ref)
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64) ###NOTE complex64 straight works here because zeros and not pylab.zeros 
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/3

    def alpha2(self,f,t=0,polarisation='plus',ref=2):
        TDIinfo = self.__TDIalpha2(ref)
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/3

    def alpha3(self,f,t=0,polarisation='plus',ref=3):
        TDIinfo = self.__TDIalpha3(ref)
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/3

    def A(self,f,t=0,polarisation='plus',ref=1):
        alpha3 = self.alpha3(f,t,polarisation,ref)
        alpha1 = self.alpha1(f,t,polarisation,ref)
        A = (alpha3 - alpha1)/sqrt(2)
        return A
    def E(self,f,t=0,polarisation='plus',ref=1):
        alpha1 = self.alpha1(f,t,polarisation,ref)
        alpha2 = self.alpha2(f,t,polarisation,ref)
        alpha3 = self.alpha3(f,t,polarisation,ref)
        E = (alpha1 - 2*alpha2 + alpha3)/sqrt(6)
        return E
    def T(self,f,t=0,polarisation='plus',ref=1):
        alpha1 = self.alpha1(f,t,polarisation,ref)
        alpha2 = self.alpha2(f,t,polarisation,ref)
        alpha3 = self.alpha3(f,t,polarisation,ref)
        T = (alpha1 + alpha2 + alpha3)/sqrt(3)
        return T

    def AA(self,f,t=0,refI=1,refJ=1):
        FIp = self.A(f,t,'plus',refI)
        FIc = self.A(f,t,'cross',refI)
        FJp = self.A(f,t,'plus',refJ)
        FJc = self.A(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        AA = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return AA
    def EE(self,f,t=0,refI=1,refJ=1):
        FIp = self.E(f,t,'plus',refI)
        FIc = self.E(f,t,'cross',refI)
        FJp = self.E(f,t,'plus',refJ)
        FJc = self.E(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        EE = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return EE
    def TT(self,f,t=0,refI=1,refJ=1):
        FIp = self.T(f,t,'plus',refI)
        FIc = self.T(f,t,'cross',refI)
        FJp = self.T(f,t,'plus',refJ)
        FJc = self.T(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        TT = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return TT
    def AE(self,f,t=0,refI=1,refJ=1):
        FIp = self.A(f,t,'plus',refI)
        FIc = self.A(f,t,'cross',refI)
        FJp = self.E(f,t,'plus',refJ)
        FJc = self.E(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        AE = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return AE
    def EA(self,f,t=0,refI=1,refJ=1):
        FIp = self.E(f,t,'plus',refI)
        FIc = self.E(f,t,'cross',refI)
        FJp = self.A(f,t,'plus',refJ)
        FJc = self.A(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        EA = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return EA
    def AT(self,f,t=0,refI=1,refJ=1):
        FIp = self.A(f,t,'plus',refI)
        FIc = self.A(f,t,'cross',refI)
        FJp = self.T(f,t,'plus',refJ)
        FJc = self.T(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        AT = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return AT
    def ET(self,f,t=0,refI=1,refJ=1):
        FIp = self.E(f,t,'plus',refI)
        FIc = self.E(f,t,'cross',refI)
        FJp = self.T(f,t,'plus',refJ)
        FJc = self.T(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        ET = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return ET



    "////////// Sagnac Gm Ldot-closed TDI observables \\\\\\\\\\"
    def __mTDIalpha1(self,ref=1):
        alpha1m = (({'lnk':'123','sgn':1,'tdi':[(-1,-2),(-1,-1),(-1,-3),(-1,3),(-1,1)]},
                    {'lnk':'312','sgn':1,'tdi':[(-1,-2),(-1,-1),(-1,-3),(-1,3)]},
                    {'lnk':'231','sgn':1,'tdi':[(-1,-2),(-1,-1),(-1,-3)]},
                    {'lnk':'1-32','sgn':1,'tdi':[(-1,-2),(-1,-1)]},
                    {'lnk':'2-13','sgn':1,'tdi':[(-1,-2)]},
                    {'lnk':'3-21','sgn':1,'tdi':[]},
                    {'lnk':'231','sgn':-1,'tdi':[]},
                    {'lnk':'312','sgn':-1,'tdi':[(-1,3)]},
                    {'lnk':'123','sgn':-1,'tdi':[(-1,3),(-1,1)]},
                    {'lnk':'3-21','sgn':-1,'tdi':[(-1,3),(-1,1),(-1,2)]},
                    {'lnk':'2-13','sgn':-1,'tdi':[(-1,3),(-1,1),(-1,2),(-1,-2)]},
                    {'lnk':'1-32','sgn':-1,'tdi':[(-1,3),(-1,1),(-1,2),(-1,-2),(-1,-1)]}),
                   ref)
        return alpha1m
    def __mTDIalpha2(self,ref=2):
        alpha2m = (({'lnk':'231','sgn':1,'tdi':[(-1,-3),(-1,-2),(-1,-1),(-1,1),(-1,2)]},
                    {'lnk':'123','sgn':1,'tdi':[(-1,-3),(-1,-2),(-1,-1),(-1,1)]},
                    {'lnk':'312','sgn':1,'tdi':[(-1,-3),(-1,-2),(-1,-1)]},
                    {'lnk':'2-13','sgn':1,'tdi':[(-1,-3),(-1,-2)]},
                    {'lnk':'3-21','sgn':1,'tdi':[(-1,-3)]},
                    {'lnk':'1-32','sgn':1,'tdi':[]},
                    {'lnk':'312','sgn':-1,'tdi':[]},
                    {'lnk':'123','sgn':-1,'tdi':[(-1,1)]},
                    {'lnk':'231','sgn':-1,'tdi':[(-1,1),(-1,2)]},
                    {'lnk':'1-32','sgn':-1,'tdi':[(-1,1),(-1,2),(-1,3)]},
                    {'lnk':'3-21','sgn':-1,'tdi':[(-1,1),(-1,2),(-1,3),(-1,-3)]},
                    {'lnk':'2-13','sgn':-1,'tdi':[(-1,1),(-1,2),(-1,3),(-1,-3),(-1,-2)]}),
                   ref)
        return alpha2m
    def __mTDIalpha3(self,ref=3):
        alpha3m = (({'lnk':'312','sgn':1,'tdi':[(-1,-1),(-1,-3),(-1,-2),(-1,2),(-1,3)]},
                    {'lnk':'231','sgn':1,'tdi':[(-1,-1),(-1,-3),(-1,-2),(-1,2)]},
                    {'lnk':'123','sgn':1,'tdi':[(-1,-1),(-1,-3),(-1,-2)]},
                    {'lnk':'3-21','sgn':1,'tdi':[(-1,-1),(-1,-3)]},
                    {'lnk':'1-32','sgn':1,'tdi':[(-1,-1)]},
                    {'lnk':'2-13','sgn':1,'tdi':[]},
                    {'lnk':'123','sgn':-1,'tdi':[]},
                    {'lnk':'231','sgn':-1,'tdi':[(-1,2)]},
                    {'lnk':'312','sgn':-1,'tdi':[(-1,2),(-1,3)]},
                    {'lnk':'2-13','sgn':-1,'tdi':[(-1,2),(-1,3),(-1,1)]},
                    {'lnk':'1-32','sgn':-1,'tdi':[(-1,2),(-1,3),(-1,1),(-1,-1)]},
                    {'lnk':'3-21','sgn':-1,'tdi':[(-1,2),(-1,3),(-1,1),(-1,-1),(-1,-3)]}),
                   ref)
        return alpha3m

    def alpha1m(self,f,t=0,polarisation='plus',ref=1):
        TDIinfo = self.__mTDIalpha1(ref)
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/6
    def alpha2m(self,f,t=0,polarisation='plus',ref=2):
        TDIinfo = self.__mTDIalpha2(ref)
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/6
    def alpha3m(self,f,t=0,polarisation='plus',ref=3):
        TDIinfo = self.__mTDIalpha3(ref)
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/6

    def Am(self,f,t=0,polarisation='plus',ref=1):
        alpha3m = self.alpha3m(f,t,polarisation,ref)
        alpha1m = self.alpha1m(f,t,polarisation,ref)
        Am = (alpha3m - alpha1m)/sqrt(2)
        return Am
    def Em(self,f,t=0,polarisation='plus',ref=1):
        alpha1m = self.alpha1m(f,t,polarisation,ref)
        alpha2m = self.alpha2m(f,t,polarisation,ref)
        alpha3m = self.alpha3m(f,t,polarisation,ref)
        Em = (alpha1m - 2*alpha2m + alpha3m)/sqrt(6)
        return Em
    def Tm(self,f,t=0,polarisation='plus',ref=1):
        alpha1m = self.alpha1m(f,t,polarisation,ref)
        alpha2m = self.alpha2m(f,t,polarisation,ref)
        alpha3m = self.alpha3m(f,t,polarisation,ref)
        Tm = (alpha1m + alpha2m + alpha3m)/sqrt(3)
        return Tm

    def AAm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Am(f,t,'plus',refI)
        FIc = self.Am(f,t,'cross',refI)
        FJp = self.Am(f,t,'plus',refJ)
        FJc = self.Am(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        AAm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return AAm
    def EEm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Em(f,t,'plus',refI)
        FIc = self.Em(f,t,'cross',refI)
        FJp = self.Em(f,t,'plus',refJ)
        FJc = self.Em(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        EEm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return EEm
    def TTm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Tm(f,t,'plus',refI)
        FIc = self.Tm(f,t,'cross',refI)
        FJp = self.Tm(f,t,'plus',refJ)
        FJc = self.Tm(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        TTm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return TTm
    def AEm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Am(f,t,'plus',refI)
        FIc = self.Am(f,t,'cross',refI)
        FJp = self.Em(f,t,'plus',refJ)
        FJc = self.Em(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        AEm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return AEm
    def EAm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Em(f,t,'plus',refI)
        FIc = self.Em(f,t,'cross',refI)
        FJp = self.Am(f,t,'plus',refJ)
        FJc = self.Am(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        EAm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return EAm
    def ATm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Am(f,t,'plus',refI)
        FIc = self.Am(f,t,'cross',refI) 
        FJp = self.Tm(f,t,'plus',refJ)
        FJc = self.Tm(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        ATm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return ATm
    def TAm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Tm(f,t,'plus',refI)
        FIc = self.Tm(f,t,'cross',refI)
        FJp = self.Am(f,t,'plus',refJ)
        FJc = self.Am(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        TAm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return TAm
    def ETm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Em(f,t,'plus',refI)
        FIc = self.Em(f,t,'cross',refI)
        FJp = self.Tm(f,t,'plus',refJ)
        FJc = self.Tm(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        ETm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return ETm
    def TEm(self,f,t=0,refI=1,refJ=1):
        FIp = self.Tm(f,t,'plus',refI)
        FIc = self.Tm(f,t,'cross',refI)
        FJp = self.Em(f,t,'plus',refJ)
        FJc = self.Em(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        TEm = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return TEm

    
    "////////// Michelson G0 closed TDI observables \\\\\\\\\\"
    def __TDIX0( self , ref=1 ):
        TDIdict = ( ( { 'lnk':'123'  , 'sgn': 1 , 'tdi':[ (-1,-2) ] } ,
                      { 'lnk':'3-21' , 'sgn': 1 , 'tdi':[] } ,
                      { 'lnk':'231'  , 'sgn':-1 , 'tdi':[] } , 
                      { 'lnk':'1-32' , 'sgn':-1 , 'tdi':[ (-1,3) ] } ) ,
                    ref )
        return TDIdict
    def __TDIY0( self , ref=2 ):
        TDIdict = ( ( { 'lnk':'231'  , 'sgn': 1 , 'tdi':[ (-1,-3) ] } ,
                      { 'lnk':'1-32' , 'sgn': 1 , 'tdi':[] } ,
                      { 'lnk':'312'  , 'sgn':-1 , 'tdi':[] } , 
                      { 'lnk':'2-13' , 'sgn':-1 , 'tdi':[ (-1,1) ] } ) ,
                    ref )
        return TDIdict    
    def __TDIZ0( self , ref=3 ):
        TDIdict = ( ( { 'lnk':'312'  , 'sgn': 1 , 'tdi':[ (-1,-1) ] } ,
                      { 'lnk':'2-13' , 'sgn': 1 , 'tdi':[] } ,
                      { 'lnk':'123'  , 'sgn':-1 , 'tdi':[] } , 
                      { 'lnk':'3-21' , 'sgn':-1 , 'tdi':[ (-1,2) ] } ) ,
                    ref )
        return TDIdict    
    "Response functions"
    def X0( self , f , t=0 , polarisation='plus' , ref=1 ):
        TDIinfo = self.__TDIX0( ref )
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/2
    def Y0( self , f , t=0 , polarisation='plus' , ref=2 ):
        TDIinfo = self.__TDIY0( ref )
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/2    
    def Z0( self , f , t=0 , polarisation='plus' , ref=3 ):
        TDIinfo = self.__TDIZ0( ref )
        lnklist = TDIinfo[0]
        refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t)
            tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk)
            kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/2
    "Overlap-reduction functions"
    def X0X0(self,f,t=0,refI=1,refJ=1):
        FIp = self.X0(f,t,'plus',refI)
        FIc = self.X0(f,t,'cross',refI)
        FJp = self.X0(f,t,'plus',refJ)
        FJc = self.X0(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def X0Y0(self,f,t=0,refI=1,refJ=2):
        FIp = self.X0(f,t,'plus',refI)
        FIc = self.X0(f,t,'cross',refI)
        FJp = self.Y0(f,t,'plus',refJ)
        FJc = self.Y0(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def X0Z0(self,f,t=0,refI=1,refJ=1):
        FIp = self.X0(f,t,'plus',refI)
        FIc = self.X0(f,t,'cross',refI)
        FJp = self.Z0(f,t,'plus',refJ)
        FJc = self.Z0(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ    
    def Y0Y0(self,f,t=0,refI=1,refJ=1):
        FIp = self.Y0(f,t,'plus',refI)
        FIc = self.Y0(f,t,'cross',refI)
        FJp = self.Y0(f,t,'plus',refJ)
        FJc = self.Y0(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ                                                                            
    def Y0Z0(self,f,t=0,refI=1,refJ=1):
        FIp = self.Y0(f,t,'plus',refI)
        FIc = self.Y0(f,t,'cross',refI)
        FJp = self.Z0(f,t,'plus',refJ)
        FJc = self.Z0(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def Z0Z0(self,f,t=0,refI=1,refJ=1):
        FIp = self.Z0(f,t,'plus',refI)
        FIc = self.Z0(f,t,'cross',refI)
        FJp = self.Z0(f,t,'plus',refJ)
        FJc = self.Z0(f,t,'cross',refJ)
        pI = self.getp(refI,t)
        pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ    


    "////////// Michelson Gm L-closed TDI observables \\\\\\\\\\"
    def __TDIXm( self , ref=1 ):
        TDIdict = ( ( { 'lnk':'1-32' , 'sgn': 1 , 'tdi':[ (-1,-2) , (-1,2) , (-1,3) ] } ,
                      { 'lnk':'231' , 'sgn': 1 , 'tdi':[ (-1,-2) , (-1,2) ] } ,
                      { 'lnk':'123' , 'sgn': 1 , 'tdi':[ (-1,-2) ] } ,
                      { 'lnk':'3-21' , 'sgn': 1 , 'tdi':[] } ,
                      { 'lnk':'231' , 'sgn': -1 , 'tdi':[] } , 
                      { 'lnk':'1-32' , 'sgn': -1 , 'tdi':[ (-1,3) ] } ,
                      { 'lnk':'3-21' , 'sgn': -1 , 'tdi':[ (-1,3) , (-1,-3) ] } ,
                      { 'lnk':'123' , 'sgn': -1 , 'tdi':[ (-1,3) , (-1,-3) , (-1,-2) ] } ) ,
                    ref )
        return TDIdict
    def __TDIYm( self , ref=1 ):
        TDIdict = ( ( { 'lnk':'2-13' , 'sgn': 1 , 'tdi':[ (-1,-3) , (-1,3) , (-1,1) ] } ,
                      { 'lnk':'312' , 'sgn': 1 , 'tdi':[ (-1,-3) , (-1,3) ] } ,
                      { 'lnk':'231' , 'sgn': 1 , 'tdi':[ (-1,-3) ] } ,
                      { 'lnk':'1-32' , 'sgn': 1 , 'tdi':[] } ,
                      { 'lnk':'312' , 'sgn': -1 , 'tdi':[] } , 
                      { 'lnk':'2-13' , 'sgn': -1 , 'tdi':[ (-1,1) ] } ,
                      { 'lnk':'1-32' , 'sgn': -1 , 'tdi':[ (-1,1) , (-1,-1) ] } ,
                      { 'lnk':'231' , 'sgn': -1 , 'tdi':[ (-1,1) , (-1,-1) , (-1,-3) ] } ) ,
                    ref )
        return TDIdict
    def __TDIZm( self , ref=1 ):
        TDIdict = ( ( { 'lnk':'3-21' , 'sgn': 1 , 'tdi':[ (-1,-1) , (-1,1) , (-1,2) ] } ,
                      { 'lnk':'123' , 'sgn': 1 , 'tdi':[ (-1,-1) , (-1,1) ] } ,
                      { 'lnk':'312' , 'sgn': 1 , 'tdi':[ (-1,-1) ] } ,
                      { 'lnk':'2-13' , 'sgn': 1 , 'tdi':[] } ,
                      { 'lnk':'123' , 'sgn': -1 , 'tdi':[] } , 
                      { 'lnk':'3-21' , 'sgn': -1 , 'tdi':[ (-1,2) ] } ,
                      { 'lnk':'2-13' , 'sgn': -1 , 'tdi':[ (-1,2) , (-1,-2) ] } ,
                      { 'lnk':'312' , 'sgn': -1 , 'tdi':[ (-1,2) , (-1,-2) , (-1,-1) ] } ) ,
                    ref )
        return TDIdict


    "////////// Michelson G2 Ldot-closed TDI observables \\\\\\\\\\"
    def __TDIX2( self , ref=1 ) :
        TDIdict = (
            ( { 'lnk':'1-32' , 'sgn': 1 , 'tdi':[ (-1,3),(-1,-3),(-1,-2),(-1,2),(-1,-2),(-1,2),(-1,3) ] } ,
              { 'lnk':'231' , 'sgn': 1 , 'tdi':[ (-1,3),(-1,-3),(-1,-2),(-1,2),(-1,-2),(-1,2) ] } ,
              { 'lnk':'123' , 'sgn': 1 , 'tdi':[ (-1,3),(-1,-3),(-1,-2),(-1,2),(-1,-2) ] } ,
              { 'lnk':'3-21' , 'sgn': 1 , 'tdi':[ (-1,3),(-1,-3),(-1,-2),(-1,2) ] } ,
              { 'lnk':'123' , 'sgn': 1 , 'tdi':[ (-1,3),(-1,-3),(-1,-2) ] } ,
              { 'lnk':'3-21' , 'sgn': 1 , 'tdi':[ (-1,3),(-1,-3) ] } ,
              { 'lnk':'1-32' , 'sgn': 1 , 'tdi':[ (-1,3) ] } ,
              { 'lnk':'231' , 'sgn': 1 , 'tdi':[] } ,
              { 'lnk':'3-21' , 'sgn': -1 , 'tdi':[] } ,
              { 'lnk':'123' , 'sgn': -1 , 'tdi':[ (-1,-2) ] } ,
              { 'lnk':'231' , 'sgn': -1 , 'tdi':[ (-1,-2),(-1,2) ] } ,
              { 'lnk':'1-32' , 'sgn': -1 , 'tdi':[ (-1,-2),(-1,2),(-1,3) ] } ,
              { 'lnk':'231' , 'sgn': -1 , 'tdi':[ (-1,-2),(-1,2),(-1,3),(-1,-3) ] } ,
              { 'lnk':'1-32' , 'sgn': -1 , 'tdi':[ (-1,-2),(-1,2),(-1,3),(-1,-3),(-1,3) ] } ,
              { 'lnk':'3-21' , 'sgn': -1 , 'tdi':[ (-1,-2),(-1,2),(-1,3),(-1,-3),(-1,3),(-1,-3) ] } ,
              { 'lnk':'123' , 'sgn': -1 , 'tdi':[ (-1,-2),(-1,2),(-1,3),(-1,-3),(-1,3),(-1,-3),(-1,-2) ] } 
              ) ,
            ref )
        return TDIdict
    def __TDIY2( self , ref=2 ) :
        TDIdict = (
            ( { 'lnk':'2-13' , 'sgn': 1 , 'tdi':[ (-1,1),(-1,-1),(-1,-3),(-1,3),(-1,-3),(-1,3),(-1,1) ] } ,
              { 'lnk':'312' , 'sgn': 1 , 'tdi':[ (-1,1),(-1,-1),(-1,-3),(-1,3),(-1,-3),(-1,3) ] } ,
              { 'lnk':'231' , 'sgn': 1 , 'tdi':[ (-1,1),(-1,-1),(-1,-3),(-1,3),(-1,-3) ] } ,
              { 'lnk':'1-32' , 'sgn': 1 , 'tdi':[ (-1,1),(-1,-1),(-1,-3),(-1,3) ] } ,
              { 'lnk':'231' , 'sgn': 1 , 'tdi':[ (-1,1),(-1,-1),(-1,-3) ] } ,
              { 'lnk':'1-32' , 'sgn': 1 , 'tdi':[ (-1,1),(-1,-1) ] } ,
              { 'lnk':'2-13' , 'sgn': 1 , 'tdi':[ (-1,1) ] } ,
              { 'lnk':'312' , 'sgn': 1 , 'tdi':[] } ,
              { 'lnk':'1-32' , 'sgn': -1 , 'tdi':[] } ,
              { 'lnk':'231' , 'sgn': -1 , 'tdi':[ (-1,-3) ] } ,
              { 'lnk':'312' , 'sgn': -1 , 'tdi':[ (-1,-3),(-1,3) ] } ,
              { 'lnk':'2-13' , 'sgn': -1 , 'tdi':[ (-1,-3),(-1,3),(-1,1) ] } ,
              { 'lnk':'312' , 'sgn': -1 , 'tdi':[ (-1,-3),(-1,3),(-1,1),(-1,-1) ] } ,
              { 'lnk':'2-13' , 'sgn': -1 , 'tdi':[ (-1,-3),(-1,3),(-1,1),(-1,-1),(-1,1) ] } ,
              { 'lnk':'1-32' , 'sgn': -1 , 'tdi':[ (-1,-3),(-1,3),(-1,1),(-1,-1),(-1,1),(-1,-1) ] } ,
              { 'lnk':'231' , 'sgn': -1 , 'tdi':[ (-1,-3),(-1,3),(-1,1),(-1,-1),(-1,1),(-1,-1),(-1,-3) ] } 
              ) ,
            ref )
        return TDIdict
    def __TDIZ2( self , ref=3 ) :
        TDIdict = (
            ( { 'lnk':'3-21' , 'sgn': 1 , 'tdi':[ (-1,2),(-1,-2),(-1,-1),(-1,1),(-1,-1),(-1,1),(-1,2) ] } ,
              { 'lnk':'123' , 'sgn': 1 , 'tdi':[ (-1,2),(-1,-2),(-1,-1),(-1,1),(-1,-1),(-1,1) ] } ,
              { 'lnk':'312' , 'sgn': 1 , 'tdi':[ (-1,2),(-1,-2),(-1,-1),(-1,1),(-1,-1) ] } ,
              { 'lnk':'2-13' , 'sgn': 1 , 'tdi':[ (-1,2),(-1,-2),(-1,-1),(-1,1) ] } ,
              { 'lnk':'312' , 'sgn': 1 , 'tdi':[ (-1,2),(-1,-2),(-1,-1) ] } ,
              { 'lnk':'2-13' , 'sgn': 1 , 'tdi':[ (-1,2),(-1,-2) ] } ,
              { 'lnk':'3-21' , 'sgn': 1 , 'tdi':[ (-1,2) ] } ,
              { 'lnk':'123' , 'sgn': 1 , 'tdi':[] } ,
              { 'lnk':'2-13' , 'sgn': -1 , 'tdi':[] } ,
              { 'lnk':'312' , 'sgn': -1 , 'tdi':[ (-1,-1) ] } ,
              { 'lnk':'123' , 'sgn': -1 , 'tdi':[ (-1,-1),(-1,1) ] } ,
              { 'lnk':'3-21' , 'sgn': -1 , 'tdi':[ (-1,-1),(-1,1),(-1,2) ] } ,
              { 'lnk':'123' , 'sgn': -1 , 'tdi':[ (-1,-1),(-1,1),(-1,2),(-1,-2) ] } ,
              { 'lnk':'3-21' , 'sgn': -1 , 'tdi':[ (-1,-1),(-1,1),(-1,2),(-1,-2),(-1,2) ] } ,
              { 'lnk':'2-13' , 'sgn': -1 , 'tdi':[ (-1,-1),(-1,1),(-1,2),(-1,-2),(-1,2),(-1,-2) ] } ,
              { 'lnk':'312' , 'sgn': -1 , 'tdi':[ (-1,-1),(-1,1),(-1,2),(-1,-2),(-1,2),(-1,-2),(-1,-1) ] } 
              ) ,
            ref )
        return TDIdict    
    "Response functions"
    def X2( self , f , t=0 , polarisation='plus' , ref=1 ):
        TDIinfo = self.__TDIX2( ref )
        lnklist = TDIinfo[0] ; refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t) ; tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk) ; kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/8
    def Y2( self , f , t=0 , polarisation='plus' , ref=1 ):
        TDIinfo = self.__TDIY2( ref )
        lnklist = TDIinfo[0] ; refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t) ; tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk) ; kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/8
    def Z2( self , f , t=0 , polarisation='plus' , ref=1 ):
        TDIinfo = self.__TDIZ2( ref )
        lnklist = TDIinfo[0] ; refp = self.getp(TDIinfo[1],t)
        F = zeros(len(f),complex64)
        for lnk in lnklist:
            td = self.timedelay(lnk['tdi'],t) ; tlnk = t + td
            lnkp = self.getp(int(lnk['lnk'][-1]),tlnk) ; kp = dot(self.getk(),lnkp-refp)
            lnkF = self.linkF(f,lnk['lnk'],tlnk,polarisation)
            F += lnk['sgn']*lnkF*exp(1j*2*pi*f*(td-kp))
        return F/8    
    def a2(self,f,t=0,polarisation='plus',ref=1) :
        X2 , Y2 , Z2 = self.X2(f,t,polarisation,ref) , self.Y2(f,t,polarisation,ref) , self.Z2(f,t,polarisation,ref)
        F = ( 2*X2 - Y2 - Z2 ) / 3
        return F
    def e2(self,f,t=0,polarisation='plus',ref=1) :
        Y2 , Z2 = self.Y2(f,t,polarisation,ref) , self.Z2(f,t,polarisation,ref)
        F = ( Z2 - Y2 )/ sqrt(3)
        return F
    def t2(self,f,t=0,polarisation='plus',ref=1) :
        X2 , Y2 , Z2 = self.X2(f,t,polarisation,ref) , self.Y2(f,t,polarisation,ref) , self.Z2(f,t,polarisation,ref)
        F = ( X2 + Y2 + Z2  ) / 3
        return F
    "Overlap-reduction functions"
    def X2X2(self,f,t=0,refI=1,refJ=1):
        FIp = self.X2(f,t,'plus',refI) ; FIc = self.X2(f,t,'cross',refI)
        FJp = self.X2(f,t,'plus',refJ) ; FJc = self.X2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def X2Y2(self,f,t=0,refI=1,refJ=1):
        FIp = self.X2(f,t,'plus',refI) ; FIc = self.X2(f,t,'cross',refI)
        FJp = self.Y2(f,t,'plus',refJ) ; FJc = self.Y2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def X2Z2(self,f,t=0,refI=1,refJ=1):
        FIp = self.X2(f,t,'plus',refI) ; FIc = self.X2(f,t,'cross',refI)
        FJp = self.Z2(f,t,'plus',refJ) ; FJc = self.Z2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ    
    def Y2Y2(self,f,t=0,refI=1,refJ=1):
        FIp = self.Y2(f,t,'plus',refI) ; FIc = self.Y2(f,t,'cross',refI)
        FJp = self.Y2(f,t,'plus',refJ) ; FJc = self.Y2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def Y2Z2(self,f,t=0,refI=1,refJ=1):
        FIp = self.Y2(f,t,'plus',refI) ; FIc = self.Y2(f,t,'cross',refI)
        FJp = self.Z2(f,t,'plus',refJ) ; FJc = self.Z2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def Z2Z2(self,f,t=0,refI=1,refJ=1):
        FIp = self.Z2(f,t,'plus',refI) ; FIc = self.Z2(f,t,'cross',refI)
        FJp = self.Z2(f,t,'plus',refJ) ; FJc = self.Z2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def a2a2(self,f,t=0,refI=1,refJ=1):
        FIp = self.a2(f,t,'plus',refI) ; FIc = self.a2(f,t,'cross',refI)
        FJp = self.a2(f,t,'plus',refJ) ; FJc = self.a2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def a2e2(self,f,t=0,refI=1,refJ=1):
        FIp = self.a2(f,t,'plus',refI) ; FIc = self.a2(f,t,'cross',refI)
        FJp = self.e2(f,t,'plus',refJ) ; FJc = self.e2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ    
    def a2t2(self,f,t=0,refI=1,refJ=1):
        FIp = self.a2(f,t,'plus',refI) ; FIc = self.a2(f,t,'cross',refI)
        FJp = self.t2(f,t,'plus',refJ) ; FJc = self.t2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ    
    def e2a2(self,f,t=0,refI=1,refJ=1):
        FIp = self.e2(f,t,'plus',refI) ; FIc = self.e2(f,t,'cross',refI)
        FJp = self.a2(f,t,'plus',refJ) ; FJc = self.a2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def e2e2(self,f,t=0,refI=1,refJ=1):
        FIp = self.e2(f,t,'plus',refI) ; FIc = self.e2(f,t,'cross',refI)
        FJp = self.e2(f,t,'plus',refJ) ; FJc = self.e2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def e2t2(self,f,t=0,refI=1,refJ=1):
        FIp = self.e2(f,t,'plus',refI) ; FIc = self.e2(f,t,'cross',refI)
        FJp = self.t2(f,t,'plus',refJ) ; FJc = self.t2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def t2a2(self,f,t=0,refI=1,refJ=1):
        FIp = self.t2(f,t,'plus',refI) ; FIc = self.t2(f,t,'cross',refI)
        FJp = self.a2(f,t,'plus',refJ) ; FJc = self.a2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def t2e2(self,f,t=0,refI=1,refJ=1):
        FIp = self.t2(f,t,'plus',refI) ; FIc = self.t2(f,t,'cross',refI)
        FJp = self.e2(f,t,'plus',refJ) ; FJc = self.e2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ
    def t2t2(self,f,t=0,refI=1,refJ=1):
        FIp = self.t2(f,t,'plus',refI) ; FIc = self.t2(f,t,'cross',refI)
        FJp = self.t2(f,t,'plus',refJ) ; FJc = self.t2(f,t,'cross',refJ)
        pI = self.getp(refI,t) ; pJ = self.getp(refJ,t)
        kp = dot(self.getk(),pI-pJ)
        IJ = (conj(FIp)*FJp + conj(FIc)*FJc)*exp(1j*2*pi*f*kp)
        return IJ    

        
class mySpharmt(Spharmt):
    def __init__(self,nlon,nlat,rsphere=1e2,gridtype='regular',legfunc='stored'):
        Spharmt.__init__(self,nlon,nlat,rsphere,gridtype,legfunc)
        if self.gridtype=='gaussian':
            raise Exception, "Can't deal with Gaussian grid at the moment."
        else:
            deltalat = 180./(nlat-1)
            self.lats = 90 - deltalat*arange(nlat)
        deltalon = 360./nlon
        self.lons = deltalon*arange(nlon)
        return


class OverlapReduction(object):
    def __init__(self,LISA,mySpharmt):
        self.Lisa = LISA
        self.sky = mySpharmt
        return

    "[[[[[[[[ Pixel basis ]]]]]]]]"
    "Sagnac-defined G1 |L|-closed optimal TDI observables"
    def AA(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.AA(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfAA = orf
        return orf
    def AE(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.AE(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfAE = orf
        return orf
    def EE(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.EE(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfEE = orf
        return orf
    def TT(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.TT(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfTT = orf
        return orf    
    "Sagnac-defined Gm Ldot-closed optimal TDI observables"
    def AAm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.AAm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfAAm = orf
        return orf
    def EEm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.EEm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfEEm = orf
        return orf
    def TTm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.TTm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfTTm = orf
        return orf
    def AEm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.AEm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfAEm = orf
        return orf
    def EAm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.EAm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfEAm = orf
        return orf
    def ATm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.ATm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfATm = orf
        return orf
    def TAm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.TAm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfTAm = orf
        return orf
    def ETm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.ETm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfETm = orf
        return orf
    def TEm(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.TEm(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfTEm = orf
        return orf
    "Michelson G0 closed TDI observables"
    def X0X0(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.X0X0(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfX0X0 = orf
        return orf
    def X0Y0(self,f,t=0,refI=1,refJ=2):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.X0Y0(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfX0Y0 = orf
        return orf
    def X0Z0(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.X0Z0(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfX0Z0 = orf
        return orf
    def Y0Y0(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.Y0Y0(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfY0Y0 = orf
        return orf
    def Y0Z0(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.Y0Z0(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfY0Z0 = orf
        return orf
    def Z0Z0(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.Z0Z0(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfZ0Z0 = orf
        return orf
    "Michelson-defined G2 Ldot-closed optimal TDI observables"
    def a2a2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.a2a2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfa2a2 = orf
        return orf
    def a2e2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.a2e2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfa2e2 = orf
        return orf
    def a2t2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.a2t2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfa2t2 = orf
        return orf    
    def e2a2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.e2a2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfe2a2 = orf
        return orf
    def e2e2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.e2e2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfe2e2 = orf
        return orf
    def e2t2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.e2t2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orfe2t2 = orf
        return orf
    def t2a2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.t2a2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orft2a2 = orf
        return orf
    def t2e2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.t2e2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orft2e2 = orf
        return orf    
    def t2t2(self,f,t=0,refI=1,refJ=1):
        nlat = self.sky.nlat
        lats = self.sky.lats
        nlon = self.sky.nlon
        lons = self.sky.lons
        orf = zeros((nlat,nlon,len(f)),complex64)
        for x,elatd in enumerate(lats):
            for y,elond in enumerate(lons):
                elat=elatd*pi/180
                elon=elond*pi/180
                response = planeGWresponse(self.Lisa,elat,elon,0)
                orf[x,y,:] = response.t2t2(f,t,refI,refJ)
        self.f = f
        self.t = t
        self.orft2t2 = orf
        return orf




    "(((((((( Spherical harmonic basis ))))))))"
    "Sagnac-defined G1 |L|-closed optimal TDI observables"
    def getspecAA(self,ntrunc=2):
        sky = self.sky
        orf = self.orfAA
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrAA = specr
        self.speciAA = speci
        return specr,speci
    def getspecAE(self,ntrunc=2):
        sky = self.sky
        orf = self.orfAE
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrAE = specr
        self.speciAE = speci
        return specr,speci
    def getspecEE(self,ntrunc=2):
        sky = self.sky
        orf = self.orfEE
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrEE = specr
        self.speciEE = speci
        return specr,speci
    def getspecTT(self,ntrunc=2):
        sky = self.sky
        orf = self.orfTT
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrTT = specr
        self.speciTT = speci
        return specr,speci    
    "Sagnac-defined Gm Ldot-closed optimal TDI observables"
    def getspecAAm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfAAm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrAAm = specr
        self.speciAAm = speci
        return specr,speci
    def getspecEEm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfEEm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrEEm = specr
        self.speciEEm = speci
        return specr,speci
    def getspecTTm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfTTm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrTTm = specr
        self.speciTTm = speci
        return specr,speci
    def getspecAEm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfAEm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrAEm = specr
        self.speciAEm = speci
        return specr,speci
    def getspecEAm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfEAm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrEAm =  specr
        self.speciEAm =  speci
        return specr,speci
    def getspecATm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfATm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrATm = specr
        self.speciATm = speci
        return specr,speci
    def getspecTAm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfTAm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrTAm =  specr
        self.speciTAm =  speci
        return specr,speci
    def getspecETm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfETm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrETm = specr
        self.speciETm = speci
        return specr,speci
    def getspecTEm(self,ntrunc=2):
        sky = self.sky
        orf = self.orfTEm
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrTEm =  specr
        self.speciTEm =  speci
        return specr,speci
    "Michelson G0 closed TDI observables"
    def getspecX0X0(self,ntrunc=2):
        sky = self.sky
        orf = self.orfX0X0
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrX0X0 = specr
        self.speciX0X0 = speci
        return specr,speci
    def getspecX0Y0(self,ntrunc=2):
        sky = self.sky
        orf = self.orfX0Y0
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrX0Y0 = specr
        self.speciX0Y0 = speci
        return specr,speci    
    def getspecX0Z0(self,ntrunc=2):
        sky = self.sky
        orf = self.orfX0Z0
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrX0Z0 = specr
        self.speciX0Z0 = speci
        return specr,speci
    def getspecY0Y0(self,ntrunc=2):
        sky = self.sky
        orf = self.orfY0Y0
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrY0Y0 = specr
        self.speciY0Y0 = speci
        return specr,speci
    def getspecY0Z0(self,ntrunc=2):
        sky = self.sky
        orf = self.orfY0Z0
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrY0Z0 = specr
        self.speciY0Z0 = speci
        return specr,speci
    def getspecZ0Z0(self,ntrunc=2):
        sky = self.sky
        orf = self.orfZ0Z0
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64)
        speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrZ0Z0 = specr
        self.speciZ0Z0 = speci
        return specr,speci
    "Michelson-defined G2 Ldot-closed optimal TDI observables"
    def getspeca2a2(self,ntrunc=2):
        sky = self.sky ; orf = self.orfa2a2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specra2a2 = specr ; self.specia2a2 = speci
        return specr,speci
    def getspeca2e2(self,ntrunc=2):
        sky = self.sky ; orf = self.orfa2e2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specra2e2 = specr ; self.specia2e2 = speci
        return specr,speci    
    def getspeca2t2(self,ntrunc=2):
        sky = self.sky ; orf = self.orfa2t2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specra2t2 = specr ; self.specia2t2 = speci
        return specr,speci
    def getspece2a2(self,ntrunc=2):
        sky = self.sky ; orf = self.orfe2a2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specre2a2 = specr ; self.specie2a2 = speci
        return specr,speci
    def getspece2e2(self,ntrunc=2):
        sky = self.sky ; orf = self.orfe2e2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specre2e2 = specr ; self.specie2e2 = speci
        return specr,speci
    def getspece2t2(self,ntrunc=2):
        sky = self.sky ; orf = self.orfe2t2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specre2t2 = specr ; self.specie2t2 = speci
        return specr,speci    
    def getspect2a2(self,ntrunc=2):
        sky = self.sky ; orf = self.orft2a2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrt2a2 = specr ; self.specit2a2 = speci
        return specr,speci
    def getspect2e2(self,ntrunc=2):
        sky = self.sky ; orf = self.orft2e2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrt2e2 = specr ; self.specit2e2 = speci
        return specr,speci
    def getspect2t2(self,ntrunc=2):
        sky = self.sky ; orf = self.orft2t2
        nlm = (ntrunc+1)*(ntrunc+2)/2
        numf = size(orf,2)
        specr = zeros((nlm,numf),complex64) ; speci = zeros((nlm,numf),complex64)
        for i in range(numf):
            specr[:,i]= sky.grdtospec(real(orf[:,:,i]),ntrunc)
            speci[:,i]= sky.grdtospec(imag(orf[:,:,i]),ntrunc)
        self.ntrunc = ntrunc
        self.specrt2t2 = specr ; self.specit2t2 = speci
        return specr,speci



                                                                                        
    def getAngularPowerAE(self):
        p = self.specrAE
        q = self.speciAE
        specindx = self.pairspecindx()
        sense = zeros((self.ntrunc+1,len(self.f)))
        for l in range(self.ntrunc+1):
            for m in range(l+1):
                i = specindx.index((m,l))
                if m==0:
                    dsigl = abs(p[i,:])**2 + abs(q[i,:])**2 + \
                        1j*(q[i,:]*conj(p[i,:])-p[i,:]*conj(q[i,:]))
                else:
                    dsigl = 2*(abs(p[i,:])**2 + abs(q[i,:])**2)
                sense[l,:] += dsigl/(2*l+1)
        apower = sqrt(abs(sense))
        self.apowerAE = apower
        return apower

    def getAngularPowerAEm(self):
        p = self.specrAEm
        q = self.speciAEm
        specindx = self.pairspecindx()
        sense = zeros((self.ntrunc+1,len(self.f)))
        for l in range(self.ntrunc+1):
            for m in range(l+1):
                i = specindx.index((m,l))
                if m==0:
                    dsigl = abs(p[i,:])**2 + abs(q[i,:])**2 + \
                        1j*(q[i,:]*conj(p[i,:])-p[i,:]*conj(q[i,:]))
                else:
                    dsigl = 2*(abs(p[i,:])**2 + abs(q[i,:])**2)
                sense[l,:] += dsigl/(2*l+1)
        apower = sqrt(abs(sense))        
        self.apowerAEm = apower
        return apower

        
    def pairspecindx(self):
        m,l = getspecindx(self.ntrunc)
        specindx = zip(m,l)
        self.specindx = specindx
        return specindx




def getMLvec(ntrunc,m='pn'):
    if ntrunc==0:
        M = [0]
        L = [0]
    elif m=='pn':    
        m,l = getspecindx(ntrunc)
        M = list(m) + list(-m[ntrunc+1:])
        L = list(l) + list(l[ntrunc+1:])
    elif m=='p':
        M,L = getspecindx(ntrunc)
    return zip(M,L)


def ORFfromFiles(realrpath, realipath, imagrpath, imagipath, sign='pn'):
    if not (sign=='p' or sign=='pn'):
        raise ValueError,"Keyword parameter sign must either be 'p' or 'pn'."

    realr = load(realrpath)
    reali = load(realipath)
    imagr = load(imagrpath)
    imagi = load(imagipath)
    try:
        assert shape(realr)==shape(reali)==shape(imagr)==shape(imagi)
    except AssertionError:
        print "The four arrays uploaded need to have the same dimensions."

    ntrunc = int( -1.5 + sqrt( 8*shape(realr)[0] + 1 ) / 2  ) 
    numf = shape(realr)[1]

    if (sign == 'p'):
        indxpn = getMLvec(ntrunc,'p')
        p = realr + 1j*reali
        q = imagr + 1j*imagi
        g = zeros( ( len(indxpn),numf ) , complex64)
        g = p + 1j*q
    else:
        indxp  = getMLvec(ntrunc,'p')
        indxpn = getMLvec(ntrunc)
        p = realr + 1j*reali
        q = imagr + 1j*imagi
        g = zeros( ( len(indxpn),numf ) , complex64)
        for i,ml in enumerate(indxpn):
            m = ml[0]
            l = ml[1]
            if m >= 0:
                k = indxp.index(ml)
                g[i,:] = p[k,:] + 1j*q[k,:]
            else:
                ml = (-m,l)
                k = indxp.index(ml)
                g[i,:] = (-1)**m*( conj(p[k,:]) + 1j*conj(q[k,:]) )

    return indxpn,g




def PSD_pm( f ) :
    """
    returns noise spectral density of standatd proof-mass noise
    INPUT:
    f --- frequency/frequencies
    OUTPUT:
    psdpm --- noise spectral density [Hz^{-1}]
    """
    psdpm = 2.54e-48 * f**-2
    return psdpm

def PSD_op( f ) :
    """
    returns noise spectral density of standard optical path noise
    INPUT:
    f --- frequency/frequencies
    OUTPUT:
    psdop --- noise spectral density [Hz^{-1}]
    """
    psdop = 1.76e-37 * f**2
    return psdop


def PSD_Xm( f ) :
    """
    returns noise spectral density of X, the L-closed, Gm Michelson TDI observable (8 links)
    INPUT:
    f --- frequency/frequencies
    OUTPUT:
    psd --- noise spectral density [Hz^{-1}]
    """
    L , w = 16.67 , 2*numpy.pi*f
    Spm , Sop = PSD_pm( f ) , PSD_op( f )
    psd = 16*numpy.sin( L*w )**2 * ( 2*( 1 + numpy.cos( L*w )**2 )*Spm + Sop )
    return psd

def PSD_Ym( f ) :
    """
    returns noise spectral density of Y, the L-closed, Gm Michelson TDI observable (8 links)
    INPUT:
    f --- frequency/frequencies
    OUPUT:
    psd --- noise spectral density [Hz^{-1}]
    """
    psd = PSD_Xm( f )
    return psd

def PSD_Zm( f ) :
    """
    returns noise spectral density of Z, the L-closed, Gm Michelson TDI observable (8 links)
    INPUT:
    f --- frequency/frequencies
    OUPUT:
    psd --- noise spectral density [Hz^{-1}]
    """
    psd = PSD_Xm( f )
    return psd    

def PSD_X1( f ) :
    """
    returns noise spectral density of X1, the Ldot-closed, G2 Michelson TDI observable (16 links)
    INPUT:
    f --- frequency/frequencies
    OUPUT:
    psd --- noise spectral density [Hz^{-1}]
    """
    L , w = 16.67 , 2 * numpy.pi * f
    psd = 4 * numpy.sin( 2*L*w )**2 * PSD_Xm( f )
    return psd

def PSD_X2( f ) :
    """
    returns noise spectral density of X2, the Ldot-closed, G2 Michelson TDI observable (16 links)
    INPUT:
    f --- frequency/frequencies
    OUPUT:
    psd --- noise spectral density [Hz^{-1}]
    """
    psd = PSD_X1( f ) 
    return psd

def PSD_X3( f ) :
    """
    returns noise spectral density of X3, the Ldot-closed, G2 Michelson TDI observable (16 links)
    INPUT:
    f --- frequency/frequencies
    OUPUT:
    psd --- noise spectral density [Hz^{-1}]
    """
    psd = PSD_X1( f ) 
    return psd


def PSD_A2( f ) :
    """
    returns noise spectral density of A, G2
    INPUT:
    f --- frequency/frequencies
    OUTPUT:
    S --- noise spectral density [Hz^{-1}]
    """
    L , w = 16.67 , 2*numpy.pi*f
    Spm , Sop = PSD_pm( f ) , PSD_op( f )
    S = 32*numpy.sin(.5*L*w)**2 * numpy.sin(1.5*L*w)**2 * ( (6 + 4*numpy.cos(L*w) + 2*numpy.cos(2*L*w))*Spm
                                                            + (2 + numpy.cos(L*w))*Sop )
    return S


def PSD_E2( f ) :
    """
    returns noise spectral density of E, G2
    INPUT:
    f --- frequency/frequencies
    OUTPUT:
    S --- noise spectral density [Hz^{-1}]
    """
    S = PSD_A2(f)
    return S


def PSD_T2( f ) :
    """
    returns noise spectral density of T, G2
    INPUT:
    f --- frequency/frequencies
    OUTPUT:
    S --- noise spectral density [Hz^{-1}]
    """
    L , w = 16.67 , 2*numpy.pi*f
    Spm , Sop = PSD_pm( f ) , PSD_op( f )
    S = 8* (1+2*numpy.cos(L*w))**2 * numpy.sin(1.5*L*w)**2 * ( 4*numpy.sin(.5*L*w)**2 * Spm + Sop )
    return S
