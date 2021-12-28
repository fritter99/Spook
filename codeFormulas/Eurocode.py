import numpy as np

'''
Unit:
    strength: Mpa
    length  : m

D: tube diameter
B: tube side length
t: tube thickness
fc: concrete strength
fs: steel strength
L: column length
flange: angle flange
angleThickness: angle thickness
fa: angle strength
fy: hooping strength
s: vertical spacing of hooping
s0: clear spacing of hooping
Es: modulus of steel
Ec: modulus of concrete
time:FRR
buckling:buckling factor
'''

class Frr:
    def __init__(self,D,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es,Ec,time,buckling):
        #CFST parameters
        self.D=D
        self.t=t
        self.fc=fc
        self.fs=fs
        self.L=L
        self.Es=Es
        self.Ec=Ec
        self.buckling=buckling
        
        #angle parameters
        self.fa=fa
        self.flange=flange
        self.thickness=angleThickness
        self.Aa=(2*self.flange-self.thickness)*self.thickness*4

        #hooping parameters
        self.s=s
        if fy!=0:
            self.ds=D-2*t-2*s0-0.01
            #print(self.ds)
        else:
            self.ds=0
        self.s0=s0
        self.fy=fy
        #print(D,t,a)
        self.A=self.ds**2/4*3.14

        #FRR
        self.time=time
        
    def _fcc(self):
        A=7.85e-5
        Ke=(1-0.5*self.s0/self.ds)**2
        f1=Ke*2*self.fy*A/self.s/self.ds
        fs,fc=self._equivalentStrength()
        fcc=fc*(-1.254+2.254*(1+7.94*f1/fc)**0.5-2*f1/fc)
        return fcc


    
class EcFrr(Frr):
    def __init__(self,D,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es,Ec,time,buckling):
        super().__init__(D,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es,Ec,time,buckling)

    def _Ecfcc(self):
        A=7.85e-5
        Ke=(1-0.5*self.s0/self.ds)**2
        f1=Ke*2*self.fy*A/self.s/self.ds
        for i in range(len(self.d_conc)):
            if self.d_conc[i]<self.ds:
                fc=self.fct[i]
                self.fct[i]=fc*(-1.254+2.254*(1+7.94*f1/fc)**0.5-2*f1/fc)
        

    def _EcEulerStress(self):
        self.Est=self.Es*self.r1
        self.r3=0.27*np.array(self.r3)
        self.Ect=0.4*self.fct/self.r3
        self.Ncr=1000*3.14**2*(self.Est*self.Is+(self.Ect*self.layer_moments).sum())/((self.buckling*self.L)**2)
        
    def _EcN0T(self):
        #fs,fc=self._equivalentStrength()
        self.fst=self.fs*self.r0
        self.fct=np.array([self.fc*i for i in self.r2])
        if self.fy!=0:
            self._Ecfcc()
            self.N0T=1000*(self.fst*self.As+self.fa*self.Aa+\
                           (self.fct*self.layer_areas).sum())
        else:
            #self.N0T=1000*(fs*self.As+fc*self.Ac+self.fa*self.Aa)
            self.N0T=1000*(self.fst*self.As+(self.fct*self.layer_areas).sum()+self.fa*self.Aa)
        #print(self.N0T)

    def _Eclambda(self):
        self._EcEulerStress()
        _lambda=(self.N0T/self.Ncr)**0.5
        return _lambda

    def _Ecphi(self):
        _lambda=self._Eclambda()
        Phi=0.5*(_lambda**2+0.49*(_lambda-0.2)+1)
        phi=1/(Phi+(Phi**2-_lambda**2)**0.5)
        return phi

    def ECNuT(self):
        self._EcN0T()
        phi=self._Ecphi()
        NuT=phi*self.N0T
        return NuT
    
class UnifiedFrr(Frr):
    
    def __init__(self,D,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es,Ec,time,buckling):
        super().__init__(D,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es,Ec,time,buckling)

    def _UnifiedN0T(self):
        fs,fc=self._equivalentStrength()
        ksi=self.As*fs/self.Ac/fc
        eta=0.5*self.ke*ksi/(1+ksi)
        if self.fy!=0:
            fcc=self._fcc()
            #print(fc,self.ds,self.Ac,self.A,fcc)
            #print(1000*fc*(self.Ac-self.A))
            #print(1000*fcc*self.A)
            self.N0T=1000*(fs*self.As+fc*(self.Ac-self.A)+self.fa*self.Aa+fcc*self.A)
            #print(self.N0T)
            self.N0T*=(1+eta)
        else:
            self.N0T=(1+eta)*1000*(fs*self.As+fc*self.Ac+self.fa*self.Aa)
            #print(1000*self.fa*self.Aa)
    
    def _equivalentTemperature(self):
        d_=((self.Ac+self.As)/3.14)**0.5-(self.Ac/3.14)**0.5
        A=1200
        B=0.337+8.5*d_
        C=0.996+14*d_
        Ts=A*(1-1/(1+(self.time/60/B)**C))+20
        L_=(self.Ac/3.14)**0.5
        A=120+1080*np.exp(-4.47*L_)
        B=0.337+8.5*d_+30*L_*(L_**2-1.46*L_+0.64)
        C=0.996+14*d_
        Tc=2*A*(1-1/(1+(self.time/60/B)**C))+20
        return Ts,Tc
    
    def _equivalentStrength(self):
        Ts,Tc=self._equivalentTemperature()
        r1=np.exp(-((Tc-20)/622)**2.5)
        #if Tc<450:r1=1
        #else:r1=2.011-2.353*(Tc-20)/1000
        if Ts<=400:
            fs=self.fs
        else:
            r0=np.exp(-(((Ts-400)/240)**1.5))
            fs=self.fs*r0
        fc=self.fc*r1
        #print(r0,r1,r2,r3)
        return fs,fc

    def _equivalentModulus(self):
        Ts,Tc=self._equivalentTemperature()
        r0=np.exp(-((Ts-20)/560)**2.5)
        r1=np.exp(-(Tc-20)/211)
        Es=self.Es*r0
        Ec=self.Ec*r1
        return Es,Ec

    def _UnifiedEulerStress(self):
        Es,Ec=self._equivalentModulus()
        #print(Es,Ec,self.Is,self.Ic)
        Ncr=1000*3.14**2*(Es*self.Is+Ec*self.Ic)/((self.buckling*self.L)**2)
        return Ncr

    def _Unifiedlambda(self):
        Ncr=self._UnifiedEulerStress()
        _lambda=(self.N0T/Ncr)**0.5
        return _lambda        

    def _Unifiedphi(self):
        K=0.25-0.09*self.ke
        _lambda=self._Unifiedlambda()
        phi=(_lambda**2+K*_lambda+1-np.sqrt((_lambda**2+K*_lambda+1)**2\
            -4*_lambda**2))/(2*_lambda**2)
        #print(phi)
        return phi

    def UNIFIEDNuT(self):
        self._UnifiedN0T()
        phi=self._Unifiedphi()
        NuT=phi*self.N0T
        return NuT
    
class CircularCFST(EcFrr,UnifiedFrr):
    def __init__(self,D,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es=2e5,Ec=3e4,time=4.5,buckling=0.5):
        super().__init__(D,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es,Ec,time,buckling)
        self.Ac=(D-2*t)**2*3.14/4
        self.As=D**2*3.14/4-self.Ac
        self.Ic=3.14*(D-2*t)**4/64
        self.Is=3.14*D**4*(1-((D-2*t)/D)**4)/64
        self.ke=1

    def EcNuT(self,layer_num,curve='ISO',a=0):
        from . import TemperatureField as TF
        tf=TF.CircularTemperatureField(self.D,self.t,layer_num,self.time,curve)
        self.r0,self.r1,self.r2,self.r3=tf.reductionCoefficient()
        self.Tc=tf.Tc
        self.Ts=tf.Ts
        if self.fa!=0:
            T=tf.temperatureInterpolation((self.D/2-2*self.t)-a)
            self.fa*=tf.steelReduction(T)
        if self.fy!=0:
            T=tf.temperatureInterpolation((self.D/2-2*self.t)-self.s0)
            self.fy*=tf.steelReduction(T)            
        def layerArea(d1,d2):
            return (d1**2-d2**2)*3.14/4
        def layerSecondmoment(d1,d2):
            alpha=d2/d1
            return 3.14*d1**4*(1-alpha**4)/64
        self.d_conc=2*tf.r_conc
        self.layer_areas,self.layer_moments=[],[]
        for i in range(len(self.d_conc)-1):
            d1,d2=self.d_conc[i],self.d_conc[i+1]
            self.layer_areas.append(layerArea(d1,d2))
            self.layer_moments.append(layerSecondmoment(d1,d2))
        self.layer_areas.append(self.d_conc[-1]**2*3.14/4)
        self.layer_moments.append(self.d_conc[-1]**4*3.14/64)
        self.layer_areas=np.array(self.layer_areas)
        self.layer_moments=np.array(self.layer_moments)
        NuT=self.ECNuT()
        return NuT

    def UnifiedNuT(self):
        NuT=self.UNIFIEDNuT()
        return NuT



class SquareCFST(EcFrr,UnifiedFrr):
    def __init__(self,B,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es=2e5,Ec=3e4,time=4.5,buckling=0.5):
        super().__init__(B,t,fc,fs,L,flange,angleThickness,fa,fy,s,s0,Es,Ec,time,buckling)
        self.Ac=(B-2*t)**2
        self.As=B**2-self.Ac
        self.Ic=(B-2*t)**4/12
        self.Is=B**4/12-self.Ic
        self.ke=1/3
        
    def EcNuT(self,layer_num,curve='ISO',a=0):
        import TemperatureField as TF
        tf=TF.SquareTemperatureField(self.D,self.t,layer_num,self.time,curve)
        self.r0,self.r1,self.r2,self.r3=tf.reductionCoefficient()
        self.Tc=tf.Tc
        self.Ts=tf.Ts
        if self.fa!=0:
            T=tf.temperatureInterpolation((self.D/2-2*self.t)-a)
            #self.fa*=tf.steelReduction(T)
        if self.fy!=0:
            T=tf.temperatureInterpolation(self.ds)
            #self.fy*=tf.steelReduction(T)             
        '''def layerArea(d1,d2):
            return (d1**2-d2**2)
        def layerSecondmoment(d1,d2):
            return d1**4/12-d2**4/12'''
        def layerArea(d1,d2):
            return (d1**2-d2**2)*3.14/4
        def layerSecondmoment(d1,d2):
            alpha=d2/d1
            return 3.14*d1**4*(1-alpha**4)/64        
        self.d_conc=2*tf.r_conc
        self.layer_areas,self.layer_moments=[],[]
        for i in range(len(self.d_conc)-1):
            d1,d2=self.d_conc[i],self.d_conc[i+1]
            self.layer_areas.append(layerArea(d1,d2))
            self.layer_moments.append(layerSecondmoment(d1,d2))
        self.layer_areas.append(self.d_conc[-1]**2*3.14/4)
        self.layer_moments.append(self.d_conc[-1]**4*3.14/64)
        self.layer_areas=np.array(self.layer_areas)
        self.layer_moments=np.array(self.layer_moments)
        NuT=self.ECNuT()
        return NuT

    def UnifiedNuT(self):
        NuT=self.UNIFIEDNuT()
        return NuT
