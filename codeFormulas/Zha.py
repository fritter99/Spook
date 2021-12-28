import numpy as np

class ZhaFrr:
    def __init__(self,D,T,L,fy,fc,t,is_circular=True,bc=1):
        self.D=D
        self.T=T#钢管厚度
        self.t=t#时间
        self.fy=fy
        self.fc=fc
        self.is_circular=is_circular
        d=D-2*T
        self.Ec=fc**0.5*4700
        self.Es=2e5
        if is_circular:
            self.Ac=d**2*3.14/4
            self.As=D**2*3.14/4-self.Ac
            self.Is=3.14/64*(D**4-d**4)
            self.Ic=3.14/64*d**4
            self.Isc=3.14/64*D**4
            if bc==1:
                self.lambda_sc=2*L/D
            elif bc==0:
                self.lambda_sc=2.8*L/D
            else:
                self.lambda_sc=4*L/D
        else:
            self.Ac=d**2
            self.As=D**2-self.Ac
            self.Is=(D**4-d**4)/12
            self.Ic=d**4/12
            self.Isc=D**4/12
            if bc==1:
                self.lambda_sc=3**0.5*L/D
            elif bc==0:
                self.lambda_sc=1.4*3**0.5*L/D
            else:
                self.lambda_sc=2*3**0.5**L/D
    def get_temperature(self):
        t=self.t
        Ac=self.Ac
        As=self.As
        ds=((As+Ac)/3.14)**0.5-(Ac/3.14)**0.5
        A=1200
        B=20.22+0.51*ds
        C=0.996+0.014*ds
        self.Ts=A*(1-1/(1+(t/B)**C))+20
        Le=(Ac/3.14)**0.5
        A=120+1080*np.exp(-.00447*Le)
        B=20.22+0.51*ds+1.8*Le*(Le**2*1e-6-0.00146*Le+0.64)
        C=0.996+0.014*ds
        self.Tc=2*A*(1-1/(1+(t/B)**C))+20
        #print(ds,Le,Ts,Tc)

    def get_property(self):
        self.get_temperature()
        T=self.Ts
        fy=self.fy
        fc=self.fc
        Es=self.Es
        Ec=self.Ec
        self.fyt=fy*np.exp(-((T-20)/601)**2.5)
        self.Est=Es*np.exp(-((T-20)/652)**3)
        T=self.Tc
        self.Ect=self.Ec*np.exp(-(T-20)/300)
        self.Esct=(self.Ect*self.Ic+self.Est*self.Is)/self.Isc
    def get_fsct(self):
        self.get_property()
        Ac=self.Ac
        As=self.As
        fy=self.fy
        fc=self.fc
        if self.Tc>938:self.fct=0
        else:self.fct=fc*(1-(self.Tc-20)/918)
        fyt=self.fyt
        fct=self.fct
        if not self.is_circular:
            B_=0.131*fy/213+0.723
            C=-0.07*fc/14.4+0.026
        else:
            B_=0.176*fy/213+0.974
            C=-0.104*fc/14.4+0.031
        theta=(As/Ac)*fy/fc
        self.fsc=(1.212+B_*theta+C*theta**2)*fc
        ksct=(Ac*fct+As*fyt)/(Ac*fc+As*fy)
        self.fsct=ksct*self.fsc

    def get_NuT(self):
        self.get_fsct()
        self.N0T=self.fsct*(self.Ac+self.As)/1000
        lambda_sct=self.lambda_sc/3.14*(self.fsct/self.Esct)**0.5
        self.phi_t=1/(2*lambda_sct**2)*(lambda_sct**2+0.25*lambda_sct+1-((lambda_sct**2+0.25*lambda_sct+1)**2-4*lambda_sct**2)**0.5)
        self.NuT=self.phi_t*self.N0T
