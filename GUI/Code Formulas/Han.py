import numpy as np

class HanFrr:

    def __init__(self,D,T,L,fy,fc,t,is_circular=True,bc=1):
        self.D=D
        self.fy=fy
        self.fc=fc
        if bc==1:L*=0.5
        elif bc==0:L*=0.7
        else:L=L
        if is_circular:
            self.C=D*3.14
            self.C0=self.C/400/3.14
            self._lambda=4*L/D
            self.Ac=(D-2*T)**2*3.14/4
            self.As=D**2*3.14/4-self.Ac
        else:
            self.C=D*4
            self.C0=self.C/1600
            self._lambda=2*3**0.5*L/D
            self.Ac=(D-2*T)**2
            self.As=D**2-self.Ac
        self.t0=t/100
        self._lambda0=self._lambda/40
        self.is_circular=is_circular

    def get_Nu(self):
        As=self.As
        Ac=self.Ac
        fy=self.fy
        fc=self.fc
        if self.is_circular:
            self.Nu=(1.14+1.02*As*fy/Ac/fc)*(As+Ac)*fc/1000
        else:
            self.Nu=(1.18+0.85*As*fy/Ac/fc)*(As+Ac)*fc/1000

    def get_kT(self):
        _lambda=self._lambda
        lambda0=self._lambda0
        C0=self.C0
        t0=self.t0
        if self.is_circular:
            a=(-0.13*lambda0**3+0.92*lambda0**2-0.39*lambda0+.74)*(-2.85*C0+19.45)
            b=(-1.59*lambda0**2+13*lambda0-3)*C0**(-.46)
            k=(-0.1*lambda0**2+1.36*lambda0+.04)*(0.0034*C0**3-0.0465*C0**2+0.21*C0-0.33)
            t1=(-0.0131*lambda0**3+0.17*lambda0**2-0.72*lambda0+1.49)*(.0072*C0**2-.02*C0+.27)
            t2=(.007*lambda0**3+0.209*lambda0**2-1.035*lambda0+1.868)*(.006*C0**2-.009*C0+.362)
            
            if t0<=t1:
                self.kT=1/(1+a*t0**2.5)
            elif t1<t0<=t2:
                self.kT=1/(1+a*t1**2.5+b*(t0-t1))
            else:
                self.kT=1/(1+a*t1**2.5+b*(t2-t1))+k*(t0-t2)

        else:
            a=(.015*lambda0**2-.025*lambda0+1.04)*(-2.56*C0+16.08)
            b=(-.19*lambda0**3+1.48*lambda0**2-.95*lambda0+.86)*(-.19*C0**2+.15*C0+9.05)
            k=.042*(lambda0**3-3.08*lambda0**2-.21*lambda0+.23)
            t1=.38*(.02*lambda0**3-.13*lambda0**2+.05*lambda0+.95)
            t2=(.03*lambda0**2-0.29*lambda0+1.21)*(.022*C0**2-.105*C0+.696)
            if t0<=t1:self.kT=1/(1+a*t0**2)
            elif t1<t0<=t2:self.kT=1/(b*t0**2+1+(a-b)*t1**2)
            else:self.kT=1/(b*t2**2+1+(a-b)*t1**2)+k*(t0-t2)

    def get_NfiRd(self):
        self.get_Nu()
        self.get_kT()
        self.NfiRd=self.Nu*self.kT
