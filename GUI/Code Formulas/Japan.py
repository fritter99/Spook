class JapanFrr:
    def __init__(self,D,fc,T,is_circular=True):
        if is_circular:
            self.Ac=D**2*3.14/4
        else:
            self.Ac=D**2
        self.fc=fc
        self.T=T#时间
        self.is_circular=is_circular

    def get_NfiRd(self):
        Ac=self.Ac
        fc=self.fc
        T=self.T
        is_circular=self.is_circular
        '''if fc<=36:
            if is_circular:
                NfiRd=Ac*fc*(5.75*1e-5*36**2.63*T+1)**(-.214)
            else:
                NfiRd=Ac*fc*(5.75*1e-5*fc**2.63*T+1)**(-.214)
        else:
            if is_circular:
                NfiRd=Ac*fc*(2.55*1e-5*36**1.735*T+1)**(-.225)
            else:
                NfiRd=Ac*fc*(2.55*1e-5*fc**1.735*T+1)**(-.225)'''
        if fc<=42:
            if is_circular:
                NfiRd=1.95*Ac*fc*(1/T)**0.313
            else:
                NfiRd=2.177*Ac*fc*(1/T)**0.367
        else:

            NfiRd=1.7*Ac*fc*(1/T)**0.318
        self.NfiRd=NfiRd/1000

