import numpy as np
import matplotlib.pyplot as plt

class TemperatureField:
    def __init__(self,D,t,layer_num,time,curve):
        self.D=D
        self.t=t
        self.layer_num=layer_num
        self.time=time
        self.r_steel=D/2
        self.w_steel=t
        self.w_conc=np.ones(layer_num)
        self.r_conc=np.ones(layer_num)
        self.curve=curve

    def heatTransfer(self):
        self.heattransfer=35
        self.emmissFire=1
        self.emmissSteel=0.2
        self.bolzman=5.67e-8
        self.steelDensity=7850


    def timeSplit(self):
        if self.t<=5e-3:self.dt=0.1
        elif 5e-3<self.t<=7e-3:self.dt=0.25
        elif 7e-3<self.t<=12e-3:self.dt=0.5
        else:self.dt=1.
        self.n=int(60*self.time/self.dt)
        self.timeSegment=np.linspace(0,60*self.time,self.n+1)
        if self.curve!='ISO':self.ambientTemperature=20+750*(1-np.exp(-3.79553*(self.timeSegment/3600)**0.5))+\
                                                          170.1*(self.timeSegment/3600)**0.5
        else:self.ambientTemperature=20+345*np.log10(8*(self.timeSegment/60)+1)

        
    def heatCalculation(self):
        self.heatTransfer()
        self.partitionLayer()
        self.timeSplit()
        water_content=5
        if 0<=water_content<1.5:cp=380*water_content+900
        elif 1.5<=water_content<3:cp=366.66667*water_content+920
        else:cp=511.4286*water_content+485.7143
        self.t_steel=20*np.ones(self.n+1)
        self.hc=np.ones(self.n)
        self.hr=np.ones(self.n)
        self.hnet=np.ones(self.n)
        self.specific_steel=np.ones(self.n)
        self.conductance_steel=np.ones(self.n)

        self.t_conc=20*np.ones([self.layer_num,self.n+1])
        self.conc_density=2300*np.ones([self.layer_num,self.n])
        self.conductance_conc=np.ones([self.layer_num,self.n])
        self.specific_conc=np.ones([self.layer_num,self.n])
        for i in range(self.n):
            self.hc[i]=self.heattransfer*(self.ambientTemperature[i+1]-self.t_steel[i])
            self.hr[i]=self.emmissSteel*self.emmissFire*self.bolzman*((self.ambientTemperature[i+1]+273)**4-\
                                                                      (self.t_steel[i]+273)**4)
            self.hnet[i]=self.hc[i]+self.hr[i]
            self.specific_steel[i]=450+0.28*self.t_steel[i]-2.91*(self.t_steel[i])**2*1e-4+1.34*((self.t_steel[i])**3)*1e-7
            self.conductance_steel[i]=14.6+1.27*self.t_steel[i]/100

            for m in range(self.layer_num):
                self.conductance_conc[m,i]=2.34-0.272*(self.t_conc[m,i]/100)+0.0112*(self.t_conc[m,i]/100)**2
                if 20<=self.t_conc[m,i]<=100:
                    self.specific_conc[m,i]=900
                elif 100<self.t_conc[m,i]<=115:
                    self.specific_conc[m,i]=(cp-900)*(self.t_conc[m,i]-100)/15+900
                elif 115<self.t_conc[m,i]<=200:
                    self.specific_conc[m,i]=(1000-cp)*(self.t_conc[m,i]-200)/85+1000
                elif 200<self.t_conc[m,i]<=400:
                    self.specific_conc[m,i]=1000+(self.t_conc[m,i]-200)/2
                else:
                    self.specific_conc[m,i]=1100
            self.t_steel[i+1]=self.t_steel[i]+(self.r_steel/((self.r_steel+self.r_steel)/2))*((1/self.w_steel)*self.hnet[i]*self.dt/(self.specific_steel[i]*self.steelDensity))-\
            (self.dt*self.conductance*self.r_conc[0]*(self.t_steel[i]-self.t_conc[0,i]))/((self.specific_steel[i]*self.steelDensity*self.w_steel)*((self.r_steel+self.r_conc[0])/2))
            
            self.t_conc[0,i+1]= self.t_conc[0,i]+(self.conductance*self.r_conc[0]*(self.t_steel[i]- self.t_conc[0,i])*self.dt-\
                        self.dt*self.r_conc[1]*((self.conductance_conc[0,i]+self.conductance_conc[1,i])/2)*(self.t_conc[0,i]-\
                        self.t_conc[1,i])/(self.w_conc[0]+self.w_conc[1]/2))/(self.specific_conc[0,i]*self.conc_density[0,i]*self.w_conc[0]*((self.r_conc[0]+self.r_conc[1])/2))

            for m in range(1,self.layer_num-1):
                self.t_conc[m,i+1]=self.t_conc[m,i]+(self.dt*self.r_conc[m]*((self.conductance_conc[m-\
                1,i]+self.conductance_conc[m,i])/2)*(self.t_conc[m-1,i]-self.t_conc[m,i])/(self.w_conc[m-1]/2+self.w_conc[m]/2)-\
                self.dt*self.r_conc[m+1]*((self.conductance_conc[m,i]+self.conductance_conc[m+1,i])/2)*(self.t_conc[m,i]-\
                self.t_conc[m+1,i])/(self.w_conc[m]/2+self.w_conc[m+1]/2))/(self.specific_conc[m,i]*self.conc_density[m,i]*\
                                                                         ((self.r_conc[m]+self.r_conc[m+1])/2)*self.w_conc[m])
            ele_num=self.layer_num-1
            self.t_conc[ele_num,i+1]=self.t_conc[ele_num,i]+(self.dt*2*((self.conductance_conc[ele_num-\
            1,i]+self.conductance_conc[ele_num,i])/2)*(self.t_conc[ele_num-1,i]- self.t_conc[ele_num,i])/(self.w_conc[ele_num-1]\
                /2+self.w_conc[ele_num]))/(self.specific_conc[ele_num,i]*self.conc_density[ele_num,i]*(self.w_conc[ele_num]))
            
    def linearInterpolation(self,xs,ys,x):
        n=len(xs)
        for i in range(n):
            if x>xs[i]:
                continue
                #print(xs[i],ys[i],xs[i-1],ys[i-1])
                
            y=ys[i-1]*(x-xs[i])/(xs[i-1]-xs[i])+ys[i]*(x-xs[i-1])/(xs[i]-xs[i-1])
                #print(y)
            return y
        
    def reductionCoefficient(self):
        self.xs=[20,100,200,300,400,500,600,700,800,900,1000,1100,1200]
        self.fys=[1,1,1,1,1,0.78,0.47,0.23,0.11,0.06,0.04,0.02,0]
        self.Ess=[1,1,0.9,0.8,0.7,0.6,0.31,0.13,0.09,0.0675,0.045,0.0225,0]
        self.Ts=self.t_steel[-1]
        self.Tc=np.ones(self.layer_num)
        
        for i in range(self.layer_num):
            self.Tc[i]=self.t_conc[i,-1]
        fy=self.linearInterpolation(self.xs,self.fys,self.Ts)
        Es=self.linearInterpolation(self.xs,self.Ess,self.Ts)
        fcs_siliceous=[1,0.95,0.9,0.85,0.75,0.6,0.45,0.3,0.15,0.08,0.04,0.01,0]
        fcs_calcerous=[1,0.97,0.94,0.91,0.85,0.74,0.6,0.43,0.27,0.15,0.06,0.02,0]
        epsilons=[1e-3*i for i in [2.5,3.5,4.5,6,7.5,9.5,12.5,14,14.5,15,15,15,15]]
        fct=[]
        epsilont=[]
        for T in self.Tc:
            fc=self.linearInterpolation(self.xs,fcs_siliceous,T)
            fct.append(fc)
            epsilon=self.linearInterpolation(self.xs,epsilons,T)
            epsilont.append(epsilon)
        return fy,Es,fct,epsilont

    def temperatureInterpolation(self,s):
        T=self.linearInterpolation(2*self.r_conc,self.Tc,s)
        return T

    def steelReduction(self,T):
        fy=self.linearInterpolation(self.xs,self.fys,T)
        return fy

    def plotField(self):
        fig=plt.figure()
        ax=fig.add_subplot(111,polar=True)
        j=0
        for i in range(len(self.r_conc)-1):
            j+=1
            ax.fill_between(np.linspace(0,2*np.pi,100),self.r_conc[i]*np.ones(100),self.r_conc[i+1]\
                            *np.ones(100),color=(j*self.r_conc[i],j*self.r_conc[i],j*self.r_conc[i]),\
                            linestyle='-')
        ax.fill(np.linspace(0,2*np.pi,100),self.r_conc[-1]*np.ones(100),\
                color=(self.r_conc[-1],self.r_conc[-1],self.r_conc[-1]))
        plt.show()



class CircularTemperatureField(TemperatureField):
    
    def __init__(self,B,t,layer_num,time,curve='ISO'):
        super().__init__(B,t,layer_num,time,curve)
        if self.D<=0.3:
            self.conductance=516*(10*self.D)**(-2.373)
        else:self.conductance=38.1
        self.heatCalculation()
    def partitionLayer(self):
        self.w_conc[0]=(self.D/2-self.t)/(self.layer_num-1)/2
        self.w_conc[self.layer_num-1]=self.w_conc[0]
        self.r_conc[0]=(self.D/2-self.t)
        for i in range(1,self.layer_num-1):
            self.w_conc[i]=((self.D/2-self.t)-self.w_conc[0]-self.w_conc[self.layer_num-1])/(self.layer_num-2)
            self.r_conc[i]=self.r_conc[i-1]-self.w_conc[i-1]
        self.r_conc[self.layer_num-1]=self.r_conc[self.layer_num-2]-self.w_conc[self.layer_num-2]

class SquareTemperatureField(TemperatureField):

    def __init__(self,B,t,layer_num,time,curve='ISO'):
        super().__init__(B,t,layer_num,time,curve='ISO')
        if self.D<=0.3:
            self.conductance=115*(10*self.D)**(-0.85)
        else:self.conductance=45.2
        self.heatCalculation()
        
    def partitionLayer(self):
        self.w_conc[0]=(self.D/2-self.t)/self.layer_num
        self.w_conc[self.layer_num-1]=self.w_conc[0]
        self.r_conc[0]=(self.D/2-self.t)
        self.b_conc=np.ones(self.layer_num)
        self.b_conc[0]=self.r_conc[0]
        for i in range(1,self.layer_num-1):
            self.w_conc[i]=((self.D/2-self.t)-self.w_conc[0]-self.w_conc[self.layer_num-1])/(self.layer_num-2)
            self.r_conc[i]=self.r_conc[i-1]-self.w_conc[i-1]
            self.b_conc[i]=self.b_conc[i-1]-self.w_conc[i-1]
        self.r_conc[self.layer_num-1]=self.r_conc[self.layer_num-2]-self.w_conc[self.layer_num-2]
        self.b_conc[self.layer_num-1]=self.b_conc[self.layer_num-2]-self.w_conc[self.layer_num-2]

