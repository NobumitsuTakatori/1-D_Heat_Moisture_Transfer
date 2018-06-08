# 温度とポテンシャルから水蒸気圧・飽和水蒸気圧を求める
import numpy as np
#import matplotlib.pyplot as plt

def Pv(Tem, Miu):
    Rv = 8316.96/18.016
    RH = np.exp(Miu/Rv/Tem)
    Pv = RH *Pvs(Tem)    
    return Pv
    
def Pvs(Tem):
    EW =-5800.22060/Tem + 1.3914993 -4.8640239E-2*Tem + 4.1764768E-5*(Tem**2.0)\
    -1.4452093E-8*(Tem**3.0)+6.5459673*np.log(Tem)
    return np.exp(EW) 

def DPvs(Tem):
    DP=10.795740*273.160/Tem/Tem-5.0280/Tem/np.log(10.0)\
    +(1.50475E-4)*8.2969/273.16*np.log(10.0)\
    *(10.0**(-8.29690*(Tem/273.160 -1.0)))\
    +(0.42873E-3)*4.769550*273.160/Tem/Tem*np.log(10.0)\
    *(10.0**(4.769550*(1.0 -273.160/Tem)))
    DPvs = Pvs(Tem) *DP*np.log(10.0)
    return DPvs

###     計算・グラフの描画      ###
#plt.xscale("log")
#plt.yscale("log")
#plt.grid(which="both")
#miu = np.array([-10000.0,-1000.0,-100.0,-10.0,-1.0])
#tem = np.array([293.15]*5)
#y  = Pv_cal(tem,miu)

###     グラフの描画    ###
#plt.plot(miu,y)
#plt.show()
###########################