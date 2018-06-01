#import matplotlib.pyplot as plt
#import numpy as np
import van_Genuchten
import chemical_potential

##########################################
###     材料情報の入力        #############
Phimax = 0.23
Ksat   = 2.0e-7
LAMDP  = 2.0E-10

###     van-Genuchten用情報    ##########
Alfa_num = 10.0/98.0
n_num = 2.0
m_num = 1.0 -(1.0/n_num)
l_num = 0.5

########################################
#   移動係数
#   熱伝導率
def LAM_cal():
    LAM = 1.2
    return LAM

#   水分化学ポテンシャル勾配に関する液相水分伝導率
def LDML_cal(Miu):
    Alfa = Alfa_num
    n = n_num
    m = m_num
    l = l_num
    LDML = van_Genuchten.Lamdml_cal_vanGenuchten(Ksat,Alfa,Miu,m,n,l)
    return LDML

#   含水率勾配に関する液相水分伝導率


#   水分化学ポテンシャル勾配に関する気相水分伝導率
def LDMG_cal(Tem,Miu):
    Phi  = Phi_cal(Miu)
    LDP  = LDP_cal(Phi)
    RH   = chemical_potential.RH_cal(Tem,Miu) 
    DPDM = chemical_potential.DPDM_cal(Tem,RH)
    LDMG = LDP*DPDM
    return LDMG

#   温度勾配に関する気相水分伝導率
def LDTG_cal(Tem,Miu):
    Phi  = Phi_cal(Miu)
    LDP  = LDP_cal(Phi)
    RH   = chemical_potential.RH_cal(Tem,Miu) 
    DPDT = chemical_potential.DPDT_cal(Tem,RH)
    LDTG = LDP*DPDT
    return LDTG

#   水蒸気圧勾配に関する気相水分伝導率
def LDP_cal(Miu):
    Phi = Phi_cal(Miu)
    Sl = Phi/Phimax
    LDP = LAMDP*(1.0-Sl*0.9)    #注意
    return LDP

#######################################
#   容量系
#   含水率 from 水分化学ポテンシャル
def Phi_cal(Miu):
    Alfa = Alfa_num
    n = n_num
    m = m_num
    Phi  = van_Genuchten.Phi_vanGenuchten(Phimax,Alfa, Miu, m, n)
    return Phi

#   含水率のポテンシャル微分
def DPhi_cal(Miu):
    Alfa = Alfa_num
    n = n_num
    m = m_num
    DPhi = van_Genuchten.DPhi_vanGenuchten(Phimax,Alfa, Miu, m, n)
    return DPhi

#   水分化学ポテンシャル from 含水率
def Miu_cal(Phi):
    Alfa = Alfa_num
    n = n_num
    m = m_num
    Miu = van_Genuchten.Miu_vanGenuchten(Phimax, Phi, Alfa, m, n)
    return Miu

#   熱容量
def crow_cal(Miu):
    Phi = Phi_cal(Miu)
    crow = 1479.25*750.0 +1000.0 *Phi *4.18605E+3
    return crow    

#   水分容量(ポテンシャル使用時)
def cmiu_cal(Miu):
    roww = 1000.0
    dphi = DPhi_cal(Miu)
    cmiu = roww * dphi
    return cmiu

#######################################
###     物性値の確認     ###
###     グラフの描画      ###
#plt.xscale("log")
#plt.yscale("log")
#plt.grid(which="both")
#x  = np.arange(-1000000,0.0,1.0)
#y = Phi_cal(x)
#plt.plot(-x,y)
#plt.show()
###########################
