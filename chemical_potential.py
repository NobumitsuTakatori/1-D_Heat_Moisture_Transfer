# -*- coding: utf-8 -*-
import numpy as np
import vapour_pressure

Rv_num = 8316.96/18.016

def Miu_cal(Tem,RH):
    Rv = Rv_num
    Miu= Rv*Tem*np.log(RH)
    return Miu 

def DPDM_cal(Tem,RH):
    Miu= Miu_cal(Tem,RH)
    Rv = Rv_num
    Pvs = vapour_pressure.Pvs_cal(Tem)
    DPDM=Pvs*np.exp(Miu/Rv/Tem)/Rv/Tem
    return DPDM

def DPDT_cal(Tem,RH):
    Miu= Miu_cal(Tem,RH)
    Rv = Rv_num
    Pvs = vapour_pressure.Pvs_cal(Tem)
    DPvs= vapour_pressure.DPvs_cal(Tem)    
    DPDT= np.exp(Miu/Rv/Tem)*(DPvs-Pvs*Miu/Rv/(Tem**2))
    return DPDT

def RH_cal(Tem,Miu):
    Rv = Rv_num
    RH = np.exp(Miu/Rv/Tem)
    return RH
