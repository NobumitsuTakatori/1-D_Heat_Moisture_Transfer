########################################################
import numpy as np

def Sl_cal_vanGenuchten(Alfa, Miu, m, n):
    Sl  = (1.0+(-Alfa *Miu)**n)**(-m)
    return Sl 

def Kl_cal_vanGenuchten(Alfa,Miu,m,n,l):
    Sl = Sl_cal_vanGenuchten(Alfa ,Miu ,m ,n)
    Kl = (Sl**l)*((1.0-(1.0-Sl**(1.0/m))**m)**2.0)
    return Kl

########################################################    

def Lamdml_cal_vanGenuchten(Ksat,Alfa,Miu,m,n,l):
    Kl = Kl_cal_vanGenuchten(Alfa,Miu,m,n,l)
    row = 1000.0        #塩溶液では密度変わるのでは
    g = 9.8
    Lamdml = Ksat*Kl*row/g
    return Lamdml

def Phi_vanGenuchten(Phimax, Alfa, Miu, m, n):
    Sl = Sl_cal_vanGenuchten(Alfa, Miu, m, n)
    Phi = Phimax *Sl
    return Phi

def DPhi_vanGenuchten(Phimax, Alfa, Miu, m, n):
    Sl = Sl_cal_vanGenuchten(Alfa, Miu, m, n)
    DPhi = -(Alfa*m*Phimax)/(1.0-m)*(Sl**(1.0/m))*((1.0-Sl**(1.0/m))**m)
    DPhi = np.abs(DPhi)
    return DPhi

def Miu_vanGenuchten(Phimax, Phi, Alfa, m, n):
    Sl = Phi/Phimax
    Miu = - (((Sl**(-1.0/m))-1.0)**(1.0/n))/Alfa
    return Miu
