import numpy as np
#################################################
#                   流量計算                    #
#################################################

grav = 9.806650
Mw = 18.02

#################################################
###     熱流     ################################
###     条件分岐
def Heat_Flux(AL,Lamda,Dtem,dx2,bcon):
    Qs = np.zeros_like(bcon)
    Qs = np.where(bcon==0, Flux_Qair(AL,Dtem), Qs)
    Qs = np.where(bcon==1, Flux_Qsol(Lamda,Dtem,dx2), Qs)
    return Qs   

###     固体内（熱伝導率）
def Flux_Qsol(Lamda,Dtem,dx2):
    Qsol = Lamda *Dtem/dx2
    return Qsol

###     材料表面（熱伝達）
def Flux_Qair(AL, Dtem):
    Qair = AL *Dtem       #流出を正
    return Qair

###     相変化熱
def Ratent_heat_cal(tem):
    r = (597.5-0.559*(tem-273.16))*4186.05
    return r


#################################################
###     水蒸気流    #############################
    
####################################
###     圧力勾配（条件分岐）
def Flux_vapour_pressure(ALD,LDP,DPv,dx2,bcon):
    Jv = np.zeros_like(bcon)
    Jv = np.where(bcon==0, Flux_Jvair_pressure(ALD,DPv), Jv)
    Jv = np.where(bcon==1, Flux_Jvsol_pressure(LDP,DPv,dx2), Jv)
    return Jv

###     材料内
def Flux_Jvsol_pressure(LDP,DPv,dx2):
    Jvsol = LDP *DPv/dx2
    return Jvsol

###     材料表面
def Flux_Jvair_pressure(ALD,DPv):
    Jvair = ALD *DPv         #流出を正
    return Jvair

####################################
###     ポテンシャル勾配
def Flux_vapour_potential(LDMG,LDTG,Dmiu,Dtem,dx2,nx):
    g = grav
    Jvsol = LDMG *(Dmiu/dx2-nx*g) +LDTG*Dtem/dx2
    return Jvsol

####################################
###     水蒸気密度勾配
def Dg_cal(pg,Tem):
    p0 = 101325.0
    Dg = 2.17e-3*p0/pg*((Tem/273.15)**1.88)
    return Dg
def Ft_cal(Phi,Sl):
    x = 1.4
    Ft = (Phi**(2*x))*((1-Sl)**(2*x+2))
    return Ft
def Flux_Jvsol_row(Ft,Dg,RowvA,RowvB,RowgA,RowgB):
    Jv = -RowgA/Mw *Ft *Dg *(RowvA-RowvB)/(RowgA-RowgB)
    return Jv


#################################################
###     液相水分流    ###########################

####################################
###     ポテンシャル勾配（条件分岐）
def Flux_liquid_potential(LDML,Dmiu,Dtem,dx2,nx,bcon):
    Jl = np.zeros_like(bcon)
    Jl = np.where(bcon==0, Flux_Jlair(), Jl)
    Jl = np.where(bcon==1, Flux_Jlsol_potential(LDML,Dmiu,Dtem,dx2,nx), Jl)
    return Jl

###     材料表面 = 0 (液相水分は流れない)
def Flux_Jlair():
    Jlair = 0.0
    return Jlair

###     材料内
###     含水率勾配
def Flux_Jlsol_moisture_content(Kl,Dphi,dx2):
    Jvsol = Kl *Dphi/dx2
    return Jvsol    

###     ポテンシャル勾配
def Flux_Jlsol_potential(LDML,Dmiu,Dtem,dx2,nx):
    g = grav
    Jlsol = LDML *(Dmiu/dx2-nx*g)
    return Jlsol

###     圧力勾配：ダルシー則
def Flux_Jsol_darcy(LDPS,PlA,PlB,dx2):
    Jsol= LDPS*(PlA-PlB)/dx2
    return Jsol







