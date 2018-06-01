import numpy as np
import Bentheimer_Sandstone

####################################################
#   各種容量系への変換
#   property番号と物性値の種類

### 熱容量
def crow_cal(tem,miu,pro):
    crow = np.zeros_like(pro)
    crow = 1000000.0
    crow = np.where(pro==1, 100000.0, crow)
    crow = np.where(pro==10, 200000.0, crow)
    crow = np.where(pro==30, Bentheimer_Sandstone.crow_cal(miu), crow)
    return crow

### 水分容量（水分化学ポテンシャルベース）
def cmiu_cal(tem,miu,pro):
    cmiu = np.zeros_like(pro)
    cmiu = 1000000.0
    cmiu = np.where(pro==1, 100000.0, cmiu)
    cmiu = np.where(pro==10, 200000.0, cmiu)
    cmiu = np.where(pro==30, Bentheimer_Sandstone.cmiu_cal(miu), cmiu)
    return cmiu

### 含水率ー水分化学ポテンシャル関係
def phi_cal(tem,miu,pro):
    phi = np.zeros_like(pro)
    phi = np.where(pro==1, 100000.0, phi)
    phi = np.where(pro==10, 200000.0, phi)
    phi = np.where(pro==30, Bentheimer_Sandstone.Phi_cal(miu), phi)
    return phi

####################################################
#   各種移動係数への変換
#   property番号と物性値の種類

### 熱伝導率
def LAM_cal(tem,miu,pro):
    LAM = np.zeros_like(pro)
#    LAM = np.where(pro==0, 0.0, LAM)
    LAM = np.where(pro==1, 1.0, LAM)
    LAM = np.where(pro==10, 2.0, LAM)
    LAM = np.where(pro==30, Bentheimer_Sandstone.LAM_cal(), LAM)
    return LAM

### 液相水分伝導率
###   （水分化学ポテンシャル勾配）
def LDML_cal(tem,miu,pro):
    LDML = np.zeros_like(pro)
    LDML = np.where(pro==1, 1.0e-7, LDML)
    LDML = np.where(pro==10, 2.0e-7, LDML)
    LDML = np.where(pro==30, Bentheimer_Sandstone.LDML_cal(miu), LDML)
    return LDML

### 気相水分伝導率
###   （水蒸気圧勾配）
def LDP_cal(tem,miu,pro):
    LDP = np.zeros_like(pro)
    LDP = np.where(pro==1, 1.0e-7, LDP)
    LDP = np.where(pro==10, 2.0e-7, LDP)
    LDP = np.where(pro==30, Bentheimer_Sandstone.LDP_cal(miu), LDP)
    return LDP


####################################################
# 調和平均計算 
def Harmonic_mean_cal(value,dx):
    Dvalue = np.diff(value*dx)
    Ddx = np.diff(dx)
    Average = (2*np.delete(value*dx,0)-Dvalue)/(2*np.delete(dx,0)-Ddx)
    return Average

####################################################
# プログラムの検証
value = np.array([10,20,25,30,32,40])
dx = np.array([10,20,15,10,5,20])
Average = Harmonic_mean_cal(value,dx)
print(Average)

tem = np.array([293.15]*3)
miu = np.array([-1000,-100,-10])
pro = np.array([10,20,30])
LAM = LAM_cal(tem,miu,pro)
LDML= LDML_cal(tem,miu,pro)
LDP = LDP_cal(tem,miu,pro)
#print(LAM,LDML,LDP)