import numpy as np
import matplotlib.pyplot as plt
import vapour_pressure
import Flux
import prop

x = []
y = []
rows    = 2872.0
vant    = 2.0
iphi    = 3.0
mw      = 18.02
ms      = 58.44
de      = 1.6e-9

####################
#     解析モデル    #
####################
# 解析モデルの大きさ
L = 20
dx  = np.array([0.1]*(L+1))

# 材料プロパティ
pro = np.array([30]*(L+1))
pro[0] = 0
pro[L] = 0

# 境界条件プロパティ
bcon= np.array([1]*(L))
bcon[0] = 0
bcon[L-1]=0 
nx = np.array([0.0]*L)

# 差分計算
ddx = np.diff(dx,n=1)
dx2 = np.delete(dx,0)-(ddx/2.0)
dx2[1]  = dx[1]+dx[2]/2.0
dx2[L-1]= dx[L-1]+dx[L-2]/2.0

####################
#     計算条件     #
####################
# 初期条件

tem0 = 283.15
tem = np.array([tem0]*(L+1))

miu0 = -1000.0
miu = np.array([miu0]*(L+1))

pv0  = vapour_pressure.Pv_cal(tem0,miu0)
pv  = np.array([pv0]*(L+1))

phi = np.array([0.22]*(L+1))

dt = 1000.0

# 境界条件の入力　
tem[0] = 283.15
tem[L] = 285.15

miu[0] = -1000.0
miu[L] = -900.0

pv[0]  = vapour_pressure.Pv_cal(tem[0],miu[0])
pv[L]  = vapour_pressure.Pv_cal(tem[L],miu[L])

alpha = (4.9+4.4)       #室内外で分ける時はリスト化
alphad= 4.9/(1005.0*1.205*(8314.41/18.02)*293.15)

######################
#   ループ計算開始   ###
######################

num = 0
while(num<3000):

############################################
# 境界の温度・湿度条件
    
############################################
# 物性値への換算
    # メッシュ毎の物性値
    lam = np.array(prop.LAM_cal(tem,miu,pro))
    ldp = np.array(prop.LDP_cal(tem,miu,pro))
    ldml= np.array(prop.LDML_cal(tem,miu,pro))

    # メッシュ間の物性値     差分型に取り直す
    dlam = prop.Harmonic_mean_cal(lam,dx)
    dldp = prop.Harmonic_mean_cal(ldp,dx)
    dldml= prop.Harmonic_mean_cal(ldml,dx)

############################################
# 流量計算
    #　差分計算
    dtem = np.diff(tem,n=1)
    dmiu = np.diff(miu,n=1)
    dpv  = np.diff(pv, n=1)

    #　熱流・水分流計算
    qf = Flux.Heat_Flux(alpha,dlam,dtem,dx2,bcon)
    jg = Flux.Flux_vapour_pressure(alphad,dldp,dpv,dx2,bcon)
    qr = jg* Flux.Ratent_heat_cal(dtem)
    jl = 0.0
    jl = Flux.Flux_liquid_potential(dldml,dmiu,dtem,dx2,nx,bcon)
    jw = jg + jl
    qs = qf + qr

############################################
# 収支計算
    # 熱収支
    dqs = np.array([0.0])
    dqs = np.append(dqs,[np.diff(qs,n=1)])
    dqs = np.append(dqs,[0.0])
    crow = prop.crow_cal(tem,miu,pro)/dt*dx
    ntem = tem+dqs/crow
    
    # 水分収支
    djw = np.array([0.0])
    djw = np.append(djw,[np.diff(jw,n=1)])
    djw = np.append(djw,[0.0])
    cmiu = prop.cmiu_cal(tem,miu,pro)/dt*dx
    nmiu = miu+djw/cmiu
    
############################################
# 値の変換
    tem  = ntem
    miu  = nmiu
    phi  = prop.phi_cal(tem,miu,pro)
    pv   = vapour_pressure.Pv_cal(tem,miu)
    rh   = pv/(vapour_pressure.Pvs_cal(tem))
    x.append(num)
    y.append(tem)
    
# 計算時間に関する計算    
    num = num+1

############################################

#plt.plot(y[1999])
#plt.yscale("log")
plt.grid(which="both")
plt.plot(x,y)
plt.show()
