
# coding: utf-8

# In[39]:


#import matplotlib.pyplot as plt
#import numpy as np
import van_Genuchten as vG


# ### Bentheimer Sandstoneの物性値  
# 参考：V.Voronina, L. Pel and K. Kopinga: The influence of osmotic pressure on poulticing treatments for heritage objects, Material and Structures, vol.46, pp221-231, 2013  
# $\phi_{max}$：空隙率[-]  
# $K_{sat}$：飽和透水係数  
# $\lambda^{'}_P$：水蒸気圧勾配に対する気相水分伝導率  
# $\rho$；材料の密度[kg/m3]  
# $C$：比熱[J/(kg・K)]  
# $\rho$：水の密度[kg/m3]  
# $r$：水の相変化熱量  
# 
# 物性値を書く際における注意事項  
# ・物性情報は必ず(tem:温度、moisture：水分状態)からなる関数する。（変数を必要としない場合でも書くこと）  
# ・
# 
# 【理想】  
# ・クラス内に"何の物性が含まれているか"、"それぞれの物性がどの水分状態の関数であるか"、"物性値を与える関数"が記載されている。  
# ・それぞれの名称を重複することなく一度定義で完結させれること。

# In[59]:


class Property():
    
##########################################
###     材料情報の入力        #############
    Phimax = 0.23
    Ksat   = 2.0e-7
    LAMDP  = 2.0E-10
    row  = 1479.25
    C    = 750.0
    
    roww = 1000.0
    r    = 4.18605E+3

###     van-Genuchten用情報    ##########
    Alfa = 10.0/98.0
    n = 2.0
    m = 1.0 -(1.0/n)
    l = 0.5
    
###     水分を表す指標（水分化学ポテンシャル） ###
    def __init__(self):
        self.proplist = {
            'crow' : 'Miu', 
            'LAM' : 'Miu', 
            'Phi' : 'Miu', 
            'Miu' : 'Phi', 
            'Dw' : 'Miu', 
            'DP' : 'Miu', 
            'DPhi' : 'Miu'
        }
        
### 熱物性 ##############################
#   熱容量
    def crow(self,tem,miu):
        return self.row*self.C +self.row *self.Phi(tem,miu) *self.r    
#   熱伝導率
    def LAM(self,tem,miu):
        return 1.2

#######################################
#   水分物性
#   含水率 from 水分化学ポテンシャル
    def Phi(self,tem,miu):
        return vG.Phi(self.Phimax,self.Alfa, miu, self.m, self.n)
#   水分化学ポテンシャル from 含水率
    def Miu(self,tem,phi):
        return vG.Miu(self.Phimax, phi, self.Alfa, self.m, self.n)                                              
#   含水率勾配に関する液相水分伝導率
    def Dw(self,tem,miu):
        return self.Ksat*vG.Kl(self.Alfa ,miu ,self.m ,self.n ,self.l)
#   水蒸気圧勾配に関する気相水分伝導率
    def DP(self,tem,miu):
        Phi = self.Phi(tem,miu)
        Sl = Phi/self.Phimax
        return self.LAMDP*(1.0-Sl*0.9)
    
#   含水率のポテンシャル微分
    def DPhi(self,tem,miu):
        return vG.DPhi(self.Phimax ,self.Alfa , miu, self.m, self.n)


# In[57]:


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


# In[60]:


m = Property()
m.proplist

