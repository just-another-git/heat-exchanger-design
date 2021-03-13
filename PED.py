#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import CoolProp as cp
import math
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


print('DESIGN OF SHELL AND TUBE EXCHANGER')
print('Enter inputs')
#taking inputs
Hot_fluid = input("Enter your Hot fluid: ")
Cold_fluid=input("Enter your cold fluid: ")
m_h=4+0.03*float(input("Enter group no.: "))
T_hin=int(input("Enter Inlet Temperature of hot fluid(C): "))+273
T_hout=int(input("Enter outlet Temperature of hot fluid(C): "))+273
T_cin=int(input("Enter Inlet Temperature of cold fluid(C): "))+273


# In[3]:


T_cout=T_hout
h_cal_temp=(T_hin+T_hout)/2
c_cal_temp=(T_cin+T_cout)/2
r=h_cal_temp/c_cal_temp
CP_c=4200
CP_h=2202.25
del_T1=T_hin-T_cout
del_T2=T_hout-T_cin
LMTD=(del_T1-del_T2)/math.log(del_T1/del_T2)


# In[4]:


m_c=m_h*CP_h*(T_hin-T_hout)/(CP_c*(T_cout-T_cin))


# In[5]:


Q=m_h*CP_h*(T_hin-T_hout)


# In[6]:


data=pd.read_excel(r'/home/saikat/Downloads/Tube_Count_1-in_254-mm_Tubes_on_Square_90-Degree_Array_Pitch_125_in_.xlsx')
df = pd.DataFrame(data, columns= ['(Internal Shell Diameter) in','(Internal Shell Diameter) mm','(Number of Passes) 4'])


# In[20]:


U_cal=0
U_assume=300
U_assume1=U_assume
while (abs((U_cal-U_assume1)/U_assume1)>0.005):
    inc=0.0254
    U_assume1=U_assume
    Ft=0.965
    mew_55=4.5*10**-4
    mew_40=6.53*10**-4
    k_w=0.62856
    A=Q/(Ft*U_assume*LMTD)
    L=20*12*inc
    od_t=1*inc
    ind_t=3/4*inc
    na=A/(22/7*od_t*L)
    j=0
    while df['(Number of Passes) 4'][j]<na:
        j=j+1
    nt=df['(Number of Passes) 4'][j]
    ind_s=df['(Internal Shell Diameter) in'][j]*inc
    Re_t=4*m_c*(4/nt)/(3.14*ind_t*mew_55)
    #tube side heat transfer co-efficient
    pr=4.2*mew_40/k_w*10**3
    Nu=0.27*(pr**0.4)*((mew_40/mew_55)**1.4)*(Re_t**0.8)
    hi=Nu*k_w/ind_t
    #shell side heat transfer co-efficient
    k_n=0.14989
    den_n=807
    mew_n=1.4579*10**-3
    P_t=5/4*inc
    B=0.5*ind_s
    C=P_t-od_t
    a_s=C*B*ind_s/P_t
    V_s=m_h/(den_n*a_s)
    De=4*(P_t**2-(3.14/4)*od_t**2)/(3.14*od_t)
    Re_s=den_n*V_s*De/mew_n
    JH = 0.352978*(Re_s**0.549474) + 0.836215
    Cp_n=2175.04
    h0=JH*k_n**(2/3)/((mew_n*Cp_n)**(-1/3)*0.6**0.14*De)
    #calculate Ucal
    Kw=16.26
    R_ind_t=0.00075
    U_cal=(1/h0+(od_t/ind_t)**2*((od_t-ind_t)/(2*Kw)+R_ind_t+1/hi))
    U_cal=1/U_cal
    U_assume=U_cal 


# In[13]:


#pressure Calculation
f=Re_t**(-0.265293)*10**(-2.50515)*12*12
del_Pf=0.5*f*1000*(15.365*4/nt)**2*L/ind_t
del_Pr=4*4*(15.365*4/nt)**2/(2*9.8)
del_P=(del_Pf+del_Pr)/6894.76
del_P


# In[22]:


#overdesign
ov=abs((na-nt)/nt)*100
if ov<10:
    print("ok")
else:
    print("check")


# In[ ]:




