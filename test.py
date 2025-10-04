from state_models.mf_52 import MF52
import numpy as np
import matplotlib.pyplot as plt
mf = MF52()
g = np.deg2rad(5)
Fx_list=[]
Fy_list=[]
kappa_set = np.linspace(-3,3,1000)
SA_set = np.linspace(-25,25,100)
for Fz in np.linspace(100,1000,5):
    for kappa in kappa_set:
        Fx = mf.Fx(Fz=Fz,Kappa=kappa,Gamma=g)
        Fx_list.append(Fx)
    for SA in SA_set:
        alpha = np.deg2rad(SA)
        Fy = mf.Fy(Fz=Fz,Alpha=alpha,Gamma=g)
        Fy_list.append(Fy)

    # plt.plot(SA_set,Fy_list)
    # plt.title('Fy vs SA')
    # plt.show()
    plt.figure(1)
    plt.plot(SA_set,Fy_list)
    plt.figure(2)
    plt.plot(kappa_set,Fx_list)
    Fx_list=[]
    Fy_list=[]
plt.figure(1)
plt.title('Fy vs SA')
plt.figure(2)
plt.title('FX vs K')
plt.show()

