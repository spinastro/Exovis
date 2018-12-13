import astropy.units as u
import pandas as pd
from scipy.stats import uniform
import numpy as np
from simplanet import *

def vis_calc(Ms, Mp, P, e, w,ob):
    RV,err,t=simplanet_MC(Ms*u.Msun, Mp*u.Mearth, P*u.day, e, w,ob)
    v = visibility_MC(P*u.day,RV,t)
    if v:
        return 1
    else:
        return 0

def MCsim(ob,Mp_min=10,Mp_max=30,P_min=3,P_max=100,Ms=1,e=0.01,Niter=100):
    uniform_dist_Mpl = uniform(loc = Mp_min, scale = Mp_max)
    uniform_dist_P = uniform(loc = P_min, scale = P_max)
    Mp=uniform_dist_Mpl.rvs(size = Niter, random_state = 123)
    P=uniform_dist_P.rvs(size = Niter, random_state = 321)
    Ms=np.zeros(Niter)+Ms
    e=np.zeros(Niter)+e
    w=np.zeros(Niter)
    data = pd.DataFrame({"Mp": Mp,"P": P,"e": e,"Ms": Ms,"w": w})
    data["Visibility"]=list(map(lambda x,y: vis_calc(1.,x,y,0.01,0.,ob), Mp,P))
    return data

def plot_MCsim(data):
    cmap = plt.get_cmap("brg")
    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(132)
    ax.set_title('Planet visibility')
    ax.scatter(data["Mp"],data["P"], c=data["Visibility"], cmap=cmap, s=100, vmin=-1, vmax=1)
    ax.set_xlabel('Mp [Mearth]')
    ax.set_ylabel('P [day]')


