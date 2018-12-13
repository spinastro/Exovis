import numpy as np
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.stats import LombScargle
from kepler import *

def simplanet(Planet,Observations):
    rnd = np.random.RandomState(seed=123)
    t = np.sort(rnd.uniform(0, Observations.time_range, Observations.n_data) * u.day)
    RV = create_Vs(Planet.Mstar,Planet.Mplanet,Planet.P,Planet.e,Planet.w,t)
    err = np.full_like(t.to(u.day)/u.day, Observations.Precision) * u.m/u.s
    RV = RV + rnd.normal(0, err.value)*err.unit
    Observations.RV=RV
    Observations.errRV=err
    Observations.t=t

def simplanet_MC(Ms,Mp,P,e,w,Observations):
    rnd = np.random.RandomState(seed=123)
    t = np.sort(rnd.uniform(0, Observations.time_range, Observations.n_data) * u.day)
    RV = create_Vs(Ms,Mp,P,e,w,t)
    err = np.full_like(t.to(u.day)/u.day, Observations.Precision) * u.m/u.s
    RV = RV + rnd.normal(0, err.value)*err.unit
    return RV, err, t

def visibility(Planet,Observations):
    LS=LombScargle(Observations.t.value, Observations.RV.value)
    frequency, power = LS.autopower()
    false_alarm = LS.false_alarm_level(0.01, method='bootstrap')
    prob = LS.power(1./Planet.P.value)
    if power.max() > false_alarm:
        return True
    else:
        return False

def visibility_MC(P,RV,t):
    LS=LombScargle(t.value, RV.value)
    frequency, power = LS.autopower()
    false_alarm = LS.false_alarm_level(0.01, method='bootstrap')
    prob = LS.power(1./P.value)
    if power.max() > false_alarm:
        return True
    else:
        return False

def plot_observations(Observations):
    fig, ax = plt.subplots(1, 1, figsize=(6,6))
    ax.errorbar(Observations.t.value, Observations.RV.value, yerr=Observations.errRV.value, fmt='o')
    ax.set_xlabel("Time [day]")
    ax.set_ylabel("RV [m/s]")

def plot_LombScargle(Planet,Observations):
    LS=LombScargle(Observations.t.value, Observations.RV.value, Observations.errRV.value)
    frequency, power = LS.autopower()
    plt.plot(frequency, power)
    plt.xlabel('frequency')
    plt.ylabel('Lomb-Scargle Power')
    plt.axvline(1./Planet.P.value, lw=2, color='red', alpha=0.4)
    plt.axhline(LS.false_alarm_level(0.1, method='bootstrap'), linestyle=':',  lw=2, color='black', alpha=0.4)
    plt.axhline(LS.false_alarm_level(0.01, method='bootstrap'), linestyle='-',  lw=2, color='black', alpha=0.4)
    plt.ylim([0.,1.])
