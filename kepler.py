import numpy as np
import astropy.units as u
from scipy.optimize import fsolve

G = 4*np.pi**2 /(356.25)**2 # G in AU^3/(day^2 Solar Mass)

def solve_kep_eq_1(l,e,P):
    """ Solve Keplers equation x - e*sin(x) = l for x"""
    """ Where x is the eccentric anomaly, e is the eccentricity and l is the mean anomaly"""
    try:
        l[0]
        res = np.zeros(l.shape)
        for i,li in enumerate(l):
            tmp,= fsolve(lambda x: x-e*np.sin(x) - (li.to(u.day)/u.day)*(2*np.pi*u.day/P.to(u.day)),li)
            res[i] = tmp
    except IndexError:
        res, = fsolve(lambda x: x - e*np.sin(x)- (l.to(u.day)/u.day)*(2*np.pi*u.day/P.to(u.day)),l)
    
    return res

def solve_kep_eq_2(E,e):
    """ Solve Keplers equation tan(x/2) = sqrt((1 + e)/(1 - e)) * tan(E/2) for x"""
    """ Where x is the true anomaly, E is the eccentric anomaly and e is the eccentricity"""
    try:
        E[0]
        res = np.zeros(E.shape)
        for i,valE in enumerate(E):
            tmp,= fsolve(lambda x: np.tan(x/2.) - np.sqrt((1 + e)/(1 - e)) * np.tan(valE/2.) ,valE)
            res[i] = tmp
    except IndexError:
        res, = fsolve(lambda x: np.tan(x/2.) - np.sqrt((1 + e)/(1 - e)) * np.tan(E/2.), E)
    
    return res

def semi_amplitude(Ms,Mp,P,e):
    K = 28.4*u.m/u.s * (Mp.to(u.Mjup)/u.Mjup) * (P.to(u.yr)/u.yr)**(-1./3.) * (Ms.to(u.Msun)/u.Msun)**(-2./3.) * (1/np.sqrt(1-e**2))
    return K

def create_Vs(Ms,Mp,P,e,w,t):
    K = semi_amplitude(Ms,Mp,P,e)
    E = solve_kep_eq_1(t,e,P)
    f = solve_kep_eq_2(E,e)
    Vr = K * (np.cos(f+w) + e*np.cos(w))
    return Vr

