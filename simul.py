import astropy.units as u

class Planet:
    """Planet objects contain information about the planet.
        """
    def __init__(self, Mstar=None, Mplanet=None, P=None, e=None, w=None):
        self.Mstar = Mstar * u.Msun
        self.Mplanet = Mplanet * u.Mearth
        self.P = P * u.day
        self.e = e
        self.w = w
        print('Planet successfully created.')

class Observations:
    """Planet objects contain information about the planet.
        """
    def __init__(self, n_data=None, time_range=None, Precision=None, RV=None, errRV=None, t=None):
        self.n_data = n_data
        self.time_range = time_range
        self.Precision = Precision
        self.RV = RV
        self.errRV = errRV
        self.t = t
        print('Observational strategy successfully created.')
