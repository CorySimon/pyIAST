###
#   Isotherm models
#   Author: CoryMSimon@gmail.com
###
__author__ = 'Cory M. Simon'
__all__ = ["LangmuirIsotherm", "QuadraticIsotherm", "BETIsotherm", "SipsIsotherm", "InterpolatorIsotherm", "plot_isotherm"]

import scipy.optimize
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class LangmuirIsotherm:
    """
    Langmuir isotherm object to store pure-component adsorption isotherm.

    The Langmuir adsorption isotherm is:

    .. math::
    
        L(P) = M\\frac{KP}{1+KP},

    where :math:`L` is the gas uptake, :math:`P` is pressure (fugacity), :math:`M` is the saturation loading, and :math:`K` is the Langmuir constant.
    """

    def __init__(self, df, loading_key=None, pressure_key=None):
        """
        Instantiation. A LangmuirIsotherm object is instantiated by passing it the pure component adsorption isotherm in the form of a Pandas DataFrame. The least squares data fitting is done here.

        :param df: DataFrame adsorption isotherm data
        :param loading_key: String key for loading column in df
        :param pressure_key: String key for pressure column in df

        :return: self
        :rtype: LangmuirIsotherm
        """
        # store isotherm data in self
        #: Pandas DataFrame on which isotherm was fit
        self.df = df
        if loading_key==None or pressure_key == None:
            raise Exception("Provide names of pressure and loading cols in dataframe")
        #: name of loading column
        self.loading_key = loading_key
        #: name of pressure column
        self.pressure_key = pressure_key

        def RSS(params):
            """
            Residual Sum of Squares between Langmuir model and data in df
            :param params: Array params = [K, M]
            M: float saturation loading (units: loading)
            K: float Langmuir constant (units: 1/pressure)
            """
            return np.sum((df[loading_key].values -
                           params[1] * params[0] * df[pressure_key].values /
                           (1.0 + params[0] * df[pressure_key].values)) ** 2)

        # for guess as starting point in minimizing RSS
        M_guess = np.max(df[pressure_key].values)  # guess saturation loading to be highest loading
        # guess K using M_guess and lowest pressure point
        idx_min = np.argmin(df[loading_key].values)
        K_guess = df[loading_key].iloc[idx_min] / df[pressure_key].iloc[idx_min] / (
            M_guess - df[pressure_key].iloc[idx_min])

        # minimize RSS
        opt_res = scipy.optimize.minimize(RSS, [K_guess, M_guess], method='Nelder-Mead')
        if opt_res.success == False:
            print(opt_res.message)
            print "M_guess = ", M_guess
            print "K_guess = ", K_guess
            raise Exception("Minimization of RSS for Langmuir isotherm fitting failed... Try a different guess.")

        # assign params
        #: Langmuir constant K (units: 1 / pressure)
        self.K = opt_res.x[0]
        #: Saturation loading (units: loading)
        self.M = opt_res.x[1]

    def loading(self, P):
        """
        Given stored Langmuir parameters, compute loading at pressure P.

        :param P: Float or Array pressure (in corresponding units as df in instantiation)
        :return: loading at pressure P (in corresponding units as df in instantiation)
        :rtype: Float or Array
        """
        return self.M * self.K * P / (1.0 + self.K * P)

    def spreading_pressure(self, P):
        """
        Calculate reduced spreading pressure at a bulk gas pressure P. (see Tarafder eqn 4)

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: spreading pressure, :math:`\\Pi`
        :rtype: Float
        """
        return self.M * np.log(1.0 + self.K * P)

    def print_params(self):
        """
        Print identified model parameters
        """
        print "Langmuir K (1/pressure) = ", self.K
        print "Saturation loading, M (loading) = ", self.M



class QuadraticIsotherm:
    """
    Quadratic isotherm object to store pure-component isotherm.

    The Quadratic adsorption isotherm is:

    .. math::
    
        L(P) = M \\frac{(K_a P + 2 K_b P)P}{1+K_aP+K_bP^2}

    where :math:`L` is the gas uptake, :math:`P` is pressure (fugacity), :math:`M` is half of the saturation loading, and constants :math:`K_a` and :math:`K_b` are model coefficients.
    """

    def __init__(self, df, loading_key=None, pressure_key=None):
        """
        Instantiation. A QuadraticIsotherm object is instantiated by passing it the pure component adsorption isotherm in the form of a Pandas DataFrame. The least squares data fitting is done here.

        :param df: DataFrame adsorption isotherm data
        :param loading_key: String key for loading column in df
        :param pressure_key: String key for pressure column in df

        :return: self
        :rtype: QuadraticIsotherm 
        """
        # store isotherm data in self
        #: Pandas DataFrame on which isotherm was fit
        self.df = df
        if loading_key==None or pressure_key == None:
            raise Exception("Provide names of pressure and loading cols in dataframe")
        #: name of loading column
        self.loading_key = loading_key
        #: name of pressure column
        self.pressure_key = pressure_key

        def RSS(params):
            """
            Residual Sum of Squares between Quadratic isotherm model and data in df
            :param params: Array params = [Ka, Kb, M]
            M: float saturation loading (units: loading)
            Ka, Kb: float Quadratic isotherm constants (units: 1/pressure, 1/pressure^2)
            """
            return np.sum((df[loading_key].values -
                           params[2] * (params[0] + 2.0 * params[1] * df[pressure_key].values) * df[
                               pressure_key].values /
                           (1.0 + params[0] * df[pressure_key].values + params[1] * df[pressure_key].values ** 2)) ** 2)

        # for guess as starting point in minimizing RSS
        M_guess = np.max(df[pressure_key].values) / 2.0  # guess saturation loading to be highest loading
        # guess K using M_guess and lowest pressure point
        idx_min = np.argmin(df[loading_key].values)
        Ka_guess = df[loading_key].iloc[idx_min] / df[pressure_key].iloc[idx_min] / (
            M_guess - df[pressure_key].iloc[idx_min])
        Kb_guess = 0.000000001  # guess as Langmuir adsorption isotherm at first

        # minimize RSS
        opt_res = scipy.optimize.minimize(RSS, [Ka_guess, Kb_guess, M_guess], method='Nelder-Mead')
        if opt_res.success == False:
            print(opt_res.message)
            print "M_guess = ", M_guess
            print "Ka_guess = ", Ka_guess
            print "Kb_guess = ", Kb_guess
            raise Exception("Minimization of RSS for Quadratic isotherm fitting failed... Try a different guess.")

        # assign params
        #: isotherm coefficient (units: 1 / pressure)
        self.Ka = opt_res.x[0]
        #: isotherm coefficient (units: 1 / pressure\ :superscript:`2`)
        self.Kb = opt_res.x[1]
        #: Half of saturation loading (units: loading)
        self.M = opt_res.x[2]

    def loading(self, P):
        """
        Given stored Quadratic isotherm parameters, compute loading at pressure P.

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: loading at pressure P (in corresponding units as df in instantiation)
        :rtype: Float or Array
        """
        return self.M * (self.Ka + 2.0 * self.Kb * P) * P / (1.0 + self.Ka * P + self.Kb * P ** 2)

    def spreading_pressure(self, P):
        """
        Calculate reduced spreading pressure at a bulk gas pressure P. (see Tarafder eqn 4)

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: spreading pressure, :math:`\\Pi`
        :rtype: Float
        """
        return self.M * np.log(1.0 + self.Ka * P + self.Kb * P ** 2)
    
    def print_params(self):
        """
        Print identified model parameters
        """
        print "Ka (1/pressure) = ", self.Ka
        print "Kb (1/pressure) = ", self.Kb
        print "M (loading) = ", self.M


class BETIsotherm:
    """
    BET isotherm object to store pure-component adsorption isotherm.

    The BET adsorption isotherm, a multi-layer model, is:

    .. math::
    
        L(P) = M\\frac{K_A P}{(1-K_B P)(1-K_B P+ K_A P)},

    where :math:`L` is the gas uptake, :math:`P` is pressure (fugacity), :math:`M` is the number of adsorption sites on the bare surface, :math:`K_A` is the Langmuir constant for the bare surface, and :math:`K_B` is the Langmuir constant for the second and higher layers.
    """

    def __init__(self, df, loading_key=None, pressure_key=None):
        """
        Instantiation. A BETIsotherm object is instantiated by passing it the pure component adsorption isotherm in the form of a Pandas DataFrame. The least squares data fitting is done here.

        :param df: DataFrame adsorption isotherm data
        :param loading_key: String key for loading column in df
        :param pressure_key: String key for pressure column in df

        :return: self
        :rtype: BETIsotherm
        """
        # store isotherm data in self
        #: Pandas DataFrame on which isotherm was fit
        self.df = df
        if loading_key==None or pressure_key == None:
            raise Exception("Provide names of pressure and loading cols in dataframe")
        #: name of loading column
        self.loading_key = loading_key
        #: name of pressure column
        self.pressure_key = pressure_key

        def RSS(params):
            """
            Residual Sum of Squares between Langmuir model and data in df
            :param params: Array params = [K_A, K_B, M]
            K_A: float Langmuir constant for first layer (units: 1/pressure)
            K_B: float Langmuir constant for > first layer (units: 1/pressure)
            M: float saturation loading of first layer (units: loading)
            """
            return np.sum((df[loading_key].values -
                           params[2] * params[0] * df[pressure_key].values /
                           (1.0 - params[1] * df[pressure_key].values) /
                           (1.0 - params[1] * df[pressure_key].values + params[0] * df[pressure_key].values)) ** 2)

        # for guess as starting point in minimizing RSS
        M_guess = np.max(df[pressure_key].values)  # guess saturation loading to be highest loading
        # guess K_A using M_guess and lowest pressure point
        idx_min = np.argmin(df[loading_key].values)
        Ka_guess = df[loading_key].iloc[idx_min] / df[pressure_key].iloc[idx_min] / (
            M_guess - df[pressure_key].iloc[idx_min])
        Kb_guess = .001 * Ka_guess  # guess small

        # minimize RSS
        opt_res = scipy.optimize.minimize(RSS, [Ka_guess, Kb_guess, M_guess], method='Nelder-Mead')
        if opt_res.success == False:
            print(opt_res.message)
            print "M_guess = ", M_guess
            print "Ka_guess = ", Ka_guess
            print "Kb_guess = ", Kb_guess
            raise Exception("Minimization of RSS for Langmuir isotherm fitting failed... Try a different guess.")

        # assign params
        #: Langmuir constant of first layer (units: 1 / pressure)
        self.Ka = opt_res.x[0]
        #: Langmuir constant of second and higher layers (units: 1 / pressure)
        self.Kb = opt_res.x[1]
        #: Saturation loading of first layer (units: loading)
        self.M = opt_res.x[2]

    def loading(self, P):
        """
        Given stored Langmuir parameters, compute loading at pressure P.

        :param P: Float or Array pressure (in corresponding units as df in instantiation)
        :return: loading at pressure P (in corresponding units as df in instantiation)
        :rtype: Float or Array
        """
        return self.M * self.Ka * P / ( (1.0 - self.Kb * P) * (1.0 - self.Kb * P + self.Ka * P) )

    def spreading_pressure(self, P):
        """
        Calculate reduced spreading pressure at a bulk gas pressure P. (see Tarafder eqn 4)

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: spreading pressure, :math:`\\Pi`
        :rtype: Float
        """
        return self.M * np.log((1.0 - self.Kb * P + self.Ka * P) / (1.0 - self.Kb * P))
    
    def print_params(self):
        """
        Print identified model parameters
        """
        print "Ka (1/pressure) = ", self.Ka
        print "Kb (1/pressure) = ", self.Kb
        print "M (loading) = ", self.M


class InterpolatorIsotherm:
    """
    Interpolator isotherm object to store pure-component adsorption isotherm.

    Here, the isotherm is characterized by linear interpolation of data.

    Loading = 0.0 at pressure = 0.0 is enforced here automatically for interpolation at low pressures.

    Default for extrapolating isotherm beyond highest pressure is to fail. Pass a value for `fill_value` in instantiation to extrapolate loading as `fill_value`.
    """

    def __init__(self, df, loading_key=None, pressure_key=None, fill_value=None):
        """
        Instantiation. InterpolatorIsotherm is instantiated by passing it the pure component adsorption isotherm in the form of a Pandas DataFrame. Contructs linear interpolator from `interp1d` function in Scipy during instantiation.

        e.g. to extrapolate loading beyond highest pressure point as 100.0, pass `fill_value=100.0`.

        :param df: DataFrame adsorption isotherm data
        :param loading_key: String key for loading column in df
        :param pressure_key: String key for pressure column in df
        :param fill_value: Float value to assume when interp1d goes outside bounds

        :return: self
        :rtype: InterpolatorIsotherm 
        """
        # if pressure = 0 not in data frame, add it for interpolation between 0 and first point.
        if 0.0 not in df[pressure_key].values:
            df = pd.concat([pd.DataFrame({pressure_key:0.0, loading_key:0.0}, index=[0]), df])

        # store isotherm data in self
        #: Pandas DataFrame on which isotherm was fit
        self.df = df.sort([pressure_key], ascending=True)
        if loading_key==None or pressure_key == None:
            raise Exception("Provide names of pressure and loading cols in dataframe")
        #: name of loading column
        self.loading_key = loading_key
        #: name of pressure column
        self.pressure_key = pressure_key
        
        if fill_value == None:
            self.interp1d = interp1d(self.df[pressure_key], self.df[loading_key])
        else:
            self.interp1d = interp1d(self.df[pressure_key], self.df[loading_key], fill_value=fill_value, bounds_error=False)

    def loading(self, P):
        """
        Interpolate isotherm to compute loading at pressure P.

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: loading at pressure P (in corresponding units as df in instantiation)
        :rtype: Float or Array
        """
        return self.interp1d(P)

    def spreading_pressure(self, P):
        """
        Calculate reduced spreading pressure at a bulk gas pressure P. (see Tarafder eqn 4)

        Use trapezoid rule on isotherm data points to compute the reduced spreading pressure via the integral:

        .. math::

            \\Pi(p) = \\int_0^p \\frac{q(\\hat{p})}{ \\hat{p}} d\\hat{p}

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: spreading pressure, :math:`\\Pi`
        :rtype: Float
        """
        # Get all data points that are at nonzero pressures
        pressures = self.df[self.pressure_key].values[self.df[self.pressure_key].values != 0.0]
        loadings = self.df[self.loading_key].values[self.df[self.pressure_key].values != 0.0]
        
        # approximate loading up to first pressure point with Henry's law
        # loading = KH * P
        # KH is the initial slope in the adsorption isotherm
        KH = loadings[0] / pressures[0]
       
        # get all data points that are less than P 
        idx = pressures < P
        if np.sum(idx) == 0:
            # if this pressure is between 0 and first pressure point...
            p_array = np.concatenate((np.array([0.0]), np.array([P])))  # array of pressures
            integrand = np.concatenate((np.array([KH]), np.array([self.loading(P) / P])))  # array of loadings
            return  np.trapz(integrand, x=p_array)
        else:
            p_array = np.concatenate((np.array([0.0]), pressures[idx], np.array([P])))  # array of pressures
            integrand = np.concatenate((np.array([KH]), loadings[idx] / pressures[idx], np.array([self.loading(P) / P])))  # array of loadings
            return  np.trapz(integrand, x=p_array)


class SipsIsotherm:
    """
    Sips isotherm object to store pure-component adsorption isotherm.

    The Sips adsorption isotherm is:

    .. math::
    
        L(P) = M\\frac{K^nP^n}{1+K^nP^n},

    where :math:`L` is the gas uptake, :math:`P` is pressure (fugacity), :math:`M` is the saturation loading, :math:`K` is the Langmuir constant, and :math:`n` is the index of heterogeneity.
    """

    def __init__(self, df, loading_key=None, pressure_key=None):
        """
        Instantiation. A SipsIsotherm object is instantiated by passing it the pure component adsorption isotherm in the form of a Pandas DataFrame. The least squares data fitting is done here.

        :param df: DataFrame adsorption isotherm data
        :param loading_key: String key for loading column in df
        :param pressure_key: String key for pressure column in df

        :return: self
        :rtype: SipsIsotherm
        """
        # store isotherm data in self
        #: Pandas DataFrame on which isotherm was fit
        self.df = df
        if loading_key==None or pressure_key == None:
            raise Exception("Provide names of pressure and loading cols in dataframe")
        #: name of loading column
        self.loading_key = loading_key
        #: name of pressure column
        self.pressure_key = pressure_key

        def RSS(params):
            """
            Residual Sum of Squares between Sips model and data in df
            :param params: Array params = [K, M, n]
            M: float saturation loading (units: loading)
            K: float Langmuir constant (units: 1/pressure)
            n: float index of heterogeneity
            """
            KPn = (params[0] * df[pressure_key].values) ** params[2]
            return np.sum((df[loading_key].values -
                           params[1] * KPn /
                           (1.0 + KPn)) ** 2)

        # for guess as starting point in minimizing RSS
        M_guess = np.max(df[pressure_key].values)  # guess saturation loading to be highest loading
        # guess K using M_guess and lowest pressure point
        idx_min = np.argmin(df[loading_key].values)
        K_guess = df[loading_key].iloc[idx_min] / df[pressure_key].iloc[idx_min] / (
            M_guess - df[pressure_key].iloc[idx_min])
        n_guess = 1.0

        # minimize RSS
        opt_res = scipy.optimize.minimize(RSS, [K_guess, M_guess, n_guess], method='Nelder-Mead')
        if opt_res.success == False:
            print(opt_res.message)
            print "M_guess = ", M_guess
            print "K_guess = ", K_guess
            print "n_guess = ", n_guess
            raise Exception("Minimization of RSS for Langmuir isotherm fitting failed... Try a different guess.")

        # assign params
        #: Langmuir constant K (units: 1 / pressure)
        self.K = opt_res.x[0]
        #: Saturation loading (units: loading)
        self.M = opt_res.x[1]
        #: index of heterogeneity (unitless)
        self.n = opt_res.x[2]

    def loading(self, P):
        """
        Given stored Sips parameters, compute loading at pressure P.

        :param P: Float or Array pressure (in corresponding units as df in instantiation)
        :return: loading at pressure P (in corresponding units as df in instantiation)
        :rtype: Float or Array
        """
        return self.M * (self.K * P) ** self.n / (1.0 + (self.K * P) ** self.n)

    def spreading_pressure(self, P):
        """
        Calculate reduced spreading pressure at a bulk gas pressure P. (see Tarafder eqn 4)

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: spreading pressure, :math:`\\Pi`
        :rtype: Float
        """
        return self.M / self.n * np.log(1.0 + (self.K * P) ** self.n)

    def print_params(self):
        """
        Print identified model parameters
        """
        print "Langmuir K (1/pressure) = ", self.K
        print "Saturation loading, M (loading) = ", self.M
        print "Index of heterogeneity n = ", self.n


def plot_isotherm(isotherm, withfit=True, xlogscale=False, ylogscale=False):
    """
    Plot isotherm data and fit using Matplotlib.
    
    :param isotherm: LangmuirIsotherm,QuadraticIsotherm,InterpolatorIsotherm the adsorption isotherm object
    :param withfit: Bool plot fit as well
    :param xlogscale: Bool log-scale on x-axis
    :param ylogscale: Bool log-scale on y-axis
    """

    fig = plt.figure()
    if withfit:
        # array of pressures to plot model
        if xlogscale:
            idx = isotherm.df[isotherm.pressure_key].values != 0.0
            min_P = np.min(isotherm.df[isotherm.pressure_key].iloc[idx])
            P = np.logspace(np.log(min_P), np.log(isotherm.df[isotherm.pressure_key].max()), 200)
            plt.plot(P, isotherm.loading(P))
        else:
            P = np.linspace(isotherm.df[isotherm.pressure_key].min(), isotherm.df[isotherm.pressure_key].max(), 200)
            plt.plot(P, isotherm.loading(P))
    plt.scatter(isotherm.df[isotherm.pressure_key], isotherm.df[isotherm.loading_key])
    if xlogscale:
        plt.xscale("log")
    if ylogscale:
        plt.yscale("log")
    plt.xlim(xmin=0.0)
    plt.ylim(ymin=0.0)
    plt.xlabel('Pressure')
    plt.ylabel('Loading')
    plt.show()
