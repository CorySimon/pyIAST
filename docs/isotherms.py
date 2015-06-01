###
#   Isotherm models
#   Author: CoryMSimon@gmail.com
###
__author__ = 'Cory M. Simon'
__all__ = ["LangmuirIsotherm", "QuadraticIsotherm", "InterpolatorIsotherm"]

import scipy.optimize
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class LangmuirIsotherm:
    """
    Langmuir isotherm object to store pure-component adsorption isotherm

    :math:`M \frac{KP}{1+KP}`
    """

    def __init__(self, df, loading_key=None, pressure_key=None):
        """
        Instantiation

        :param df: DataFrame adsorption isotherm data
        :param loading_key: String key for loading column in df
        :param pressure_key: String key for pressure column in df

        :return: self
        :rtype: LangmuirIsotherm
        """
        # store isotherm data in self
        self.df = df
        if loading_key==None or pressure_key == None:
            raise Exception("Provide names of pressure and loading cols in dataframe")
        self.loading_key = loading_key
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
        self.K = opt_res.x[0]
        self.M = opt_res.x[1]

    def loading(self, P):
        """
        Given stored Langmuir parameters, compute loading at pressure P

        :param P: Float or Array pressure (in corresponding units as df in instantiation)
        :return: loading at pressure P (in corresponding units as df in instantiation)
        :rtype: Float or Array
        """
        return self.M * (self.K * P / (1.0 + self.K * P))

    def spreading_pressure(self, P):
        """
        Calculate spreading pressure at a bulk gas pressure P

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: spreading pressure
        :rtype: Float
        """
        return self.M * np.log(1.0 + self.K * P)



class QuadraticIsotherm:
    """
    Quadratic isotherm object
    
    :math:`M \frac{(K_a P + 2 K_b P)P}{1+K_aP+K_bP^2}`
    """

    def __init__(self, df, loading_key=None, pressure_key=None):
        """
        Instantiation

        :param df: DataFrame adsorption isotherm data
        :param loading_key: String key for loading column in df
        :param pressure_key: String key for pressure column in df

        :return: self
        :rtype: QuadraticIsotherm 
        """
        # store isotherm data in self
        self.df = df
        if loading_key==None or pressure_key == None:
            raise Exception("Provide names of pressure and loading cols in dataframe")
        self.loading_key = loading_key
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
        self.Ka = opt_res.x[0]
        self.Kb = opt_res.x[1]
        self.M = opt_res.x[2]

    def loading(self, P):
        """
        Given stored Quadratic isotherm parameters, compute loading at pressure P

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: loading at pressure P (in corresponding units as df in instantiation)
        :rtype: Float or Array
        """
        return self.M * (self.Ka + 2.0 * self.Kb * P) * P / (1.0 + self.Ka * P + self.Kb * P ** 2)

    def spreading_pressure(self, P):
        """
        Calculate spreading pressure at a bulk gas pressure P

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: pi: float spreading pressure
        """
        return self.M * np.log(1.0 + self.Ka * P + self.Kb * P ** 2)


class InterpolatorIsotherm:
    """
    Interpolator isotherm object
    """

    def __init__(self, df, loading_key=None, pressure_key=None, fill_value=None):
        """
        Instantiation

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
        self.df = df.sort([pressure_key], ascending=True)
        if loading_key==None or pressure_key == None:
            raise Exception("Provide names of pressure and loading cols in dataframe")
        self.loading_key = loading_key
        self.pressure_key = pressure_key
        
        if fill_value == None:
            self.interp1d = interp1d(df[pressure_key], df[loading_key])
        else:
            self.interp1d = interp1d(df[pressure_key], df[loading_key], fill_value=fill_value, bounds_error=False)

    def loading(self, P):
        """
        Interpolate isotherm to compute loading at pressure P

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: loading at pressure P (in corresponding units as df in instantiation)
        :rtype: Float or Array
        """
        return self.interp1d(P)

    def spreading_pressure(self, P):
        """
        Calculate spreading pressure at a bulk gas pressure P

        Use trapezoid rule on isotherm data points

        Spreading pressure :math`\Pi = \int_0^p q(\hat{p})/ \hat{p} d\hat{p}` 

        :param P: float pressure (in corresponding units as df in instantiation)
        :return: spreading pressure
        :rtype: Float
        """
        # get all data points less than this P. Do not use P = 0.0, this will give nan.
        idx = (self.df[self.pressure_key].values < P) & (self.df[self.pressure_key].values != 0.0)
        if np.sum(idx) == 0:
            # if this pressure is between 0 and first pressure point
            l_array = float(self.loading(P))
            return l_array / 2.0 # area of triangle starting at P = 0
        else:
            # (1/2 * p1 * [height=loading/p1]
            # area of first triangle in integral, from 0 to first pressure point
            A_first_triangle = self.df[self.loading_key].iloc[idx][0] / 2.0
            p_array = np.concatenate((self.df[self.pressure_key].iloc[idx], np.array([P])))  # array of pressures
            l_array = np.concatenate(
                (self.df[self.loading_key].iloc[idx], np.array([self.loading(P)])))  # array of loadings
            return np.trapz(l_array / p_array, x=p_array) + A_first_triangle


def plot_isotherm(isotherm, withfit=True, xlogscale=False, ylogscale=False):
    """
    Plot isotherm data and fit
    
    :param isotherm: LangmuirIsotherm,QuadraticIsotherm,InterpolatorIsotherm the adsorption isotherm object
    :param withfit: Bool plot fit as well
    :param xlogscale: Bool log-scale on x-axis
    :param ylogscale: Bool log-scale on y-axis
    """

    fig = plt.figure()
    if withfit:
        # array of pressures to plot model
        if xlogscale:
            P = np.logspace(np.log(isotherm.df[isotherm.pressure_key].min()), np.log(isotherm.df[isotherm.pressure_key].max()), 200)
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
