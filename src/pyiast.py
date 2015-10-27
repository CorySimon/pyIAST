"""
This module performs the heart of the IAST calculations, given the
pure-component adsorption isotherm models from the `isotherms` module.
"""
__author__ = 'Cory M. Simon'

from isotherms import _MODELS, _MODEL_PARAMS, ModelIsotherm,\
    InterpolatorIsotherm, plot_isotherm
import scipy.optimize
import numpy as np


def iast(partial_pressures, isotherms, verboseflag=False, warningoff=False):
    """
    Perform IAST calculation to predict multi-component adsorption isotherm from
    pure component adsorption isotherms.

    The material is now in equilibrium with a mixture of gases with partial
    pressures in the array `partial_pressures` in units corresponding to those
    passed in the list of isotherms.

    Pass a list of pure-component adsorption isotherms `isotherms`.

    :param partial_pressures: Array or list partial pressures of gas components,
        e.g. [5.0, 10.0] (bar)
    :param isotherms: list pure-component adsorption isotherms.
        e.g. [methane_isotherm, ethane_isotherm]
    :param verboseflag: Bool print off a lot of information
    :param warningoff: Bool when False, warnings will print when the IAST
        calculation result required extrapolation of the pure-component 
        adsorption isotherm beyond the highest pressure in the data

    :return: loadings: predicted uptakes of each component
    :rtype: Array
    """
    partial_pressures = np.array(partial_pressures)
    n_components = len(isotherms)  # number of components in the mixture
    if n_components == 1:
        raise Exception("Pass list of pure component isotherms...")

    if np.size(partial_pressures) != n_components:
        print """Example use:\n
              IAST([0.5,0.5], [xe_isotherm, kr_isotherm], verboseflag=true)"""
        raise Exception("Length of partial pressures != length of array of"
                        " isotherms...")

    if verboseflag:
        print "%d components." % n_components
        for i in range(n_components):
            print "\tPartial pressure component %d = %f" % (i,
                partial_pressures[i])

    # assert that the spreading pressures of each component are equal
    def spreading_pressure_differences(adsorbed_mole_fractions):
        r"""
        Assert that spreading pressures of each component at fictitious pressure
        are equal.

        :param adsorbed_mole_fractions: array mole fractions in the adsorbed
            phase; np.size(adsorbed_mole_fractions) = n_components - 1 because
            \sum z_i = 1 asserted here automatically.
        :returns: spreading_pressure_diff: array spreading pressure difference
            between component i and i+1
        """
        spreading_pressure_diff = np.zeros((n_components - 1,))
        for i in range(n_components - 1):
            if i == n_components - 2:
                # automatically assert \sum z_i = 1
                adsorbed_mole_fraction_n = 1.0 - np.sum(adsorbed_mole_fractions)
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(
                    partial_pressures[i] / adsorbed_mole_fractions[i]) -\
                    isotherms[i + 1].spreading_pressure(
                    partial_pressures[i + 1] / adsorbed_mole_fraction_n)
            else:
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(
                    partial_pressures[i] / adsorbed_mole_fractions[i]) -\
                    isotherms[i + 1].spreading_pressure(
                    partial_pressures[i + 1] / adsorbed_mole_fractions[i + 1])
        return spreading_pressure_diff

    # solve for mole fractions in adsorbed phase by equating spreading pressures
    guess = partial_pressures / np.sum(partial_pressures)
    n_tries = 0  # try with many different guesses until result found
    while n_tries < 100:
        adsorbed_mole_fractions = scipy.optimize.fsolve(
            spreading_pressure_differences, guess[:-1])

        # mole fraction in adsorbed phase
        adsorbed_mole_fractions = np.concatenate((adsorbed_mole_fractions,
            np.array([1 - np.sum(adsorbed_mole_fractions)])))
        n_tries += 1
        # check if feasible soln is found
        if ((np.sum(adsorbed_mole_fractions >= 0.0) == n_components) &
                (np.sum(adsorbed_mole_fractions <= 1.0) == n_components)):
            break
        guess = np.random.uniform(size=n_components)
        guess = guess / np.sum(guess)

        #     z = np.concatenate((res.z, np.array([1 - np.sum(res.z)])))
    if (np.sum(adsorbed_mole_fractions < 0.0) != 0) | (
            np.sum(adsorbed_mole_fractions > 1.0) != 0):
        print "Tried %d times" % n_tries
        for i in range(n_components):
            print "\tadsorbed mole fraction [%d] = %f" % (i,
                adsorbed_mole_fractions[i])
            print "\tGuess: ", guess[i]
        raise Exception("adsorbed mole fraction not in [0,1],"
                        " solution infeasible...")

    pressure0 = partial_pressures / adsorbed_mole_fractions

    # solve for the total gas adsorbed
    inverse_loading = 0.0
    for i in range(n_components):
        inverse_loading += adsorbed_mole_fractions[i] / isotherms[i].loading(
            pressure0[i])
    loading_total = 1.0 / inverse_loading

    # get loading of each component by multiplying by mole fractions
    loadings = adsorbed_mole_fractions * loading_total
    if verboseflag:
        # print IAST loadings and corresponding pure-component loadings
        for i in range(n_components):
            print "Component ", i
            print "\tp = ", partial_pressures[i]
            print "\tp^0 = ", pressure0[i]
            print "\tLoading: ", loadings[i]
            print "\tx = ", adsorbed_mole_fractions[i]
            print "\tSpreading pressure = ", isotherms[i].spreading_pressure(
                pressure0[i])
    # print warning if had to extrapolate isotherm in spreading pressure
    if not warningoff:
        for i in range(n_components):
            if pressure0[i] > isotherms[i].df[isotherms[i].pressure_key].max():
                print """WARNING:
                  Component %d: p^0 = %f > %f, the highest pressure
                  exhibited in the pure-component isotherm data. Thus,
                  pyIAST had to extrapolate the isotherm data to achieve
                  this IAST result.""" % (i, pressure0[i],
                    isotherms[i].df[isotherms[i].pressure_key].max())

    # return loadings [component 1, component 2, ...]. same units as in data
    return loadings


def reverse_iast(adsorbed_mole_fractions, total_pressure, isotherms,
                 verboseflag=False, warningoff=False):
    """
    Perform reverse IAST to predict gas phase composition at total pressure
    `total_pressure` that will yield adsorbed mole fractions 
    `adsorbed_mole_fractions`.

    Pass a list of pure-component adsorption isotherms `isotherms`.

    :param adsorbed_mole_fractions: Array desired adsorbed mole fractions,
        e.g. [.5, .5]
    :param total_pressure: Float total bulk gas pressure
    :param isotherms: list of pure-component adsorption isotherms.
        e.g. [ethane_isotherm, methane_isotherm]
    :param verboseflag: Bool print stuff
    :param warningoff: Bool when False, warnings will print when the IAST
        calculation result required extrapolation of the pure-component
        adsorption isotherm beyond the highest pressure in the data

    :return: gas_mole_fractions, loadings: bulk gas mole fractions that yield 
    desired adsorbed mole fractions `adsorbed_mole_fractions` at 
    `total_pressure`, adsorbed component loadings according to reverse IAST
    :rtype: Array, Array
    """
    n_components = len(isotherms)  # number of components in the mixture
    adsorbed_mole_fractions = np.array(adsorbed_mole_fractions)
    if n_components == 1:
        raise Exception("Pass list of pure component isotherms...")

    if np.size(adsorbed_mole_fractions) != n_components:
        print """Example use:\n
              reverse_IAST([0.5,0.5], 1.0, [xe_isotherm, kr_isotherm],
              verboseflag=true)"""
        raise Exception("Length of desired adsorbed mole fractions != length of"
                        " array of isotherms...")

    if np.sum(adsorbed_mole_fractions) != 1.0:
        raise Exception("Desired adsorbed mole fractions should sum to 1.0...")

    if verboseflag:
        print "%d components." % n_components
        for i in range(n_components):
            print "\tDesired adsorbed phase mole fraction of component %d = %f"\
                % (i, adsorbed_mole_fractions[i])

    # assert that the spreading pressures of each component are equal
    def spreading_pressure_differences(gas_mole_fractions):
        r"""
        Assert that spreading pressures of each component at fictitious pressure
        are equal.

        :param gas_mole_fractions: array mole fractions in bulk gas phase
            np.size(y) = n_components - 1 because \sum y_i = 1 asserted here
            automatically.
        :returns: spreading_pressure_diff: array spreading pressure difference
            between component i and i+1
        """
        spreading_pressure_diff = np.zeros((n_components - 1,))
        for i in range(n_components - 1):
            if i == n_components - 2:
                # automatically assert \sum y_i = 1
                gas_mole_fraction_n = 1.0 - np.sum(gas_mole_fractions)
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(
                    total_pressure * gas_mole_fractions[i] /\
                    adsorbed_mole_fractions[i]) -\
                    isotherms[i + 1].spreading_pressure(
                    total_pressure * gas_mole_fraction_n /\
                    adsorbed_mole_fractions[i + 1])
            else:
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(
                    total_pressure * gas_mole_fractions[i] /\
                    adsorbed_mole_fractions[i]) -\
                    isotherms[i + 1].spreading_pressure(
                    total_pressure * gas_mole_fractions[i + 1] /\
                    adsorbed_mole_fractions[i + 1])
        return spreading_pressure_diff

    # solve for mole fractions in adsorbed phase by equating spreading pressures
    # guess to be the same as the adsorbed phase mole fractions
    guess = adsorbed_mole_fractions
    n_tries = 0  # try with many different guesses until result found
    while n_tries < 10:
        gas_mole_fractions = scipy.optimize.fsolve(
                spreading_pressure_differences, guess[:-1])

        # mole fraction in bulk gas phase
        gas_mole_fractions = np.concatenate((gas_mole_fractions,
                                    np.array([1 - np.sum(gas_mole_fractions)])))
        n_tries += 1
        # check if feasible soln is found
        if (np.sum(gas_mole_fractions >= 0.0) == n_components) &\
           (np.sum(gas_mole_fractions <= 1.0) == n_components):
            break
        guess = np.random.uniform(size=(n_components,))
        guess = guess / np.sum(guess)

    if (np.sum(gas_mole_fractions < 0.0) != 0) |\
       (np.sum(gas_mole_fractions > 1.0) != 0):
        print "Tried %d times" % n_tries
        for i in range(n_components):
            print "\ty[%d] = %f\n" % (i, gas_mole_fractions[i])
            print "\tGuess: ", guess
        raise Exception("y not in [0,1], solution infeasible...")

    pressure0 = total_pressure * gas_mole_fractions / adsorbed_mole_fractions

    # solve for the total gas adsorbed
    inverse_loading = 0.0
    for i in range(n_components):
        inverse_loading += adsorbed_mole_fractions[i] / isotherms[i].loading(
                            pressure0[i])
    loading_total = 1.0 / inverse_loading

    # get loading of each component by multiplying by mole fractions
    loadings = adsorbed_mole_fractions * loading_total

    if verboseflag:
        # print off IAST loadings and corresponding pure component loadings
        for i in range(n_components):
            print "Component ", i
            print "\tDesired mole fraction in adsorbed phase, x = ",\
                adsorbed_mole_fractions[i]
            print "\tBulk gas mole fraction that gives this, y = ",\
                gas_mole_fractions[i]
            print "\tSpreading pressure = ",\
                isotherms[i].spreading_pressure(pressure0[i])
            print "\tp^0 = ", pressure0[i]
            print "\tLoading: ", loadings[i]

    # print warning if had to extrapolate isotherm in spreading pressure
    if not warningoff:
        for i in range(n_components):
            if pressure0[i] > isotherms[i].df[isotherms[i].pressure_key].max():
                print """WARNING:
                  Component %d: p0 = %f > %f, the highest pressure
                  exhibited in the pure-component isotherm data. Thus,
                  pyIAST had to extrapolate the isotherm data to achieve
                  this IAST result.""" % (i, pressure0[i],
                            isotherms[i].df[isotherms[i].pressure_key].max())

    # return mole fractions in gas phase, component loadings
    return gas_mole_fractions, loadings


def print_selectivity(component_loadings, partial_pressures):
    """
    Calculate selectivity as a function of component loadings and bulk gas
    pressures

    :param component_loadings: numpy array of component loadings
    :param partial_pressures: partial pressures of components
    """
    n = np.size(component_loadings)
    for i in range(n):
        for j in range(i + 1, n):
            print "Selectivity for component %d over %d = %f" % (i, j,
                    component_loadings[i] / component_loadings[j] /\
                    (partial_pressures[i] / partial_pressures[j]))
