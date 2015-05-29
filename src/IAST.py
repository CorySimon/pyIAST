__author__ = 'Cory M. Simon'

from isotherms import *
import scipy.optimize


def IAST(p, isotherms, verboseflag=False):
    """
    Perform IAST to predict multi-component adsorption isotherm from pure component adsorption isotherms.
    
    The material is now in equilibrium with a mixture of gases with partial pressures in the array `p`. 
    
    :param p: Array partial pressures of gas components, e.g. [.5, .5]
    :param isotherms: list of pure-component adsorption isotherms. e.g. [xe_isotherm, kr_isotherm]
    :param verboseflag: Bool print stuff

    :return: predicted uptakes of each component
    :rtype: Array
    """
    n_components = len(isotherms)  # number of components in the mixture
    if n_components == 1:
        raise Exception("Pass list of pure component isotherms...")

    if (np.size(p) != n_components):
        print """Example use:\n
              IAST([0.5,0.5], [xe_isotherm, kr_isotherm], verboseflag=true)"""
        raise Exception("Length of partial pressures != length of array of isotherms...")

    if verboseflag:
        print "%d components." % n_components
        for i in range(n_components):
            print "\tPartial pressure component %d = %f" % (i, p[i])

    # assert that the spreading pressures of each component are equal
    def spreading_pressure_differences(x):
        """
        Assert that spreading pressures of each component are equal.

        :param x: array mole fractions in the adsorbed phase
        np.size(x) = n_components - 1 because \sum x_i = 1 asserted here automatically.
        :returns: spreading_pressure_diff: array spreading pressure difference between component i and i+1
        """
        spreading_pressure_diff = np.zeros((n_components - 1,))
        for i in range(n_components - 1):
            if i == n_components - 2:
                xn = 1.0 - np.sum(x)  # automatically assert \sum x_i = 1
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(p[i] / x[i]) - isotherms[
                    i + 1].spreading_pressure(p[i + 1] / xn)
            else:
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(p[i] / x[i]) - isotherms[
                    i + 1].spreading_pressure(p[i + 1] / x[i + 1])
        return spreading_pressure_diff
    
    # solve for mole fractions in adsorbed phase by equating spreading pressures
    guess = p / np.sum(p)
    n_tries = 0  # try with many different guesses until result found
    while n_tries < 10:
     #     res = scipy.optimize.root(spreading_pressure_differences, guess[:-1])
        x = scipy.optimize.fsolve(spreading_pressure_differences, guess[:-1])

        # mole fraction in adsorbed phase
        x = np.concatenate((x, np.array([1 - np.sum(x)])))
        n_tries += 1
        # check if feasible soln is found
        if ((np.sum(x >= 0.0) == n_components) & (np.sum(x <= 1.0) == n_components)):
            break
        guess = np.random.uniform(size=n_components)
        guess = guess / np.sum(guess)

 #     x = np.concatenate((res.x, np.array([1 - np.sum(res.x)])))
    if ((np.sum(x < 0.0) != 0) | (np.sum(x > 1.0) != 0)):
        print "Tried %d times" % n_tries
        for i in range(n_components):
            print "\tx[%d] = %f\n" % (i, x[i])
            print "\tGuess: ", guess
        raise Exception("x not in [0,1], solution infeasible...")

    p0 = p / x

    # solve for the total gas adsorbed
    denom = 0.0
    for i in range(n_components):
        denom += x[i] / isotherms[i].loading(p0[i])
    q_total = 1.0 / denom

    # get loading of each component by multiplying by mole fractions
    q = x * q_total
    if verboseflag:
        # print off loadings according to IAST and corresponding pure component loadings
        for i in range(n_components):
            print "Component ", i
            print "\tp = ", p[i]
            print "\tLoading: ", q[i]
            print "\tx = ", x[i]
            print "\tSpreading pressure = ", isotherms[i].spreading_pressure(p0[i])
            print "\tPure component loading at same partial pressure = ", isotherms[i].loading(p[i])

    return q  # loadings [component 1, component 2, ...]. same units as in data
