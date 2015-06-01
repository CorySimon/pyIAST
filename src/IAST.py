__author__ = 'Cory M. Simon'

from isotherms import *
import scipy.optimize


def IAST(p, isotherms, verboseflag=False):
    """
    Perform IAST calculation to predict multi-component adsorption isotherm from pure component adsorption isotherms.
    
    The material is now in equilibrium with a mixture of gases with partial pressures in the array `p`.

    Pass a list of pure-component adsorption isotherms `isotherms`.
    
    :param p: Array or list partial pressures of gas components, e.g. [.5, .5]
    :param isotherms: list pure-component adsorption isotherms. e.g. [xe_isotherm, kr_isotherm]
    :param verboseflag: Bool print stuff

    :return: q: predicted uptakes of each component
    :rtype: Array
    """
    p = np.array(p)
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
    def spreading_pressure_differences(z):
        """
        Assert that spreading pressures of each component at fictitious pressure are equal.

        :param z: array mole fractions in the adsorbed phase
        np.size(z) = n_components - 1 because \sum z_i = 1 asserted here automatically.
        :returns: spreading_pressure_diff: array spreading pressure difference between component i and i+1
        """
        spreading_pressure_diff = np.zeros((n_components - 1,))
        for i in range(n_components - 1):
            if i == n_components - 2:
                zn = 1.0 - np.sum(z)  # automatically assert \sum z_i = 1
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(p[i] / z[i]) - isotherms[
                    i + 1].spreading_pressure(p[i + 1] / zn)
            else:
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(p[i] / z[i]) - isotherms[
                    i + 1].spreading_pressure(p[i + 1] / z[i + 1])
        return spreading_pressure_diff
    
    # solve for mole fractions in adsorbed phase by equating spreading pressures
    guess = p / np.sum(p)
    n_tries = 0  # try with many different guesses until result found
    while n_tries < 10:
     #     res = scipy.optimize.root(spreading_pressure_differences, guess[:-1])
        z = scipy.optimize.fsolve(spreading_pressure_differences, guess[:-1])

        # mole fraction in adsorbed phase
        z = np.concatenate((z, np.array([1 - np.sum(z)])))
        n_tries += 1
        # check if feasible soln is found
        if ((np.sum(z >= 0.0) == n_components) & (np.sum(z <= 1.0) == n_components)):
            break
        guess = np.random.uniform(size=n_components)
        guess = guess / np.sum(guess)

 #     z = np.concatenate((res.z, np.array([1 - np.sum(res.z)])))
    if ((np.sum(z < 0.0) != 0) | (np.sum(z > 1.0) != 0)):
        print "Tried %d times" % n_tries
        for i in range(n_components):
            print "\tz[%d] = %f\n" % (i, z[i])
            print "\tGuess: ", guess
        raise Ezception("z not in [0,1], solution infeasible...")

    p0 = p / z

    # solve for the total gas adsorbed
    denom = 0.0
    for i in range(n_components):
        denom += z[i] / isotherms[i].loading(p0[i])
    q_total = 1.0 / denom

    # get loading of each component by multiplying by mole fractions
    q = z * q_total
    if verboseflag:
        # print off loadings according to IAST and corresponding pure component loadings
        for i in range(n_components):
            print "Component ", i
            print "\tp = ", p[i]
            print "\tLoading: ", q[i]
            print "\tz = ", z[i]
            print "\tSpreading pressure = ", isotherms[i].spreading_pressure(p0[i])
            print "\tPure component loading at same partial pressure = ", isotherms[i].loading(p[i])

    return q  # loadings [component 1, component 2, ...]. same units as in data

def reverse_IAST(z, P_total, isotherms, verboseflag=False):
    """
    Perform reverse IAST to predict gas phase composition at total pressure `P_total` that will yield adsorbed mole fractions `z`.

    Pass a list of pure-component adsorption isotherms `isotherms`.
    
    :param z: Array desired adsorbed mole fractions, e.g. [.5, .5]
    :param P_total: Float total bulk gas pressure
    :param isotherms: list of pure-component adsorption isotherms. e.g. [xe_isotherm, kr_isotherm]
    :param verboseflag: Bool print stuff

    :return: y, q: bulk gas mole fractions that yield desired adsorbed mole fractions `z` at `P_total`, component loadings
    :rtype: Array, Array
    """
    n_components = len(isotherms)  # number of components in the mixture
    z = np.array(z)
    if n_components == 1:
        raise Exception("Pass list of pure component isotherms...")

    if (np.size(z) != n_components):
        print """Example use:\n
              reverse_IAST([0.5,0.5], 1.0, [xe_isotherm, kr_isotherm], verboseflag=true)"""
        raise Exception("Length of desired adsorbed mole fractions != length of array of isotherms...")

    if np.sum(z) != 1.0:
        raise Exception("Desired adsorbed mole fractions should sum to 1.0...")

    if verboseflag:
        print "%d components." % n_components
        for i in range(n_components):
            print "\tDesired adsorbed phase mole fraction of component %d = %f" % (i, z[i])

    # assert that the spreading pressures of each component are equal
    def spreading_pressure_differences(y):
        """
        Assert that spreading pressures of each component at fictitious pressure are equal.

        :param y: array mole fractions in bulk gas phase
        np.size(y) = n_components - 1 because \sum y_i = 1 asserted here automatically.
        :returns: spreading_pressure_diff: array spreading pressure difference between component i and i+1
        """
        spreading_pressure_diff = np.zeros((n_components - 1,))
        for i in range(n_components - 1):
            if i == n_components - 2:
                yn = 1.0 - np.sum(y)  # automatically assert \sum y_i = 1
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(P_total * y[i] / z[i]) - isotherms[
                    i + 1].spreading_pressure(P_total * yn / z[i + 1])
            else:
                spreading_pressure_diff[i] = isotherms[i].spreading_pressure(P_total * y[i] / z[i]) - isotherms[
                    i + 1].spreading_pressure(P_total * y[i + 1] / z[i + 1])
        return spreading_pressure_diff
    
    # solve for mole fractions in adsorbed phase by equating spreading pressures
    guess = z
    n_tries = 0  # try with many different guesses until result found
    while n_tries < 10:
     #     res = scipy.optimize.root(spreading_pressure_differences, guess[:-1])
        y = scipy.optimize.fsolve(spreading_pressure_differences, guess[:-1])

        # mole fraction in bulk gas phase
        y = np.concatenate((y, np.array([1 - np.sum(y)])))
        n_tries += 1
        # check if feasible soln is found
        if ((np.sum(y >= 0.0) == n_components) & (np.sum(y <= 1.0) == n_components)):
            break
        guess = np.random.uniform(size=(n_components, ))
        guess = guess / np.sum(guess)

    if ((np.sum(y < 0.0) != 0) | (np.sum(y > 1.0) != 0)):
        print "Tried %d times" % n_tries
        for i in range(n_components):
            print "\ty[%d] = %f\n" % (i, y[i])
            print "\tGuess: ", guess
        raise Exception("y not in [0,1], solution infeasible...")

    p0 = P_total * y / z
    
    # solve for the total gas adsorbed
    denom = 0.0
    for i in range(n_components):
        denom += z[i] / isotherms[i].loading(p0[i])
    q_total = 1.0 / denom

    # get loading of each component by multiplying by mole fractions
    q = z * q_total

    if verboseflag:
        # print off loadings according to IAST and corresponding pure component loadings
        for i in range(n_components):
            print "Component ", i
            print "\tDesired mole fraction in adsorbed phase, z = ", z[i]
            print "\tBulk gas mole fraction that gives this, y = ", y[i]
            print "\tSpreading pressure = ", isotherms[i].spreading_pressure(p0[i])
            print "\tLoading: ", q[i]
            print "\tPure component loading at same partial pressure = ", isotherms[i].loading(P_total * y[i])

    return y, q  # mole fractions in gas phase, component loadings
