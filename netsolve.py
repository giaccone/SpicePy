from scipy.sparse.linalg import spsolve
import numpy as np

def dc_solve(net):
    """
    "dc_solve" solves DC network

    :return:
        * x, solution
    """

    net.conductance_matrix()
    net.rhs_matrix()

    # linear system definition
    net.x = spsolve(net.G, net.rhs)

def ac_solve(net):
    """

    :param net:
    :return:
    """

    net.conductance_matrix()
    net.dynamic_matrix()
    net.rhs_matrix()

    # frequency
    f = float(net.analysis[-1])

    # linear system definition
    net.x = spsolve(net.G + 1j * 2 * np.pi * f* net.C, net.rhs)

def net_solve(net):
    """
    
    :param net:
    :return:
    """
    if net.analysis[0] == '.op':
        dc_solve(net)
    elif net.analysis[0] == '.ac':
        ac_solve(net)