# ===========================================================
# About the code
# ===========================================================
# This code is part of the project 'SpicePy'.
# See README.md for more details
#
# Licensed under the MIT license (see LICENCE)
# Copyright (c) 2017 Luca Giaccone (luca.giaccone@polito.it)
# ===========================================================

# ==================
# imported modules
# ==================
from scipy.sparse.linalg import spsolve
import numpy as np
import spicepy.transient_sources as tsr


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

    # matrix and frequency
    net.conductance_matrix()
    net.dynamic_matrix()
    net.rhs_matrix()
    net.frequency_span()

    if np.isscalar(net.f) == 1:
        # single frequency solution
        net.x = spsolve(net.G + 1j * 2 * np.pi * net.f * net.C, net.rhs)
    else:
        net.x = np.zeros((net.G.shape[0], net.f.size), dtype=np.complex)
        for k, f in enumerate(net.f):
            net.x[:, k] = spsolve(net.G + 1j * 2 * np.pi * f * net.C, net.rhs)


def transient_solve(net):

    # ========================
    # compute OP related to IC
    # ========================

    # (deep) copy of the network
    from copy import deepcopy
    net_op = deepcopy(net)

    # dictionary to track changes (L <-> I and C <-> V)
    track_change = {}

    # reorder original network and get index
    net_op.reorder()
    indexL = sorted(net_op.isort[1])
    indexC = sorted(net_op.isort[2])
    indexV = sorted(net_op.isort[3])
    indexI = sorted(net_op.isort[4])

    nv = 1  # get max IDnumber for voltage sources
    for iv in indexV:
        Vnum = int(net.names[iv][1:])

        if isinstance(net.values[iv],list):
            tsr_fun = getattr(tsr, net.source_type[net.names[iv]])
            net_op.values[iv] = tsr_fun(*net.values[iv], t=0)

        if Vnum >= nv:
            nv = Vnum + 1
    ni = 1  # get max IDnumber for current sources
    for ii in indexI:
        Inum = int(net.names[ii][1:])

        if isinstance(net.values[ii],list):
            tsr_fun = getattr(tsr, net.source_type[net.names[ii]])
            net_op.values[ii] = tsr_fun(*net.values[ii], t=0)

        if Inum >= ni:
            ni = Inum + 1

    # transform inductors (to current sources)
    for k, il in enumerate(indexL):
        new_name = 'I' + str(ni + k)
        track_change[new_name] = net_op.names[il]
        net_op.values[il] = net_op.IC[net_op.names[il]]
        net_op.names[il] = new_name

    # transform capacitors (to voltage sources)
    for k, ic in enumerate(indexC):
        new_name = 'V' + str(nv + k)
        track_change[new_name] = net_op.names[ic]
        net_op.values[ic] = net_op.IC[net_op.names[ic]]
        net_op.names[ic] = new_name

    # reorder new network (to avoid confusion)
    net_op.reorder()
    # change type of analysis and solve (to get IC)
    net_op.analysis = ['.op']
    net_solve(net_op)

    # ==================
    # transient analysis
    # ==================

    # get time step and tend
    h = float(net.convert_unit(net.analysis[1]))
    tend = float(net.convert_unit(net.analysis[2]))
    # create time array
    net.t = np.arange(0, tend, h)  # if tend is not multiple of h --> net.t[-1] < tend

    # build required matrices
    net.conductance_matrix()
    net.dynamic_matrix()
    rhs_fun = net.rhs_matrix()

    # initialize solution space
    net.x = np.zeros((net.G.shape[0], net.t.size))

    # fill with initial conditions
    NV = len(net.isort[3])
    NE = len(net.isort[5])
    NH = len(net.isort[8])
    net.x[:, 0] = np.concatenate((net_op.x[:net_op.node_num],
                                  np.array(net_op.values)[sorted(net.isort[1])],
                                  net_op.x[net_op.node_num:(net_op.node_num + NV)],
                                  net_op.x[(net_op.node_num + NV):(net_op.node_num + NV + NE)],
                                  net_op.x[(net_op.node_num + NV + NE):(net_op.node_num + NV + NE + NH)]))

    # Solution (Integration using trepezoidal rule. Ref: Vlach, eq 9.4.6, pag. 277)
    K = net.C + 0.5 * h * net.G
    for k in range(1, net.t.size):
        rhs = (net.C - 0.5 * h * net.G) * net.x[:, k - 1] + 0.5 * h * (rhs_fun(net.t[k - 1]) + rhs_fun(net.t[k]))
        net.x[:, k] = spsolve(K, rhs)


def net_solve(net):
    """

    :param net:
    :return:
    """
    if net.analysis[0].lower() == '.op':
        dc_solve(net)
    elif net.analysis[0].lower() == '.ac':
        ac_solve(net)
    elif net.analysis[0].lower() == '.tran':
        transient_solve(net)