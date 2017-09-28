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
from scipy.sparse import csr_matrix
import numpy as np

# ==================
# constants
# ==================
pi = np.pi

# ==================
# Class
# ==================
class Network:
    """
    Class that defines the network.

    Basilar attributes (always assigned):
        * self.names: component names
        * self.values: component values
        * self.IC: initial conditions for dynamic components (stored as dict)
        * self.nodes: component nodes (only two-port right now)
        * self.node_num: number of nodes in the network
        * self.analysis: type of analysis
        * self.plot_cmd: plot directive for transient analysis
    Other attributes (default is None):
        * self.A:
        * self.G:
        * self.C:
        * self.rhs:
        * self.isort:
        * self.t:
        * self.x:
        * self.vb:
        * self.ib:
        * self.pb:
    Methods:
        * read_netlist(self): reads a SPICE netlist (used by '__init__')
        * incidence_matrix(self):
    """

    def __repr__(self):
        return "SpicePy.Network: {} analysis".format(self.analysis[0])

    def __str__(self):
        msg = '------------------------\n'
        msg += '    SpicePy.Network:\n'
        msg += '------------------------\n'
        for ele, nodes, val in zip(self.names, self.nodes, self.values):
            if np.iscomplex(val):
                msg += "{} {} {} {} {}\n".format(ele, *nodes, np.abs(val), np.angle(val) * 180/np.pi)
            elif ele[0].upper() == 'C' or ele[0].upper() == 'L':
                if ele in self.IC:
                    msg += "{} {} {} {} ic={}\n".format(ele, *nodes, val, self.IC[ele])
                else:
                    msg += "{} {} {} {}\n".format(ele, *nodes, val)
            else:
                msg += "{} {} {} {}\n".format(ele, *nodes, val)

        msg += " ".join(self.analysis) + '\n'

        if self.plot_cmd is not None:
            msg += self.plot_cmd + '\n'

        msg += '------------------------\n'
        msg += '* number of nodes {}\n'.format(self.node_num)
        msg += '* number of branches {}\n'.format(len(self.names))
        msg += '------------------------\n'

        return msg


    def __init__(self, filename):
        """
        __init__ initializes common attributes using read_netlist method

        :param filename:
        """

        # common attributes
        (self.names,
         self.values,
         self.IC,
         self.nodes,
         self.node_num,
         self.analysis,
         self.plot_cmd) = self.read_netlist(filename)

        # initialization of other possible attributes
        self.A = None
        self.G = None
        self.C = None
        self.rhs = None
        self.isort = None
        self.t = None
        self.x = None
        self.vb = None
        self.ib = None
        self.pb = None

    def read_netlist(self, filename):
        """
        'readNetlist' reads a SPICE netlist

        :param filename: file name with the netlist
        :return:
            * names: element names
            * values: element values
            * nodes: element nodes
            * Nn: number of nodes in the netlist
            * analysis: type of analysis
        """

        # initialize counter for number of nodes
        Nn = 0

        # initialize variable names and values
        names = []
        values = []
        nodes = []
        IC = {}
        plot_cmd = None
        analysis = None

        # get the analysis type
        with open(filename) as f:
            # cycle on lines
            for b, line in enumerate(f):

                if line[0][0] == '.':  # analysis/command identifier
                    if (line.lower().find('.plot') != -1):
                        plot_cmd = line
                    else:
                        # split into a list
                        sline = line.split()
                        analysis = sline
        f.close()

        # read the netlist (according to analysis type)
        with open(filename) as f:
            # cycle on lines
            for b, line in enumerate(f):
                # split into a list
                sline = line.split()

                if sline[0][0] != '.':  # analysis identifier
                    # get branch nodes
                    N1 = int(sline[1])
                    N2 = int(sline[2])

                # detect element type
                if sline[0][0].upper() == 'R':  # resistance
                    # add name and value
                    names.append(sline[0])
                    values.append(float(sline[3]))
                    nodes.append([N1, N2])

                    # update counter
                    Nn = max([N1, N2, Nn])

                # inductor
                elif sline[0][0].upper() == 'L':
                    # add name and value
                    names.append(sline[0])
                    values.append(float(sline[3]))
                    nodes.append([N1, N2])

                    if analysis[0] == '.tran':
                        if len(sline) == 5:
                            if sline[4].lower().find('ic') != -1:
                                IC[sline[0]] = float(sline[4].split('=')[1])
                            else:
                                IC[sline[0]] = 'Please check this --> ' + sline[-1]
                                print("Warning: wrong definition of IC for {} --> ".format(sline[0]) + IC[sline[0]])
                        else:
                            IC[sline[0]] = 0

                    # update counter
                    Nn = max([N1, N2, Nn])

                # capacitor
                elif sline[0][0].upper() == 'C':
                    # add name and value
                    names.append(sline[0])
                    values.append(float(sline[3]))
                    nodes.append([N1, N2])

                    if analysis[0] == '.tran':
                        if len(sline) == 5:
                            if sline[4].lower().find('ic') != -1:
                                IC[sline[0]] = float(sline[4].split('=')[1])
                            else:
                                IC[sline[0]] = 'Please check this --> ' + sline[-1]
                                print("Warning: wrong definition of IC for {} --> ".format(sline[0]) + IC[sline[0]])
                        else:
                            IC[sline[0]] = 0

                    # update counter
                    Nn = max([N1, N2, Nn])

                # independent current source
                elif sline[0][0].upper() == 'I':
                    # add name and value
                    names.append(sline[0])
                    if (analysis[0] == '.ac') & (len(sline) == 5):
                        values.append(float(sline[3]) * (np.cos(float(sline[4]) * pi / 180) + np.sin(float(sline[4]) * pi / 180) * 1j))
                    else:
                        values.append(float(sline[3]))
                    nodes.append([N1, N2])

                    # update counter
                    Nn = max([N1, N2, Nn])

                # independent voltage sources
                elif sline[0][0].upper() == 'V':  # independent voltage sources
                    # add name and value
                    names.append(sline[0])
                    if (analysis[0] == '.ac') & (len(sline) == 5):
                        values.append(float(sline[3]) * (np.cos(float(sline[4]) * pi / 180) + np.sin(float(sline[4]) * pi / 180) * 1j))
                    else:
                        values.append(float(sline[3]))
                    nodes.append([N1, N2])

                    # update counter
                    Nn = max([N1, N2, Nn])

        # return network structure
        return names, values, IC, nodes, Nn, analysis, plot_cmd

    def incidence_matrix(self):
        """
        'incidence_matrix' creates the branch-2-node incidence matrix

        :return: update self with self.A
        """

        # initialize incidence matrix terms
        a = []
        a_row = []
        a_col = []

        # cycle on branches (N1 and N2)
        for b, nodes in enumerate(self.nodes):
            # get nodes
            N1, N2 = nodes

            # detect connection to ground
            if N1 == 0:
                a.append(-1)
                a_row.append(N2 - 1)
                a_col.append(b)
            elif N2 == 0:
                a.append(1)
                a_row.append(N1 - 1)
                a_col.append(b)
            else:
                a.append(1)
                a_row.append(N1 - 1)
                a_col.append(b)
                a.append(-1)
                a_row.append(N2 - 1)
                a_col.append(b)

        # create conductance matrix
        self.A = csr_matrix((a, (a_row, a_col)))

    def conductance_matrix(self):
        """
        'conductance_matrix' creates the conductance matrix

        :return: G, conductance matrix (including constant terms related to inductors and independent voltage sources)
        """
        # initialize conductance terms
        g = []
        g_row = []
        g_col = []

        # reorder if necessary
        if self.isort is None:
            self.reorder()

        # get index
        indexR = self.isort[0]
        indexL = sorted(self.isort[1])
        indexV = sorted(self.isort[3])

        # cycle on resistances
        for ir in indexR:
            # get nores
            N1, N2 = self.nodes[ir]

            # detect connection
            if (N1 == 0) or (N2 == 0): # if grounded...
                # diagonal term
                g.append(1.0 / self.values[ir])
                g_row.append(max([N1, N2]) - 1)
                g_col.append(max([N1, N2]) - 1)

            else:                      # if not grounded...
                # diagonal term
                g.append(1.0 / self.values[ir])
                g_row.append(N1 - 1)
                g_col.append(N1 - 1)

                # diagonal term
                g.append(1.0 / self.values[ir])
                g_row.append(N2 - 1)
                g_col.append(N2 - 1)

                # N1-N2 term
                g.append(-1.0 / self.values[ir])
                g_row.append(N1 - 1)
                g_col.append(N2 - 1)

                # N2-N1 term
                g.append(-1.0 / self.values[ir])
                g_row.append(N2 - 1)
                g_col.append(N1 - 1)

        # cycle on inductors
        for k, il in enumerate(indexL):
            # get nodes
            N1, N2 = self.nodes[il]

            # detect connection
            if N1 == 0:  # if grounded to N1 ...
                # negative terminal
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(self.node_num + k)

                # negative terminal
                g.append(-1)
                g_row.append(self.node_num + k)
                g_col.append(N2 - 1)

            elif N2 == 0:  # if grounded to N2 ...
                # positive terminal
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(self.node_num + k)

                # positive terminal
                g.append(1)
                g_row.append(self.node_num + k)
                g_col.append(N1 - 1)

            else:  # if not grounded ...
                # positive terminal
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(self.node_num + k)

                # positive terminal
                g.append(1)
                g_row.append(self.node_num + k)
                g_col.append(N1 - 1)

                # negative terminal
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(self.node_num + k)

                # negative terminal
                g.append(-1)
                g_row.append(self.node_num + k)
                g_col.append(N2 - 1)

        # cycle on independent voltage sources
        for k, iv in enumerate(indexV):
            # get nodes
            N1, N2 = self.nodes[iv]

            # detect connection
            if N1 == 0:               # if grounded to N1 ...
                # negative terminal
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(self.node_num + len(indexL) + k)

                # negative terminal
                g.append(-1)
                g_row.append(self.node_num + len(indexL) + k)
                g_col.append(N2 - 1)

            elif N2 == 0:              # if grounded to N2 ...
                # positive terminal
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(self.node_num + len(indexL) + k)

                # positive terminal
                g.append(1)
                g_row.append(self.node_num + len(indexL) + k)
                g_col.append(N1 - 1)

            else:                      # if not grounded ...
                # positive terminal
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(self.node_num + len(indexL) + k)

                # positive terminal
                g.append(1)
                g_row.append(self.node_num + len(indexL) + k)
                g_col.append(N1 - 1)

                # negative terminal
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(self.node_num + len(indexL) + k)

                # negative terminal
                g.append(-1)
                g_row.append(self.node_num + len(indexL) + k)
                g_col.append(N2 - 1)

        # create conductance matrix
        self.G = csr_matrix((g,(g_row,g_col)))

    def dynamic_matrix(self):
        """
        'dynamic_matrix' creates the dynamic matrix

        :return: C, dynamic matrix (inductors and capacitors)
        """

        # initialize conductance terms
        c = []
        c_row = []
        c_col = []

        # reorder if necessary
        if self.isort is None:
            self.reorder()

        # get index
        indexL = sorted(self.isort[1])
        indexC = self.isort[2]

        # cycle on inductors
        for k, il in enumerate(indexL):
            c.append(-self.values[il])
            c_row.append(self.node_num + k)
            c_col.append(self.node_num + k)

        # cycle on capacitors
        for ic in indexC:
            # get nores
            N1, N2 = self.nodes[ic]

            # detect connection
            if (N1 == 0) or (N2 == 0):  # if grounded...
                # diagonal term
                c.append(self.values[ic])
                c_row.append(max([N1, N2]) - 1)
                c_col.append(max([N1, N2]) - 1)

            else:  # if not grounded...
                # diagonal term
                c.append(self.values[ic])
                c_row.append(N1 - 1)
                c_col.append(N1 - 1)

                # diagonal term
                c.append(self.values[ic])
                c_row.append(N2 - 1)
                c_col.append(N2 - 1)

                # N1-N2 term
                c.append(-self.values[ic])
                c_row.append(N1 - 1)
                c_col.append(N2 - 1)

                # N2-N1 term
                c.append(-self.values[ic])
                c_row.append(N2 - 1)
                c_col.append(N1 - 1)

        # create dynamic matrix
        self.C = csr_matrix((c, (c_row, c_col)), shape=self.G.shape)

    def rhs_matrix(self):
        """
        'rhs_matrix' creates the right hand side matrix

        :return: rhs, right hand side
        """
        # reorder if necessary
        if self.isort is None:
            self.reorder()

        # initialize rhs
        rhs = [0] * (self.node_num + len(self.isort[1]) + len(self.isort[3]))

        # get index
        NL = len(self.isort[1])
        indexV = sorted(self.isort[3])
        indexI = self.isort[4]

        # cycle on independent voltage sources
        for k, iv in enumerate(indexV):
            # update rhs
            rhs[self.node_num + NL + k] += self.values[iv]

        # cycle on independent current sources
        for ii in indexI:
            # get nodes
            N1, N2 = self.nodes[ii]

            if N1 == 0:
                # update rhs
                rhs[N2 - 1] += self.values[ii]

            elif N2 == 0:
                # update rhs
                rhs[N1 - 1] -= self.values[ii]

            else:
                # update rhs
                rhs[N1 - 1] -= self.values[ii]
                rhs[N2 - 1] += self.values[ii]

        self.rhs = np.array(rhs)

    def branch_voltage(self):
        """
        "branch_voltage"  computes the branch voltages

        :return:
            * self.vb
        """
        # check if the incidence matrix is available
        if self.A is None:
            self.incidence_matrix()

        # check if the solution is available
        if self.x is None:
            print("No solution available")
            return None

        # branch voltages
        self.vb = self.A.transpose() * self.x[:self.node_num, ...]

    def branch_current(self):
        """
        "branch_current"  computes the branch currents

        :return:
            * self.ib
        """
        # check is branch voltages are available
        if self.vb is None:
            self.branch_voltage()

        ibranch = np.zeros_like(self.vb)
        cnt_l = 0
        cnt_v = 0
        for k, (name, val, voltage) in enumerate(zip(self.names, self.values, self.vb)):
            if name[0].upper() == 'R':
                ibranch[k, ...] = voltage / val
            elif name[0].upper() == 'L':
                ibranch[k, ...] = self.x[self.node_num + cnt_l, ...]
                cnt_l += 1
            elif name[0].upper() == 'C':
                if self.analysis[0] == '.op':    # the current is zero, hence, do nothing
                    pass
                elif self.analysis[0] == '.ac':
                    f = float(self.analysis[-1])
                    Xc = -1.0 / (2 * np.pi * f * val)
                    ibranch[k] = voltage / (Xc * 1j)
                elif self.analysis[0] == '.tran':
                    from scipy.interpolate import CubicSpline
                    cs = CubicSpline(self.t, voltage)
                    csd = cs.derivative()
                    ibranch[k, ...] = val * csd(self.t)

            elif name[0].upper() == 'V':
                ibranch[k, ...] = self.x[self.node_num + len(self.isort[1]) + cnt_v, ...]
                cnt_v += 1
            elif name[0].upper() == 'I':
                ibranch[k, ...] = val

        self.ib = ibranch

    def branch_power(self):
        """
        "branch_power"  computes the branch power
        (passive sign convention)

        :return:
            * self.pb: real power for '.op' and complex power for '.ac'
        """

        # check is branch voltages are available
        if self.vb is None:
            self.branch_voltage()

        # check is branch current are available
        if self.ib is None:
            self.branch_current()

        if self.analysis[0] == '.ac':
            self.pb = self.vb * np.conj(self.ib)
        else:
            self.pb = self.vb * self.ib

    def get_voltage(self, arg):
        """
        "get_voltage" computes a voltage across components of between nodes

        :param arg:
            * can be a string in the form 'R1 C1 (2,3) (3,0) (2)'
            * can be a list of node-pair in the form [[2,3],[3,0]]
        :return: voltages (numpy.array)
        """

        if isinstance(arg, str):    # in the input is a string
            # make all uppercase and split
            arg = arg.upper()
            voltage_list = arg.split()

            # initialize output
            if self.analysis[0] == '.tran':
                v = np.zeros((len(voltage_list), self.t.size), dtype=self.x.dtype)
            else:
                v = np.zeros(len(voltage_list), dtype=self.x.dtype)

            # cycle on voltages
            for k, variable in enumerate(voltage_list):

                # remove unused symbols
                remove_char = ('(', ')')
                for char in remove_char:
                    variable = variable.replace(char, '')

                # check variable/node-pair
                if variable in self.names:
                    id = self.names.index(variable)
                    nodes = [n - 1 for n in self.nodes[id] if n != 0]
                    if len(nodes) == 2:
                        v[k,...] = self.x[nodes[0], ...] - self.x[nodes[1], ...]
                    else:
                        v[k,...] = self.x[nodes[0], ...]
                else:
                    nodes = [int(k) - 1 for k in variable.split(',') if k != '0']
                    if len(nodes) == 2:
                        v[k,...] = self.x[nodes[0], ...] - self.x[nodes[1], ...]
                    else:
                        v[k,...] = self.x[nodes[0], ...]

        else: # in the input is a node-pair list

            # check is a single node-pair is given
            if not isinstance(arg[0], list):
                arg = [arg]

            # initialize output
            if self.analysis[0] == '.tran':
                v = np.zeros((len(arg), self.t.size))
            else:
                v = np.zeros(len(arg))

            # cycle on voltages
            for k, index in enumerate(arg):
                nodes = [n - 1 for n in index if n != 0]
                if len(nodes) == 2:
                    v[k, ...] = self.x[nodes[0], ...] - self.x[nodes[1], ...]
                else:
                    v[k, ...] = self.x[nodes[0], ...]


        # remove one dimension for single voltage in .tran
        if len(v.shape) == 2 and v.shape[0] == 1:
            v = v.flatten()

        return v

    def get_current(self, arg):
        """
        "get_current" computes a current in components/branches

        :param arg:
            * can be a string in the form 'R1 C1 (1) (0)'
            * can be a list of branch-index in the form [1, 3, 0]
        :return: currents (numpy.array)
        """

        if isinstance(arg, str):    # in the input is a string
            # make all uppercase and split
            arg = arg.upper()
            current_list = arg.split()
        else:
            current_list = []
            for index in arg:
                current_list.append(self.names[index])


        # initialize output
        if self.analysis[0] == '.tran':
            i = np.zeros((len(current_list), self.t.size), dtype=self.x.dtype)
        else:
            i = np.zeros(len(current_list), dtype=self.x.dtype)

        # cycle on voltages
        for k, variable in enumerate(current_list):

            # remove unused symbols
            remove_char = ('(', ')')
            for char in remove_char:
                variable = variable.replace(char, '')

            if variable not in self.names:
                variable = self.names[int(variable)]

            if (variable[0] == 'R') or (variable[0] == 'C'):

                v = self.get_voltage(variable)
                id = self.names.index(variable)
                if (variable[0] == 'R'):
                    i[k, ...] = v / self.values[id]
                elif (variable[0] == 'C'):
                    from scipy.interpolate import CubicSpline
                    cs = CubicSpline(self.t, v)
                    csd = cs.derivative()
                    i[k, ...] = self.values[id] * csd(self.t)

            elif (variable[0] == 'L'):
                indexL = sorted(self.isort[1])
                for h, il in enumerate(indexL):
                    if variable == self.names[il]:
                        n = self.node_num + h

                i[k, ...] = self.x[n, ...]

            elif (variable[0] == 'V'):
                indexV = sorted(self.isort[3])
                for h, iv in enumerate(indexV):
                    if variable == self.names[iv]:
                        n = self.node_num + len(self.isort[1]) + h

                i[k, ...] = self.x[n, ...]

            elif (variable[0] == 'I'):
                id = self.names.index(variable)
                i[k, ...] = self.values[id]


        # remove one dimension for single voltage in .tran
        if len(i.shape) == 2 and i.shape[0] == 1:
            i = i.flatten()

        return i


    def reorder(self):
        """
        'reorder' reorder the netlist. Order: R, L, C, V and I
        Elements of the same type are ordered in ascending order (eg. R1, R2, R3,...)

        :return:
            * self.isort (index for reordering the network)
        """
        ires = []
        iind = []
        icap = []
        ivolt = []
        icur = []
        for k,ele in enumerate(self.names):
            if ele[0].upper() == 'R':
                ires.append([int(ele[1]),k])
            elif ele[0].upper() == 'L':
                iind.append([int(ele[1]),k])
            elif ele[0].upper() == 'C':
                icap.append([int(ele[1]),k])
            elif ele[0].upper() == 'V':
                ivolt.append([int(ele[1]),k])
            elif ele[0].upper() == 'I':
                icur.append([int(ele[1]),k])

        self.isort = []
        self.isort.append([k for foo, k in sorted(ires)])
        self.isort.append([k for foo, k in sorted(iind)])
        self.isort.append([k for foo, k in sorted(icap)])
        self.isort.append([k for foo, k in sorted(ivolt)])
        self.isort.append([k for foo, k in sorted(icur)])

    def print(self, variable='all', polar=False, message=False):

        # common formatter
        voltage_fmt = 'v({}) = {:10.4f} V\n'
        voltage_fmt_polar = 'v({}) = {:10.4f} V < {:10.4f}°\n'
        current_fmt = 'i({}) = {:10.4f} A\n'
        current_fmt_polar = 'i({}) = {:10.4f} A < {:10.4f}°\n'
        power_fmt = 'p({}) = {:10.4f} {}\n'
        power_fmt_polar = 'p({}) = {:10.4f} {} < {:10.4f}°\n'
        
        # if necessary reorder
        if self.isort is None:
            self.reorder()

        if variable.lower() == 'voltage':
            msg = '==============================================\n'
            msg += '             branch voltages\n'
            msg += '==============================================\n'

            for k, index in enumerate(self.isort):
                if polar:
                    if k == 0:  # resistors
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 1:  # voltage sources
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 2:  # voltage sources
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 3:  # voltage sources
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 4:  # current sources
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += '----------------------------------------------'
                else:
                    if k == 0:  # resistors
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += '----------------------------------------------\n'
                    elif k == 1:  # voltage sources
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += '----------------------------------------------\n'
                    elif k == 2:  # voltage sources
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += '----------------------------------------------\n'
                    elif k == 3:  # voltage sources
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += '----------------------------------------------\n'
                    elif k == 4:  # current sources
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += '----------------------------------------------'

        elif variable.lower() == 'current':
            msg = '=============================================='
            msg += '             branch currents'
            msg += '=============================================='

            for k, index in enumerate(self.isort):
                if polar:
                    if k == 0:  # resistors
                        for h in index:
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 1:  # inductors
                        for h in index:
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 2:  # capacitors
                        for h in index:
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 3:  # voltage sources
                        for h in index:
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 4:  # current sources
                        for h in index:
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += '----------------------------------------------'
                else:
                    if k == 0:  # resistors
                        for h in index:
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += '----------------------------------------------\n'
                    elif k == 1:  # inductors
                        for h in index:
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += '----------------------------------------------\n'
                    elif k == 2:  # capacitors
                        for h in index:
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += '----------------------------------------------\n'
                    elif k == 3:  # voltage sources
                        for h in index:
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += '----------------------------------------------\n'
                    elif k == 4:  # current sources
                        for h in index:
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += '----------------------------------------------'

        if variable.lower() == 'power':

            if self.analysis[0] == '.op':
                unitsR = 'W'
                unitsL = 'W'
                unitsC = 'W'
                unitsV = 'W'
                unitsI = 'W'
            elif self.analysis[0] == '.ac':
                unitsR = 'W'
                unitsL = 'var'
                unitsC = 'var'
                unitsV = 'VA'
                unitsI = 'VA'
            elif self.analysis[0] == '.tran':
                print("Function not supported for transient")
                return -1

            msg = '==============================================\n'
            msg += '             branch powers\n'
            msg += '==============================================\n'

            for k, index in enumerate(self.isort):
                if polar:
                    if k == 0:  # resistors
                        for h in index:
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsR, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 1:  # inductors
                        for h in index:
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsL, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 2:  # capacitors
                        for h in index:
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsC, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 3:  # voltage sources
                        for h in index:
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsV, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 4:  # current sources
                        for h in index:
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsI, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------'
                else:
                    if k == 0:  # resistors
                        for h in index:
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsR)
                            msg += '----------------------------------------------\n'
                    elif k == 1:  # inductors
                        for h in index:
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsL)
                            msg += '----------------------------------------------\n'
                    elif k == 2:  # capacitors
                        for h in index:
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsC)
                            msg += '----------------------------------------------\n'
                    elif k == 3:  # voltage sources
                        for h in index:
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsV)
                            msg += '----------------------------------------------\n'
                    elif k == 4:  # current sources
                        for h in index:
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsI)
                            msg += '----------------------------------------------'

        elif variable.lower() == 'all':

            if self.analysis[0] == '.op':
                unitsR = 'W'
                unitsL = 'W'
                unitsC = 'W'
                unitsV = 'W'
                unitsI = 'W'
            elif self.analysis[0] == '.ac':
                unitsR = 'W'
                unitsL = 'var'
                unitsC = 'var'
                unitsV = 'VA'
                unitsI = 'VA'
            elif self.analysis[0] == '.tran':
                print("Function not supported for transient")
                return -1

            msg = '==============================================\n'
            msg += '               branch quantities\n'
            msg += '==============================================\n'

            for k, index in enumerate(self.isort):
                if polar:
                    if k == 0:  # resistors
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsR , np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 1:  # inductors
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsL , np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 2:  # capacitors
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsC , np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 3:  # voltage sources
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsV, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif k == 4:  # current sources
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsI, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------'

                else:
                    if k == 0:  # resistors
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsR)
                            msg += '----------------------------------------------\n'
                    elif k == 1:  # inductors
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsL)
                            msg += '----------------------------------------------\n'
                    elif k == 2:  # capacitors
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsC)
                            msg += '----------------------------------------------\n'
                    elif k == 3:  # voltage sources
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsV)
                            msg += '----------------------------------------------\n'
                    elif k == 4:  # current sources
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsI)
                            msg += '----------------------------------------------'

        if message:
            return msg
        else:
            print(msg)
            return None

    def plot(self, to_file=False, filename=None, dpi_value=150):

        # check if the analysis is '.tran'
        if self.analysis[0] != '.tran':
            print("plot not supported of analysis: '{}".format(self.analysis[0]))
            return -1
        else:
            # import required library
            import matplotlib.pyplot as plt

            # ===============
            # analyze command
            # ===============
            # set all to uppercase
            plot_cmd = self.plot_cmd.upper()

            # check at least one voltage/current has to be plotted
            plotV = plot_cmd.find('V(')
            plotI = plot_cmd.find('I(')
            if (plotV == -1) and (plotI == -1):
                print("no variables has been provided in this command: {}".format(plot_cmd))
                return -1
            elif (plotV != -1) and (plotI == -1):
                makesubplot = False
                Ylbl = 'voltage (V)'
            elif (plotV == -1) and (plotI != -1):
                makesubplot = False
                Ylbl = 'current (A)'
            elif (plotV != -1) and (plotI != -1):
                makesubplot = True

            # get variables to be plotted
            plot_list = plot_cmd.split()[1:]
            legend_entries = plot_cmd.split()[1:]

            # initialize figure
            if makesubplot:
                hf, axs = plt.subplots(2, 1)
            else:
                hf = plt.figure()

            # cycle on variables
            for k, variable in enumerate(plot_list):
                if variable[0] == 'V':
                    remove_char = ('V','(',')')

                    for char in remove_char:
                        variable = variable.replace(char, '')

                    v = self.get_voltage(variable)

                    if makesubplot:
                        plt.axes(axs[0])
                    plt.plot(self.t, v, label=legend_entries[k])


                elif variable[0] == 'I':
                    remove_char = ('I', '(', ')')

                    for char in remove_char:
                        variable = variable.replace(char, '')

                    i = self.get_current(variable)

                    if makesubplot:
                        plt.axes(axs[1])
                    plt.plot(self.t, i, label=legend_entries[k])



            if makesubplot:
                plt.axes(axs[0])
                plt.ylabel('voltage (V)', fontsize=16)
                plt.grid()
                plt.legend()
                plt.tight_layout()

                plt.axes(axs[1])
                plt.xlabel('time (s)', fontsize=16)
                plt.ylabel('current (A)', fontsize=16)
                plt.grid()
                plt.legend()
                plt.tight_layout()
            else:
                plt.xlabel('time (s)', fontsize=16)
                plt.ylabel(Ylbl, fontsize=16)
                plt.grid()
                plt.legend()
                plt.tight_layout()

            if to_file:
                if filename is None:
                    filename='transient_plot.png'

                hf.savefig(filename, dpi=dpi_value)

        return hf