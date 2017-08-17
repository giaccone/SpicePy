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
        * self.nodes: component nodes (only two-port right now)
        * self.node_num: number of nodes in the network
    Other attributes (default is None):
        * self.A:
        * self.G:
        * self.rhs:
        * self.isort:
        * self.x:
        * self.vb:
        * self.ib:
    Methods:
        * read_netlist(self): reads a SPICE netlist (used by '__init__')
        * incidence_matrix(self):
    """

    def __init__(self, filename):
        """
        TO DO...
        """
        # basilar attributes
        (self.names,
         self.values,
         self.nodes,
         self.node_num,
         self.analysis) = self.read_netlist(filename)

        # initialization of other possible attributes
        self.A = None
        self.G = None
        self.C = None
        self.rhs = None
        self.isort = None
        self.x = None
        self.vb = None
        self.ib = None

    def read_netlist(self, filename):
        """
        'readNetlist' reads a SPICE netlist

        :param filename: file name with the netlist
        :return: tuple including:
            * A incidence matrix, node-2-branch (in sparse form)
            * G conductance matrix (in sparse form)
            * rhs right end side (in sparse form)
        """

        # initialize counter for number of nodes
        Nn = 0

        # initialize variable names and values
        names = []
        values = []
        nodes = []

        # read netlist
        with open(filename) as f:
            # cycle on lines
            for b, line in enumerate(f):
                # split into a list
                sline = line.split()

                if sline[0][0] == '.':  # analysis identifier
                    analysis = sline
        f.close()

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

                    # update counter
                    Nn = max([N1, N2, Nn])

                # capacitor
                elif sline[0][0].upper() == 'C':
                    # add name and value
                    names.append(sline[0])
                    values.append(float(sline[3]))
                    nodes.append([N1, N2])

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
        return names, values, nodes, Nn, analysis

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

        :return: G, conductance matrix (including terms related to independent voltage sources)
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
        indexL = self.isort[1]
        indexV = self.isort[3]

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
        TO DO...

        :return:
        """

        # initialize conductance terms
        c = []
        c_row = []
        c_col = []

        # reorder if necessary
        if self.isort is None:
            self.reorder()

        # get index
        indexL = self.isort[1]
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
        indexV = self.isort[3]
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

        self.rhs = rhs

    def reorder(self):
        """
        'reorder' reorder the netlist

        :return:
            * self.isort
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
