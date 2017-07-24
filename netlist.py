# ===========================================================
# About the code
# ===========================================================
# This code is part of the project 'spicethon'.
# See README.md for more details
#
# Licensed under the MIT license (see LICENCE)
# Copyright (c) 2017 Luca Giaccone (luca.giaccone@polito.it)
# ===========================================================


# ==================
# imported modules
# ==================
from scipy.sparse import csr_matrix


# ==================
# Class
# ==================
class Network():
    """
    Class that defines the network.

    Attributes:
        * self.A: incidence matrix, node-2-branch
        * self.G: conductance matrix
        * self.rhs: right hand side
    Methods:
        * read_netlist: reads a SPICE netlist (used by '__init__')
    """

    def __init__(self, filename):
        (self.A,
         self.G,
         self.rhs,
         self.node_num,
         self.names,
         self.values) = self.read_netlist(filename)


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

        # initialize incidence matrix terms
        a = []
        a_row = []
        a_col = []

        # initialize conductance terms
        g = []
        g_row = []
        g_col = []

        # initialize rhs terms
        rhs = []
        rhs_row = []
        rhs_col = []

        # initialize source terms
        Vsource = []

        # read netlist
        with open(filename) as f:
            # cycle on lines
            for b, line in enumerate(f):
                # split into a list
                sline = line.split()

                # get branch nodes
                N1 = int(sline[1])
                N2 = int(sline[2])

                # ===================
                # incidence matrix
                # ===================
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

                # ===================
                # conductance matrix
                # ===================
                #
                # detect element type
                if sline[0][0].upper() == 'R': # resistance
                    # add name and value
                    names.append(sline[0])
                    values.append(float(sline[3]))

                    # detect connection
                    if (N1 == 0) or (N2 == 0): # if grounded...
                        # diagonal term
                        g.append(1.0 / float(sline[3]))
                        g_row.append(max([N1, N2]) - 1)
                        g_col.append(max([N1, N2]) - 1)

                        # update counter
                        Nn  = max([N1, N2, Nn])

                    else:                     # if not grounded...
                        # diagonal term
                        g.append(1.0 / float(sline[3]))
                        g_row.append(N1 - 1)
                        g_col.append(N1 - 1)

                        # diagonal term
                        g.append(1.0 / float(sline[3]))
                        g_row.append(N2 - 1)
                        g_col.append(N2 - 1)

                        # N1-N2 term
                        g.append(-1.0 / float(sline[3]))
                        g_row.append(N1 - 1)
                        g_col.append(N2 - 1)

                        # N2-N1 term
                        g.append(-1.0 / float(sline[3]))
                        g_row.append(N2 - 1)
                        g_col.append(N1 - 1)

                        # update counter
                        Nn  = max([N1, N2, Nn])

                # ideal current source
                elif sline[0][0].upper() == 'I':
                    # add name and value
                    names.append(sline[0])
                    values.append(float(sline[3]))

                    if N1 == 0:
                        rhs.append(float(sline[3]))
                        rhs_row.append(N2 - 1)
                        rhs_col.append(0)

                        # update counter
                        Nn  = max([N1, N2, Nn])

                    elif N2 == 0:
                        rhs.append(-float(sline[3]))
                        rhs_row.append(N1 - 1)
                        rhs_col.append(0)

                        # update counter
                        Nn  = max([N1, N2, Nn])
                    else:
                        rhs.append(-float(sline[3]))
                        rhs_row.append(N1 - 1)
                        rhs_col.append(0)

                        rhs.append(float(sline[3]))
                        rhs_row.append(N2 - 1)
                        rhs_col.append(0)

                        # update counter
                        Nn  = max([N1, N2, Nn])


                # temporary storage of voltage sources
                elif sline[0][0].upper() == 'V': # independent voltage sources
                    # add name and value
                    names.append(sline[0])
                    values.append(float(sline[3]))

                    # create temporary list
                    Vsource.append(sline)

                    # update counter
                    Nn  = max([N1, N2, Nn])


        # ideal voltage sources assignment
        for k, vs in enumerate(Vsource):
            # get nodes
            N1 = int(vs[1])
            N2 = int(vs[2])

            # detect connection
            if N1 == 0:               # if grounded to N1 ...
                # negative terminal
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(Nn + k)

                # negative terminal
                g.append(-1)
                g_row.append(Nn + k)
                g_col.append(N2 - 1)

                # update rhs
                rhs.append(float(vs[3]))
                rhs_row.append(Nn + k)
                rhs_col.append(0)

            elif N2 == 0:               # if grounded to N2 ...
                # positive terminal
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(Nn + k)

                # positive terminal
                g.append(1)
                g_row.append(Nn + k)
                g_col.append(N1 - 1)

                # update rhs
                rhs.append(float(vs[3]))
                rhs_row.append(Nn + k)
                rhs_col.append(0)

            else:                      # if not grounded ...
                # positive terminal
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(Nn + k)

                # positive terminal
                g.append(1)
                g_row.append(Nn + k)
                g_col.append(N1 - 1)

                # negative terminal
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(Nn + k)

                # negative terminal
                g.append(-1)
                g_row.append(Nn + k)
                g_col.append(N2 - 1)

                # update rhs
                rhs.append(float(vs[3]))
                rhs_row.append(Nn + k)
                rhs_col.append(0)

        # create conductance matrix
        A = csr_matrix((a, (a_row, a_col)))

        # create conductance matrix
        G = csr_matrix((g,(g_row,g_col)))

        # create conductance matrix
        rhs = csr_matrix((rhs,(rhs_row,rhs_col)))

        return A, G, rhs, Nn, names, values


    def dc_solve(self):
        """
        "dc_solve" computes the DC operating point

        :return:
            * self.x
        """
        from scipy.sparse.linalg import spsolve
        self.x = spsolve(self.G, self.rhs)


    def branch_voltage(self):
        """
        "branch_voltage"  computes the branch voltages

        :return:
            * self.vb
        """
        # branch voltages
        self.vb = self.A.transpose() * self.x[:self.node_num]


    def branch_current(self, verbose='y'):
        """
        "branch_current"  computes the branch currents

        :param verbose: set to 'y/n' to print results (default 'y')
        :return:
            * self.ib
        """
        # check if branch voltages are already computed
        if not hasattr(self, 'vb'):
            self.branch_voltage()

        ib = []
        cnt_v = 0
        for name, val, voltage in zip(self.names, self.values, self.vb):
            if name[0].upper() == 'R':
                ib.append(voltage/val)
            elif name[0].upper() == 'V':
                ib.append(self.x[self.node_num + cnt_v])
                cnt_v += 1
            elif name[0].upper() == 'I':
                ib.append(val)

        self.ib = ib


    def reorder(self):
        """
        'reorder' reorder the netlist

        :return:
            * self.isort
        """
        ires = []
        ivolt = []
        icur = []
        for k,ele in enumerate(self.names):
            if ele[0].upper() == 'R':
                ires.append([int(ele[1]),k])
            elif ele[0].upper() == 'V':
                ivolt.append([int(ele[1]),k])
            elif ele[0].upper() == 'I':
                icur.append([int(ele[1]),k])

        self.isort = []
        self.isort.append([k for foo, k in sorted(ires)])
        self.isort.append([k for foo, k in sorted(ivolt)])
        self.isort.append([k for foo, k in sorted(icur)])
