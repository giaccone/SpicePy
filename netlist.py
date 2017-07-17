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
        * self.G: conductance matrix
        * self.rhs: right hand side
    Methods:
        * read_netlist: reads a SPICE netlist (used by '__init__')
    """

    def __init__(self, filename):
        (self.G,
         self.rhs) = self.read_netlist(filename)

    def read_netlist(self, filename):
        """
        'readNetlist' reads a SPICE netlist

        :param filename: file name with the netlist
        :return: tuple including:
            * G conductance matrix (in sparse form)
            * rhs right end side (returned as a list)
        """

        # initialize counter for number of nodes
        Nn = 0

        # initialize conductance terms
        g = []
        g_row = []
        g_col = []

        # initialize source terms
        Vsource = []

        # read netlist
        with open(filename) as f:
            # cycle on lines
            for line in f:
                # split into a list
                sline = line.split()

                # detect element type
                if sline[0][0] == 'R': # resistance
                    # get nodes
                    N1 = int(sline[1])
                    N2 = int(sline[2])

                    # detect connection
                    if (N1 == 0) or (N2 == 0): # if grounded...
                        # diagonal term
                        g.append(1.0 / float(sline[3]))
                        g_row.append(max([N1, N2]) - 1)
                        g_col.append(max([N1, N2]) - 1)

                        # update counter
                        Nn += 1
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
                        Nn += 1

                elif sline[0][0] == 'V': # independent voltage sources
                    # create temporary list
                    Vsource.append(sline)

        # initialization of rhs
        rhs = [0] * Nn

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


        # create conductance matrix
        G = csr_matrix((g,(g_row,g_col)))

        return G, rhs
