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
import spicepy.transient_sources as tsr

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
        * self.source_type: type of transient sources (stored as dict)
        * self.nodes: component nodes (only two-port right now)
        * self.node_label2num: dictionary to convert node-labels to local node-number
        * self.node_num: number of nodes in the network
        * self.analysis: type of analysis
        * self.plot_cmd: plot directive for transient analysis
        * self.tf_cmd: transfer function definition (for .ac multi-freq)
    Other attributes (default is None):
        * self.A:
        * self.G:
        * self.C:
        * self.rhs:
        * self.isort:
        * self.t:
        * self.f:
        * self.x:
        * self.vb:
        * self.ib:
        * self.pb:
    Methods:
        * read_netlist(self): reads a SPICE netlist (used by '__init__')
        * incidence_matrix(self, filename):
        * conductance_matrix(self):
        * dynamic_matrix(self):
        * rhs_matrix(self):
        * branch_voltage(self):
        * branch_current(self):
        * branch_power(self):
        * get_voltage(self, arg):
        * get_current(self, arg):
        * reorder(self):
        * print(self, variable='all', polar=False, message=False):
        * plot(self, to_file=False, filename=None, dpi_value=150):
    """

    def __repr__(self):
        return "SpicePy.Network: {} analysis".format(self.analysis[0])

    def __str__(self):
        # create local dictionary to convert node-numbers to node-labels
        num2node_label = {num: name for name, num in self.node_label2num.items()}

        # build message to print
        msg = '------------------------\n'
        msg += '    SpicePy.Network:\n'
        msg += '------------------------\n'
        for ele, nodes, val in zip(self.names, self.nodes, self.values):
            # if val is a list --> ele is a transient source
            if isinstance(val, list):
                if self.source_type[ele] == 'pwl':
                    fmt = "{} {} {} {}(" + "{} " * (len(val[0]) - 1) + "{})\n"
                    msg += fmt.format(ele, num2node_label[nodes[0]], num2node_label[nodes[1]], self.source_type[ele], *val[0])
                else:
                    fmt = "{} {} {} {}(" + "{} " * (len(val) - 1) + "{})\n"
                    msg += fmt.format(ele, num2node_label[nodes[0]], num2node_label[nodes[1]], self.source_type[ele], *val)
            # controlled sources
            elif (ele[0].upper() == 'E') | (ele[0].upper() == 'G'):
                msg += "{} {} {} {} {} {}\n".format(ele, num2node_label[nodes[0]], num2node_label[nodes[1]],
                                                    self.control_source[ele][0],
                                                    self.control_source[ele][1],
                                                    val)
            elif (ele[0].upper() == 'F') | (ele[0].upper() == 'H'):
                msg += "{} {} {} {} {}\n".format(ele, num2node_label[nodes[0]], num2node_label[nodes[1]],
                                                 self.control_source[ele],
                                                 val)
            # if val is complex --> ele is a phasor
            elif np.iscomplex(val):
                msg += "{} {} {} {} {}\n".format(ele, num2node_label[nodes[0]], num2node_label[nodes[1]], np.abs(val), np.angle(val) * 180/np.pi)
            # if ele is C or L
            elif ele[0].upper() == 'C' or ele[0].upper() == 'L':
                # check if an i.c. is present and print it
                if ele in self.IC:
                    msg += "{} {} {} {} ic={}\n".format(ele, num2node_label[nodes[0]], num2node_label[nodes[1]], val, self.IC[ele])
                # otherwise...
                else:
                    msg += "{} {} {} {}\n".format(ele, num2node_label[nodes[0]], num2node_label[nodes[1]], val)
            # otherwise...general case -->  ele n+ n- val
            else:
                msg += "{} {} {} {}\n".format(ele, num2node_label[nodes[0]], num2node_label[nodes[1]], val)

        # add analysis
        msg += " ".join(self.analysis) + '\n'

        # if a plot command is present, add it
        if self.plot_cmd is not None:
            msg += self.plot_cmd + '\n'

        # if a transfer function definition is present, add it
        if self.tf_cmd is not None:
            msg += self.tf_cmd + '\n'

        # add number of nodes (reference node is included) and number of branches
        msg += '------------------------\n'
        msg += '* number of nodes {}\n'.format(self.node_num + 1)
        msg += '* number of branches {}\n'.format(len(self.names))
        msg += '------------------------\n'

        return msg

    def __init__(self, filename):
        """
        __init__ initializes common attributes using read_netlist method

        :param filename:
        """

        # initialization of other possible attributes
        self.A = None
        self.G = None
        self.C = None
        self.rhs = None
        self.isort = None
        self.t = None
        self.f = None
        self.x = None
        self.vb = None
        self.ib = None
        self.pb = None
        self.unit_prefix = {'meg': 'e6', 'f': 'e-15', 'p': 'e-12', 'n': 'e-9', 'u': 'e-6', 'm': 'e-3', 'k': 'e3', 'g': 'e9', 't': 'e12'}

        # common attributes
        (self.names,
         self.values,
         self.IC,
         self.source_type,
         self.control_source,
         self.nodes,
         self.node_label2num,
         self.node_num,
         self.analysis,
         self.plot_cmd,
         self.tf_cmd) = self.read_netlist(filename)

    def read_netlist(self, filename):
        """
        'readNetlist' reads a SPICE netlist

        :param filename: file name with the netlist
        :return:
            * names: element names
            * values: element values
            * IC: initial conditions
            * source_type: type of transient source
            * nodes: element nodes
            * node_labels2num: dictionary to convert node-labels to local node-number
            * Nn: number of nodes in the netlist
            * analysis: type of analysis
            * plot_cmd: plot command
        """

        # initialize counter for number of nodes
        Nn = 0

        # initialize variable names and values
        names = []
        values = []
        node_labels = []
        IC = {}
        source_type = {}
        control_source = {}
        plot_cmd = None
        tf_cmd = None
        analysis = None

        # initial letter of all available components
        initials = ['V', 'I', 'R', 'C', 'L', 'E', 'F', 'G', 'H']
        components = []

        # 1) get the analysis type
        # 2) catch plot command (if any)
        # 3) filter out comments
        with open(filename) as f:
            # cycle on lines
            for b, line in enumerate(f):

                # look for inline comments
                if ';' in line:
                    # remove comment
                    line = line[:line.index(';')]

                # if it is not a line-comment
                if line[0] != '*':

                    # check if line describes a component
                    if line[0].upper() in initials:
                        # remove carriage return
                        line = line.replace('\n', '')
                        # add to component list
                        components.append(line)

                    # check if line describes a command
                    elif line[0] == '.':
                        # remove carriage return
                        line = line.replace('\n', '')

                        if line.lower().find('.end') != -1:    # if .end is reached exit
                            break

                        elif line.lower().find('.plot') != -1:    # if .plot is reached save it
                            plot_cmd = line

                        elif line.lower().find('.tf') != -1:    # if .tf is reached save it
                            tf_cmd = line

                        elif line.lower().find('.backanno') != -1:    # skip .backanno if present
                            pass

                        else:    # save analysis type
                            # split into a list
                            analysis = line.split()

                else:
                    pass

        # cycle on component list
        for line in components:

            # look for transient sources (if analysis is .tran)
            if analysis[0] == '.tran':
                # loop on all possible transient sources
                time_sources = ['pwl', 'pulse', 'sin', 'exp']
                for source in time_sources:
                    # when one is found
                    if source in line.lower():
                        # get index of the related string
                        index = line.lower().index(source)
                        # split string before its name
                        sline = line[:index].split()
                        # 1) remove '(' and ')' in the string after its name 2) and split
                        param = line[index:].replace('(',' ').replace(')',' ').split()
                        # append transient-source name
                        sline.append(param[0])
                        # append parameters
                        sline.append(param[1:])
                        break
                # if component is not a transient source, catch it normally
                else:
                    # split into a list
                    sline = line.split()

            # if not '.tran', catch them all normally
            else:
                # split into a list
                sline = line.split()

            # detect element type
            if sline[0][0].upper() == 'R':  # resistance
                # add name and value
                names.append(sline[0])
                values.append(float(self.convert_unit(sline[3])))
                node_labels.append(sline[1:3])

            # inductor
            elif sline[0][0].upper() == 'L':
                # add name, value and nodes
                names.append(sline[0])
                values.append(float(self.convert_unit(sline[3])))
                node_labels.append(sline[1:3])

                # for '.tran'
                if analysis[0] == '.tran':
                    # check presence of i.c.
                    if len(sline) == 5:
                        if sline[4].lower().find('ic') != -1:
                            IC[sline[0]] = float(self.convert_unit(sline[4].split('=')[1]))
                        else:
                            #IC[sline[0]] = 'Please check this --> ' + sline[-1]
                            # print("Warning: wrong definition of IC for {} --> ".format(sline[0]) + IC[sline[0]])
                            raise ValueError("Warning: wrong definition of IC for {} --> ".format(sline[0]) + IC[sline[0]])

                    # add ic=0 if i.c. non provided by the user
                    else:
                        IC[sline[0]] = 0

            # capacitor
            elif sline[0][0].upper() == 'C':
                # add name and value
                names.append(sline[0])
                values.append(float(self.convert_unit(sline[3])))
                node_labels.append(sline[1:3])

                if analysis[0] == '.tran':
                    if len(sline) == 5:
                        if sline[4].lower().find('ic') != -1:
                            IC[sline[0]] = float(self.convert_unit(sline[4].split('=')[1]))
                        else:
                            #IC[sline[0]] = 'Please check this --> ' + sline[-1]
                            #print("Warning: wrong definition of IC for {} --> ".format(sline[0]) + IC[sline[0]])
                            raise ValueError("Warning: wrong definition of IC for {} --> ".format(sline[0]) + IC[sline[0]])

                    # add ic=0 if i.c. non provided by the user
                    else:
                        IC[sline[0]] = 0

            # independent current source
            elif sline[0][0].upper() == 'I':
                # add name and nodes
                names.append(sline[0])
                node_labels.append(sline[1:3])

                # if '.ac' and phase is present:
                if (analysis[0] == '.ac') & (len(sline) == 5):
                    values.append(float(self.convert_unit(sline[3])) * (
                    np.cos(float(sline[4]) * pi / 180) + np.sin(float(sline[4]) * pi / 180) * 1j))
                # if '.tran' ...
                elif analysis[0] == '.tran':
                    # if is a transient source
                    if isinstance(sline[-1], list):
                        source_type[sline[0]] = sline[-2]
                        if source_type[sline[0]] == 'pwl':
                            values.append([[float(self.convert_unit(k)) for k in sline[-1]]])
                        else:
                            values.append([float(self.convert_unit(k)) for k in sline[-1]])
                    # otherwise...
                    else:
                        values.append(float(self.convert_unit(sline[3])))
                # otherwise...
                else:
                    values.append(float(self.convert_unit(sline[3])))

            # independent voltage sources
            elif sline[0][0].upper() == 'V':
                # add name and nodes
                names.append(sline[0])
                node_labels.append(sline[1:3])

                # if '.ac' and phase is present:
                if (analysis[0] == '.ac') & (len(sline) == 5):
                    values.append(float(self.convert_unit(sline[3])) * (
                    np.cos(float(sline[4]) * pi / 180) + np.sin(float(sline[4]) * pi / 180) * 1j))
                # if '.tran'
                elif analysis[0] == '.tran':
                    # if is a transient source
                    if isinstance(sline[-1], list):
                        source_type[sline[0]] = sline[-2]
                        if source_type[sline[0]] == 'pwl':
                            values.append([[float(self.convert_unit(k)) for k in sline[-1]]])
                        else:
                            values.append([float(self.convert_unit(k)) for k in sline[-1]])
                    # otherwise...
                    else:
                        values.append(float(self.convert_unit(sline[3])))
                # otherwise...
                else:
                    values.append(float(self.convert_unit(sline[3])))

            # VCVS or VCCS
            elif (sline[0][0].upper() == 'E') | (sline[0][0].upper() == 'G'):
                # add name and nodes
                names.append(sline[0])
                node_labels.append(sline[1:3])
                # get control nodes
                control_source[sline[0]] = sline[3:5]
                # get gain
                values.append(float(self.convert_unit(sline[5])))

            # CCCS or CCVS
            elif (sline[0][0].upper() == 'F') | (sline[0][0].upper() == 'H'):
                # add name and nodes
                names.append(sline[0])
                node_labels.append(sline[1:3])
                # get control Vsens
                control_source[sline[0]] = sline[3]
                # get gain
                values.append(float(self.convert_unit(sline[4])))

        # reordering nodes
        unique_names, ii = np.unique(node_labels, return_inverse=True)
        if '0' not in unique_names:
            raise ValueError("Error: the network does not include node '0'")

        nodes = np.reshape(ii, (len(node_labels),2))
        # link name-2-number
        node_labels2num = {}
        for k , label in enumerate(np.unique(node_labels)):
            node_labels2num[label] = k

        Nn = nodes.max()

        # return network structure
        return names, values, IC, source_type, control_source, nodes, node_labels2num, Nn, analysis, plot_cmd, tf_cmd

    def convert_unit(self, string_value):
        """
        'convert_unit' convert a unit-prefix in a string with the releted value

        :param string_value: a string with a unit prefic (e.g. '10.5k')
        :return: the same sting with the numerical value of the unit-prefix (e.g. '10.5e3')
        """
        for prefix, prefix_value in self.unit_prefix.items():
            if prefix in string_value.lower():
                string_value = string_value.lower().replace(prefix, prefix_value)
                break

        return string_value

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
        indexE = sorted(self.isort[5])
        indexF = sorted(self.isort[6])
        indexG = sorted(self.isort[7])
        indexH = sorted(self.isort[8])

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


        # cycle on VCVS
        for k, ie in enumerate(indexE):
            # get nodes
            N1, N2 = self.nodes[ie]

            # detect connection
            if N1 == 0:  # if grounded to N1 ...
                # negative terminal
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(self.node_num + len(indexL) + len(indexV) + k)

                # negative terminal
                g.append(-1)
                g_row.append(self.node_num + len(indexL) + len(indexV) + k)
                g_col.append(N2 - 1)

            elif N2 == 0:  # if grounded to N2 ...
                # positive terminal
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(self.node_num + len(indexL) + len(indexV) + k)

                # positive terminal
                g.append(1)
                g_row.append(self.node_num + len(indexL) + len(indexV) + k)
                g_col.append(N1 - 1)

            else:  # if not grounded ...
                # positive terminal
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(self.node_num + len(indexL) + len(indexV) + k)

                # positive terminal
                g.append(1)
                g_row.append(self.node_num + len(indexL) + len(indexV) + k)
                g_col.append(N1 - 1)

                # negative terminal
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(self.node_num + len(indexL) + len(indexV) + k)

                # negative terminal
                g.append(-1)
                g_row.append(self.node_num + len(indexL) + len(indexV) + k)
                g_col.append(N2 - 1)

            # get control nodes
            N1, N2 = [self.node_label2num[k] for k in self.control_source[self.names[ie]]]

            # detect connection
            if N1 == 0:  # if grounded to N1 ...
                # negative terminal
                g.append(self.values[ie])
                g_row.append(self.node_num + len(indexL) + len(indexV) + k)
                g_col.append(N2 - 1)

            elif N2 == 0:  # if grounded to N2 ...
                # positive terminal
                g.append(-self.values[ie])
                g_row.append(self.node_num + len(indexL) + len(indexV) + k)
                g_col.append(N1 - 1)

            else:  # if not grounded ...
                # positive terminal
                g.append(-self.values[ie])
                g_row.append(self.node_num + len(indexL) + len(indexV) + k)
                g_col.append(N1 - 1)

                # negative terminal
                g.append(self.values[ie])
                g_row.append(self.node_num + len(indexL) + len(indexV) + k)
                g_col.append(N2 - 1)

        # cycle on CCCSs
        for k, indF in enumerate(indexF):
            # get nodes
            N1, N2 = self.nodes[indF]
            # get Vsens
            Vsens = self.control_source[self.names[indF]]
            # get index of Vsens
            if Vsens[0].upper() == 'V':
                h = sorted(self.isort[3]).index(self.names.index(Vsens))
                n = self.node_num + len(self.isort[1]) + h
            elif Vsens[0].upper() == 'E':
                h = sorted(self.isort[5]).index(self.names.index(Vsens))
                n = self.node_num + len(self.isort[1]) + len(self.isort[3]) + h
            elif Vsens[0].upper() == 'H':
                h = sorted(self.isort[8]).index(self.names.index(Vsens))
                n = self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + h

            if N1 == 0: # if grounded to N1 ...
                g.append(-self.values[indF])
                g_row.append(N2 - 1)
                g_col.append(n)
            elif N2 == 0: # if grounded to N2 ...
                g.append(self.values[indF])
                g_row.append(N1 - 1)
                g_col.append(n)
            else: # if not grounded ...
                g.append(self.values[indF])
                g_row.append(N1 - 1)
                g_col.append(n)
                #
                g.append(-self.values[indF])
                g_row.append(N2 - 1)
                g_col.append(n)

        # cycle on VCCSs
        for k, iG in enumerate(indexG):
            # get nodes
            N1, N2 = self.nodes[iG]
            # get control nodes
            Nc1, Nc2 = [self.node_label2num[k] for k in self.control_source[self.names[iG]]]

            if (N1 != 0) & (Nc1 != 0):
                g.append(self.values[iG])
                g_row.append(N1 - 1)
                g_col.append(Nc1 - 1)
            if (N2 != 0) & (Nc2 != 0):
                g.append(self.values[iG])
                g_row.append(N2 - 1)
                g_col.append(Nc2 - 1)
            if (N1 != 0) & (Nc2 != 0):
                g.append(-self.values[iG])
                g_row.append(N1 - 1)
                g_col.append(Nc2 - 1)
            if (N2 != 0) & (Nc1 != 0):
                g.append(-self.values[iG])
                g_row.append(N2 - 1)
                g_col.append(Nc1 - 1)

        # cycle on CCVSs
        for k, iH in enumerate(indexH):
            # get nodes
            N1, N2 = self.nodes[iH]
            # get Vsens
            Vsens = self.control_source[self.names[iH]]
            # get index of Vsens
            if Vsens[0].upper() == 'V':
                h = sorted(self.isort[3]).index(self.names.index(Vsens))
                n = self.node_num + len(self.isort[1]) + h
            elif Vsens[0].upper() == 'E':
                h = sorted(self.isort[5]).index(self.names.index(Vsens))
                n = self.node_num + len(self.isort[1]) + len(self.isort[3]) + h
            elif Vsens[0].upper() == 'H':
                h = sorted(self.isort[8]).index(self.names.index(Vsens))
                n = self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + h

            if N1 != 0:
                g.append(1)
                g_row.append(N1 - 1)
                g_col.append(self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + k)
                #
                g.append(1)
                g_row.append(self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + k)
                g_col.append(N1 - 1)
            if N2 != 0:
                g.append(-1)
                g_row.append(N2 - 1)
                g_col.append(self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + k)
                #
                g.append(-1)
                g_row.append(self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + k)
                g_col.append(N2 - 1)

            g.append(-self.values[iH])
            g_row.append(self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + k)
            g_col.append(n)

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

        if self.analysis[0] == '.tran':
            def fun(t):
                # initialize rhs
                rhs = [0] * (self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + len(self.isort[8]))

                # get index
                NL = len(self.isort[1])
                indexV = sorted(self.isort[3])
                indexI = self.isort[4]

                # cycle on independent voltage sources
                for k, iv in enumerate(indexV):
                    if isinstance(self.values[iv], list):
                        tsr_fun = getattr(tsr, self.source_type[self.names[iv]])
                        rhs[self.node_num + NL + k] += tsr_fun(*self.values[iv], t=t)
                    else:
                        # update rhs
                        rhs[self.node_num + NL + k] += self.values[iv]

                # cycle on independent current sources
                for ii in indexI:
                    # get nodes
                    N1, N2 = self.nodes[ii]

                    if isinstance(self.values[ii], list):
                        tsr_fun = getattr(tsr, self.source_type[self.names[ii]])
                        if N1 == 0:
                            # update rhs
                            rhs[N2 - 1] += tsr_fun(*self.values[ii], t=t)

                        elif N2 == 0:
                            # update rhs
                            rhs[N1 - 1] -= tsr_fun(*self.values[ii], t=t)

                        else:
                            # update rhs
                            rhs[N1 - 1] -= tsr_fun(*self.values[ii], t=t)
                            rhs[N2 - 1] += tsr_fun(*self.values[ii], t=t)

                    else:
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

                return np.array(rhs)

            return fun

        else:

            # initialize rhs
            rhs = [0] * (self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + len(self.isort[8]))

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

        self.ib = np.zeros_like(self.vb)

        for k, name in enumerate(self.names):
            self.ib[k, ...] = self.get_current(name)

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

        if self.analysis[0].lower() == '.ac':
            self.pb = self.vb * np.conj(self.ib)
        else:
            self.pb = self.vb * self.ib

    def frequency_span(self):
        """
        "frequency_span" generates the frequency span for .ac analysys

        :return:
            * self.f: frequency array
        """
        if self.analysis[0].lower() != '.ac':
            raise ValueError("frequency_span works only for .ac analyses")

        if self.analysis[1].lower() == 'lin':
            npt = int(self.analysis[2])
            fs = float(self.convert_unit(self.analysis[3]))
            fe = float(self.convert_unit(self.analysis[4]))
            self.f = np.linspace(fs, fe, npt)

        elif self.analysis[1].lower() == 'dec':
            npt_d = float(self.analysis[2])
            fs = np.log10(float(self.convert_unit(self.analysis[3])))
            fe = np.log10(float(self.convert_unit(self.analysis[4])))
            self.f = np.logspace(fs, fe, int(np.ceil(npt_d * (fe - fs)).item()))

        elif self.analysis[1].lower() == 'oct':
            npt_d = float(self.analysis[2])
            fs = np.log2(float(self.convert_unit(self.analysis[3])))
            fe = np.log2(float(self.convert_unit(self.analysis[4])))
            self.f = np.logspace(fs, fe, int(np.ceil(npt_d * (fe - fs)).item()), base=2)

        if self.f.size == 1:
            self.f = np.asscalar(self.f)

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
            arg = arg
            voltage_list = arg.split()

            # initialize output
            if self.analysis[0].lower() == '.tran':
                v = np.zeros((len(voltage_list), self.t.size), dtype=self.x.dtype)
            elif (self.analysis[0].lower() == '.ac') & (not np.isscalar(self.f)):
                v = np.zeros((len(voltage_list), self.f.size), dtype=self.x.dtype)
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
                        if self.nodes[id][0] != 0:
                            sign = 1
                        else:
                            sign = -1
                        
                        v[k,...] = self.x[nodes[0], ...] * sign
                else:
                    node_labels = variable.split(',')
                    node_number = [self.node_label2num[k] for k in node_labels]
                    nodes = [int(k) - 1 for k in node_number if k != 0]
                    if len(nodes) == 2:
                        v[k,...] = self.x[nodes[0], ...] - self.x[nodes[1], ...]
                    else:
                        if node_number[0] != 0:
                            sign = 1
                        else:
                            sign = -1

                        v[k,...] = self.x[nodes[0], ...] * sign

        else:  # if the input is a node-pair list

            # check is a single node-pair is given
            if not isinstance(arg[0], list):
                arg = [arg]

            # initialize output
            if self.analysis[0].lower() == '.tran':
                v = np.zeros((len(arg), self.t.size), dtype=self.x.dtype)
            else:
                v = np.zeros(len(arg), dtype=self.x.dtype)

            # cycle on voltages
            for k, node_labels in enumerate(arg):
                node_number = [self.node_label2num[str(k)] for k in node_labels]
                nodes = [int(k) - 1 for k in node_number if k != 0]

                if len(nodes) == 2:
                    v[k, ...] = self.x[nodes[0], ...] - self.x[nodes[1], ...]
                else:
                    if node_number[0] != 0:
                        sign = 1
                    else:
                        sign = -1
                    
                    v[k, ...] = self.x[nodes[0], ...] * sign

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
        if self.analysis[0].lower() == '.tran':
            i = np.zeros((len(current_list), self.t.size), dtype=self.x.dtype)
        elif (self.analysis[0].lower() == '.ac') & (not np.isscalar(self.f)):
            i = np.zeros((len(current_list), self.f.size), dtype=self.x.dtype)
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

            # resistor or capacitor
            if (variable[0] == 'R') or (variable[0] == 'C'):

                # get voltage and index
                v = self.get_voltage(variable)
                id = self.names.index(variable)
                # resistor
                if variable[0] == 'R':
                    i[k, ...] = v / self.values[id]
                # capacitor
                elif variable[0] == 'C':
                    # time derivative for .tran
                    if self.analysis[0].lower() == '.tran':
                        from scipy.interpolate import CubicSpline
                        cs = CubicSpline(self.t, v)
                        csd = cs.derivative()
                        i[k, ...] = self.values[id] * csd(self.t)
                    # symbolic method in .ac
                    elif self.analysis[0].lower() == '.ac':
                        Xc = -1.0 / (2 * np.pi * self.f * self.values[id])
                        i[k, ...] = v / (Xc * 1j)
                    # .OP: do nothing --> ic already zero

            # inductor
            elif variable[0] == 'L':
                # get sub-index of 'L'
                h = sorted(self.isort[1]).index(self.names.index(variable))
                # index of the related current in the solution
                n = self.node_num + h
                # get current
                i[k, ...] = self.x[n, ...]

            # voltage generator
            elif variable[0] == 'V':
                # get sub-index of 'V'
                h = sorted(self.isort[3]).index(self.names.index(variable))
                # index of the related current in the solution
                n = self.node_num + len(self.isort[1]) + h
                # get current
                i[k, ...] = self.x[n, ...]

            # current generator
            elif variable[0] == 'I':
                id = self.names.index(variable)
                if isinstance(self.values[id], list):
                    tsr_fun = getattr(tsr, self.source_type[self.names[id]])
                    i[k, ...] = tsr_fun(*self.values[id], self.t)
                else:
                    i[k, ...] = self.values[id]

            # VCVS
            elif variable[0] == 'E':
                # get sub-index of 'E'
                h = sorted(self.isort[5]).index(self.names.index(variable))
                # index of the related current in the solution
                n = self.node_num + len(self.isort[1]) + len(self.isort[3]) + h
                # get current
                i[k, ...] = self.x[n, ...]

            # CCCS
            elif variable[0] == 'F':
                # get Vsens
                Vsens = self.control_source[variable]
                # get current
                i[k, ...] = self.get_current(Vsens) * self.values[self.names.index(variable)]

            # VCCS
            elif variable[0] == 'G':
                i[k, ...] = self.values[self.names.index(variable)] * self.get_voltage(self.control_source[variable])

            # CCVS
            elif variable[0] == 'H':
                # get sub-index of 'H'
                h = sorted(self.isort[8]).index(self.names.index(variable))
                # index of the related current in the solution
                n = self.node_num + len(self.isort[1]) + len(self.isort[3]) + len(self.isort[5]) + h
                # get current
                i[k, ...] = self.x[n, ...]

        # remove one dimension for single voltage in .tran
        if len(i.shape) == 2 and i.shape[0] == 1:
            i = i.flatten()

        return i


    def reorder(self):
        """
        'reorder' reorder the netlist. Order: R, L, C, V , I, E
        Elements of the same type are ordered in ascending order (eg. R1, R2, R3,...)

        :return:
            * self.isort (index for reordering the network)
        """
        ires = []
        iind = []
        icap = []
        ivolt = []
        icur = []
        ivcvs = []
        icccs = []
        ivccs = []
        iccvs = []
        for k, ele in enumerate(self.names):
            if ele[0].upper() == 'R':
                ires.append([int(ele[1:]), k])
            elif ele[0].upper() == 'L':
                iind.append([int(ele[1:]), k])
            elif ele[0].upper() == 'C':
                icap.append([int(ele[1:]), k])
            elif ele[0].upper() == 'V':
                ivolt.append([int(ele[1:]), k])
            elif ele[0].upper() == 'I':
                icur.append([int(ele[1:]), k])
            elif ele[0].upper() == 'E':
                ivcvs.append([int(ele[1:]), k])
            elif ele[0].upper() == 'F':
                icccs.append([int(ele[1:]), k])
            elif ele[0].upper() == 'G':
                ivccs.append([int(ele[1:]), k])
            elif ele[0].upper() == 'H':
                iccvs.append([int(ele[1:]), k])

        self.isort = []
        self.isort.append([k for foo, k in sorted(ires)])
        self.isort.append([k for foo, k in sorted(iind)])
        self.isort.append([k for foo, k in sorted(icap)])
        self.isort.append([k for foo, k in sorted(ivolt)])
        self.isort.append([k for foo, k in sorted(icur)])
        self.isort.append([k for foo, k in sorted(ivcvs)])
        self.isort.append([k for foo, k in sorted(icccs)])
        self.isort.append([k for foo, k in sorted(ivccs)])
        self.isort.append([k for foo, k in sorted(iccvs)])

    def print(self, variable='all', polar=False, message=False):

        # common formatter
        voltage_fmt = 'v({}) = {:10.4g} V\n'
        voltage_fmt_polar = 'v({}) = {:10.4g} V < {:10.4g}°\n'
        current_fmt = 'i({}) = {:10.4g} A\n'
        current_fmt_polar = 'i({}) = {:10.4g} A < {:10.4g}°\n'
        power_fmt = 'p({}) = {:10.4g} {}\n'
        power_fmt_polar = 'p({}) = {:10.4g} {} < {:10.4g}°\n'
        
        # if necessary reorder
        if self.isort is None:
            self.reorder()

        if variable.lower() == 'voltage':
            msg = '==============================================\n'
            msg += '             branch voltages\n'
            msg += '==============================================\n'

            # check presence of branch voltages
            if self.vb is None:
                self.branch_voltage()

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
                    elif (k == 3) | (k == 5) | (k ==8):  # voltage sources or VCVS or CCVS
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif (k == 4) | (k == 6) | (k == 7):  # current sources or CCCS or VCCS
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += '----------------------------------------------\n'
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
                    elif (k == 3) | (k == 5) | (k ==8):  # voltage sources or VCVS or CCVS
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += '----------------------------------------------\n'
                    elif (k == 4) | (k == 6) | (k == 7):  # current sources or CCCS or VCCS
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += '----------------------------------------------\n'

        elif variable.lower() == 'current':
            msg = '==============================================\n'
            msg += '             branch currents\n'
            msg += '==============================================\n'

            # check presence of branch currents
            if self.ib is None:
                self.branch_current()

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
                    elif (k == 3) | (k == 5) | (k ==8):  # voltage sources or VCVS or CCVS
                        for h in index:
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif (k == 4) | (k == 6) | (k == 7):  # current sources or CCCS or VCCS
                        for h in index:
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += '----------------------------------------------\n'
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
                    elif (k == 3) | (k == 5) | (k ==8):  # voltage sources or VCVS or CCVS
                        for h in index:
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += '----------------------------------------------\n'
                    elif (k == 4) | (k == 6) | (k == 7):  # current sources or CCCS or VCCS
                        for h in index:
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += '----------------------------------------------\n'

        if variable.lower() == 'power':

            # check presence of branch powers
            if self.pb is None:
                self.branch_power()

            if self.analysis[0].lower() == '.op':
                unitsR = 'W'
                unitsL = 'W'
                unitsC = 'W'
                unitsV = 'W'
                unitsI = 'W'
            elif self.analysis[0].lower() == '.ac':
                unitsR = 'W'
                unitsL = 'var'
                unitsC = 'var'
                unitsV = 'VA'
                unitsI = 'VA'
            elif self.analysis[0].lower() == '.tran':
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
                    elif (k == 3) | (k == 5) | (k ==8):  # voltage sources or VCVS or CCVS
                        for h in index:
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsV, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif (k == 4) | (k == 6) | (k == 7):  # current sources or CCCS or VCCS
                        for h in index:
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsI, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
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
                    elif (k == 3) | (k == 5) | (k ==8):  # voltage sources or VCVS or CCVS
                        for h in index:
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsV)
                            msg += '----------------------------------------------\n'
                    elif (k == 4) | (k == 6) | (k == 7):  # current sources or CCCS or VCCS
                        for h in index:
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsI)
                            msg += '----------------------------------------------\n'

        elif variable.lower() == 'all':

            # check presence of branch powers
            if self.pb is None:
                self.branch_power()

            if self.analysis[0].lower() == '.op':
                unitsR = 'W'
                unitsL = 'W'
                unitsC = 'W'
                unitsV = 'W'
                unitsI = 'W'
            elif self.analysis[0].lower() == '.ac':
                unitsR = 'W'
                unitsL = 'var'
                unitsC = 'var'
                unitsV = 'VA'
                unitsI = 'VA'
            elif self.analysis[0].lower() == '.tran':
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
                    elif (k == 3) | (k == 5) | (k ==8):  # voltage sources or VCVS or CCVS
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsV, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'
                    elif (k == 4) | (k == 6) | (k == 7):  # current sources or CCCS or VCCS
                        for h in index:
                            msg += voltage_fmt_polar.format(self.names[h], np.abs(self.vb[h]), np.angle(self.vb[h], deg=True))
                            msg += current_fmt_polar.format(self.names[h], np.abs(self.ib[h]), np.angle(self.ib[h], deg=True))
                            msg += power_fmt_polar.format(self.names[h], np.abs(self.pb[h]), unitsI, np.angle(self.pb[h], deg=True))
                            msg += '----------------------------------------------\n'

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
                    elif (k == 3) | (k == 5) | (k ==8):  # voltage sources or VCVS or CCVS
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsV)
                            msg += '----------------------------------------------\n'
                    elif (k == 4) | (k == 6) | (k == 7):  # current sources or CCCS or VCCS
                        for h in index:
                            msg += voltage_fmt.format(self.names[h], self.vb[h])
                            msg += current_fmt.format(self.names[h], self.ib[h])
                            msg += power_fmt.format(self.names[h], self.pb[h], unitsI)
                            msg += '----------------------------------------------\n'

        if message:
            return msg
        else:
            print(msg)
            return None

    def plot(self, to_file=False, filename=None, dpi_value=150):

        # check if the analysis is '.tran'
        if self.analysis[0].lower() != '.tran':
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
                    remove_char = ('V(',')')

                    for char in remove_char:
                        variable = variable.replace(char, '')

                    v = self.get_voltage(variable)

                    if makesubplot:
                        plt.sca(axs[0])
                    plt.plot(self.t, v, label=legend_entries[k])


                elif variable[0] == 'I':
                    remove_char = ('I(', ')')

                    for char in remove_char:
                        variable = variable.replace(char, '')

                    i = self.get_current(variable)

                    if makesubplot:
                        plt.sca(axs[1])
                    plt.plot(self.t, i, label=legend_entries[k])



            if makesubplot:
                plt.sca(axs[0])
                plt.ylabel('voltage (V)', fontsize=16)
                plt.grid()
                plt.legend()
                plt.tight_layout()

                plt.sca(axs[1])
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

    def bode(self, decibel=False, to_file=False, filename=None, dpi_value=150):
        """
        bode print the bode diagram for the defined transfer functions (.ac multi-freq analysis only)

        :param
            * decibel: set to True to use decibel (default is False)
            * to_file: set to True to save on file (default is False)
            * filename: set the filename (default is bode_plot.png)
            * dpi_value: set the dpi (default is 150)

        :return:
        """

        # check if the analysis is '.ac'
        if self.analysis[0].lower() != '.ac':
            raise ValueError("bode() method available only for .ac analyses")
        # check if is a multi-freq analysis
        if np.isscalar(self.f):
            raise ValueError("bode() method useful for multi-frequency analyses. Use print() for single-frequency.")

        # convert the N tf in a Nx2 array
        tf_array = np.array(self.tf_cmd.upper().split()[1:]).reshape(-1,2)

        # compute tf
        hf = []
        for tf in tf_array:
            # output variable
            if 'V(' in tf[0]:
                out = self.get_voltage(tf[0].replace('V(','').replace(')',''))
            elif 'I(' in tf[0]:
                out = self.get_current(tf[0].replace('I(', '').replace(')', ''))
            # input variable
            if 'V(' in tf[1]:
                input = self.get_voltage(tf[1].replace('V(','').replace(')',''))
            elif 'I(' in tf[1]:
                input = self.get_current(tf[1].replace('I(', '').replace(')', ''))
            # tf
            H = out / input

            # plot
            import matplotlib.pyplot as plt
            fig, axs = plt.subplots(2, 1)
            plt.sca(axs[0])
            plt.title('tf: ' + tf[0] + '/' + tf[1], fontsize=14)
            if decibel:
                plt.semilogx(self.f, 20 * np.log10(np.abs(H)))
                plt.ylabel('magnitude (dB)', fontsize=14)
            else:
                plt.semilogx(self.f, np.abs(H))
                plt.ylabel('magnitude', fontsize=14)
            plt.grid()

            plt.sca(axs[1])
            plt.semilogx(self.f, np.angle(H) * 180 / np.pi)
            plt.xlabel('frequency (Hz)', fontsize=14)
            plt.ylabel('phase (deg)', fontsize=14)
            plt.grid()
            plt.tight_layout()

            hf.append(fig)

        # save to file
        if to_file:
            if filename is None:
                filename = 'bode_plot.png'
            else:
                if filename[-4:].lower() != '.png':
                    filename += '.png'

            if len(hf) == 0:
                pass
            elif len(hf) == 1:
                hf = hf[0]
                hf.savefig(filename, dpi=dpi_value)
            else:
                for k, fig in enumerate(hf):
                    fig.savefig(filename.replace('.png', '_' + str(k) + '.png'), dpi=dpi_value)
        else:
            if len(hf) == 1:
                hf = hf[0]

        return hf
