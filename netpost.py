def branch_voltage(net):
    """
    "branch_voltage"  computes the branch voltages

    :return:
        * net.vb
    """
    # check if the incedence matrix is available
    if net.A is None:
        net.incidence_matrix()

    # check if the solituon is available
    if net.x is None:
        printf("No solution available")
        return None

    # branch voltages
    net.vb = net.A.transpose() * net.x[:net.node_num]


def branch_current(net):
    """
    "branch_current"  computes the branch currents

    :return:
        * net.net.ib
    """
    # check is branch voltages are available
    if net.vb is None:
        branch_voltage(net)

    net.ib = []
    cnt_v = 0
    for name, val, voltage in zip(net.names, net.values, net.vb):
        if name[0].upper() == 'R':
            net.ib.append(voltage/val)
        elif name[0].upper() == 'V':
            net.ib.append(net.x[net.node_num + cnt_v])
            cnt_v += 1
        elif name[0].upper() == 'I':
            net.ib.append(val)
