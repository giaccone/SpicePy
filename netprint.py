def print_branch_voltage(net):
    """
    'print_branch_voltage' print the branch voltages on terminal

    """

    # if necessary reorder
    if net.isort is None:
        net.reorder()

    print('=====================')
    print('   branch voltages   ')
    print('=====================')

    for k, index in enumerate(net.isort):
        if k == 0:  # resistors
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
        elif k == 1:  # voltage sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
        elif k == 2:  # voltage sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
        elif k == 3:  # voltage sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
        elif k == 4:  # current sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                print('=====================')


def print_branch_current(net):
    """
    'print_branch_current' print the branch voltages on terminal

    """

    # if necessary reorder
    if net.isort is None:
        net.reorder()

    print('=====================')
    print('   branch currents   ')
    print('=====================')

    for k, index in enumerate(net.isort):
        if k == 0: # resistors
            for h in index:
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
        elif k == 1:  # inductors
            for h in index:
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
        elif k == 2:  # capacitors
            for h in index:
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
        elif k == 3:  # voltage sources
            for h in index:
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
        elif k == 4:  # current sources
            for h in index:
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                print('=====================')


def print_branch_quantity(net):
    """
    'print_branch_quantity' print the branch voltages and currents on terminal

    """

    # if necessary reorder
    if net.isort is None:
        net.reorder()

    print('==============================================')
    print('               branch quantities              ')
    print('==============================================')

    for k, index in enumerate(net.isort):
        if k == 0:  # resistors
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]), end='    ')
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
        elif k == 1:  # inductors
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]), end='    ')
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
        elif k == 2:  # capacitors
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]), end='    ')
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
        elif k == 3:  # voltage sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]), end='    ')
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
        elif k == 4:  # current sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]), end='    ')
                print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                print('==============================================')
