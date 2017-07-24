def print_branch_voltage(network):
    """
    'print_branch_voltage' print the branch voltages on terminal

    """

    # if necessary reorder
    if not hasattr(network, 'isort'):
        network.reorder()

    print('=====================')
    print('   branch voltages   ')
    print('=====================')

    for k, index in enumerate(network.isort):
        if k == 0: # resistors
            for h in index:
                print('v({}) = {:10.4f} V'.format(network.names[h], network.vb[h]))
        elif k == 1: # voltage sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(network.names[h], network.vb[h]))
        elif k == 2: # current sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(network.names[h], network.vb[h]))
                print('=====================')


def print_branch_current(network):
    """
    'print_branch_current' print the branch voltages on terminal

    """

    # if necessary reorder
    if not hasattr(network, 'isort'):
        network.reorder()

    print('=====================')
    print('   branch currents   ')
    print('=====================')

    for k, index in enumerate(network.isort):
        if k == 0: # resistors
            for h in index:
                print('i({}) = {:10.4f} A'.format(network.names[h], network.ib[h]))
        elif k == 1: # voltage sources
            for h in index:
                print('i({}) = {:10.4f} A'.format(network.names[h], network.ib[h]))
        elif k == 2: # current sources
            for h in index:
                print('i({}) = {:10.4f} A'.format(network.names[h], network.ib[h]))
                print('=====================')


def print_branch_quantity(network):
    """
    'print_branch_quantity' print the branch voltages and currents on terminal

    """

    # if necessary reorder
    if not hasattr(network, 'isort'):
        network.reorder()

    print('==============================================')
    print('               branch quantities              ')
    print('==============================================')

    for k, index in enumerate(network.isort):
        if k == 0: # resistors
            for h in index:
                print('v({}) = {:10.4f} V'.format(network.names[h], network.vb[h]), end='    ')
                print('i({}) = {:10.4f} A'.format(network.names[h], network.ib[h]))
        elif k == 1: # voltage sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(network.names[h], network.vb[h]), end='    ')
                print('i({}) = {:10.4f} A'.format(network.names[h], network.ib[h]))
        elif k == 2: # current sources
            for h in index:
                print('v({}) = {:10.4f} V'.format(network.names[h], network.vb[h]), end='    ')
                print('i({}) = {:10.4f} A'.format(network.names[h], network.ib[h]))
                print('==============================================')
