import numpy as np

def print_branch_voltage(net, polar=False):
    """
    'print_branch_voltage' print the branch voltages on terminal

    """

    # if necessary reorder
    if net.isort is None:
        net.reorder()

    print('==============================================')
    print('             branch voltages')
    print('==============================================')

    for k, index in enumerate(net.isort):
        if polar:
            if k == 0:  # resistors
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h], deg=True)))
                    print('----------------------------------------------')
            elif k == 1:  # voltage sources
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h], deg=True)))
                    print('----------------------------------------------')
            elif k == 2:  # voltage sources
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h], deg=True)))
                    print('----------------------------------------------')
            elif k == 3:  # voltage sources
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h], deg=True)))
                    print('----------------------------------------------')
            elif k == 4:  # current sources
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h], deg=True)))
                    print('----------------------------------------------')
        else:
            if k == 0:  # resistors
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('----------------------------------------------')
            elif k == 1:  # voltage sources
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('----------------------------------------------')
            elif k == 2:  # voltage sources
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('----------------------------------------------')
            elif k == 3:  # voltage sources
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('----------------------------------------------')
            elif k == 4:  # current sources
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('----------------------------------------------')


def print_branch_current(net, polar=False):
    """
    'print_branch_current' print the branch voltages on terminal

    """

    # if necessary reorder
    if net.isort is None:
        net.reorder()

    print('==============================================')
    print('             branch currents')
    print('==============================================')

    for k, index in enumerate(net.isort):
        if polar:
            if k == 0: # resistors
                for h in index:
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h], deg=True)))
                    print('----------------------------------------------')
            elif k == 1:  # inductors
                for h in index:
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h], deg=True)))
                    print('----------------------------------------------')
            elif k == 2:  # capacitors
                for h in index:
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h], deg=True)))
                    print('----------------------------------------------')
            elif k == 3:  # voltage sources
                for h in index:
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h], deg=True)))
                    print('----------------------------------------------')
            elif k == 4:  # current sources
                for h in index:
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h], deg=True)))
                    print('----------------------------------------------')
        else:
            if k == 0: # resistors
                for h in index:
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')
            elif k == 1:  # inductors
                for h in index:
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')
            elif k == 2:  # capacitors
                for h in index:
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')
            elif k == 3:  # voltage sources
                for h in index:
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')
            elif k == 4:  # current sources
                for h in index:
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')


def print_branch_quantity(net, polar=False):
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
        if polar:
            if k == 0:  # resistors
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h],deg=True)))
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h],deg=True)))
                    print('----------------------------------------------')
            elif k == 1:  # inductors
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h],deg=True)))
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h],deg=True)))
                    print('----------------------------------------------')
            elif k == 2:  # capacitors
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h],deg=True)))
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h],deg=True)))
                    print('----------------------------------------------')
            elif k == 3:  # voltage sources
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h],deg=True)))
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h],deg=True)))
                    print('----------------------------------------------')
            elif k == 4:  # current sources
                for h in index:
                    print('v({}) = {:10.4f} V < {:10.4f}°'.format(net.names[h], np.abs(net.vb[h]), np.angle(net.vb[h],deg=True)))
                    print('i({}) = {:10.4f} A < {:10.4f}°'.format(net.names[h], np.abs(net.ib[h]), np.angle(net.ib[h],deg=True)))
                    print('----------------------------------------------')

        else:
            if k == 0:  # resistors
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')
            elif k == 1:  # inductors
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')
            elif k == 2:  # capacitors
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')
            elif k == 3:  # voltage sources
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')
            elif k == 4:  # current sources
                for h in index:
                    print('v({}) = {:10.4f} V'.format(net.names[h], net.vb[h]))
                    print('i({}) = {:10.4f} A'.format(net.names[h], net.ib[h]))
                    print('----------------------------------------------')

