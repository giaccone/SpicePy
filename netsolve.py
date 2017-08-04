from scipy.sparse.linalg import spsolve

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
