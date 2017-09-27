# import SciPy modules
import netlist as ntl
from netsolve import net_solve
import matplotlib.pyplot as plt
plt.ion()

# read netlist
net = ntl.Network('tran_network3.net')

# compute the circuit solution
net_solve(net)

# plot results
net.plot()