# import SciPy modules
import spicepy.netlist as ntl
from spicepy.netsolve import net_solve
import matplotlib.pyplot as plt
plt.ion()

# read netlist
net = ntl.Network('ac_low_high_pass_filter.net')

# compute the circuit solution
net_solve(net)

# plot bode diagrams
net.bode()