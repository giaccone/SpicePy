# import SpicePy modules
import spicepy.netlist as ntl
from spicepy.netsolve import net_solve

# read netlist
net = ntl.Network('op_network.net')

# compute the circuit solution
net_solve(net)

# compute branch quantities
net.branch_voltage()
net.branch_current()
net.branch_power()

# print results
net.print()