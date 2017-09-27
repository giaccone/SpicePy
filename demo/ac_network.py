# import SciPy modules
import netlist as ntl
from netsolve import net_solve

# read netlist
net = ntl.Network('ac_network.net')

# compute the circuit solution
net_solve(net)

# compute branch quantities
net.branch_voltage()
net.branch_current()
net.branch_power()

# print results
net.print()

# print results in polar format
net.print(polar=True)