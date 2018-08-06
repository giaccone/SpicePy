# import SpicePy modules
import netlist as ntl
from netsolve import net_solve

# read netlist
net = ntl.Network('VCCS_and_CCVS.net')

# compute the circuit solution
net_solve(net)

# compute branch quantities
net.branch_voltage()
net.branch_current()
net.branch_power()

# print results
net.print()