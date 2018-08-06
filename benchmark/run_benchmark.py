import netlist as ntl
from netsolve import net_solve
import numpy as np


# check rms
def check(rms, toll):
    if rms < toll:
        print('OK')
    else:
        print('*FAILED*')


print('Start benchmark:')

# DC network
print(' * Testing DC network ...', end=' ',flush=True)
net = ntl.Network('../demo/dc_network.net')
net_solve(net)
ref = np.loadtxt('./dc_network/solution.txt')
rms = np.sqrt(np.sum((ref - net.x) ** 2)) / ref.size
check(rms, 1e-10)

# OP network
print(' * Testing OP network ...', end=' ',flush=True)
net = ntl.Network('../demo/op_network.net')
net_solve(net)
ref = np.loadtxt('./op_network/solution.txt')
rms = np.sqrt(np.sum((ref - net.x) ** 2)) / ref.size
check(rms, 1e-5)

# AC network (single frequency)
print(' * Testing AC network (single frequency) ...', end=' ',flush=True)
net = ntl.Network('../demo/ac_single_frequency.net')
net_solve(net)
ref = np.loadtxt('./ac_single_frequency/solution.txt')
ref = ref[:,0] * np.exp(ref[:,1] * 1j * np.pi / 180)
rms_r = np.sqrt(np.sum((ref.real - net.x.real) ** 2))
rms_i = np.sqrt(np.sum((ref.imag - net.x.imag) ** 2))
rms = np.sqrt(rms_r ** 2 + rms_i ** 2) / ref.size
check(rms, 1e-5)

# AC network (low-high pass filter)
print(' * Testing AC network (low-high pass filter) ...', end=' ',flush=True)
net = ntl.Network('../demo/ac_low_high_pass_filter.net')
net_solve(net)
ref = np.loadtxt('./ac_low_high_pass_filter/solution.txt')
ref1 = 10 ** (ref[:,1]/20) * np.exp(ref[:,2] * 1j * np.pi / 180)
ref2 = 10 ** (ref[:,3]/20) * np.exp(ref[:,4] * 1j * np.pi / 180)
vc1 = net.get_voltage('C1')
vr1 = net.get_voltage('R1')
rms_r1 = np.sqrt(np.sum((ref1.real - vc1.real) ** 2))
rms_i1 = np.sqrt(np.sum((ref1.imag - vc1.imag) ** 2))
rms_r2 = np.sqrt(np.sum((ref2.real - vr1.real) ** 2))
rms_i2 = np.sqrt(np.sum((ref2.imag - vr1.imag) ** 2))
rms = np.sqrt(rms_r1 ** 2 + rms_i1 ** 2 + rms_r2 ** 2 + rms_i2 ** 2) / ref.shape[0]
check(rms, 1e-5)

# AC network (band pass filter)
print(' * Testing AC network (band pass filter) ...', end=' ',flush=True)
net = ntl.Network('../demo/ac_band_pass_filter.net')
net_solve(net)
ref = np.loadtxt('./ac_band_pass_filter/solution.txt')
ref = 10 ** (ref[:,1]/20) * np.exp(ref[:,2] * 1j * np.pi / 180)
vr1 = net.get_voltage('R1')
rms_r = np.sqrt(np.sum((ref.real - vr1.real) ** 2))
rms_i = np.sqrt(np.sum((ref.imag - vr1.imag) ** 2))
rms = np.sqrt(rms_r1 ** 2 + rms_i1 ** 2) / ref.shape[0]
check(rms, 1e-5)

# tran network1
print(' * Testing TRAN network1 ...', end=' ',flush=True)
net = ntl.Network('../demo/tran_network1.net')
net_solve(net)
val = net.get_current('L1')
ref = np.loadtxt('./tran_network1/solution.txt', skiprows=1)
ref = np.interp(net.t, ref[:,0],ref[:,1])
rms = np.sqrt(np.sum((ref - val) ** 2)) / ref.size / np.abs(ref).max()
check(rms, 1e-5)

# tran network2
print(' * Testing TRAN network2 ...', end=' ',flush=True)
net = ntl.Network('../demo/tran_network2.net')
net_solve(net)
val = net.get_current('L1')
ref = np.loadtxt('./tran_network2/solution.txt', skiprows=1)
ref = np.interp(net.t, ref[:,0],ref[:,1])
rms = np.sqrt(np.sum((ref[1:] - val[1:]) ** 2)) / ref[1:].size / np.abs(ref).max()
check(rms, 1e-5)

# tran network3
print(' * Testing TRAN network3 ...', end=' ',flush=True)
net = ntl.Network('../demo/tran_network3.net')
net_solve(net)
val1 = net.get_voltage('C1')
val2 = net.get_current('L1')
ref = np.loadtxt('./tran_network3/solution.txt', skiprows=1)
ref1 = np.interp(net.t, ref[:,0], ref[:,1])
ref2 = np.interp(net.t, ref[:,0], ref[:,2])
rms1 = np.sqrt(np.sum((ref1 - val1) ** 2)) / ref1.size / np.abs(ref1).max()
rms2 = np.sqrt(np.sum((ref2 - val2) ** 2)) / ref2.size / np.abs(ref2).max()
rms = np.sqrt(rms1 ** 2 + rms2 ** 2)
check(rms, 1e-4)

# tran network4
print(' * Testing TRAN network4 ...', end=' ',flush=True)
net = ntl.Network('../demo/tran_network4.net')
net_solve(net)
val1 = net.get_voltage('C1')
val2 = net.get_current('L1')
ref = np.loadtxt('./tran_network4/solution.txt', skiprows=1)
ref1 = np.interp(net.t[1:], ref[:,0], ref[:,1])
ref2 = np.interp(net.t[1:], ref[:,0], ref[:,2])
rms1 = np.sqrt(np.sum((ref1 - val1[1:]) ** 2)) / ref1.size / np.abs(ref1).max()
rms2 = np.sqrt(np.sum((ref2 - val2[1:]) ** 2)) / ref2.size / np.abs(ref2).max()
rms = np.sqrt(rms1 ** 2 + rms2 ** 2)
check(rms, 1e-4)

# tran network5
print(' * Testing TRAN network5 ...', end=' ',flush=True)
net = ntl.Network('../demo/tran_network5.net')
net_solve(net)
val = net.get_current('L1')
ref = np.loadtxt('./tran_network5/solution.txt', skiprows=1)
ref = np.interp(net.t, ref[:,0],ref[:,1])
rms = np.sqrt(np.sum((ref - val) ** 2)) / ref.size / np.abs(ref).max()
check(rms, 1e-5)

# tran network7
print(' * Testing TRAN network6 ...', end=' ',flush=True)
net = ntl.Network('../demo/tran_network6.net')
net_solve(net)
val = net.get_voltage('C1')
ref = np.loadtxt('./tran_network6/solution.txt', skiprows=1)
ref = np.interp(net.t, ref[:,0],ref[:,1])
rms = np.sqrt(np.sum((ref - val) ** 2)) / ref.size / np.abs(ref).max()
check(rms, 1e-5)

# tran network6
print(' * Testing TRAN network7 ...', end=' ',flush=True)
net = ntl.Network('../demo/tran_network7.net')
net_solve(net)
val = net.get_current('L1')
ref = np.loadtxt('./tran_network7/solution.txt', skiprows=1)
ref = np.interp(net.t, ref[:,0],ref[:,1])
rms = np.sqrt(np.sum((ref - val) ** 2)) / ref.size / np.abs(ref).max()
check(rms, 1e-5)

# VCVS and CCCS
print(' * Testing VCVS & CCCS ...', end=' ',flush=True)
net = ntl.Network('../demo/VCVS_and_CCCS.net')
net_solve(net)
ref = np.loadtxt('./VCVS_and_CCCS/solution.txt')
rms = np.sqrt(np.sum((ref - net.x) ** 2)) / ref.size
check(rms, 1e-9)
