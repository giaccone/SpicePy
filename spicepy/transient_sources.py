# ===========================================================
# About the code
# ===========================================================
# This code is part of the project 'SpicePy'.
# See README.md for more details
#
# Licensed under the MIT license (see LICENCE)
# Copyright (c) 2017 Luca Giaccone (luca.giaccone@polito.it)
# ===========================================================

# ==================
# imported modules
# ==================
import numpy as np
from scipy.interpolate import interp1d


def pwl(pairs, t):
    """
    PWL describes a piecewise linear waveform

    :param pairs: a 2D list. Each rom is a pair [time, value]
    :param t: array with times where the function has to be evaluated
    :return: the function values at times defined in t
    """

    # convert pairs to np.array is necessary
    if isinstance(pairs, list):
        pairs = np.array(pairs)

    # check pairs format
    if pairs.ndim == 1:
        pairs.shape = (-1,2)

    # get pairs
    x = []
    y = []
    for xk, yk in pairs:
        x.append(xk)
        y.append(yk)

    # create linear interpolator
    fun = interp1d(x, y, bounds_error=False, fill_value=(y[0], y[-1]), assume_sorted=True)
    # interpolation
    out = fun(t)

    return out


def pulse(V1, V2, Td=0, Tr=None, Tf=None, Pw=None, Period=None, t=None):
    """
    PULSE provides a pulse according to this scheme:
        * t < Td ==> V1
        * Td < t < Td + Tr ==>  linear ramp from V1 to V2
        * Td + Tr < t < Td + Tr + Pw ==>  V2
        * Td + Tr + Pw < t < Td + Tr + Pw + Tf ==>  linear ramp from V2 to V1
        * Td + Tr + Pw + Tf < t < (Period - Td + Tr + Pw + Tf) ==> V2

    :param V1: low value
    :param V2: high value
    :param Td: delay time (s)
    :param Tr: rise time (s)
    :param Tf: fall time (s)
    :param Pw: pulse width (s)
    :param Period: period (s)
    :param t: array with times where the function has to be evaluated
    :return: the function values at times defined in t
    """

    # check presence of time input
    if t is None:
        raise TypeError('Missing time input')

    # check if t is scalar
    if isinstance(t, (int, float)):
        # convert to numpy array
        t = np.array([t])

        # check presence of Tr, and Tf (raise error)
        if (Tr == None) or (Tf == None):
            raise TypeError('Missing Tr and/or Tf')
    else:

        # check presence of Tr, and Tf (assign default)
        time_step = t[1] - t[0]
        if Tr == None:
            Tr = np.copy(time_step)
        if Tf == None:
            Tf = np.copy(time_step)

    # check presence of Pw and Period fot t > Td
    if t[-1] > Td:

        if Pw == None:
            Pw = t[-1] - Td

        if Period == None:
            Period = t[-1] - Td

    # get number of finite periods
    N = np.ceil((t[-1] - Td) / Period) + 1


    # initialize pairs (with or without delay)
    if Td != 0:
        x = [0, Td]
        y = [V1, V1]
        time = np.copy(Td)
    else:
        x = [0]
        y = [V1]
        time = 0

    # update pairs for N periods
    for k in range(int(N)):
        x.append(time + Tr)
        y.append(V2)
        if Pw != 0:
            x.append(time + Tr + Pw)
            y.append(V2)
        x.append(time + Tr + Pw + Tf)
        y.append(V1)
        x.append(time + Period)
        y.append(V1)
        time = time + Period


    # create linear interpolator
    fun = interp1d(x, y,bounds_error=False, fill_value=(y[0], y[-1]), assume_sorted=True)
    # interpolation
    out = fun(t)

    # if input is scalar convert out to scalar too
    if out.size == 1:
        out = np.asscalar(out)

    return out


def sin(Vo, Va, Freq=None, Td=0, Df=0, Phase=0, t=None):
    """
    SIN provides a damped sinusoidal waveform in the form
    Vo + Va * np.sin(2 * np.pi * Freq * t + Phase * (np.pi / 180))

    The waveforms is:
        * t < Td ==> Vo + Va * np.sin(Phase * (np.pi / 180))
        * t > Td ==> Vo + Va * np.sin(2 * np.pi * Freq * (t - Td) + Phase * (np.pi / 180)) * np.exp(-(t - Td) * Df)

    :param Vo: offset
    :param Va: amplitude (peak) of the waveform
    :param Freq: frequency (Hz)
    :param Td: delay time (s)
    :param Df: damping factor (1/s)
    :param Phase: voltage phase (deg)
    :param t: array with times where the function has to be evaluated
    :return: the function values at times defined in t
    """

    # check presence of time array
    if t is None:
        raise TypeError('Missing time array')

    # check if t is scalar
    if isinstance(t, (int, float)):
        t = np.array([t])

    # check presence of Freq
    if Freq is None:
        Freq = 1 / t[-1]

    out = np.zeros_like(t)
    out[t <= Td] = Vo + Va * np.sin(Phase * (np.pi / 180))
    out[t > Td] = Vo + Va * np.sin(2 * np.pi * Freq * (t[t > Td] - Td) + Phase * (np.pi / 180)) * np.exp(-(t[t > Td] - Td) * Df)

    # if input is scalar convert out to scalar too
    if out.size == 1:
        out = np.asscalar(out)

    return out


def exp(V1, V2, Td1=0, tau1=None, Td2=None, tau2=None, t=None):
    """
    EXP provides an exponential waveform according to:
        * t < Td1 ==> V1
        * Td1 < t < Td2 ==> from V1 to V2 with exponential decay (time constant tau1)
        (if Td2 < 5*T1 the voltage reaches V(Td2) < V2 )
        * t > Td2 ==> from V2 to V1 with exponential decay (time constant tau2)
        (if V(Td2) < V2 the exponential decay is from V(Td2) to V1)

    :param V1: initial value
    :param V2: peak
    :param Td1: rise time delay (s)
    :param tau1: rise time constant (s)
    :param Td2: fall time delay (s)
    :param tau2: fall time constant (s)
    :param t: array with times where the function has to be evaluated
    :return: the function values at times defined in t
    """

    # check presence of time array
    if t is None:
        raise TypeError('Missing time array')

    if isinstance(t, (int, float)):
        t = np.array([t])

        # check presence of T1, Td2 and T2
        if (tau1 is None) or (Td2 is None) or (tau2 is None):
            raise TypeError('Missing tau1 and/or Td2 and/or tau2')
    else:

        # check presence of T1, Td2 and T2
        time_step = t[1] - t[0]
        if tau1 is None:
            tau1 = np.copy(time_step)
        if Td2 is None:
            Td2 = t[t > Td1][1]
        if tau2 is None:
            tau2 = np.copy(time_step)


    out = np.zeros_like(t)

    out[t <= Td1] = V1

    id = (t > Td1) & (t <= Td2)
    out[id] = V1 + (V2 - V1) * (1 - np.exp(-(t[id] - t[id][0]) / tau1))

    Vend = out[id][-1]
    id = (t > Td2)
    out[id] = V1 + (Vend - V1) * np.exp(-(t[id] - t[id][0]) / tau2)

    # if input is scalar convert out to scalar too
    if out.size == 1:
        out = np.asscalar(out)

    return out