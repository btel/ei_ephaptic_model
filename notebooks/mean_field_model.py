#!/usr/bin/env python
#coding=utf-8

import numpy as np
from scipy import special
from scipy import optimize

#SI base units
s = 1
kg = 1
m = 1
A = 1

#derived units
S = s**3*A**2/(kg*m**2)
V = kg*m**2*s**-3*A**-1
F = s**4 * A**2 * m**-2 * kg ** -1
Hz = 1/s

#with prefixes
nS = 1e-9 * S
uS = 1e-6 * S
mV = 1e-3 * V
pF = 1e-12 * F
ms = 1e-3 * s
nA = 1e-9 * A
pA = 1e-12 * A

# constants
dt = 10 * ms
E_e = 0*mV
E_i = -75*mV
E_l = -70*mV
g_l = 1./60*uS
C = 250*pF
v_reset = -60*mV
threshold = -50*mV
threshold_inh = -53 * mV
tau_ref = 2*ms

tau_e = 0.2 * ms  # width of excitatory PSC (ms)
tau_i = 2 * ms    # width of inhibitory PSC (ms)
B_e = 7.1 * nS    # peak excitatory conductance (nS)
B_i = 3.7 * nS    # peak inhibitory conductance (nS)
fr_e = 9655 * Hz    # total firing rate of excitatory population (Hz)
fr_i = 4473 * Hz     # total firing rate of excitatory population (Hz)
f_ext = 5000 * Hz # etracellular_firing_rate

n_e = 350. # number of excitatory neurons
n_i = n_e / 4. # number of inhibitory neurons

def calc_membrane_stats(fr_e=fr_e, fr_i=fr_i):

    mu_ge = fr_e * B_e * tau_e * np.exp(1)
    mu_gi = fr_i * B_i * tau_i * np.exp(1)

    gtot = g_l + mu_ge + mu_gi

    mu_u = (E_l * g_l + E_e * mu_ge + E_i * mu_gi) / gtot

    tau_eff = C / gtot

    epsp_int = (E_e - mu_u) * B_e * tau_e * np.exp(1) * tau_eff / C
    ipsp_int = (E_i - mu_u) * (B_i * tau_i * np.exp(1) * tau_eff / C)                          
    epsp_sq =  epsp_int ** 2 * (2 * tau_eff + tau_e) /(4 * (tau_eff + tau_e)**2)
    ipsp_sq =  ipsp_int ** 2 * (2 * tau_eff + tau_i) /(4 * (tau_eff + tau_i)**2)
    sigma_sq_u = fr_e * epsp_sq + fr_i * ipsp_sq
    
    return gtot, mu_u, tau_eff, sigma_sq_u

def kuhn_transfer_function(tau_eff, threshold, mu_u, sigma_sq_u):
    return 1. /(2 * tau_eff) * special.erfc((threshold - mu_u) / (np.sqrt(2) * np.sqrt(sigma_sq_u)))

def derivative_kuhn_transfer_function(tau_eff, threshold, mu_u, sigma_sq_u):
    return 1. /(2 * tau_eff * np.sqrt(2 * sigma_sq_u)) * (2/ np.sqrt(np.pi)) * np.exp(-(threshold - mu_u) ** 2 / (2 * sigma_sq_u))

def calc_output_rate_inh_exc(fr_e, fr_i, n_e, n_i, threshold_exc=-50 * mV, threshold_inh=-50 * mV):
    gtot_exc, mu_u_exc, tau_eff_exc, sigma_sq_u_exc = calc_membrane_stats(fr_e=fr_e, fr_i=fr_i)
    gtot_inh, mu_u_inh, tau_eff_inh, sigma_sq_u_inh = calc_membrane_stats(fr_e=fr_e, fr_i=fr_i)
    
    rexc = kuhn_transfer_function(tau_eff_exc, threshold_exc, mu_u_exc, sigma_sq_u_exc)
    rinh = kuhn_transfer_function(tau_eff_inh, threshold_inh, mu_u_inh, sigma_sq_u_inh)

    return rexc * n_e, rinh * n_i

