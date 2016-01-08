#====== CAMB Power Spectrum =============
"""
List of functions:
    EoS: plots the equation of state for quintessence DE
    hubble: plot the dimensionless hubble of L/q CDM
    density: plot the density parameter of L/q CDM
    reltchange_hubble: plot the relative difference of q/L CDM
"""
from __future__ import division
from numpy import *
from matplotlib import *
from matplotlib import cm
from pylab import *
import var_def as var
import constants as cc
from model import model_check
from itertools import cycle
from mpl_toolkits.mplot3d.axes3d import Axes3D
#from matplotlib.collections import PolyCollection
from matplotlib.collections import *
#from matplotlib.colors import colorConverter

import os.path

from params_ini import power_primordial

#__________ Plotting setup _______
rc("font", size="15"); rc("axes",labelsize="15")
lines = ["-","--","-.",":"]; linecycler = cycle(lines)

#======================================================: NOTES !!
#--: Quintessence !!
def w_x(model):
    w_x = var.back_vars(model)[3]
    semilogx(cc.a, w_x, next(linecycler), linewidth = 2, label = model)
    legend(loc = 'best')
    xlabel('a')
    ylabel('$w_x(a)$')

#def quintess_vars(model):
#    phi = var.quin_vars(model)[0]
#    Lambda = var.quin_vars(model)[1]
#    print phi[-1], Lambda[-1]

#--: Background !!
def omega_m(model):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma
    
    omega_m = var.back_vars(model)[0]; omega_x = var.back_vars(model)[1]
    semilogx(cc.a, omega_m, next(linecycler), linewidth = 2, label = label1 + ' :: $\Omega_m$' )
    semilogx(cc.a, omega_x, next(linecycler),linewidth = 2, label =  label1 + ' :: $\Omega_x$')
    
    legend(loc = 'best', prop = {'size':9}) #, ncol = 3
    xlabel('a')
    xlim((4 * 10**-2, 1.0))
    ylim((-0.05 , 1.05))
    #print model, '::'
    #print '$\Omega_m: $', omega_m[0][-1], '&&', '$\Omega_x: $',omega_x[0][-1]


def h(model):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma

    h = var.back_vars(model)[2]
    loglog(cc.a, h, next(linecycler), linewidth = 2, label = label1)
    legend(loc = 'best', prop = {'size':9})
    xlim((10**-1, 1.0))
    ylim(( 1.0, 0.2 * 10**2))
    xlabel('a')
    ylabel('h(a)')
    #print model, '::'
    #print 'h: ',h[0][-1]


#perturb_vars = Delta_m, Delta_x, u_m, u_x, Phi
#perturb_vars_ini = Delta_m_ini, Delta_x_ini, u_m_ini, u_x_ini, Phi_ini
#--: Perturbation !!
def Delta_m(model):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma

    Delta_m = var.perturb_vars(model)[0]
    semilogx(cc.k, Delta_m, next(linecycler), linewidth = 2, label = label1)
    legend(loc = 'best', prop = {'size':9})
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$\Delta_m$(k, a = 1)')

def Delta_m_ini(model, **cosmo):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma

    Delta_m_ini = var.perturb_vars_ini(model, **cosmo)[0]
    conv = 1/(cc.CH_0**(3/2))
    semilogx(cc.k, Delta_m_ini, next(linecycler), linewidth = 2, label = label1)
    legend(loc = 'best', prop = {'size':9})
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$\Delta_m$(k, a = 1)')

def u_m(model):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma

    u_m = var.perturb_vars(model)[2]
    semilogx(cc.k, u_m, next(linecycler), linewidth = 2, label = label1)
    legend(loc = 'best', prop = {'size':9})
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$u_m$(k, a = 1)')

def u_m_ini(model, **cosmo):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma

    u_m_ini = var.perturb_vars_ini(model, **cosmo)[2]
    conv = 1/(cc.CH_0**(3/2))
    semilogx(cc.k, u_m_ini, next(linecycler), linewidth = 2, label = label1)
    legend(loc = 'best', prop = {'size':9})
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$u_m$(k, a = 1)')

def Phi(model):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma

    Phi = var.perturb_vars(model)[4]
    semilogx(cc.k, Phi, next(linecycler), linewidth = 2, label = label1)
    legend(loc = 'best', prop = {'size':9})
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$\Phi$(k, a = 1)')

def Phi_ini(model, **cosmo):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma

    Phi_ini = var.perturb_vars_ini(model, **cosmo)[4]
    conv = 1/(cc.CH_0**(3/2))
    semilogx(cc.k, Phi_ini, next(linecycler), linewidth = 2, label = label1)
    legend(loc = 'best', prop = {'size':9})
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$\Phi$(k, a = 1)')

#--: Growth Rate functions !!
def g_Delta_m(model, BP, **cosmo):
    g_Delta_m = var.grwoth_factor(model, BP,**cosmo)[0]
    loglog(cc.k, g_Delta_m, next(linecycler), linewidth = 2, label = model)
    legend(loc = 'best', prop = {'size':9})
    xlim((0.3 * 10**-3, 0.3))
    xlabel('k')
    ylabel('$g_m$(k, a = 1)')


#--: Primordial power spectrum !!
def P_Phi_p(fv, **cosmo):
    P_Phi_p = power_primordial(fv, **cosmo)
    loglog(cc.k, P_Phi_p, next(linecycler), linewidth = 2, label = '$k^{-2}$-factor: %s' %fv)
    legend(loc = 'best', prop = {'size':9})
    axvline(x = 0.17 * 10**-2, color='r', ls = '--') # vertical line at k_eq
    axvline(x = 0.14 * 10**-1, color='b', ls = '-.') # vertical line at k_eq
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$P^p_{Phi}$(k)')

#--: Linear matter Power Spectra !!
def P_Phi(model, fv, **cosmo):
    P_Phi_p = power_primordial(fv, **cosmo);  g_Delta_m = var.grwoth_factor(model,**cosmo)[0]
    g_Phi = var.grwoth_factor(model,**cosmo)[1]
    P_Delta = (g_Delta_m * g_Phi)**2 * P_Phi_p
    loglog(cc.k, P_Delta, next(linecycler), linewidth = 2, label = '$k^{-2}$-factor: %s' %fv)
    legend(loc = 'best', prop = {'size':9})
    axvline(x = 0.17 * 10**-2, color='r', ls = '--') # vertical line at k_eq
    axvline(x = 0.14 * 10**-1, color='b', ls = '-.') # vertical line at k_eq
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$P_m$(k, a = 1)')

#--: Linear matter Power Spectra !!
def P_Delta_m(model, lcdm_model, **cosmo):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    label1 = model_name + 'z = %s,' %z + ' $w_x$ = %s,' %w_x + ' $c^2_x$ = %s,' %c2_x  + ' $\Gamma$ =  %s' %gamma

    P_Delta_m = var.power_spectra_linear(model, lcdm_model, **cosmo)
    loglog(cc.k, P_Delta_m, next(linecycler), linewidth = 2, label = label1)
    legend(loc = 'best', prop = {'size':9})
    axvline(x = 0.17 * 10**-2, color='r', ls = '--') # vertical line at k_eq
    axvline(x = 0.14 * 10**-1, color='b', ls = '-.') # vertical line at k_eq
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$P_m$(k, a = 1)')

def diff_P_Delta_m(BP, **cosmo):
    diff_P_Delta_m_0 = (var.power_spectra_linear("wcdm_%s" %cc.w_x, BP, **cosmo) - var.power_spectra_linear("wcdm_%s" %cc.w_x, BP, **cosmo))/var.power_spectra_linear("wcdm_%s" %cc.w_x, BP, **cosmo)
    diff_P_Delta_m = (var.power_spectra_linear("gwcdm_%s" %cc.gamma, BP, **cosmo) - var.power_spectra_linear("wcdm_%s" %cc.w_x, BP, **cosmo))/var.power_spectra_linear("wcdm_%s" %cc.w_x, BP, **cosmo)
    #diff_P_Delta_m = (var.power_spectra_linear("gwcdm_%s" %cc.gamma_0, BP, **cosmo) - var.power_spectra_linear("lcdm", BP, **cosmo))/var.power_spectra_linear("lcdm", BP, **cosmo)
    semilogx(cc.k, diff_P_Delta_m_0 * 100, next(linecycler), linewidth = 2)
    semilogx(cc.k, diff_P_Delta_m * 100, next(linecycler), linewidth = 2, label = "$\Gamma$: %s" %cc.gamma)
    legend(loc = 'best')
    axvline(x = 0.17 * 10**-2, color='r', ls = '--') # vertical line at k_eq
    axvline(x = 0.14 * 10**-1, color='b', ls = '-.') # vertical line at k_eq
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$P^{\Gamma}_m$(k, a = 1) - $P^{\Lambda}_m$(k, a = 1)/$P^{\Lambda}_m$(k, a = 1) %')

#--: Bias(linear & non-gaussian) !!
def bias(model, f_NL, **cosmo):
    bias = var.galaxy_bias(model, f_NL, **cosmo)[0]
    loglog(cc.k, bias, next(linecycler), linewidth = 2, label = model+ "-$f_{NL}$_%s" %f_NL)
    legend(loc = 'best')
    xlim((0.3 * 10**-3, 0.3* 10**-2))
    xlabel('k')
    ylabel('$b_g$(k, a = 1)')

def delta_bias(model, f_NL, **cosmo):
    delta_bias = var.galaxy_bias(model, f_NL, **cosmo )[1]
    semilogx(cc.k, delta_bias, next(linecycler), linewidth = 2, label = model + "-$f_{NL}$_%s" %f_NL)
    legend(loc = 'best')
    xlim((0.3 * 10**-3, 0.3* 10**-2))
    xlabel('k')
    ylabel('$\Delta b_g$(k, a = 1)')

#--:  Galaxy Power Specra !!
def Pg_Delta_m(model, lcdm_model, f_NL, **cosmo):
    Pg_Delta_m = var.galaxy_power(model, lcdm_model, f_NL, **cosmo)
    loglog(cc.k, Pg_Delta_m, next(linecycler), linewidth = 2, label = model + "_%s" %f_NL)
    legend(loc = 'best')
    axvline(x = 0.33 * 10**-3, color='r', ls = '--') # vertical line at k = H_0
    axvline(x = 0.14 * 10**-1, color='b', ls = '-.') # vertical line at k_eq
    xlim((cc.k_min, cc.k_max))
    xlabel('k')
    ylabel('$Pg_m$(k, a = 1)')

def gfnl():
    gamma = cc.gamma;
    for j in range(len(cc.w_x)):
        for i in range(len(cc.z)):
            fNL = var.f_NL_eff(cc.z[i], cc.w_x[j]); alpha = 1.21
            plot(fNL, gamma, '-.o', linewidth = 2, label = "DXT - z: %s, $w_x$: %s" %(cc.z[i], cc.w_x[j]))
            #plot(fNL, gamma, next(linecycler), linewidth = 2, label = "DXT - z: %s, $w_x$: %s" %(cc.z[i], cc.w_x[j]))
    #plot(gamma, alpha * gamma, '--', linewidth = 2, label = '$f^{eff}_{NL}$ = 1.2 $\Gamma/H_0$')
    legend(loc = 'best', prop = {'size':7})
    xlabel('$f_{NL}$')
    ylabel('$\Gamma/H_0$')


#--: Angular Power (Linear) !!
def AP_m(model, BP, **cosmo):
    AP_Delta_m = var.angular_power(model, BP, **cosmo)
    semilogy(cc.l, AP_Delta_m, next(linecycler), linewidth = 2, label = model)
    legend(loc = 'best')
    xlabel('$l$')
    ylabel('$l(l+1)/(2\pi) C_l$(a=1)')

#--: Angular Power (Galaxy)!!
def AP_g(model, BP, f_NL, **cosmo):
    AP_Delta_m = var.galaxy_angular_power(model, BP, f_NL, **cosmo)
    semilogy(cc.l, AP_Delta_m, next(linecycler), linewidth = 2, label = model + " - $f_{NL}: %s$" %f_NL)
    legend(loc = 'best')
    xlabel('$l$')
    ylabel('$l(l+1)/(2\pi) C_l$(a=1)')



#---------------------- END ----------------------