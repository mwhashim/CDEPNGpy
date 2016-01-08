from __future__ import division
from numpy import *
from scipy import *
import os, sys
import constants as cc
import var_def as var
from model import model_check
from time import *

from scipy.interpolate import Rbf, interp1d

#---------- Save file path/name --------------
tf1 = strftime("%Y-%m-%d"); folder_name = "Data_%s" %tf1
save_path_ini = folder_name + '/data_ini'

#______ Transfer Function _______________________
def transfer_func(root):
    #------ cosmo parameters ----------------------
    cosmo = {'omega_M_0': 0.315, 'omega_b_0': 0.045, 'omega_lambda_0': 1 - 0.315}
    omega_M_0 = cosmo['omega_M_0']; omega_b_0 = cosmo['omega_b_0']; omega_lambda_0 = cosmo['omega_lambda_0']
    
    if root == "camb":
        k2, Tk2_Delta_m = loadtxt('camb_transfer/camb2010_transfer_out.dat', unpack=True,usecols = [0, 6])
    
        new_length = cc.nn
        new_k2 = linspace(k2.min(), k2.max(), new_length)
        new_Tk2_Delta_m = interp1d(k2, Tk2_Delta_m, kind='cubic')(new_k2)
        return new_Tk2_Delta_m
    elif root == "BBKS":
        #---: transfer k_eq
        Gamma =   omega_M_0 * 0.7 * exp(-omega_b_0 * (1 - sqrt(2 * 0.7)/omega_M_0)) # shape function
        kk_eq = cc.k/Gamma  # constant ??? # with baryons
        x = kk_eq
        Tk = log(1 + 2.34 * x)/(2.34 * x) * (1 + 3.89 * x + (16.19 * x)**2 + (5.46 * x)**3 + (16.71 * x)**4)**(-1/4)
        return Tk
    else:
        print "Please Enter Vaild model name"

def power_primordial(fv, **cosmo):
    #------ cosmo parameters ----------------------
    omega_M_0 = cosmo['omega_M_0']; omega_b_0 = cosmo['omega_b_0']; omega_lambda_0 = cosmo['omega_lambda_0']
    
    AA = 50/9 * pi**2 * cc.delta_H**2 * cc.A**2
    Phi_p = - sqrt(AA) * (cc.k)**(-3/2) * (cc.zk)**((cc.n-1)/2) * omega_M_0
    T = 9/10 * transfer_func("BBKS"); b = 10**3
    return ((b + (b - 1) * fv/cc.k**2)* cc.k**2)**2 * T**2 * Phi_p**2

#===================================================

def vars_ini(model, BP, i, **cosmo):
    
    #---: file name !!
    name_of_file_pert_ini = model + "_PERT_ini"
    completeName_pert_ini = os.path.join(save_path_ini, name_of_file_pert_ini + ".npy")
    
    #------ cosmo parameters ----------------------
    omega_M_0 = cosmo['omega_M_0']; omega_b_0 = cosmo['omega_b_0']; omega_lambda_0 = cosmo['omega_lambda_0']
    
    #------ Primordial Conditions ------------------
    AA = 50/9 * pi**2 * cc.delta_H**2 * cc.A**2 # scaling constant & A is added
    Phi_p = sqrt(AA) * (cc.k)**(-3/2) * (cc.zk)**((cc.n-1)/2) * omega_M_0 # change the sign of the init with respect to the convension (v), it could take +/- since sqrt(P^inf_phi).

    #--: Background Initial conditions
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    
    #--: Setteing Initial Conditions !!
    if BP == 0:
        h_ini = 1; omega_m_ini = omega_M_0; omega_x_ini = 1 - omega_M_0; ww_x_ini = -1.0
        phi_ini = 10**-2; phi_dash_ini = 10**-2
    elif BP == 1:
        h_ini = var.back_vars_ini(model)[2]; omega_m_ini = var.back_vars_ini(model)[0]
        omega_x_ini = var.back_vars_ini(model)[1]; ww_x_ini = 0.0
        phi_ini = 10**-2; phi_dash_ini = 10**-2
    elif BP == 2:
        h_ini = var.back_vars_ini(model)[2]; omega_m_ini = var.back_vars_ini(model)[0]
        omega_x_ini = var.back_vars_ini(model)[1]; ww_x_ini = 0.0
        phi_ini = 10**-2; phi_dash_ini = 10**-2
    else:
        print "Invalid Input !!"

    #--: Interaction model !! # change sign + --> -
    if cc.Inrct == "DMT":
        Q_m_ini = gamma * omega_m_ini; Q_x_ini = - Q_m_ini
    elif cc.Inrct == "DXT":
        Q_m_ini = - gamma * omega_x_ini; Q_x_ini = - Q_m_ini # change sign here
    elif cc.Inrct == "VARSPHI":
        Q_m_ini = cc.beta0 * exp(cc.beta11 * phi_ini); Q_x_ini = - Q_m_ini
    else:
        print "Enter Valid Interaction Term"

    #--: Perturb Initial conditions
    w_m_eff_ini = - Q_m_ini/(3 * h_ini * omega_m_ini)
    w_x_eff_ini = w_x - Q_x_ini/(3 * h_ini * omega_x_ini)
    
    rhodash_mx_ini = (omega_x_ini/omega_m_ini) * (1 + w_x_eff_ini)/(1 + w_m_eff_ini) # check(v)
    mu_delta = rhodash_mx_ini/(1 - rhodash_mx_ini) # check(v)

    Phi_ini = 9/10 * Phi_p * transfer_func("BBKS")
    # change sign for metric convention !!
    Delta_m_ini = 2/3 * (cc.zk/(cc.a_init * h_ini))**2 * Phi_ini/omega_m_ini * (1 + mu_delta)
    Delta_x_ini = 2/3 * (cc.zk/(cc.a_init * h_ini))**2 * Phi_ini/omega_x_ini * mu_delta * PP

    #u_m_ini = cc.a_init * h_ini * omega_m_ini/(cc.zk**2) * (1 + rhodash_mx_ini)/(1 + w_x * omega_x_ini) * Delta_m_ini
    u_m_ini = 2/3 * 1/((1 + w_x * omega_x_ini) * h_ini) * Phi_ini
    u_x_ini = PP * u_m_ini
    
    perturb_vars_ini = Delta_m_ini, Delta_x_ini, u_m_ini, u_x_ini, Phi_ini
    
    #--: Saving initial conditions !!
    save(completeName_pert_ini, perturb_vars_ini)

    return array([ omega_m_ini, omega_x_ini, h_ini, ww_x_ini,
                  Delta_m_ini[i] * BP, Delta_x_ini[i] * BP, u_m_ini[i] * BP, u_x_ini[i] * BP, Phi_ini[i] * BP,
                  phi_ini, phi_dash_ini ])
