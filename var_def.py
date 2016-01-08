from __future__ import division
from scipy import *
from numpy import *
import os, sys

import constants as cc
import params_ini as pini
from scipy.special import *
from scipy.integrate import *
from model import model_check
from time import *
#from params_ini import power_primordial

#---------- Save file path/name --------------
#save_path !!
tf1 = strftime("%Y-%m-%d"); folder_name = "Data_%s" %tf1
save_path = folder_name; save_path_ini = folder_name +'/data_ini'
#save_path = 'data_output'; save_path_ini = 'data_output/data_ini'

##------ Equation of state --------------
#def eqos(model):
#    w_x = load("output_3/"+model+"_w_x.npy")
#    return w_x
##------ Quintenssence vars -------------
#def quin_vars(model):
#    phi = load("output_3/"+model+"_phi.npy")
#    Lambda = load("output_3/"+model+"_Lambda.npy")
#    return phi, Lambda

#------ Output variables ---------------
#----: BACK VARS !!
#--: Initial Conditions !!
def back_vars_ini(model):
    name_of_file_back_ini = model + "_BACK_ini"
    completeName_back_ini = os.path.join(save_path_ini, name_of_file_back_ini + ".npy")
    return load(completeName_back_ini)

def back_vars(model):
    name_of_file_back = model + "_BACK"
    completeName_back = os.path.join(save_path, name_of_file_back + ".npy")
    return load(completeName_back)

#----: PERT VARS !!
#--: Initial Conditions !!
def perturb_vars_ini(model, **cosmo):
    name_of_file_pert_ini = model + "_PERT_ini"
    completeName_pert_ini = os.path.join(save_path_ini, name_of_file_pert_ini + ".npy")
    return  load(completeName_pert_ini)

def perturb_vars(model):
    name_of_file_pert = model + "_PERT"
    completeName_pert = os.path.join(save_path, name_of_file_pert + ".npy")
    return load(completeName_pert)

#---: Effective f_NL !!
def f_NL_eff(z, w_x):
    f_NL_eff = load("data_output/gamma_f_NL/f_NL_eff_%s_%s.npy" %(z, w_x))
    return f_NL_eff

#------ Growth factor ------------------
def grwoth_factor(model, **cosmo):
    Delta_m, Delta_x, u_m, u_x, Phi = perturb_vars(model)
    Delta_m_ini, Delta_x_ini, u_m_ini, u_x_ini, Phi_ini = perturb_vars_ini(model, **cosmo)
    g_Delta_m = Delta_m/Delta_m_ini; g_Phi = Phi/Phi_ini
    return g_Delta_m, g_Phi

#------ Power Spectra ------------------
#def power_spectra_linear_phi(model, fv, **cosmo):
#    P_Phi_p = power_primordial(fv, **cosmo)
#    P_Delta_m = g_Delta_m**2 * P_Phi_p
#    return P_Delta_m

#------ Power Spectra ------------------
def power_spectra_linear(model, ref_model, **cosmo):
    conv = 1.0/grwoth_factor(model, **cosmo)[1]
    P_Delta_m = conv**2 * (perturb_vars(model)[0])**2
    return P_Delta_m

#------ Galaxy Bias --------------------
def galaxy_bias(model, f_NL, **cosmo):
    #------ cosmo parameters ----------------------
    omega_M_0 = cosmo['omega_M_0']; omega_b_0 = cosmo['omega_b_0']; omega_lambda_0 = cosmo['omega_lambda_0']
    #----: redshift !!
    z = model_check(model)[0]; b_g = cc.b_0 * sqrt(1 + z)

    TT =  3/5 * cc.zk**2 * pini.transfer_func("BBKS") * grwoth_factor(model, **cosmo)[0] 
    bias = b_g + 3 * (b_g - 1) * f_NL * cc.delta_c * omega_M_0/TT
    delta_bias = 3 * (b_g - 1) * f_NL * cc.delta_c * omega_M_0/TT
    return bias, delta_bias

#------ Galaxy power spectrum ----------
def galaxy_power(model, lcdm_model, f_NL, **cosmo):
    Pg_Delta_m = galaxy_bias(model, f_NL, **cosmo)[0]**2 * power_spectra_linear(model, lcdm_model, **cosmo)
    return Pg_Delta_m

#------ Angular Power Spectra ----------
def angular_power(model, BP, **cosmo):
    #romb(y, dx=1.0, axis=-1, show=False)
    chi = trapz(1.0/(cc.a * back_vars(model)[2][0])) #* cc.CH_0

    Cm_l = []
    for j in range(len(cc.l)):
        sph_jn_k = array([sph_jn(cc.l[j], chi * cc.k[i])[0][cc.l[j]] for i in range(len(cc.k))])
        Cm_l.append(2.0/(9.0 * pi) * trapz(cc.k**2 * sph_jn_k**2 * power_spectra_linear(model, BP, **cosmo)) )
    
    Cm_l = array(Cm_l)
    c_l = array([cc.l[i] * (cc.l[i] + 1)/(2 * pi) * Cm_l[i]  for i in range(len(cc.l))])

    return c_l

#--: Galaxy angular power !!
def galaxy_angular_power(model, BP, f_NL, **cosmo):
    #romb(y, dx=1.0, axis=-1, show=False)
    chi = trapz(1.0/(cc.a * back_vars(model)[2][0])) #* cc.CH_0
    
    Cm_l = []
    for j in range(len(cc.l)):
        sph_jn_k = array([sph_jn(cc.l[j], chi * cc.k[i])[0][cc.l[j]] for i in range(len(cc.k))])
        Cm_l.append(2.0/(9.0 * pi) * trapz(cc.k**2 * sph_jn_k**2 * galaxy_power(model, BP, f_NL, **cosmo)) )
    
    Cm_l = array(Cm_l)
    c_l = array([cc.l[i] * (cc.l[i] + 1)/(2 * pi) * Cm_l[i]  for i in range(len(cc.l))])
    
    return c_l

#--: Effective f_NL-Gamma !!
def f_NL_effect():
    #---------- Cosmo parameters -----------------
    cosmo = {'omega_M_0': 0.315, 'omega_b_0': 0.045, 'omega_lambda_0': 1 - 0.315}
    
    print "Calculating The Effective f_NL....."
    f_NL = linspace(cc.f_NL[0], cc.f_NL[-1], cc.nn)
    for jj in range(len(cc.w_x)):
        for kk in range(len(cc.z)):
            f_Nl_eff = []
            for j in range(len(cc.gamma)):
                ref_g_model = "gwcdm_%s_%s_%s_%s_%s" %(cc.z[kk], cc.w_x[jj], cc.c2_s[0], cc.gamma[j], cc.Inrct)
                ref_l_model = "lcdm_%s" %cc.z[kk]
                PPg_Delta_m01 = galaxy_power(ref_g_model, ref_l_model, cc.f_NL[0], **cosmo)[0]
                for i in range(len(f_NL)):
                    ref_w_model = "wcdm_%s_%s_%s" %(cc.z[kk], cc.w_x[jj], cc.c2_s[0])
                    PPg_Delta_m1 = galaxy_power(ref_w_model, ref_l_model, f_NL[i], **cosmo)[0]
                    if PPg_Delta_m01 >= PPg_Delta_m1:
                        f_NL01 = f_NL[i]
                    else:
                        exit
                f_Nl_eff.append(f_NL01)
            
            save("data_output/gamma_f_NL/f_NL_eff_%s_%s.npy" %(cc.z[kk], cc.w_x[jj]), f_Nl_eff)
    return

#______ 3D var_def ______________________
#A, K = meshgrid(cc.a, cc.k)