"""
    Cosmological background model, like LCDM, qCDM, etc.
"""
from __future__ import division
from numpy import *
import constants as cc
import re

#----------- General Module for different DE models -------------
"""
Perturbation variables:
"""
#--: Model Check
def model_check(model):
    
    model_name = re.split(r"[-+]?\d*\.\d+|\d+", model)[0]
    model_params = re.findall(r"[-+]?\d*\.\d+|\d+", model)
    params = [float(model_params[i]) for i in range(len(model_params))] # params (z, w_x, c2_s, gamma, Inrct)
    z = params[0]
    
    if model_name == "lcdm_":
        w_x = -1.0; c2_x = 0.0; PP = 0.0; gamma = 0.0
    elif model_name == "wcdm_":
        w_x = params[1]; c2_x = params[2]; PP = 1.0; gamma = 0.0
    elif model_name == "gwcdm_":
        w_x = params[1]; c2_x = params[2]; PP = 1.0; gamma = params[3]
    elif model_name == "qcdm_":
        w_x = -1.0; c2_x = 1.0; PP = 1.0; gamma = 0.15
    else:
        w_x = 0; c2_x = 0; PP = 0; gamma = 0.0

    return z, w_x, c2_x, PP, gamma, model_name


#---: Module of differential system of equations !!
def module(vars, x, zk, BP, model):
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    
    Omega_m, Omega_x, h, ww_x = vars[0], vars[1], vars[2], vars[3] # back_vars !!
    
    if model_name == "qcdm_":
        w_x = ww_x
    else:
        w_x = w_x

    Delta_m, Delta_x, u_m, u_x, Phi = vars[4], vars[5], vars[6], vars[7], vars[8] # perturb_vars !!
    
    phi, phi_dash = vars[9], vars[10] # quintessence_vars !!
    
    #---: Decay model selection !! # change sign + --> -
    if cc.Inrct == "DMT":
        Q_m = gamma * Omega_m; Q_x = - Q_m
    elif cc.Inrct == "DXT":
        Q_m = -gamma * Omega_x; Q_x = - Q_m # change sign here
    elif cc.Inrct == "VARSPHI":
        Q_m = cc.beta0 * exp(cc.beta11 * phi); Q_x = - Q_m
    else:
        print "Enter Valid Interaction Term"

    #---: Quintessence_potential !!
    def dVdphi(potential):
        if potential == "exp":
            return - cc.alpha * cc.AP * exp(-cc.alpha * phi)
        elif potential == "Dexp":
            return -(cc.alpha1 * cc.AP1 * exp(-cc.alpha1 * phi) + cc.beta1 * cc.AP2 * exp(-cc.beta1 * phi))
        elif potential == "sugra":
            return cc.AP3 * exp(phi**2)/phi**(cc.alpha2) * (2 * phi - cc.alpha2 * 1/phi)
        elif potential == "Reta-Peebles":
            return -(cc.alpha3/phi) * cc.AP4/phi**(cc.alpha3)
        else:
            return 0.0

    eq_Omega_m = 1/x * (3 * w_x * Omega_m * Omega_x + Q_m/h) # check (v)
    eq_Omega_x = -1/x * (3 * w_x * Omega_m * Omega_x - Q_x/h) # check (v) 
    eq_h = - 3.0/(2.0 * x) * (1 + w_x * Omega_x) * h # check (v)

    
    w_eff = w_x * Omega_x
    u = 1/(1 + w_eff) * (Omega_m * u_m + (1 + w_x) * Omega_x * PP * u_x)
    
    w_m_eff = - Q_m/(3 * h * Omega_m)
    w_x_eff = w_x - Q_x/(3 * h * Omega_x)
    
    delta_Q_m = Q_m * (Delta_x + 3 * x * h * (1 + w_x_eff) * u_x) #check (v)
    delta_Q_x = - delta_Q_m
    
    if cc.Inrct_v == "DXMT":
        #f_m = Q_m * (u - u_m); f_x = - f_m
        f_m = Q_m * (u_m - u); f_x = - f_m
    elif cc.Inrct_v == "DXXT":
        #f_m = Q_m * (u - u_x); f_x = - f_m
        f_m = Q_m * (u_x - u); f_x = - f_m
    else:
        print "Enter Valid Interaction Term"



    #---: change when changin interacting models !!
    if cc.Inrct == "DMT":
        Omega_dash_m =  0.0; Omega_dash_x =  -3 * (w_m_eff - w_x_eff)
    elif cc.Inrct == "DXT":
        Omega_dash_m = 3 * (w_m_eff - w_x_eff); Omega_dash_x =  0.0
    else:
        print "Enter Valid Interaction Term"

#    RHS_m = (- x * Q_m/Omega_m * Omega_dash_m * u_m # check()
#             - x * Q_m/Omega_m * (3 + Q_m/(h * Omega_m)) * (u - u_m)
#             - x * Q_m/Omega_m * (3 + 1/(h * Omega_m)) * f_m
#             - Q_m/(Omega_m * h) * Delta_m
#             + 2 * Q_m/(Omega_m * h) * Phi
#              - x * Q_m/Omega_m * (3 + Q_m/(h * Omega_m)) * u_m
#             + delta_Q_m/(Omega_m * h)
#             ) #check(v)

    #--: change sign of Q_A = - Q_A !!
    RHS_m = (  x * Q_m/Omega_m * Omega_dash_m * u_m # check()
             - x * Q_m/Omega_m * (3 + Q_m/(h * Omega_m)) * (u - u_m)
             - x/Omega_m * (3 + Q_m/(h * Omega_m)) * f_m
             - Q_m/(Omega_m * h) * Delta_m
             + 2 * Q_m/(Omega_m * h) * Phi
             + x * Q_m/Omega_m * (3 + Q_m/(h * Omega_m)) * u_m
             + delta_Q_m/(Omega_m * h)
             ) #check(v)

    eq_Phi = Phi/x + 3/2 * h * (1 + w_eff) * u # check(v) # change sign (- --> +) # check very carfully !!!
    
    RHS_m_u = 1/(Omega_m * x * h) * (Q_m * (u - u_m) + f_m) #check(v)
    eq_u_m = (- u_m/x - Phi/(x**2 * h) + RHS_m_u)
              
    eq_Delta_m = zk**2/(x**2 * h) * u_m + 9/2 * h * (1 + w_eff) * (u_m - u) + RHS_m # check(v) # check sign (v) # + 6 * Phi/x
   
    if model_name == "lcdm_":
        eq_ww_x = 0.0; eq_u_x = 0.0; eq_Delta_x = 0.0; eq_phi = 0.0; eq_phi_dash = 0.0

    elif model_name =="qcdm_":
        eq_w_x = -3/x * (1 - w_x**2) - 2/(x**2 * h**2) * (1 + w_x) * dVdphi("exp")/phi_dash
        
        RHS_x = (- x * Q_x/Omega_x * Omega_dash_x * u_x # check()
                 - x * Q_x/Omega_x * (3 + Q_x/((1 + w_x) * h * Omega_x)) * (u - (1 + 2 * c2_x) * u_x)
                 - x/Omega_x * (3 + Q_x/((1 + w_x) * h * Omega_x)) * f_x
                 + Q_x/(Omega_x * h) * (c2_x/(1 + w_x) - 1) * Delta_x
                 + 2 * Q_x/(Omega_x * h) * Phi
                 - x * Q_x/Omega_x * (3 * (1 + w_x) + Q_x/(h * Omega_x)) * u_x
                 + delta_Q_x/(Omega_x * h)
                 ) #check(v)
        RHS_x_u = 1/((1 + w_x) * Omega_x * x * h) * (Q_x * (u - (1 + 2 * c2_x) * u_x) + f_x) #check(v)
        eq_u_x =  PP * (- u_x/x - c2_x/((1 + w_x) * x**2 * h) * Delta_x - Phi/(x**2 * h) + RHS_x_u)
                        
                                 
        eq_Delta_x = PP * (3 * w_x/x * Delta_x + zk**2/(x**2 * h) * (1 + w_x) * u_x + 9/2 * h * (1 + w_x) * (1 + w_eff) * (u_x - u) + RHS_x) #check(v)
        eq_phi = phi_dash
        eq_phi_dash = -3/x * phi_dash -  1/(x**2 * h**2) * dVdphi("exp")
                        
    else:
        eq_ww_x = 0.0
#        RHS_x = (- x * Q_x/Omega_x * Omega_dash_x * u_x # check()
#                 - x * Q_x/Omega_x * (3 + Q_x/((1 + w_x) * h * Omega_x)) * (u - (1 + 2 * c2_x) * u_x)
#                 - x * Q_x/Omega_x * (3 + 1/((1 + w_x) * h * Omega_x)) * f_x
#                 + Q_x/(Omega_x * h) * (c2_x/(1 + w_x) - 1) * Delta_x
#                 + 2 * Q_x/(Omega_x * h) * Phi
#                 - x * Q_x/Omega_x * (3 * (1 + w_x) + Q_x/(h * Omega_x)) * u_x
#                 + delta_Q_x/(Omega_x * h)
#                 ) #check(v)

        #--: change sign Q_A = -Q_A !!
        RHS_x = (  x * Q_x/Omega_x * Omega_dash_x * u_x # check()
                  - x * Q_x/Omega_x * (3 + Q_x/((1 + w_x) * h * Omega_x)) * (u - u_x)
                  - x/Omega_x * (3 + Q_x/((1 + w_x) * h * Omega_x)) * f_x
                  - Q_x/(Omega_x * h) * (c2_x/(1 + w_x) + 1) * Delta_x
                  + 2 * Q_x/(Omega_x * h) * Phi
                  + x * Q_x/Omega_x * (3 * (1 + w_x) + Q_x/(h * Omega_x)) * u_x
                  + delta_Q_x/(Omega_x * h)
                  ) #check(v)
        
        RHS_x_u = 1/((1 + w_x) * Omega_x * x * h) * (Q_x * (u - u_x) + f_x) #check(v)
        eq_u_x =  PP * (- u_x/x - c2_x/((1 + w_x) * x**2 * h) * Delta_x - Phi/(x**2 * h) + RHS_x_u)
                        
        eq_Delta_x = PP * (3 * w_x/x * Delta_x + zk**2/(x**2 * h) * (1 + w_x) * u_x + 9/2 * h * (1 + w_x) * (1 + w_eff) * (u_x - u) + RHS_x) #check(v)
    
        eq_phi = 0.0
        eq_phi_dash = 0.0

    diff_sys = array([eq_Omega_m, eq_Omega_x, eq_h, eq_ww_x,
                      eq_Delta_m * BP, eq_Delta_x * BP, eq_u_m * BP, eq_u_x * BP, eq_Phi * BP,
                      eq_phi, eq_phi_dash ])
    return diff_sys

