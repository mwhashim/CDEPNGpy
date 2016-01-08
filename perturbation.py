"""
Perturbations calculations for cosmology !!    
"""
from __future__ import division
from numpy import *
from scipy import *
from scipy.integrate import *
from scipy.special import *
import os, sys

import constants as cc
from model import module
from model import model_check
from params_ini import vars_ini

from time import *

#---------- Save file path/name --------------
#save_path !!
tf1 = strftime("%Y-%m-%d"); folder_name = "Data_%s" %tf1
save_path = folder_name; save_path_ini = folder_name +'/data_ini'

#----------- General Module for different DE models -------------
"""
    Perturbation variables:
"""
#---: Module_v1
def modulle_sol(model, BP, **cosmo):
    """
    """
    #---: file name !!
    name_of_file_back = model + "_BACK"
    name_of_file_back_ini = model + "_BACK_ini"
    name_of_file_pert = model + "_PERT"
    
    completeName_back_ini = os.path.join(save_path_ini, name_of_file_back_ini + ".npy")
    completeName_back = os.path.join(save_path, name_of_file_back + ".npy")
    completeName_pert = os.path.join(save_path, name_of_file_pert + ".npy")

    #--: Model Parameters defs
    z, w_x, c2_x, PP, gamma, model_name = model_check(model)
    
    #--: Model def printing !!
    print "Model:", model
    print "Background:", cosmo

    #____ Timing: Start ________
    t_init = time()
    
    #------- Redshift -----------------------------
    a_q = 1.0/(1.0 + z)
    aa_q = linspace(cc.a_init, a_q, cc.nn)
    print 'Redshift:', z, '-', 'Scalar Factor:', a_q
    
    #------ Solving diff_sys ----------------------
    if BP == 0:
        vars = odeint(module, vars_ini(model, BP, 0, **cosmo), cc.a_revs, args = (cc.zk[0], BP, model), mxstep = cc.n_step)
        #------ Variable definitions ------------------
        Omega_m = vars[:,0]; Omega_x = vars[:,1]; h = vars[:,2]; ww_x = vars[:,3]
        back_vars_ini = Omega_m[-1], Omega_x[-1], h[-1], ww_x[-1]
        #------ Variable Save ------------------
        save(completeName_back_ini, back_vars_ini)
    
    elif BP == 1:
        vars = odeint(module, vars_ini(model, BP, 0, **cosmo), cc.a, args = (cc.zk[0], BP, model), mxstep = cc.n_step)
        #------ Variable definitions ------------------
        Omega_m = vars[:,0]; Omega_x = vars[:,1]; h = vars[:,2]; ww_x = vars[:,3]
        phi = vars[:,9]; Lambda = vars[:,10]
        back_vars = Omega_m, Omega_x, h, ww_x
        #------ Variable Save ------------------
        save(completeName_back, back_vars)
    
    elif BP == 2:
        vars = [odeint(module, vars_ini(model, BP, i, **cosmo),  aa_q, args = (cc.zk[i], BP, model), mxstep = cc.n_step) for i in range(cc.nn)] # i for k
        #------ Variable definitions ------------------
        Delta_m = array([vars[j][:,4][-1] for j  in range(cc.nn)])
        Delta_x = array([vars[j][:,5][-1] for j  in range(cc.nn)])
        u_m = array([vars[j][:,6][-1] for j  in range(cc.nn)])
        u_x = array([vars[j][:,7][-1] for j  in range(cc.nn)])
        Phi = array([vars[j][:,8][-1] for j  in range(cc.nn)])

        perturb_vars = Delta_m, Delta_x, u_m, u_x, Phi
        #------ Variable Save ------------------
        save(completeName_pert, perturb_vars)

    else:
        print "Invalid Input !!"

    #_______ Timing: Final _________
    print 'Time Elapsed: ', (time() - t_init)/60.0, 'Mins'
    return
