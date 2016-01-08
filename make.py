#!/bin/sh

#============= ORDER INPUT ====================
#order_list = ["SIC-lcdm", "BACK-lcdm", "SIC-wcdm", "BACK-wcdm","SIC-gwcdm", "BACK-gwcdm"]
order_list = [#"PRUN-lcdm", "PRUN-wcdm",
                            "PRUN-gwcdm"]
#order_list = ["PLOT"]

#--: Effective Gamma-f_NL
#order_list = ["f_NL-eff"]

#--: Deleting data and figures
#order_list = ["RMVF"]
#order_list = ["RMVD"]
#order_list = ["MKFLDR"]
#==============================================

#_________________________ MAKE NGLSSpy SETTING ___________________________
import constants as cc
from constants import models
import perturbation as cp
import var_def as var
import os, sys
import params_ini as pini
from time import *

#---------- Cosmo parameters -----------------
cosmo = {'omega_M_0': 0.32, 'omega_b_0': 0.045, 'omega_lambda_0': 1 - 0.315}

#--- Models definitions !!
models, lcdm_models, wcdm_models, gwcdm_models, qcdm_models = models()

#____ Timing: Start ________
t_init = time()
#save_path !!
tf1 = strftime("%Y-%m-%d")
folder_name = "Data_%s" %tf1
save_path = folder_name

def JOB(order):
    if order == "SIC-lcdm":
        #---------- Setting Initial Conditions -------
        print "Setting Initial Conditions ....."
        print "MODELS: ", lcdm_models
        for i in range(len(lcdm_models)):
            cp.modulle_sol(lcdm_models[i] , 0, **cosmo)
    elif order == "SIC-wcdm":
        print "MODELS: ", wcdm_models
        for i in range(len(wcdm_models)):
            cp.modulle_sol(wcdm_models[i] , 0, **cosmo)
    elif order == "SIC-gwcdm":
        print "MODELS: ", gwcdm_models
        for i in range(len(gwcdm_models)):
            cp.modulle_sol(gwcdm_models[i] , 0, **cosmo)
    elif order == "SIC-qcdm":
        print "MODELS: ", qcdm_models
        for i in range(len(qwcdm_models)):
            cp.modulle_sol(qcdm_models[i] , 0, **cosmo)

    elif order == "BACK-lcdm":
        #---------- Background RUN -------
        print "Solve the Background ....."
        print "MODELS: ", lcdm_models
        for i in range(len(lcdm_models)):
            cp.modulle_sol(lcdm_models[i] , 1, **cosmo)
    elif order == "BACK-wcdm":
        print "MODELS: ", wcdm_models
        for i in range(len(wcdm_models)):
            cp.modulle_sol(wcdm_models[i] , 1, **cosmo)
    elif order == "BACK-gwcdm":
        print "MODELS: ", gwcdm_models
        for i in range(len(gwcdm_models)):
            cp.modulle_sol(gwcdm_models[i] , 1, **cosmo)
    elif order == "BACK-qcdm":
        print "MODELS: ", qcdm_models
        for i in range(len(qcdm_models)):
            cp.modulle_sol(qcdm_models[i] , 1, **cosmo)

    elif order == "PRUN-lcdm":
        #---------- Perturbations run ----------------
        print "Solve the Perturbation System ....."
        print "MODELS: ", lcdm_models
        for i in range(len(lcdm_models)):
            cp.modulle_sol(lcdm_models[i] , 2, **cosmo)
    elif order == "PRUN-wcdm":
        print "MODELS: ", wcdm_models
        for i in range(len(wcdm_models)):
            cp.modulle_sol(wcdm_models[i] , 2, **cosmo)
    elif order == "PRUN-gwcdm":
        print "MODELS: ", gwcdm_models
        for i in range(len(gwcdm_models)):
            cp.modulle_sol(gwcdm_models[i] , 2, **cosmo)
    elif order == "PRUN-qcdm":
        print "MODELS: ", qcdm_models
        for i in range(len(qcdm_models)):
            cp.modulle_sol(qcdm_models[i] , 2, **cosmo)
    
    elif order == "f_NL-eff":
        var.f_NL_effect()

    #----------- Fills Orders -----------------------

    elif order == "PLOT":
        #---------- Plot -----------------------------
        print "Plotting ..........."
        import make_plot

    elif order == "RMVD":
        #---------- Removing unneccessary data files -
        print "WARNING: Removing Data Files !!"
        os.system("rm %s/*.npy" %folder_name)
        os.system("rm %s" %folder_name + "/data_ini/*.npy")
        os.system("rm %s" %folder_name + "/gamma_f_NL/*.npy")

    elif order == "VIWF":
        #---------- Removing unneccessary data files -
        print "VIEWING Figures !!"
        os.system("open figures/*.pdf")

    elif order == "RMVF":
        #---------- Removing unneccessary data files -
        print "WARNING: Removing Figures !!"
        os.system("rm figures/*.pdf")
        
    elif order == "MKFLDR":
        #---------- Removing unneccessary data files -
        print "PROCESSING: Making Folder!!"
        os.system("mkdir %s" %folder_name)
        os.system("mkdir %s/data_ini" %folder_name)
        os.system("mkdir %s/gamma_f_NL" %folder_name)

    else:
        print "Please Enter a Valid Order?"
    return

#---------------------------------
print "TASK STARTED !!", "Time Started: ", strftime("%H:%M:%S")
for i in range(len(order_list)):
    JOB(order_list[i])
    
#_______ Timing: Final _________
print 'Total Time Elapsed: ', (time() - t_init)/60.0, 'Mins'
print "TASK FINISHED !!", "Time Ended: ", strftime("%H:%M:%S")