
#============= ORDER INPUT ==================================================================
#plot_order_list = ["BACK","PERT", "LPOW", "VPLOT"]                                         #
#plot_order_list = ["BACK-OMEG", "BACK-h", "PERT-Delta", "PERT-Delta_ini", "PERT-u_m",      #
#                   "PERT-u_m_ini", "PERT-Phi", "PERT-Phi_ini", "LPOW-Delta", "GPOW"]       #
#plot_order_list = ["PERT-Delta", "PERT-Delta_ini"]                                         #
plot_order_list = ["BACK-OMEG", "BACK-h"]#, "LPOW-Delta"]#
#plot_order_list =["EOS"]
#plot_order_list =["LPOW-Phi_p", "LPOW-Phi"]                                                            #
plot_order_list =["GPOW"]
#plot_order_list =["APOW-Delta"]
#--: Effective f_NL
#plot_order_list =["Gamma-FNL"]
#--: variable plot !!
#plot_order_list =["VPLOT"]
#plot_order_list =["AVPLOT"]
#plot_order_list =["fPLOT"] #
#plot_order_list =["Gamma-FNL"]
#plot_order_list =["APOW-galaxy"]
#--: Final Plot !!
#plot_order_list =["LPOW-Delta", "VPLOT", "Gamma-FNL", "APOW-Delta", "APOW-galaxy"]
                  
#==========================================================================================

#_________________________ MAKE_PLOT NGLSSpy SETTING ___________________________
#--: Plotting
import plot as plt
import var_plot as plt2
import pylab  as pplt
import constants as cc
from constants import models
import os, sys
import time

#---------- Save file path/name --------------
save_path = 'figures'

#--- Models definitions !!
models, lcdm_models, wcdm_models, gwcdm_models, qcdm_models = models()
#---------- Cosmo parameters -----------------
cosmo = {'omega_M_0': 0.315, 'omega_b_0': 0.045, 'omega_lambda_0': 1 - 0.315}

def PLOT_JOB(order):
    #---: file name !!
    tf = time.strftime("%H:%M:%S") # %I
    name_of_file = order
    completeName = os.path.join(save_path, name_of_file + "_%s" %tf + ".pdf")
    
    if order == "EOS":
        #---: Equation of state
        fig = pplt.figure()
        #for i in range(len(qcdm_models)):
        fig.show(( plt.w_x(qcdm_models[0]) ))
        fig.savefig(completeName)

    elif order == "QUIN":
        #---: quintessence
        plt.quintess_vars("qcdm")

    elif order == "BACK-OMEG":
        #---: Back - omega's
        fig = pplt.figure()
        for i in range(len(qcdm_models)):
            fig.show(( plt.omega_m(qcdm_models[i]) ))
        fig.savefig(completeName)

    elif order == "BACK-h":
        #---: Back - h
        fig = pplt.figure()
        for i in range(len(qcdm_models)):
            fig.show(( plt.h(qcdm_models[i]) ))
        fig.savefig(completeName)
                     
    elif order == "PERT-Delta":
        #---: Perturb - Delta_m
        fig = pplt.figure()
        for i in range(len(lcdm_models)):
            fig.show(( plt.Delta_m(lcdm_models[i]) ))
        fig.savefig(completeName)
        
        #---:
    elif order == "PERT-Delta_ini":
        fig = pplt.figure()
        for i in range(len(lcdm_models)):
            fig.show(( plt.Delta_m_ini(lcdm_models[i], **cosmo) ))
        fig.savefig(completeName)
    
        #---: Perturb - u_m
    elif order == "PERT-u_m":
        fig = pplt.figure()
        for i in range(len(models)):
            fig.show(( plt.u_m(models[i]) ))
        fig.savefig(completeName)
        
        #---: Perturb - u_m_ini
    elif order == "PERT-u_m_ini":
        fig = pplt.figure()
        for i in range(len(models)):
            fig.show(( plt.u_m_ini(models[i], **cosmo) ))
        fig.savefig(completeName)

        #---: Perturb - Phi
    elif order == "PERT-Phi":
        fig= pplt.figure()
        for i in range(len(models)):
            fig.show(( plt.Phi(models[i]) ))
        fig.savefig(completeName)

    elif order == "PERT-Phi_ini":
        fig = pplt.figure()
        for i in range(len(models)):
            fig.show(( plt.Phi_ini(models[i], **cosmo) ))
        fig.savefig(completeName)

    elif order == "LPOW-Phi_p":
        #---: Power Spectra: linear
        #-: P_Delta_m_p
        fig = pplt.figure()
        for i in range(len(cc.fv)):
            fig.show((plt.P_Phi_p(cc.fv[i], **cosmo)))
        fig.savefig(completeName)

    elif order == "LPOW-Phi":
        #---: Power Spectra: linear
        #-: P_Delta_m_p
        fig = pplt.figure()
        for i in range(len(cc.fv)):
            fig.show((plt.P_Phi("lcdm_0.0", cc.fv[i], **cosmo)))
        fig.savefig(completeName)

    elif order == "LPOW-Delta":
        #---: Power Spectra: linear
        #-: P_Delta_m
        fig = pplt.figure()
            #for i in range(len(models)):
            #for j in range(len(wcdm_models)):
        fig.show((  plt.P_Delta_m(qcdm_models[0], lcdm_models[0],**cosmo)  )) # plt.P_Delta_m(lcdm_models[0], lcdm_models[0],**cosmo),
        # plt.P_Delta_m(lcdm_models[1], lcdm_models[1],**cosmo), plt.P_Delta_m(wcdm_models[1], lcdm_models[1],**cosmo) ))
        fig.savefig(completeName)

    elif order == "BIAS":
        #---: galaxy_bias:
        #-: bias
        pplt.show(( plt.bias("lcdm-%s" %cc.z, cc.f_NL[0], **cosmo), plt.bias("wcdm_-0.8", cc.f_NL[1], **cosmo), plt.bias("wcdm_-0.8", cc.f_NL[2], **cosmo )))

        #-: delta_bias
        #pplt.show(( plt.delta_bias("wcdm_-0.8", cc.f_NL[0], **cosmo), plt.delta_bias("wcdm_-0.8", cc.f_NL[1], **cosmo), plt.delta_bias("wcdm_-0.8", cc.f_NL[2], **cosmo )))
        
    elif order == "GPOW":
        #---: galaxy power spectrum
        fig = pplt.figure()
        fig.show((#plt.Pg_Delta_m(lcdm_models[0], lcdm_models[0], cc.f_NL_p[0], **cosmo),
                  #plt.Pg_Delta_m(lcdm_models[0], lcdm_models[0], cc.f_NL_p[-1], **cosmo),
                  plt.Pg_Delta_m(wcdm_models[0], lcdm_models[0], cc.f_NL_p[0], **cosmo),
                  plt.Pg_Delta_m(wcdm_models[0], lcdm_models[0], cc.f_NL_p[-1] , **cosmo), 
                  plt.Pg_Delta_m(gwcdm_models[0], lcdm_models[0], cc.f_NL_p[0], **cosmo),
                  ))
        fig.savefig(completeName)

    elif order == "APOW-Delta":
        #---: Power Spectra: Angular !!
        #-: Angular_P_Delta_m
        fig = pplt.figure()
        fig.show((plt.AP_m("lcdm-%s" %cc.z, 1,**cosmo), plt.AP_m("wcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.z), 1,**cosmo), plt.AP_m("gwcdm_%s-%s-%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma, cc.Inrct, cc.z), 1, **cosmo) ))
        fig.savefig(completeName)
    
    elif order == "APOW-galaxy":
        #---: Power Spectra: Angular !!
        #-: Angular_P_Delta_m
        fig = pplt.figure()
#        fig.show((plt.AP_g("lcdm-%s" %cc.z, 1, cc.f_NL[1],**cosmo), plt.AP_g("wcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.z), 1, cc.f_NL[1],**cosmo), plt.AP_g("gwcdm_%s-%s-%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma, cc.Inrct, cc.z), 1, cc.f_NL[1], **cosmo) ))

        fig.show((plt.AP_g("lcdm-%s" %cc.z, 1, cc.f_NL[0],**cosmo), plt.AP_g("lcdm-%s" %cc.z, 1, cc.f_NL[1],**cosmo) ))

#        fig.show(( plt.AP_g("wcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.z), 1, cc.f_NL[0],**cosmo),  plt.AP_g("wcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.z), 1, cc.f_NL[1],**cosmo) ))
#
#        fig.show(( plt.AP_g("gwcdm_%s-%s-%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma, cc.Inrct, cc.z), 1, cc.f_NL[0], **cosmo), plt.AP_g("gwcdm_%s-%s-%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma, cc.Inrct, cc.z), 1, cc.f_NL[1], **cosmo) ))

        fig.savefig(completeName)

    elif order == "d-GPOW":
        #---: galaxy power spectrum
    #    pplt.show(( plt.Pg_Delta_m("wcdm_%s-%s" %(cc.w_x, cc.c2_s), 1,cc.f_NL[0], **cosmo), plt.Pg_Delta_m("wcdm_%s-%s" %(cc.w_x, cc.c2_s),1, cc.f_NL[1], **cosmo) ))
    #    pplt.show(( plt.Pg_Delta_m("gwcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma), 1,cc.f_NL[0], **cosmo), plt.Pg_Delta_m("gwcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma),1, cc.f_NL[1], **cosmo) ))
        fig = pplt.figure()
        fig.show(( plt.Pg_Delta_m("wcdm_%s-%s" %(cc.w_x, cc.c2_s), 1,cc.f_NL[0], **cosmo), plt.Pg_Delta_m("wcdm_%s-%s" %(cc.w_x, cc.c2_s),1, cc.f_NL[1], **cosmo), plt.Pg_Delta_m("gwcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma), 1,cc.f_NL[0], **cosmo) ))
        fig.savefig('figures_3/gamma_fNLwcdm_Pg_5.pdf')
    
    elif order == "Gamma-FNL":
        fig = pplt.figure()
        fig.show(( plt.gfnl() ))
        fig.savefig(completeName)

    #---: Vraible plotting !!
    elif order == "VPLOT":
        fig = pplt.figure()
        fig.show((plt2.IN_Pg_var_plot()))
        fig.savefig(completeName)
    
    #---: Effective f_NL - variable plots !!
    elif order == "fPLOT":
        fig = pplt.figure()
        fig.show((plt2.IN_gfnl_var_plot() ))
        fig.savefig(completeName)
    #---: Angular power spectrum variable plotting !!
    elif order == "AVPLOT":
        fig = pplt.figure()
        fig.show((plt2.IN_APg_var_plot()))
        fig.savefig(completeName)

    else:
        print "Please Enter a Valid Order?"
    return

#----------------------------
for i in range(len(plot_order_list)):
    PLOT_JOB(plot_order_list[i])


