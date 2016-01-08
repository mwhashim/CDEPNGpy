"""
    Various constants used by NGLSSPy code.
    
    Unit abreviations are appended to the name, but powers are not
    specified. For instance, the gravitational constant has units "Mpc^3
    msun^-1 s^-2", but is called "G_const_Mpc_Msun_s".
    
    Most of these numbers are from google calculator.
    
    Constants are:
    ::
    
"""

from __future__ import division
from numpy import *


doc = ""

#=============== MODEL PARAMETRS ==================
#-: Unit convert !!
ln = 35
c = 3.0 * 10**5  # speed of sound !
CH_0 = 100/c  #[h/M_pc]
#-: k-limits
zk_min = 10**(-1); zk_max = 10**(1)
#zk_min = 10**(-1); zk_max = 10**(4) # change to start for k/H_0 = 0.1 # 1
z0 = 0.0 # normalization redshift !!
pk = -1 # pivot scale !!
k_min =  CH_0 * zk_min; k_max =  CH_0 * zk_max

#-: REDSHIFT
#z = 1.0
z = array([0.1]) #array([0.0, 0.5, 1.0])#, 1.5, 2.0])
#-: EOS
w_x = array([-0.9, -1.1])
#-: Speed of Sound
c2_s = array([1.0])

#---: COUPLING
Inrct = "DXT"; Inrct_v = "DXXT" #"VARSPHI" #"DXT" #"DMT", "DXT"
#--: Coupling parameters !!
beta0 = 3; beta11 = 0.1
#-: gamma
gg = 10**(-2)
gamma_p = gg * linspace(0.0, 3.0, 4)
#gamma_p = gg * array([3.0])
gamma_n = gg * linspace(-3.0, 0.0, 4)
#gamma_n = gg * array([-3.0])
#==================================================

#--- RUNNING Models --------------------------------------------
def models():
    models = []; lcdm_models = []; wcdm_models = []; gwcdm_models = []; qcdm_models = []
    for i in range(len(z)):
        lcdm_models.append("lcdm_%s" %z[i])
        qcdm_models.append("qcdm_%s" %z[i])
        models.append( "lcdm_%s" %z[i]); models.append( "qcdm_%s" %z[i])
        for i1 in range(len(w_x)):
            if w_x[i1] > -1.0 : # reverse the sign for check
                gamma  = gamma_p
            elif w_x[i1] < -1.0 :
                gamma  = gamma_n
            else:
                print "Enter avalid value"
            
            for i2 in range(len(c2_s)):
                wcdm_models.append( "wcdm_%s_%s_%s" %(z[i], w_x[i1], c2_s[i2]))
                models.append( "wcdm_%s_%s_%s" %(z[i], w_x[i1], c2_s[i2]))
                for i3 in range(len(gamma)):
                    gwcdm_models.append( "gwcdm_%s_%s_%s_%s_%s" %(z[i], w_x[i1], c2_s[i2], gamma[i3], Inrct) )
                    models.append( "gwcdm_%s_%s_%s_%s_%s" %(z[i], w_x[i1], c2_s[i2], gamma[i3], Inrct) )
    return models, lcdm_models, wcdm_models, gwcdm_models, qcdm_models
#---------------------------------------------------------------

#------ Quintessence Potentials --------------
#---: Exponential potential !!
AP = 0.0219
alpha = 0.08
#---: Douple Exponential potential !!
AP1 = 0.4; AP2 =AP1
alpha1 = 0.1; beta1 = -20
#---: Sugra poential !!
alpha2 = 0.65705469
AP3 = 0.45 * 0.45**(alpha2/4)
##---: Reeta-Peebles potential !!
alpha3 = 0.1
AP4 = 0.45

#----------------------------------------------------------------
doc += "current value of scale: a_0"
a_0 = 1.0

doc += "initial value of scale: a_init"
a_init = 10**(-3)

doc += "length of a,k array: nn"
nn = 10**3

doc += "number of steps: n_step"
n_step = 10**3

doc += "A: A"
A = 1.06

doc += "delta_H: delta_H"
delta_H = 1.9 * 10**(-5)

doc += "spectral index: n"
n = 0.96

doc += "scalar parameter: a"
a = linspace(a_init, a_0, nn)
a_revs = a[::-1]

doc += "normalized k: zk"
zk = linspace(zk_min, zk_max, nn)

doc += "unit convert: CH_0 "
k = CH_0 * zk

l = range(0,ln)

#-: Non-Gaussianity
doc += "non-gaussianity local amplitude: f_NL"
f_NL_p = linspace(0.0, 6.0, nn)
f_NL_n = linspace(-0.6, 0.0, nn)

#-------------- Bias -----------------
doc += "gaussian bais: b_g"
b_0  = 2.0

doc += "critical density: delta_c"
delta_c =1.68

#------- Angular Power Spectrum ----

#------- Primordial power spectrum ---
# find a relation between fv and f_NL and Gamma for Non-gaussianity and coupling constant
# for wcdm aproper definition of Delta could rise it out
amp = 10**-7
fv = array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0]) * amp

#------------------------------------------------
__doc__ += "\n".join(sorted(doc.split("\n")))

