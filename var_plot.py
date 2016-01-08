from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import constants as cc
import var_def as var
import plot as plt

from numpy import save

from decimal import Decimal, ROUND_DOWN

#-----: Discrete Slider !!
class DiscreteSlider(Slider):
    """A matplotlib slider widget with discrete steps."""
    def __init__(self, *args, **kwargs):
        """
            Identical to Slider.__init__, except for the new keyword 'allowed_vals'.
            This keyword specifies the allowed positions of the slider
            """
        self.allowed_vals = kwargs.pop('allowed_vals',None)
        self.previous_val = kwargs['valinit']
        Slider.__init__(self, *args, **kwargs)
        if self.allowed_vals==None:
            self.allowed_vals = [self.valmin,self.valmax]
    
    def set_val(self, val):
        discrete_val = self.allowed_vals[abs(val-self.allowed_vals).argmin()]
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % discrete_val)
        if self.drawon:
            self.ax.figure.canvas.draw()
        self.val = val
        if self.previous_val!=discrete_val:
            self.previous_val = discrete_val
            if not self.eventson:
                return
            for cid, func in self.observers.iteritems():
                func(discrete_val)

#---: Interacting-Non-Gaussianity matching Galaxy Power Spectrum !!
def IN_Pg_var_plot():
    #---------- Cosmo parameters -----------------
    cosmo = {'omega_M_0': 0.315, 'omega_b_0': 0.045, 'omega_lambda_0': 1 - 0.315}
    
    ax = subplot(111)
    subplots_adjust(left=0.25, bottom=0.25)
    
    f_NL0 = cc.f_NL[0]
    f_NL1  = cc.f_NL[1]
    
    #galaxy_power(model, lcdm_model, f_NL, **cosmo)
    ref_l_model = "lcdm_%s" %cc.z[0]
    ref_w_model = "wcdm_%s_%s_%s" %(cc.z[0], cc.w_x[0], cc.c2_s[0])
    ref_g_model = "gwcdm_%s_%s_%s_%s_%s" %(cc.z[0], cc.w_x[0], cc.c2_s[0], cc.gamma[2], cc.Inrct)
    
    Pg_Delta_m0 = var.galaxy_power(ref_w_model, ref_l_model, f_NL0, **cosmo)
    Pg_Delta_m1 = var.galaxy_power(ref_w_model, ref_l_model, f_NL1, **cosmo)
    Pg_Delta_m01 = var.galaxy_power(ref_g_model, ref_l_model, f_NL0, **cosmo)
    
    PPg_Delta_m01 = var.galaxy_power(ref_g_model, ref_l_model, f_NL0, **cosmo)[0]

    
    f_NL = linspace(f_NL0, f_NL1, 1000)
    for i in range(len(f_NL)):
        PPg_Delta_m1 = var.galaxy_power(ref_w_model, ref_l_model, f_NL[i], **cosmo)[0]
        if PPg_Delta_m01 >= PPg_Delta_m1:
            f_NL01  = f_NL[i]
        else:
            exit
    print "Effective f_NL: ", f_NL01

    text(0.05, 0.9,'$f^{eff}_{NL}$ = %s' %round(f_NL01, 3) , ha ='left', va ='top', transform = ax.transAxes)

    l, = loglog(cc.k, Pg_Delta_m0, '-', linewidth = 2)
    l1, = loglog(cc.k, Pg_Delta_m1, '--', linewidth = 2, label = "wcdm - $f_{NL}$")
    l2, = loglog(cc.k, Pg_Delta_m01, '-.', linewidth = 2, label = "$\Gamma$-wcdm")
    l3, = loglog(cc.k, Pg_Delta_m01, ':', linewidth = 2, label = "$\Gamma$-wcdm + $f_{NL}$" )
    
    axvline(x = 0.17 * 10**-2, color='r', ls = '--') # vertical line at k_eq
    axvline(x = 0.14 * 10**-1, color='b', ls = '-.') # vertical line at k_eq
    
    xlim((cc.k_min, cc.k_max))
    #xlim((10**-3, cc.k_max))
    axvline(x = 10**-3, color='g', ls = '--') # vertical line at k_eq
    legend(loc = 'best')
    
    xlabel('k')
    ylabel('$Pg_m$(k, a = 1)')
    
    axcolor = 'lightgoldenrodyellow'

    
    axf_NL = axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
    axgamma = axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axw_x = axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    axz = axes([0.025, 0.2, 0.1, 0.04], axisbg=axcolor)

    sf_NL = Slider(axf_NL, '$f_{NL}$', f_NL0, f_NL1, valfmt='%1.3f', valinit = f_NL01)
    sgamma = DiscreteSlider(axgamma, '$\Gamma$', 0.0, 1.0, allowed_vals = cc.gamma, valfmt='%1.1f', valinit = cc.gamma[2])
    sz = DiscreteSlider(axz,'z', cc.z[0], cc.z[-1], valfmt='%1.1f' , allowed_vals = cc.z, valinit = cc.z[0])
    sw_x = DiscreteSlider(axw_x, '$w_x$', cc.w_x[0], cc.w_x[-1], allowed_vals = cc.w_x , valfmt='%1.1f', valinit = cc.w_x[0])
    
    def update(val):
        f_NL = sf_NL.val
        gamma = Decimal(sgamma.val).quantize(Decimal('0.1'),rounding=ROUND_DOWN)
        z = Decimal(int(sz.val)).quantize(Decimal('0.0'),rounding=ROUND_DOWN)
        w_x = Decimal(sw_x.val).quantize(Decimal('0.1'),rounding=ROUND_DOWN)

        ref_l_model1 = "lcdm_%s" %z
        ref_w_model1 = "wcdm_%s_%s_%s" %(z, w_x, cc.c2_s[0])
        ref_g_model1 = "gwcdm_%s_%s_%s_%s_%s" %(z, w_x, cc.c2_s[0], gamma, cc.Inrct)

        Pg_Delta_m = var.galaxy_power(ref_w_model1, ref_l_model1, f_NL, **cosmo)
        Pg_Delta_m11 = var.galaxy_power(ref_g_model1, ref_l_model1, f_NL0, **cosmo)
        Pg_Delta_mt = var.galaxy_power(ref_g_model1, ref_l_model1, f_NL, **cosmo)
        
        l.set_ydata(Pg_Delta_m0)
        l1.set_ydata(Pg_Delta_m)
        l2.set_ydata(Pg_Delta_m11)
        l3.set_ydata(Pg_Delta_mt)
        
        draw()
    
    sf_NL.on_changed(update)
    sgamma.on_changed(update)
    sz.on_changed(update)
    sw_x.on_changed(update)
    
    resetax = axes([0.8, 0.0, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    
    def reset(event):
        sf_NL.reset()
        sgamma.reset()
        sz.reset()
        sw_x.reset()

    button.on_clicked(reset)
    show()
    
    return

#--: Interacting-Non-Gaussianity matching Angular Power Spectrum !!
def IN_APg_var_plot():
    #---------- Cosmo parameters -----------------
    cosmo = {'omega_M_0': 0.315, 'omega_b_0': 0.045, 'omega_lambda_0': 1 - 0.315}
    
    ax = subplot(111)
    subplots_adjust(left=0.25, bottom=0.25)
    
    f_NL0 = cc.f_NL[0]
    f_NL1  = cc.f_NL[1]
    
    APg_Delta_m0 = var.galaxy_angular_power("wcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.z), 1, f_NL0, **cosmo)
    APg_Delta_m1 = var.galaxy_angular_power("wcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.z), 1, f_NL1, **cosmo)
    APg_Delta_m01 = var.galaxy_angular_power("gwcdm_%s-%s-%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma, cc.Inrct, cc.z), 1, f_NL0, **cosmo)
    
    PPg_Delta_m01 = var.galaxy_angular_power("gwcdm_%s-%s-%s-%s-%s" %(cc.w_x, cc.c2_s, cc.gamma, cc.Inrct, cc.z), 1, f_NL0, **cosmo)[0]
    
    
#    f_NL = linspace(f_NL0, f_NL1, 1000)
#    for i in range(len(f_NL)):
#        PPg_Delta_m1 = var.galaxy_angular_power("wcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, cc.z), 1, f_NL[i], **cosmo)[0]
#        if PPg_Delta_m01 >= PPg_Delta_m1:
#            f_NL01  = f_NL[i]
#        else:
#            exit
#    print "Effective f_NL: ", f_NL01

    l, = semilogy(cc.l, APg_Delta_m0, '-', linewidth = 2)#, label = "wcdm - $f_{NL}$: %s" % cc.f_NL[0] )
    l1, = semilogy(cc.l, APg_Delta_m1, '--', linewidth = 2, label = "wcdm - $f_{NL}$")
    l2, = semilogy(cc.l, APg_Delta_m01, '-.', linewidth = 2, label = "$\Gamma$-wcdm")
    l3, = semilogy(cc.l, APg_Delta_m01, ':', linewidth = 2, label = "$\Gamma$-wcdm + $f_{NL}$" )
    
    #axvline(x = 0.17 * 10**-2, color='r', ls = '--') # vertical line at k_eq
    #axvline(x = 0.14 * 10**-1, color='b', ls = '-.') # vertical line at k_eq
    
    #xlim((cc.k_min, cc.k_max))
    #xlim((10**-3, cc.k_max))
    #axvline(x = 10**-3, color='g', ls = '--') # vertical line at k_eq
    legend(loc = 'best')
    
    xlabel('l')
    ylabel('$l(l+1)/(2\pi) C_l$(a=1)')
    
    axcolor = 'lightgoldenrodyellow'
    
    
    axf_NL = axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
    axgamma = axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axz = axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    
    sf_NL = Slider(axf_NL, '$f_{NL}$', f_NL0, 10.0, valfmt='%1.3f', valinit = f_NL1)
    sgamma = Slider(axgamma, '$\Gamma$', 0.0, 1.0, valfmt='%1.1f', valinit = cc.gamma)
    sz = Slider(axz, '$z$', 0.0, 3.0, valfmt='%1.0f', valinit = cc.z)
    
    def update(val):
        f_NL = sf_NL.val
        gamma = Decimal(sgamma.val).quantize(Decimal('0.1'),rounding=ROUND_DOWN)
        z = Decimal(int(sz.val)).quantize(Decimal('0.0'),rounding=ROUND_DOWN)
        
        APg_Delta_m = var.galaxy_angular_power("wcdm_%s-%s-%s" %(cc.w_x, cc.c2_s, z), 1, f_NL, **cosmo)
        APg_Delta_m11 = var.galaxy_angular_power("gwcdm_%s-%s-%s-%s-%s" %(cc.w_x, cc.c2_s, gamma, cc.Inrct, z), 1, f_NL0, **cosmo)
        APg_Delta_mt = var.galaxy_angular_power("gwcdm_%s-%s-%s-%s-%s" %(cc.w_x, cc.c2_s, gamma, cc.Inrct, z), 1, f_NL, **cosmo)
        
        l.set_ydata(APg_Delta_m0)
        l1.set_ydata(APg_Delta_m)
        l2.set_ydata(APg_Delta_m11)
        l3.set_ydata(APg_Delta_mt)
        
        draw()
    
    sf_NL.on_changed(update)
    sgamma.on_changed(update)
    sz.on_changed(update)
    
    resetax = axes([0.8, 0.0, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    
    def reset(event):
        sf_NL.reset()
        sgamma.reset()
        sz.reset()
    
    button.on_clicked(reset)
    show()
    
    return

#--: Effective Non-Gaussainity !!
def IN_gfnl_var_plot():
    #---------- Cosmo parameters -----------------
    cosmo = {'omega_M_0': 0.315, 'omega_b_0': 0.045, 'omega_lambda_0': 1 - 0.315}
    
    ax = subplot(111)
    subplots_adjust(left=0.25, bottom=0.25)
    
    eff_f_NL0 = var.f_NL_eff(cc.z[0], cc.w_x[0])
    
    
    l, = plot(eff_f_NL0, cc.gamma, '-.o', linewidth = 2)
    
    xlabel('$f^{eff}_{NL}$')
    ylabel('$\Gamma$/$H_0$')
    
    axcolor = 'lightgoldenrodyellow'
    
    
    axz = axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
    axw_x = axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    
    sz = DiscreteSlider(axz,'z', cc.z[0], cc.z[-1], valfmt='%1.1f' , allowed_vals = cc.z, valinit = cc.z[0])
    sw_x = DiscreteSlider(axw_x, '$w_x$', cc.w_x[0], cc.w_x[-1], allowed_vals = cc.w_x , valfmt='%1.1f', valinit = cc.w_x[0])
    
    def update(val):
        
        w_x = Decimal(sw_x.val).quantize(Decimal('0.1'),rounding=ROUND_DOWN)
        z = Decimal(int(sz.val)).quantize(Decimal('0.0'),rounding=ROUND_DOWN)

        eff_f_NL1 = var.f_NL_eff(z, w_x)
        l.set_xdata(eff_f_NL1)
        
        draw()
    
    sw_x.on_changed(update)
    sz.on_changed(update)
    
    
    resetax = axes([0.8, 0.0, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    
    def reset(event):
        sw_x.reset()
        sz.reset()
    
    button.on_clicked(reset)
    show()
    
    return

