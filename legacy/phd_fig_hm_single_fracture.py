import numpy as np
import surfmesh as sm
import surfgen as sg
import surfflow as sf
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import surfplot as splt
import phd_dirtools as phd
import sys, glob, itertools
import pickle as pc
import surfcont as sc
import surfstiff as st
import surfunits as su
import scipy.stats as scst
import scipy.ndimage as spim


su.length_unit_is_micrometer()
splt.plot_format_is_pdf()

def get_hrms(L, pref=0.01):
    L = su.to_meter(L)
    return su.from_meter(pref*L**0.8)

def exp_cont_clusters():
    rc_target = 0.08
    threshold, rc, cont_code, free_code = 100, 0., 0., 255
    while abs(rc-rc_target)>0.01:
        threshold += 1
        cont_img = spim.imread('nemoto_3mm_10MPa_film.png', flatten=True)
        cont_img[cont_img<threshold] = cont_code
        cont_img[cont_img>=threshold] = free_code
        cont_count = (cont_img == cont_code).sum()
        total_count = float(cont_img.size)
        rc = cont_count/total_count
    plt.imshow(cont_img)
    plt.show()
    cont_p = np.zeros_like(cont_img)
    cont_p[cont_img==0] = 1.
    dxy = 50.0E3/np.max(cont_img.shape)
    contact = sc.Results(dxy, cont_p, cont_p)
    sc.save(contact, 'cont_nemoto_3mm_10MPa')
    sc.plot_contact_cluster_areas(contact)


def comp_contacts():
    cnum = sc.load('cont_num_3mm_10MPa')
    cexp = sc.load('cont_nemoto_3mm_10MPa')
    fname = 'cont_cluster_comp'
    sc.plot_contact_cluster_areas([cexp,cnum], label=['Nemoto et al. (2009)','Numerical'], figsize=(8,4), save_as=fname)


def comp_mult_contacts():
    colors, alphas = ['blue'], [1.0]
    cexps = [sc.load('cont_nemoto_3mm_10MPa')]
    cfnames, cnums = glob.glob('cont_num_3mm_10MPa_seed*'), []
    for cfname in cfnames:
        cnums.append(sc.load(cfname))  
        colors.append('red')
        alphas.append(0.3)
    fname = 'cont_cluster_comp_mult'   
    print 'contacts:', len(cnums)
    def hook(ax):
        blue_line = mpl.lines.Line2D([], [], color='blue', label='Nemoto et al., 2009')
        red_line = mpl.lines.Line2D([], [], color='red', alpha=0.3, label='Numerical')
        plt.legend(handles=[blue_line,red_line])
    sc.plot_contact_cluster_areas(itertools.chain(cexps,cnums), figsize=(8,4), save_as=fname, colors=colors, alphas=alphas, hook=hook)


def contact_ratios():
    pk_nem = np.loadtxt('3mm_offset_p_vs_rc.txt')
    p_nem_3 = pk_nem[:,0]/1.0E6 # MPa
    rc_nem_3 = pk_nem[:,1]      # -
    pk_nem = np.loadtxt('2mm_offset_p_vs_rc.txt')
    p_nem_2 = pk_nem[:,0]/1.0E6 # MPa
    rc_nem_2 = pk_nem[:,1]      # -

    f = open('results_rc', 'rb')
    P = pc.load(f)
    AHs = pc.load(f)
    AH_perps = pc.load(f)
    AH_paras = pc.load(f)
    RCs = pc.load(f)
    f.close() 

    fig, ax = plt.subplots()  
    p = np.array(P)/1.0E6  
    for i in range(len(RCs)):
        rc = np.array(RCs[i])
        ax.plot(p, rc, marker='x', lw=0, mec='k', markersize=7)
    ax.plot(p_nem_3, rc_nem_3, marker='o', lw=0)
    ax.plot(p_nem_2, rc_nem_2, marker='o', lw=0)
    ax.set_xlabel('$\mathit{p}$ (MPa)')
    ax.set_ylabel('$\mathit{R_c}$ (-)')
    ax.set_xlim(0,60)
    plt.tight_layout()    
    fname = 'rc_vs_p'
    splt.file_name(fname)
    plt.savefig(fname)
    plt.show()


def comp_apertures():
    f = open('results', 'rb')
    P = pc.load(f)
    AHs = pc.load(f)
    AH_perps = pc.load(f)
    AH_paras = pc.load(f)
    RCs = pc.load(f)
    f.close()

    pk_nem = np.loadtxt('zero_offset_p_vs_k.txt')
    p_nem_0 = pk_nem[:,0]/1.0E6 # MPa
    k_nem_0 = pk_nem[:,1]       # m2
    pk_nem = np.loadtxt('1mm_offset_p_vs_k.txt')
    p_nem_1 = pk_nem[:,0]/1.0E6 # MPa
    k_nem_1 = pk_nem[:,1]       # m2
    pk_nem = np.loadtxt('3mm_offset_p_vs_k.txt')
    p_nem_3 = pk_nem[:,0]/1.0E6 # MPa
    k_nem_3 = pk_nem[:,1]       # m2

    p = np.array(P)/1.0E6

    # k vs p
    pmax_exp = 55
    fig, ax = plt.subplots()    
    colors = ['darkred', 'blue', 'forestgreen']
    for i in range(len(AHs)):
        k = np.array(AHs[i])**2/12
        k_perp = np.array(AH_perps[i])**2/12
        k_para = np.array(AH_paras[i])**2/12
        color = colors[i]
        ax.plot(p, k, marker='x', lw=0, markersize=8, zorder=1, mec=color)
        ax.plot(p, k_perp, marker=r'$\perp$', lw=0, mec='none', markersize=11, zorder=1, mfc=color)
        ax.plot(p, k_para, marker=r'$\parallel$', lw=0, mec='none', markersize=11, zorder=1, mfc=color)
        #ax.errorbar(p, k, yerr=k_para)
    ax.plot(p_nem_0[p_nem_0<pmax_exp], k_nem_0[p_nem_0<pmax_exp], marker='v', lw=0, ms=9, mfc='none', mew=2, mec='k', alpha=0.4, zorder=0)
    ax.plot(p_nem_1[p_nem_1<pmax_exp], k_nem_1[p_nem_1<pmax_exp], marker='^', lw=0, ms=9, mfc='none', mew=2, mec='k', alpha=0.4, zorder=0)
    ax.plot(p_nem_3[p_nem_3<pmax_exp], k_nem_3[p_nem_3<pmax_exp], marker='s', lw=0, ms=9, mfc='none', mew=2, mec='k', alpha=0.4, zorder=0)

    # shear annotations
    ax.text(0.9, 0.94, '$\mathit{\delta_T}$', transform=ax.transAxes, fontsize=22)
    ax.text(0.87, 0.84, '3 mm', transform=ax.transAxes)
    ax.text(0.87, 0.72, '1 mm', transform=ax.transAxes)
    ax.text(0.87, 0.18, '0 mm', transform=ax.transAxes)

    # legend
    x, y = 0.05, 0.05
    ax.plot([x+0.10], [y], marker='v', lw=0, ms=9, mfc='none', mew=2, mec='k', alpha=0.4, zorder=0, transform=ax.transAxes)
    ax.plot([x+0.05], [y], marker='^', lw=0, ms=9, mfc='none', mew=2, mec='k', alpha=0.4, zorder=0, transform=ax.transAxes)
    ax.plot([x], [y], marker='s', lw=0, ms=9, mfc='none', mew=2, mec='k', alpha=0.4, zorder=0, transform=ax.transAxes)
    ax.text(x+0.15, y, '$\mathit{k_{\perp}}$ Watanabe et al. (2008)', va='center',  transform=ax.transAxes, fontsize=14)
    y = y+0.075
    ax.plot([x], [y], marker=r'$\perp$', lw=0, mec='none', markersize=11, mfc='k', transform=ax.transAxes)
    ax.plot([x+0.05], [y], marker=r'$\parallel$', lw=0, mec='none', markersize=11, mfc='k', transform=ax.transAxes)
    ax.plot([x+0.10], [y], marker='x', lw=0, markersize=8, mec='k', transform=ax.transAxes)
    ax.text(x+0.15, y, r'$\mathit{k_{\perp}, k_{\parallel}}$, and harmonic mean, numerical', va='center', transform=ax.transAxes, fontsize=14)
    b = mpatch.Rectangle((0.01,0.01), 0.75, 0.17, color='none', ec='k', lw=1, transform=ax.transAxes)
    ax.add_artist(b)

    ax.set_xlabel('$\mathit{p}$ (MPa)')
    ax.set_ylabel('$\mathit{k}$ ($m^2$)')
    ax.set_yscale('log')
    ax.set_ylim(10**-13, 10**-8)
    ax.set_xlim(0,60)
    #ax.text(0.05, 0.05, 'E = %d GPa, $L_{mm}$ = %d $\mu$m, N = %d, pref = %.4f'%(int(E/1E9), int(L_mm), N, pref), transform=ax.transAxes)
    plt.tight_layout()    
    fname = 'E%d_Lmm%d_N%d_pref%d'%(int(E/1E9), int(L_mm), N, int(1000000)*pref)
    fname = 'k_vs_p_watamabe_comparison'
    fname = splt.file_name(fname)
    plt.savefig(fname)
    plt.show()


def num_hm(seed, fname_contains_seed=False):
    lo = sg.self_affine_psd_based_ext(L, hrms, 0.8, N, seed=seed)
    up = sg.self_affine_psd_based_ext(L, hrms, 0.8, N, seed=seed, lambda_L_over_lambda_mm=lambda_L_over_lambda_mm)
    
    # compute numerical results and save to binary
    O = [0., 1000., 3000.] 
    P = [10.0E6,20.0E6,30.0E6,40.0E6,50.0E6]
    #O = [3000.] 
    #P = [10.0E6]
    AHs, AH_perps, AH_paras, RCs = [], [], [], []
    for offset in O:
        sco = sm.composite_surface_in_x(up.h, lo.h, lo.a, offset)
        #sm.plot2D(sco, ticks=False)            
        AH, AH_perp, AH_para, RC = [], [], [], []
        for p in P:
            contact = sc.contact_FFT(sco, p, E, nu, store_aperture_field=True, verbose=True)
            if offset == 3000. and p == 10.0E6:
                cfname = 'cont_num_3mm_10MPa'
                if fname_contains_seed:
                    cfname += '_seed'+str(seed)
                sc.save(contact, cfname)
            #sc.plot_aperture_field(contact.aperture_field, contact.dxy, save_as='am_3mm_10MPa')
            #comp_contacts()
            #sys.exit()
            #sc.plot_clusters(contact)
            #sc.plot_contact_cluster_areas(contact)
            flow_x, flow_y = sf.hydraulic_aperture(contact.aperture_field, contact.dxy, verbose=True)
            ah_iso = scst.hmean([flow_x.a, flow_y.a])
            AH += [su.to_meter(ah_iso)]
            AH_perp += [su.to_meter(flow_y.a)]
            AH_para += [su.to_meter(flow_x.a)]
            RC += [contact.contact_ratio()]
        AHs += [AH]
        AH_perps += [AH_perp]
        AH_paras += [AH_para]
        RCs += [RC]

    bfname = 'results'
    if fname_contains_seed:
        bfname += '_seed'+str(seed)
    f = open(bfname, 'wb')
    pc.dump(P, f)
    pc.dump(AHs, f)
    pc.dump(AH_perps, f)
    pc.dump(AH_paras, f)
    pc.dump(RCs, f)
    f.close()
        

# C:\Users\Philipp\Box Sync\Thesis\image_data\c3_hm_single_fracture
if __name__ == '__main__':

    L, N = 100.0E3, 8
    pref = 0.0023
    pref = 0.0023
    hrms = get_hrms(L, pref)
    L_mm = 55.0
    lambda_L_over_lambda_mm = L/L_mm
    E, nu = 30.0E9, 0.15
    seed = 2
    
    if 0: # compute nemoto contact distribution and save to file
        exp_cont_clusters()


    # single fracture
    # ===============

    if 0: # compute numerical results and save to binary
        num_hm(seed)

    if 0: # load nemoto and numerical contact and plot comparison
        comp_contacts()

    if 0: # ploat numerical vs nemoto Rc vs load
        contact_ratios()

    if 0: # load nemoto and numerical apertures and plot comparison 
        comp_apertures()


    
    # multiple fracture realizations
    # ==============================

    if 0: # compute numerical results and save to binary
        for i in range(50):
            num_hm(i, True)

    if 1: # load nemoto and numerical contact and plot comparison
        comp_mult_contacts()






    

    

    
    
