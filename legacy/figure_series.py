import os
import sys
import keff
import json
import figure_hemisphere
import figure_mohr
#import figure_sphere
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pltpatches
import surfplot as splt
import scipy.stats as scst
import phd_dirtools as phd


splt.plot_format_is_pdf()


def hist_and_match(valarrs_list, marker='x', xlabel='', save_as='', colors=None, legend=None, ymax=None, dists=None, xmax=None):   
    no_bins = 12
    fig, ax = plt.subplots()

    if colors is None:
        colors = phd.color_list()
    if legend is None:
        legend = ['NA' for i in valarrs_list]
    if dists is None:
        dists = [scst.lognorm for i in valarrs_list]
    for i in range(len(valarrs_list)):
        vals = valarrs_list[i]
        color = colors[i]
        # fit
        dist = dists[i]
        param = dist.fit(vals)
        vrange = np.linspace(0., vals.max(), 500)
        pdf_fitted = dist.pdf(vrange, *param[:-2], loc=param[-2], scale=param[-1])
        if dist is scst.expon:
            vrange = vrange[5:]
            pdf_fitted = pdf_fitted[5:]
        # plot
        hist, bin_edges = np.histogram(vals, no_bins, normed=True)
        bin_centers = (bin_edges[:-1]+bin_edges[1:])/2
        ax.plot(vrange, pdf_fitted, color=color, lw=1, label=legend[i]) 
        ax.plot(bin_centers, hist, marker=marker, ls='None', color=color)

    if len(valarrs_list) > 1:
        plt.legend()
    ax.set_ylabel('Relative frequency')
    ax.set_xlabel(xlabel)
    ax.text(0.9, 0.9, '$\mathit{N}$ = %d'%len(vals), va='top', ha='right', transform=ax.transAxes)
    if ymax is not None:
        ax.set_ylim(top=ymax)
    if xmax is not None:
        ax.set_xlim(right=xmax)
    plt.tight_layout()
    if save_as:
        plt.savefig(splt.file_name(save_as))
    plt.show()


# TODO this should call a method that reads all runs into a pandas df and stores to file/reads to file
if __name__ == '__main__':
    cprops = dict()
    for cprop in ['ah_perp', 'ah_para', 'sheard', 'radius', 'ah', 'idx', 'rc']:
        cprops[cprop] = []
    maxcount = 270
    theta, phi = [], []
    ntract, stract = [], []
    theta_kmax, phi_kmax = [], []
    theta_kmed, phi_kmed = [], []
    theta_kmin, phi_kmin = [], []
    kmax, kmed, kmin = [], [], []
    bdir = os.getcwd()
    fracture_summary = dict()
    successful_runs = 0
    kres_avg = keff.permeability_results(3)
    for root, dirs, files in os.walk('.'):
        for dir in dirs:
            if successful_runs == maxcount:
                break
            os.chdir(dir)
            try: # opening files
                s = keff.read_stress_matrix() 
                s = s/s[2,2] # Pa to MPa
                with open('results.json') as f:
                    gresults = json.load(f)           # global results
                fresults = gresults['fractures']      # fracture results
                kres = keff.read_permeability_result_json(gresults)
            except:
                print ('failed', dir)
                os.chdir(bdir)
                continue

            successful_runs += 1
            
            fracture_summary[dir] = fresults
            kres_avg.add(kres)
            theta_kmax_m, phi_f_kmax_m = keff.keff_theta_phi(gresults, 'kmax') 
            theta_kmed_m, phi_f_kmed_m = keff.keff_theta_phi(gresults, 'kmed')            
            theta_kmin_m, phi_f_kmin_m = keff.keff_theta_phi(gresults, 'kmin')            
            theta_f, phi_f = keff.funorm_theta_phi(fresults)
            ntract_f, stract_f = keff.projected_tractions(fresults, s)
            kmax += [kres.k_eigen[0]]
            kmed += [kres.k_eigen[1]]
            kmin += [kres.k_eigen[2]]
            theta_kmax += [theta_kmax_m]
            phi_kmax += [phi_f_kmax_m]
            theta_kmed += [theta_kmed_m]
            phi_kmed += [phi_f_kmax_m]
            theta_kmin += [theta_kmin_m]
            phi_kmin += [phi_f_kmin_m]
            theta += theta_f.tolist()
            phi += phi_f.tolist()
            ntract += ntract_f
            stract += stract_f
            for cprop in cprops:
                cprops[cprop] += keff.fprops(fresults, cprop)
            os.chdir(bdir)
    kres_avg.normalize(successful_runs)
    radii = keff.phi_to_r(phi)
    theta_kmax, phi_kmax = keff.lower_hemisphere_only_theta_phi(theta_kmax, phi_kmax)
    theta_kmed, phi_kmed = keff.lower_hemisphere_only_theta_phi(theta_kmed, phi_kmed)
    theta_kmin, phi_kmin = keff.lower_hemisphere_only_theta_phi(theta_kmin, phi_kmin)
    radii_kmax = keff.phi_to_r(phi_kmax)
    radii_kmed = keff.phi_to_r(phi_kmed)
    radii_kmin = keff.phi_to_r(phi_kmin)
    print ('{:d} successful runs'.format(successful_runs))
    print ('{:d} fractures'.format(len(radii)))

    

    as_text = False
    tformat = '{:d}'
    norm = None
    kmatrix = 5.0E-12

    data = {'kmatrix': kmatrix, 'theta': theta, 'radii': radii.tolist(), 'cprops': cprops,
            'theta_kmax': theta_kmax.tolist(), 'radii_kmax': radii_kmax.tolist(), 'theta_kmed': theta_kmed.tolist(),
            'radii_kmed': radii_kmed.tolist(),  'theta_kmin': theta_kmin.tolist(), 
            'theta_kmin': theta_kmin.tolist(), 'radii_kmin': radii_kmin.tolist(), 'kmax': kmax,
            'kmed': kmed, 'kmin': kmin}
    with open('hmruns.json', 'w') as f:
        json.dump(data, f, indent=4)
    sys.exit()

    kmax_n = np.array(kmax)/kmatrix
    kmed_n = np.array(kmed)/kmatrix
    kmin_n = np.array(kmin)/kmatrix

    cbformat = '%.2f'
    if 0: # contact ratio
        cprop = cprops['rc']
        save_postfix = 'rc'
        label = r'$\mathit{R_c}$'   
    if 0: # radius
        cprop = cprops['radius']
        save_postfix = 'radius'
        label = r'$\mathit{L/2}$'       
    if 1: # normalized average hydraulic aperture
        save_postfix = 'ah_prime'
        cprop = np.array(cprops['ah'])/(2*np.array(cprops['radius']))
        label = r'$\mathit{a_{h}^\prime}$'
        cbformat = '%.4f'
    
    settings_postfix = '_3p1'#'_3p1_diss'#'_omega1p5_lc'#'_3p1'
    #keff.matlab_ellipsoid(kres_avg, fname='keff_ellipsoid'+settings_postfix+'.m', view=(124,14), disc_representation=False)
    #sys.exit()
    #figure_sphere.sphere_cprop(theta_kmax, phi_kmax, np.log10(kmax_n))
    np.savetxt(str('kmax'+settings_postfix+'.txt'), kmax_n)
    with open(str('fracture_summary'+settings_postfix+'.json'), 'w') as f:
        f.write(json.dumps(fracture_summary, indent=2, sort_keys=True))
    hist_and_match([np.log10(kmax_n)], xlabel=r'log$_{10}\mathit{k_{max}^\prime}$', save_as='distribution_kmax'+settings_postfix)
    #hist_and_match(np.log10(kmin_n), xlabel=r'log$_{10}\mathit{k_{min}^\prime}$', save_as='distribution_kmin'+settings_postfix)
    #sys.exit()
    alpha = 1.0
    figure_hemisphere.hemisphere_cprop(theta, radii, cprop, label=label, cbformat=cbformat, save_as='hemisphere_'+save_postfix+settings_postfix, alpha=alpha)
    #sys.exit()
    figure_hemisphere.hemisphere_cprop(theta_kmax, radii_kmax, np.log10(kmax_n), label=r'log$_{10}\mathit{k_{max}^\prime}$', cbformat='%.1f', save_as='hemisphere_kmax'+settings_postfix)
    #figure_hemisphere.hemisphere_cprop(theta_kmed, radii_kmed, np.log10(kmed_n), label=r'log$_{10}\mathit{k_{med}^\prime}$', cbformat='%.1f', save_as='hemisphere_kmed'+settings_postfix)
    #figure_hemisphere.hemisphere_cprop(theta_kmin, radii_kmin, np.log10(kmin_n), label=r'log$_{10}\mathit{k_{min}^\prime}$', cbformat='%.1f', save_as='hemisphere_kmin'+settings_postfix)
    figure_mohr.mohr_cprop(ntract, stract, cprop, s, label, save_as='mohr_'+save_postfix+settings_postfix, as_text=as_text, cbformat=cbformat, tformat=tformat, norm=norm, alpha=alpha)
