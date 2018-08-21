import surfgen as sg
import surfmesh as sm
import surfcont as sc
import matplotlib.pyplot as plt
import surfflow as sf
import sys
import surfplot as splt

plt.style.use('sansserif_paper')

# C:\Users\Philipp\Box Sync\Thesis\image_data\c2_surfcont
if __name__ == '__main__':
    if 0:      
        # generate
        L, H, N, l0 = 1., 0.8, 10, 2
        surfi = sg.self_affine_psd_based_ext(L, L/100., H, N, seed=1, lambda_L_over_lambda_0=l0)
        surfa = sg.self_affine_psd_based_ext_anisotropic(L, L/100., H, N, seed=1, lambda_L_over_lambda_0=l0, anisotropy=8)
        sm.save(surfi, 'surfi')
        sm.save(surfa, 'surfa')
        surfir = sm.interpolate_to_new_N(surfi, 2**(N-2))
        surfar = sm.interpolate_to_new_N(surfa, 2**(N-2))
        sm.save(surfir, 'surfir')
        sm.save(surfar, 'surfar')
        
        # contact
        E, nu, p = 10.0E+9, 0.3, 45.0E6
        conti = sc.contact_FFT(surfi, p, E, nu, err_lim=1.0E-6, verbose=True, store_aperture_field=True)
        conta = sc.contact_FFT(surfa, p, E, nu, err_lim=1.0E-6, verbose=True, store_aperture_field=True)
        sc.save(conti, 'conti')
        sc.save(conta, 'conta')
        contir = sc.contact_FFT(surfir, p, E, nu, err_lim=1.0E-6, verbose=True, store_aperture_field=True)
        contar = sc.contact_FFT(surfar, p, E, nu, err_lim=1.0E-6, verbose=True, store_aperture_field=True)
        sc.save(contir, 'contir')
        sc.save(contar, 'contar')

        # solve flow
        hi_x, hi_y = sf.hydraulic_aperture(contir.aperture_field, contir.dxy)
        ha_x, ha_y = sf.hydraulic_aperture(contar.aperture_field, contar.dxy)
        sf.save(hi_x, 'hi_x')
        sf.save(hi_y, 'hi_y')
        sf.save(ha_x, 'ha_x')
        sf.save(ha_y, 'ha_y')
    else:
        conti = sc.load('conti')
        conta = sc.load('conta')
        surfi = sm.load('surfi')
        surfa = sm.load('surfa')        
        hi_x, hi_y = sf.load('hi_x'), sf.load('hi_y')
        ha_x, ha_y = sf.load('ha_x'), sf.load('ha_y')
        contir = sc.load('contir')
        contar = sc.load('contar')

    # rigid aperture
    apertureri = sm.aperture_surface(surfi.h, surfi.a)
    aperturera = sm.aperture_surface(surfa.h, surfa.a)

    show = False

    # plot hm
    sf.plot_contact_velocity_field_dualmesh(conti.P, conti.dxy, hi_x.v_x, hi_x.v_y, hi_x.dxy, show=show, save_as='hm_i', skinny=True, every=3)
    sf.plot_contact_velocity_field_dualmesh(conta.P, conta.dxy, ha_x.v_x, ha_x.v_y, ha_x.dxy, show=show, save_as='hm_a', skinny=True, every=3)
    sys.exit()
    # auxilliaries
    '''
    figsize = (6,6)  
    sm.plot2D(apertureri, show=show, save_as='aperture_cbar', figsize=figsize, ticks=False,
              forcecbar=True, cbarlabel='$\mathit{a}$ (m)', cbticks=[0., 0.025, 0.05])  
    sm.plot2D(surfi, ticks=False, figsize=figsize, show=show, save_as='surf_cbar', 
              forcecbar=True, cbarlabel='$\mathit{h}$ (m)', cbticks=[-0.025,-0.012,0,0.012,0.025])    

    # plot aperture fields
    figsize = (6,6)
    sm.plot2D(apertureri, show=show, save_as='aperturer_i', figsize=figsize, ticks=False)   
    sm.plot2D(sm.generate_surface(conti.aperture_field, conti.dxy), show=show, save_as='aperturem_i', figsize=figsize, ticks=False)
    sm.plot2D(aperturera, show=show, save_as='aperturer_a', figsize=figsize, ticks=False)
    sm.plot2D(sm.generate_surface(conta.aperture_field, conta.dxy), show=show, save_as='aperturem_a', figsize=figsize, ticks=False)

    # plot surface 
    sm.plot2D(surfi, ticks=False, figsize=figsize, show=show, save_as='surf_i')
    sm.plot2D(surfa, ticks=False, figsize=figsize, show=show, save_as='surf_a') 

    
    '''

    splt.plot_format_is_pdf()
    # plot aperture histograms
    figsize = (8,4)
    sc.plot_aperture_histograms([apertureri.h, conti.aperture_field, aperturera.h, conta.aperture_field], 
                                ['Isotropic, $\mathit{a}$', 'Isotropic, $\mathit{a_m}$', 'Anisotropic, $\mathit{a}$', 'Anisotropic, $\mathit{a_m}$'], 
                                show=show, save_as='aperture_histograms', bin_count=25, figsize=figsize, ylim=(0,60)) 
    # contact cluster size    
    sc.plot_contact_cluster_areas([conti, conta], show=show, 
                                  label=['Isotropic, $\mathit{R_c = %.2f}$' % conti.contact_ratio(),
                                         'Anisotropic, $\mathit{R_c = %.2f}$' % conta.contact_ratio()],
                                  save_as='cont_clusters', figsize=figsize)


    




