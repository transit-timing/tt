import sys
from datetime import datetime
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from mandelagol import occultquad
from scipy import optimize
from scipy.optimize import minimize_scalar
import lmfit as lmfit
from utilities import *
from matplotlib.ticker import MaxNLocator
 
# This script contains code for processing 2, 10, and 30-minute TESS light curves


plt.rcParams['figure.figsize'] = (16.0, 8.0)
plt.rcParams['lines.markersize'] = np.sqrt(30)
plt.rcParams['font.size'] = 18
plt.rcParams['font.serif'] = "Times New Roman"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
 
 
direct = os.getcwd()

def analyze_system(params_system, flux, time):
    '''
    Process a TESS light curve (2-minute, 10-minute, or 30-minute cadence) and estimate mid-transit times
    '''

    # Unpack the parameters
    name, period_ref, t0_ref, duration, depth, cadence = params_system
    t0_ref = t0_ref - 2457000.0 # convert to TESS time
    
    # initial guesses and other parameters
    radratio = np.sqrt(depth) # planet/star radius ratio
    impactparam = 0.5 # impact parameter
    limbdark1 = 0.5 # first limb-darkening coefficient
    limbdark2 = 0 # second limb-darkening coefficient
    a_over_r = period_ref/(duration*np.pi)*np.sqrt(1.-impactparam**2) # ratio of planet's semi-major axis and host star's radius
    
    n_durations = 4 # number of durations 
    oot_factor = 1.5 # out-of-transit factor

    nsigma_clip = 5    
    method='leastsq'
    n_links = 1000 # number of links in MCMC
    niter_sigmaclip = 10 # number of iterations in sigma-clipping procedure
    ndeg_max = 3 # maximum degree of polynomial for detrending

    # set the expected minimum number of data points in-transit and out-of-transit based on cadence
    if int(cadence) == 120:
        min_npts_tra = np.minimum(5, 0.75*duration*24.*60 / 2) # 75% of the transit should be covered
        min_npts_oot = 0.9*duration*24*60/2       
    elif int(cadence) == 600:
        min_npts_tra = np.minimum(5, 0.75*duration*24.*60 / 10) # 75% of the transit should be covered
        min_npts_oot = 0.9*duration*24*60/10    
    elif int(cadence) == 1800:
        min_npts_tra = np.minimum(5, 0.75*duration*24.*60 / 30) # 75% of the transit should be covered
        min_npts_oot = 0.9*duration*24*60/30     

    # Make a directory for the output
    dirname = direct + '/output/' + name + '/'

    if not os.path.exists(direct + '/output/'):
        os.mkdir(direct + '/output/')

    if not os.path.exists(dirname):
        os.mkdir(dirname)
        
    # Open the logfile
    original_stdout = sys.stdout
    sys.stdout = open(dirname + name + '_log.txt','w')
 
    print(name)
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print('Starting at ', date_time)

    # open the results file
    results_file = open(dirname + name + '_results.txt', 'w')

    # remove NaNs
    ok = np.isfinite(flux)
    time=time[ok]
    flux=flux[ok]
 
    time_folded, orbit_number = fold_time(time, t0_ref, period_ref)
 
    orbits = np.unique(orbit_number) # a list of the orbits with valid timestamps

    #############################################################################################
    # select orbits that have more data points than the expected minimum number calculated above.
    #############################################################################################

    oibeo = -1.0 + 0.0*time
    #  oibeo tells us:
    #   -1 is not analyzed, too far from transit
    #    0 is in-transit
    #    1 is pre-transit
    #    2 is post-transit
    #   >0 is OOT
    rejected_orbits = 0
    for n in orbits:

        this_orbit = (orbit_number == n) & (np.abs(time_folded) <= 0.5*n_durations*duration)
        pre = this_orbit & (time_folded <= -0.5*oot_factor*duration)
        post = this_orbit & (time_folded >=  0.5*oot_factor*duration)
        oot = pre | post
        tra = this_orbit & (np.abs(time_folded) <= 0.5*duration)
        n_pre = np.sum(pre)
        n_post = np.sum(post)
        n_tra = np.sum(tra)
        if (n_pre >= min_npts_oot) & (n_tra >= min_npts_tra) & (n_post >= min_npts_oot):
            oibeo[this_orbit] = 0
            oibeo[pre] = 1
            oibeo[post] = 2
        else:
            rejected_orbits +=1
            print('Rejecting data from orbit ', n)
            print('  n_tra, min_npts_tra = ', n_tra, min_npts_tra, ' and n_pre, n_post, min_npts_oot = ', n_pre, n_post, min_npts_oot)
            
    keep = (oibeo >= 0)
    if (np.sum(keep) == 0):
        print('No transits occurred during the timespan of TESS observations.\n')
        plt.figure(figsize=(10,10))
        plt.plot([0,0],[1,1],'w.')
        plt.text(0.5,0.5, name + ': no transits this sector')
        plt.savefig(dirname + name + '_NoData.png',bbox_inches='tight')
        plt.close()
        
        sys.stdout = original_stdout
         
    
    t = time[keep]
    f = flux[keep]
    oibeo = oibeo[keep]
    orbit_number = orbit_number[keep]
    orbits = np.unique(orbit_number)

    # make a copy of flux and time array that is read only (we will not modify it throughout the script)
    flux_data = np.copy(f)
    time_data = np.copy(t) 
 
            
 
    ##################################################################################
    # Detrend each individual transit, fold the data and find a preliminary transit model 
    # (before we omit any orbits with high MAD)
    ##################################################################################

    f_unc = 0.*t # we will set this equal to the std of oot after detrending
    initial_detrended_events = [] # array to store time series detrended at this stage.
    initial_flux_uncertainties = [] # array to store flux uncertainties coming from detrending at this stage.
    n_cols = 3
    n_rows = np.int(np.round(len(orbits)/n_cols))+1
    fig = plt.figure(figsize=(15,3*n_rows))
    j=0
    for n in orbits:

        this_orbit = (orbit_number == n)
        tp = t[this_orbit]
        fp = f[this_orbit]
        oot = (oibeo[this_orbit] > 0.)

        j=j+1
        plt.subplot(n_rows,n_cols,j)
        plt.plot(tp, fp, 'k.')

        tmid = np.mean(tp[oot])
        coeff = detrend(tp[oot]-tmid, fp[oot], ndeg_max)
        f_baseline = np.polyval(coeff,tp-tmid)
        plt.plot(tp[oot], f_baseline[oot], 'r.')

        fp = fp/f_baseline
        f[this_orbit] = fp # !! de-trended flux
        scatter = np.std(fp[oot])
        f_unc[this_orbit] = scatter

        initial_detrended_events.append(fp)
        initial_flux_uncertainties.append(f_unc[this_orbit])

        print('   Detrended orbit ', n, ' with polynomial of order ', len(coeff)-1, ', scatter = ', scatter)

    plt.close(fig)

    initial_detrended_events_flat = [item for sublist in initial_detrended_events for item in sublist]
    initial_detrended_flux = np.array(initial_detrended_events_flat)
    initial_flux_uncertainties = [item for sublist in initial_flux_uncertainties for item in sublist]
    initial_flux_uncertainties = np.array(initial_flux_uncertainties)

    
    # find preliminary fit to folded light curve, before orbits with high MAD are excluded
  
    t_folded, orbit_number = fold_time(t, t0_ref, period_ref)
    
    params = lmfit.Parameters()
    params.add('period', value=period_ref, vary=False)
    params.add('t0', value=0.0, vary=True, min= -duration, max = duration)
    params.add('radratio', value=radratio, vary=True, min=0.0)
    params.add('a_over_r', value=a_over_r, vary=True)
    params.add('impactparam', value=impactparam, vary=True, min=0.0, max=1.0)
    params.add('limbdark1', value=limbdark1, vary=True, min=0.0, max=1.0)
    params.add('limbdark2', value=0, vary=False)
 
        
    print('\nFitting the folded light curve.\n')

     
    mini = lmfit.Minimizer(transit_lightcurve_residuals, params,
                           fcn_args=(t_folded, f, f_unc, cadence))
    
    out = mini.minimize(method=method)
    popt = out.params

    t0_revised = t0_ref + popt['t0']
    print(lmfit.fit_report(out))

    ##################################################################################
    # Timing each transit
    ##################################################################################
    # first, we exclude transit events that have MAD > 1.5 * overall MAD averaged over all events
    # We apply this filtering criteria only to light curves that have at least 3 transits present. 
 
    print('Before applying MAD mask')
    orbits = np.unique(orbit_number)
    no = len(orbits)
    t_obs = np.zeros((3,no))
    t_unc = np.zeros((3,no))

    n_cols = 3
    n_rows = np.int(np.round(no/n_cols))+1
    
    plt.axis('off')

    transit_residuals = []

    if no > 2: # check that there are at least 3 transits in the light curve
        i = 0
        for n in orbits:

            print('\nWorking on transit serial number ', i, ', orbit number ', n) 
            this_orbit = (orbit_number==n)
            tp = time_data[this_orbit]
            fp = flux_data[this_orbit]
            oot = (oibeo[this_orbit] > 0.)
  
            tmid = np.mean(tp[oot])
            coeff = detrend(tp[oot]-tmid, fp[oot], ndeg_max)
            f_baseline = np.polyval(coeff,tp-tmid)
            print('Intial de-trending coefficients: ', coeff)

            f_detrended = fp/f_baseline
            fp_uncertainty = np.std(f_detrended[oot])
            fp_unc = np.full(f_detrended.shape, fp_uncertainty)
 
            # estimate t0 and its uncertainty using lmfit minimizer
            t_n_guess = t0_revised + n*period_ref # initial guess for t0
            
            
            params = lmfit.Parameters()
            params.add('t0', value=t_n_guess, vary=True, min=t_n_guess-duration, max=t_n_guess+duration)
            params.add('period', value=popt['period'], vary=False)
            params.add('radratio', value=popt['radratio'].value, min=0.0, vary=False)
            params.add('a_over_r', value=popt['a_over_r'].value, vary=False)
            params.add('impactparam', value=popt['impactparam'].value, vary=False)
            params.add('limbdark1', value=popt['limbdark1'].value, vary=False)
            params.add('limbdark2', value=popt['limbdark2'].value, vary=False)
  
            
            mini = lmfit.Minimizer(transit_lightcurve_residuals, params,
                           fcn_args=(tp, f_detrended, fp_unc, cadence))

            out = mini.minimize(method=method)
            popt = out.params
            print('lmfit outputs t0: ', popt['t0'])

            t0_guess =  popt['t0'].value

            # estimate t0 and its uncertainty using MCMC
            t0_0, coefficients_0 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)
            t0_1, coefficients_1 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)
            t0_2, coefficients_2 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)
            t0_3, coefficients_3 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)
            t0_4, coefficients_4 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)

            t0_concatenated_chain = np.concatenate((t0_0, t0_1, t0_2, t0_3, t0_4))
            t0_chains = [t0_0, t0_1, t0_2, t0_3, t0_4]

            # coeffs from MCMC
            coeff_concatenated_chain = np.concatenate((coefficients_0, coefficients_1, coefficients_2, coefficients_3, coefficients_4), axis=0)
            coeff_mcmc = np.mean(coeff_concatenated_chain, axis = 0)
            print('Coeffs estimated with MCMC: ', coeff_mcmc)

            # t0 median in each chain
            medians = []
            for k in range(5):
                medians.append(np.median(t0_chains[k]))

            dispersion = np.std(np.array(medians))
            

            # calculate mean t0 and 1-sigma uncertaintiy in t0 from the concatenated chain
            t0_unc_mcmc = np.std(t0_concatenated_chain)  #t0_86th_percentile - t0_50th_percentile
            t0_mcmc = np.mean(t0_concatenated_chain)

            popt.add('t0', value = t0_mcmc, vary = True)
            
            print('Dispersion/t0_unc: ', dispersion/t0_unc_mcmc)
            print('Dispersion: ', dispersion)
            
            if dispersion < 0.5 * t0_unc_mcmc:
                print('MCMC converged!')
            else:
                print('MCMC did not converge!')

            print('MCMC results for t0 = ', t0_mcmc, t0_unc_mcmc)
            
            f_base = np.polyval(coeff_mcmc,tp-tmid)
            f_detrend= fp/f_base
            if int(cadence) == 120 or int(cadence) == 600:
                fp_calc = transit_lightcurve(popt, tp)
            else:
                fp_calc = transit_lightcurve_30min(popt, tp)
            transit_residuals.append(f_detrend-fp_calc)
 
            i=i+1


   
 
    ##################################################################################
    # Exclude orbits with high MAD
    ##################################################################################

    # flattern array of residuals to calculate MAD
    transit_residuals_flat = [item for sublist in transit_residuals for item in sublist]
    transit_residuals_flat = np.array(transit_residuals_flat)
    mad = np.median(np.absolute(transit_residuals_flat - np.median(transit_residuals_flat)))

    # check if MAD of a given transit > 1.5 overall mad defined on the previous line - this is likely bad data 
    # and we should not consider those orbits
    idx2delete = []
    for tr_num in range(len(transit_residuals)):
        tr_mad =  np.median(np.absolute(transit_residuals[tr_num] - np.median(transit_residuals[tr_num])))
        if tr_mad > 1.5 * mad:
            idx2delete.append(tr_num)
    orbits = np.delete(orbits, idx2delete)
    print('Indices of excluded orbits: ', idx2delete)


    #############################################################################################
    # Fit phase folded light curve that includes only the data that passes the filtering based on
    # 1) the minimum number of data points required
    # 2) having relatively low MAD
    #############################################################################################

    
    detrended_flux_stage2_ = [] # array for fluxes after bad events were removed
    time_stage2_ = []
    flux_uncertainty_stage2_ = []
    oibeo_stage2_ = []
    flux_data_stage2_ = []
 
    # orbits now include only events that passed the two filtering criteria
    for n in orbits:
        this_orbit = (orbit_number == n)
        detrended_flux_stage2_.append(initial_detrended_flux[this_orbit])
        time_stage2_.append(t[this_orbit]) 
        flux_uncertainty_stage2_.append(initial_flux_uncertainties[this_orbit]) 
        oibeo_stage2_.append(oibeo[this_orbit])
        flux_data_stage2_.append(flux_data[this_orbit])

    # flatten arrays
    detrended_flux_stage2_flat = [item for sublist in detrended_flux_stage2_ for item in sublist]
    detrended_flux_stage2 = np.array(detrended_flux_stage2_flat)
    time_stage2_flat = [item for sublist in time_stage2_ for item in sublist]
    time_stage2 = np.array(time_stage2_flat)
    flux_uncertainty_stage2_flat = [item for sublist in flux_uncertainty_stage2_ for item in sublist]
    flux_uncertainty_stage2 = np.array(flux_uncertainty_stage2_flat)
    oibeo_stage2_flat = [item for sublist in oibeo_stage2_ for item in sublist]
    oibeo_stage2 = np.array(oibeo_stage2_flat)
    flux_data_stage2_flat = [item for sublist in flux_data_stage2_ for item in sublist]
    flux_data_stage2 = np.array(flux_data_stage2_flat)
    
 
    # fold light curve and fit transit model
    t_folded2, orbit_number_stage2 = fold_time(time_stage2, t0_ref, period_ref)
    orbits_stage2 = np.unique(orbit_number_stage2)
 
    # transit model fit with quadratic limb-darkening law to folded time-series
    params = lmfit.Parameters()
    params.add('period', value=period_ref, vary=False)
    params.add('t0', value=0.0, vary=True, min= -duration, max = duration)
    params.add('radratio', value=radratio, vary=True, min=0.0)
    params.add('a_over_r', value=a_over_r, vary=True)
    params.add('impactparam', value=impactparam, vary=True, min=0.0, max=1.0)
    params.add('limbdark1', value=limbdark1, vary=True, min=0.0, max=1.0)
    params.add('limbdark2', value=limbdark2, vary=True, min=0.0, max=1)
    params.add('delta', expr='limbdark2+limbdark1', min = 0, max=1)
        
    print('\nFitting the folded light curve.\n')

    outlier_quadratic_ld = np.full(len(detrended_flux_stage2), False, dtype=bool)

    for i in range(niter_sigmaclip):

        n_outliers_before = np.sum(outlier_quadratic_ld)
        mini2 = lmfit.Minimizer(transit_lightcurve_residuals, params,
                           fcn_args=(t_folded2[~outlier_quadratic_ld], detrended_flux_stage2[~outlier_quadratic_ld], flux_uncertainty_stage2[~outlier_quadratic_ld], cadence))
    
        out2 = mini2.minimize(method=method)
        popt_quadratic_ld = out2.params
        f_calc = transit_lightcurve(popt_quadratic_ld, t_folded2)
        f_res =  detrended_flux_stage2 - f_calc
        oot = (~outlier_quadratic_ld) & (oibeo_stage2 > 1)
        scatter = np.std(f_res[oot])
        outlier_quadratic_ld = (np.abs(f_res/scatter) > nsigma_clip)
        print('   Sigma-clipping iteration ', i, ': total number of outliers clipped = ', np.sum(outlier_quadratic_ld))
        if (np.sum(outlier_quadratic_ld) == n_outliers_before): break

    print('transit model fit with quadratic limb-darkening law to folded time-series')
    print(lmfit.fit_report(out2))
    # transit model fit with linear limb-darkening law to folded time-series
    params = lmfit.Parameters()
    params.add('period', value=period_ref, vary=False)
    params.add('t0', value=0.0, vary=True, min= -duration, max = duration)
    params.add('radratio', value=radratio, vary=True, min=0.0)
    params.add('a_over_r', value=a_over_r, vary=True)
    params.add('impactparam', value=impactparam, vary=True, min=0.0, max=1.0)
    params.add('limbdark1', value=limbdark1, vary=True, min=0.0, max=1.0)
    params.add('limbdark2', value=0, vary=False)
 
        
    print('\nFitting the folded light curve.\n')

    outlier_linear_ld = np.full(len(detrended_flux_stage2), False, dtype=bool)

    for i in range(niter_sigmaclip):

        n_outliers_before = np.sum(outlier_linear_ld)
        mini2 = lmfit.Minimizer(transit_lightcurve_residuals, params,
                           fcn_args=(t_folded2[~outlier_linear_ld], detrended_flux_stage2[~outlier_linear_ld], flux_uncertainty_stage2[~outlier_linear_ld], cadence))
    
        out2 = mini2.minimize(method=method)
        popt_linear_ld = out2.params
        f_calc = transit_lightcurve(popt_linear_ld, t_folded2)
        f_res =  detrended_flux_stage2 - f_calc
        oot = (~outlier_linear_ld) & (oibeo_stage2 > 1)
        scatter = np.std(f_res[oot])
        outlier_linear_ld = (np.abs(f_res/scatter) > nsigma_clip)
        print('   Sigma-clipping iteration ', i, ': total number of outliers clipped = ', np.sum(outlier_linear_ld))
        if (np.sum(outlier_linear_ld) == n_outliers_before): break

    print('transit model fit with quadratic limb-darkening law to folded time-series')
    print(lmfit.fit_report(out2))
    # compare fits with linear vs linear limb-darkening law and select the one with lower chi-square
    # linear ld model
    f_calc_linear_ld = transit_lightcurve(popt_linear_ld, t_folded2)
    f_res_linear_ld =  detrended_flux_stage2 - f_calc_linear_ld
    oot_linear_ld = (~outlier_linear_ld) & (oibeo_stage2 > 1)
    scatter_linear_ld = np.std(f_res_linear_ld[oot_linear_ld])

    chi_square_linear_ld = np.sum(f_res_linear_ld[~outlier_linear_ld]**2)/scatter_linear_ld**2
    ndata = f_res_linear_ld[~outlier_linear_ld].shape[0]
    bic_linear_ld = chi_square_linear_ld + 5 * np.log(ndata)

    # quadratic ld model
    f_calc_quadratic_ld = transit_lightcurve(popt_quadratic_ld, t_folded2)
    f_res_quadratic_ld =  detrended_flux_stage2 - f_calc_quadratic_ld
    oot_quadratic_ld = (~outlier_quadratic_ld) & (oibeo_stage2 > 1)
    scatter_quadratic_ld = np.std(f_res_quadratic_ld[oot_linear_ld])
    chi_square_quadratic_ld = np.sum(f_res_quadratic_ld[~outlier_quadratic_ld]**2)/scatter_quadratic_ld**2
    ndata = f_res_quadratic_ld[~outlier_quadratic_ld].shape[0]
    bic_quadratic_ld = chi_square_quadratic_ld + 6 * np.log(ndata)
    # select model that gives lower BIC
    model_choice = 0 
    if  bic_linear_ld < bic_quadratic_ld:
        popt = popt_linear_ld
        outlier = outlier_linear_ld
        model_choice = 0 
    else:
        popt = popt_quadratic_ld
        outlier = outlier_quadratic_ld
        model_choice = 1
 
    print('Model chosen: ', model_choice)

    tp_fit2 = np.linspace(t_folded2.min(), t_folded2.max(), num=500)
    tr_model = transit_lightcurve(popt, tp_fit2)
    f_model = transit_lightcurve(popt, t_folded2) 

    if cadence == 120:
        brightness = 0.5
    else:
        brightness = 0.95
    plt.figure(figsize=(15,10))

    ax = plt.subplot(2,1,1)
    plt.plot(t_folded2, detrended_flux_stage2, 'k.',alpha=0.3)
    plt.plot(t_folded2[outlier], detrended_flux_stage2[outlier], '.', color='magenta', ms=15, alpha=brightness)
    plt.plot(tp_fit2, tr_model, 'r', linewidth=2)
    plt.ylabel('Normalized flux')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
 
    ax = plt.subplot(2,1,2)
    plt.plot(t_folded2, detrended_flux_stage2-f_model, 'k.', alpha=0.3)
    plt.plot(t_folded2[outlier], detrended_flux_stage2[outlier]-f_model[outlier], '.', color='magenta',ms=10, alpha=brightness)
    plt.xlabel('Folded time [days]')
    plt.ylabel('Residual flux')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.savefig(dirname + name + '_b_FoldedLightCurve.png',bbox_inches='tight')
    plt.close()




    ##################################################################################
    # Time each transit event (t0, t0_unc)
    ##################################################################################
    no = len(orbits)
    t_obs = np.zeros((3,no))
    t_unc = np.zeros((3,no))

    n_cols = 3
    n_rows = np.int(np.round(no/n_cols))+1
    
 

    print('Timing the transits after applying MAD mask')
    fig = plt.figure(figsize=(15,3*n_rows))
    plt.axis('off')
    outer = gridspec.GridSpec(n_rows, n_cols, figure=fig)
 
    i = 0
    for n in orbits_stage2:

        print('\nWorking on transit serial number ', i, ', orbit number ', n) 
        this_orbit = (orbit_number_stage2==n)
        tp = time_stage2[this_orbit & ~outlier]
        fp = flux_data_stage2[this_orbit & ~outlier]
        oot = (oibeo_stage2[this_orbit & ~outlier] > 0.)

        tmid = np.mean(tp[oot])
        coeff = detrend(tp[oot]-tmid, fp[oot], ndeg_max)
        f_baseline = np.polyval(coeff,tp-tmid)
        f_detrended = fp/f_baseline
        fp_uncertainty = np.std(f_detrended[oot])
        fp_unc = np.full(f_detrended.shape, fp_uncertainty)

        dt_unc_est = estimate_timing_uncertainty(params, len(tp), fp_unc[0])
        print('Theoretically estimated timing uncertainty [days, min]   = ', dt_unc_est, 24.*60.*dt_unc_est)
        
        t_n_guess = t0_revised + n*period_ref 

  
        params = lmfit.Parameters()
        params.add('t0', value=t_n_guess, vary=True, min=t_n_guess-duration, max=t_n_guess+duration)
        params.add('period', value=popt['period'], vary=False)
        params.add('radratio', value=popt['radratio'].value, min=0.0, vary=False)
        params.add('a_over_r', value=popt['a_over_r'].value, vary=False)
        params.add('impactparam', value=popt['impactparam'].value, vary=False)
        params.add('limbdark1', value=popt['limbdark1'].value, vary=False)
        params.add('limbdark2', value=popt['limbdark2'].value, vary=False)
 
        
        mini = lmfit.Minimizer(transit_lightcurve_residuals, params,
                       fcn_args=(tp, f_detrended, fp_unc, cadence))

        out = mini.minimize(method=method)
        #print(lmfit.fit_report(out))
        popt = out.params
        t_obs[0,i] = popt['t0'].value
        t_unc[0,i] = popt['t0'].stderr
        print('lmfit outputs t0: ', popt['t0'].value, popt['t0'].stderr)

 
        t0_guess =  popt['t0'].value
        print('1/40 * duration ', 1/40 * duration)

        # estimate t0 and its uncertainty using MCMC
        #print('fluxes inputted to MCMC', fp)
        t0_0, coefficients_0 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)
        t0_1, coefficients_1 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)
        t0_2, coefficients_2 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)
        t0_3, coefficients_3 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)
        t0_4, coefficients_4 = mcmc(popt, t0_guess + np.random.normal(0, 1/40 * duration), tp, fp, f_detrended, fp_unc, coeff, oot, tmid, n_links, cadence, duration)

        t0_concatenated_chain = np.concatenate((t0_0, t0_1, t0_2, t0_3, t0_4))
        t0_chains = [t0_0, t0_1, t0_2, t0_3, t0_4]

        # coeffs from MCMC
        coeff_concatenated_chain = np.concatenate((coefficients_0, coefficients_1, coefficients_2, coefficients_3, coefficients_4), axis=0)
        coeff_mcmc = np.mean(coeff_concatenated_chain, axis = 0)

 
        # t0 median in each chain
        medians = []
        for k in range(5):
            medians.append(np.median(t0_chains[k]))
            #print('median ', np.median(t0_chains[i]))

        dispersion = np.std(np.array(medians))
        

        # calculate mean t0 and 1-sigma uncertaintiy in t0 from the concatenated chain

        t0_50th_percentile = np.percentile(t0_concatenated_chain, 50)
        t0_86th_percentile = np.percentile(t0_concatenated_chain, 86)

        t0_unc_mcmc = np.std(t0_concatenated_chain)  #t0_86th_percentile - t0_50th_percentile
        t0_mcmc = np.mean(t0_concatenated_chain)

        print('t0_86th_percentile - t0_50th_percentile = ', t0_86th_percentile - t0_50th_percentile)
        print('t0_unc_mcmc = ', t0_unc_mcmc)

        #value1 = np.append(coeff_mcmc, t0_mcmc)
        print('Coeffs estimated with MCMC: ', coeff_mcmc)
 


        popt.add('t0', value = t0_mcmc, vary = True)
        
        print('Dispersion/t0_unc: ', dispersion/t0_unc_mcmc)
        print('Dispersion: ', dispersion)

     
        if dispersion < 0.5 * t0_unc_mcmc:
            mcmc_failed = 0
            print('MCMC converged')
        else:
            mcmc_failed = 1
            print('MCMC did not converge!!!')

        print('MCMC results for t0 = ', t0_mcmc, t0_unc_mcmc)
        t_obs[1,i] = t0_mcmc
        t_unc[1,i] = t0_unc_mcmc

         

        # Plot each transit
        tp_fit = np.linspace(tp.min(), tp.max(), num=500)

        if int(cadence) == 120 or int(cadence) == 600:
            fp_calc = transit_lightcurve(popt, tp_fit)
        else:
            fp_calc = transit_lightcurve_30min(popt, tp_fit)

        f_base = np.polyval(coeff_mcmc,tp-tmid)
        f_detrend = fp/f_base
 
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i])


        time_outlier = time_stage2[this_orbit & outlier]
        if time_outlier.shape[0] > 0:
            flux_outlier = flux_data_stage2[this_orbit & outlier]
            f_base = np.polyval(coeff_mcmc,time_outlier-tmid)
            flux_outlier_detrended = flux_outlier/f_base
            
        ax = plt.Subplot(fig, inner[0])
        ax.plot(tp-t0_mcmc, f_detrend,'k.')
        ax.plot(tp_fit-t0_mcmc, fp_calc,'r')
        if time_outlier.shape[0] > 0:
            ax.plot(time_outlier-t0_mcmc,flux_outlier_detrended,'.', color='magenta')

        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Orbit ' + str(np.int(n)))
        ax.set_ylim(min(f_detrend)-0.003*min(f_detrend), max(f_detrend)+0.003*max(f_detrend))
        fig.add_subplot(ax)

        if int(cadence) == 120 or int(cadence) == 600:
            fp_calc = transit_lightcurve(popt, tp)
            if time_outlier.shape[0] > 0:
                fp_calc_outlier = transit_lightcurve(popt, time_outlier)
        else:
            fp_calc = transit_lightcurve_30min(popt, tp)
            if time_outlier.shape[0] > 0:
                fp_calc_outlier = transit_lightcurve_30min(popt, time_outlier)

        ax = plt.Subplot(fig, inner[1])
        ax.plot(tp-t0_mcmc, f_detrend-fp_calc, 'k.')

        if time_outlier.shape[0] > 0:
            ax.plot(time_outlier-t0_mcmc, flux_outlier_detrended-fp_calc_outlier, '.', color='magenta')

        ax.set_ylim(-1.5*max(np.abs(f_detrend-fp_calc)), 1.5*max(np.abs(f_detrend-fp_calc)))
        ax.axhline(0,color='red')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])

        fig.add_subplot(ax) 
        i=i+1
        
    plt.savefig(dirname + name + '_b_IndividualTransitsWithFit.png')#,#bbox_inches='tight')
    plt.close()
    ##################################################################################



    ###################################################################################################
    # Plot time-series with transit events included in our analysis shown in red
    # the filtering for which events to include consists of 1) the minimum number of points requrement; 
    # 2) making sure MAD of a given transit event isn't too high.
    ###################################################################################################
 
    
    print('\nPlotting the time series.')
    if cadence == 120:
        brightness = 0.2
    else:
        brightness = 0.75
    
    fig, ax = plt.subplots()
    plt.plot(time,flux,'k.',ms=5,alpha=brightness)
    plt.plot(time_stage2,flux_data_stage2,'r.',ms=5,alpha=brightness)
    time_label = r'Time [BJD$_{\rm TDB} - 2457000$]'
    plt.xlabel(time_label)
    plt.ylabel(r'Flux [e$^{-}$/s]')
    plt.title(name.replace("_", " "))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.savefig(dirname + name + '_a_TimeSeries.png',bbox_inches='tight')
    plt.close()  
 


    ##################################################################################           
    # Transit timing deviations
    ##################################################################################
    print('\nCalculating transit timing deviations.')

    color = ['red','black','blue']
     
    i = 1 # t0 and its uncertainty estimated from MCMC

    ok = np.isfinite(t_unc[i,:])
    nok = np.sum(ok)
    tn = t_obs[i,ok]
    tn_unc = t_unc[i,ok]
    n = orbits[ok]

    bjdtdb = np.full(tn.shape[0], 'BJD_TDB')
    tesssrc = np.full(tn.shape[0], 'TESS')
    num_pts =  np.full(tn.shape[0], 1)

    # save transit times and their uncertainties in a .CSV file.
    df = pd.DataFrame({'Mid-point': tn + 2457000, 'Uncertainty': tn_unc, 'Time System': bjdtdb, 'Source': tesssrc, '#': num_pts})
    df.to_csv(dirname + name + '.csv', index = None)

   
    print('\nTimes with MCMC ', file=results_file)
    for j in range(len(tn)):
        print(tn[j] + 2457000.0, tn_unc[j], file=results_file)
    
    if (nok < 2.):
        print('\nNot enough timings to calculate ephemeris.')
        sys.stdout = original_stdout
        return mcmc_failed, rejected_orbits, len(idx2delete), model_choice, popt['limbdark1'].value, popt['limbdark2'].value

    
    print('  Uncertainty method ', i, ': number of valid transit times = ', nok)
    ephemeris = np.polyfit(n, tn, 1, w=1./tn_unc)
    tn_calc = np.polyval(ephemeris, n)
    chisqr = np.sum( ((tn-tn_calc)/tn_unc)**2 )
    print('Chisqr, Ndof, Ndata, Np = ', chisqr, nok-2, nok, 2)
    print('\nResult of linear fit to transit times using uncertainty method ', i, file=results_file)
    print('  Period = ' , ephemeris[0], file=results_file)
    print('  t0     = ' , ephemeris[1] + 2457000, file=results_file)
    print('  Chisqr, Ndof, Ndata, Nparam = ', chisqr, nok-2, nok, 2, file=results_file)
    
    #plt.figure()
    #ax = plt.figure()#.gca()
    fig, ax = plt.subplots()
    ax.errorbar(n+0.1*i,24.*60.*(tn-tn_calc),yerr=24.*60.*tn_unc,fmt='.',color=color[i],ms=10)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('Orbit Number')
    ax.set_ylabel(r'$\Delta t$ [min]')
    plt.savefig(dirname + name + '_c_TimingResiduals.png',bbox_inches='tight')
    plt.tight_layout()
    plt.close()

    i = 0 # t0 and its uncertainty estimated from LMFIT
    print(' Uncertainty estimated from LMFIT ', t_unc[i,:])
    ok = np.isfinite(t_unc[i,:])
    nok = np.sum(ok)
    tn = t_obs[i,ok]
    tn_unc = t_unc[i,ok]
    n = orbits[ok]

    print('\nTimes with lmfit optimizer ', file=results_file)
    for j in range(len(tn)):
        print(tn[j] + 2457000.0, tn_unc[j], file=results_file)
    
    if (nok < 2.):
        print('\nNot enough timings to calculate ephemeris.')
        sys.stdout = original_stdout
        return mcmc_failed, rejected_orbits, len(idx2delete), model_choice, popt['limbdark1'].value, popt['limbdark2'].value

    
    print('  Uncertainty method ', i, ': number of valid transit times = ', nok)
    ephemeris = np.polyfit(n, tn, 1, w=1./tn_unc)
    tn_calc = np.polyval(ephemeris, n)
    chisqr = np.sum( ((tn-tn_calc)/tn_unc)**2 )
    print('Chisqr, Ndof, Ndata, Np = ', chisqr, nok-2, nok, 2)
    print('\nResult of linear fit to transit times using uncertainty method ', i, file=results_file)
    print('  Period = ' , ephemeris[0], file=results_file)
    print('  t0     = ' , ephemeris[1] + 2457000, file=results_file)
    print('  Chisqr, Ndof, Ndata, Nparam = ', chisqr, nok-2, nok, 2, file=results_file)
    

    sys.stdout = original_stdout
    return mcmc_failed, rejected_orbits, len(idx2delete), model_choice, popt['limbdark1'].value, popt['limbdark2'].value

