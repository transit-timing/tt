import numpy as np
from mandelagol import occultquad
 
 
def fold_time(t, t0, period):
    """    
    Returns folded time in between -period/2 and +period/2
    as well as an orbit number, assigning t0 = orbit zero
    """  
    print(np.divmod(t-(t0-0.5*period), period))

    orbit_number, t_fold = np.divmod( t-(t0-0.5*period), period )
    t_fold = t_fold - 0.5*period
    return(t_fold, orbit_number.astype(np.int))


def detrend(t, f, deg_max):
    """
    De-trend time-series with a polynomial.
    The degree of the polynomial is either 
    1, 2, or 3, with the decision based on 
    minimizing the Bayesian Information 
    Criterion (BIC).
    """ 
    ndata = len(t)
    coeff = np.polyfit(t,f,1)
    f_calc = np.polyval(coeff,t)
    f_res = f-f_calc
    sigma = np.std(f_res)
    bic = np.sum((f_res/sigma)**2) + 1.0*np.log(ndata)
    
    if (deg_max > 1):
        for n in range(2,deg_max+1):
            coeff_trial = np.polyfit(t,f,n)
            f_calc = np.polyval(coeff_trial,t)
            f_res = f-f_calc
            bic_trial = np.sum((f_res/sigma)**2) + n*np.log(ndata)
            if (bic_trial < bic): coeff = coeff_trial
        
    return coeff
        
def transit_lightcurve(params, t):
    """
    Compute transit model with parameters params 
    and evaluate it at times t. 
    """ 
    t0 = params['t0']
    period = params['period']
    u1 = params['limbdark1']
    u2 = params['limbdark2']
    k = params['radratio']
    b = params['impactparam']
    a = params['a_over_r']
 

    phi = 2.*np.pi*(t-t0)/period
    x = a*np.sin(phi)
    y = b*np.cos(phi)
    s = np.sqrt(x*x + y*y)
    f_calc, f_no_ld = occultquad(s,u1,u2,k)
    
    return f_calc

def transit_lightcurve_30min(params, t):
    """
    Compute transit model with parameters params. 
    The model is averaged over 30-minute bins and 
    evaluated at times t. 
    """    
    ts = np.linspace(t.min() - 0.015, t.max() + 0.015, num = t.shape[0] * 100)
    f_calc = transit_lightcurve(params, ts)
    avg_f = np.zeros(t.shape[0])

    for i in range(t.shape[0]):
        bound1 = t[i] + 0.01
        bound2 = t[i] - 0.01
        m1 = ts < bound1
        m2 = ts > bound2
        mask = m1 & m2
        avg_f[i] = np.mean(f_calc[mask])
    return avg_f


def transit_lightcurve_residuals(params, t, f, f_unc, cadence):
    """
    Compute transit light curve residuals
    """
    if int(cadence) == 120:       
        f_calc = transit_lightcurve(params, t)
        res = (f-f_calc)/f_unc
    else:
        avg_f = transit_lightcurve_30min(params, t)
        res = (f-avg_f)/f_unc        
    return res


def transit_lightcurve_chisq(params, t, f, f_unc, cadence):
    """
    Compute chi-square between light curve time-series and transit model 
    """
    if int(cadence) == 120:  
        f_calc = transit_lightcurve(params, t)
        res = (f-f_calc)/f_unc
    else:
        avg_f = transit_lightcurve_30min(params, t)
        res = (f-avg_f)/f_unc            
    return np.sum(res*res)


 
def estimate_timing_uncertainty(params, n, sigma):
    """
    Estimate timing uncertainty based on theoretical model.
    """
    k = params['radratio']
    b = params['impactparam']
    if (b>0.9): b=0.9
    a = params['a_over_r']
    per = params['period']
    durtot = per/a/np.pi*np.sqrt(1.0-b*b)
    Q = k*k/sigma*np.sqrt(n)
    durpar = 2.*k*per/a/np.pi/np.sqrt(1.0-b*b)
    theta = 2.*k/(1.0 - b*b)
    dt_unc_est = durtot/Q * np.sqrt(theta/2.)

    return dt_unc_est


def mcmc(params, t0_guess, t, f, f_detrended, f_unc, coeffs, oot, tmid, n_links, cadence, duration):
    """
    Run Markov Chain Monte Carlo methods
    """   

    print('Starting MCMC with n_links = ', n_links)
    dt_unc_est = estimate_timing_uncertainty(params, len(t), f_unc[0])
    print('Theoretically  estimated error = ', dt_unc_est)
    
    t0 = np.zeros(n_links)
     
    coefficients = np.zeros((n_links, coeffs.shape[0]))
    for j in range(coeffs.shape[0]):
        coefficients[0, j] = coeffs[j]
    t0[0] = t0_guess
    print('starting MCMC with t0 = ', t0_guess)
    chisq = np.zeros(n_links)
    chisq[0] = transit_lightcurve_chisq(params, t, f_detrended, f_unc, cadence)
    n_accept=0
    for i in range(1,n_links):
        t0_trial = t0[i-1] + np.random.normal(0, 1/40 * duration)
        coeffs_trial = coefficients[i-1, :] + np.random.normal(0,0.00025, size = coeffs.shape[0])

        f_baseline = np.polyval(coeffs_trial, t-tmid)
        fp = f/f_baseline
        scatter = np.std(fp[oot])
        f_unc = scatter
        params.add('t0', value = t0_trial, vary = True)
 
        chisq_trial = transit_lightcurve_chisq(params, t, fp, f_unc, cadence)
        del_chisq = chisq_trial-chisq[i-1]
        if (del_chisq < 0):
            n_accept = n_accept + 1
            chisq[i] = chisq_trial
            t0[i] = t0_trial
            for j in range(coeffs.shape[0]):
                coefficients[i, j] = coeffs_trial[j]
        else:
            prob = np.exp(-0.5*del_chisq)
            random_number = np.random.uniform(0,1)
            if (random_number <= prob):
                n_accept = n_accept + 1
                t0[i] = t0_trial
                chisq[i] = chisq_trial
                for j in range(coeffs.shape[0]):
                    coefficients[i, j] = coeffs_trial[j]
            else:
                chisq[i] = chisq[i-1]
                t0[i] = t0[i-1]
                for j in range(coeffs.shape[0]):
                    coefficients[i, j] = coefficients[i-1, j]


    print('Done, acceptance rate 0 = ', n_accept/n_links)
    
    return t0, coefficients 


 