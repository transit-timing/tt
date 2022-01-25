import sys
import os
import pandas as pd
import numpy as np
from astropy.io import fits
from tt import *

direct = os.path.dirname(os.getcwd())

def main():
    # specify planet name and TESS sector number
    planet_name = 'WASP-012' 
    sector = 43

    # open the DataFrame with the sample information
    table = pd.read_csv(direct + '/3_database/1_target_list.csv')
    ephemerides =  pd.read_csv(direct + '/3_database/table4.csv')
    df = pd.read_csv(direct + '/3_database/table2.csv')

    cadence =  df[(df['System'] == planet_name) & (df['Sector'] == sector)]['Cadence'].iloc[0] 
    sample = table[(table['System'].str.upper() == planet_name.upper())]

    t0 = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (BJD TDB)'].iloc[0]
    period = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['Period (days)'].iloc[0]

    hdul = fits.open(direct + f'/4_data/{planet_name}/{planet_name}_{sector}.fits')
    data = hdul[1].data
    t = data['TIME']
    flux = data['FLUX']
                                            
    params = (planet_name + f'_Sector_{sector}',
              period,
              t0,
              sample['duration'].iloc[0],
              sample['depth'].iloc[0]/100,
              cadence)
 
    mcmc_check, orbits_with_few_pts, orbits_with_high_mad, model_choice, ld1, ld2  = analyze_system(params, flux, t)
            
if __name__ == "__main__":
    main()
