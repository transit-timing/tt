import lightkurve as lk
import matplotlib.pyplot as plt
import os 


# specify planet's name and which light curve to download 
target = 'NGTS-11'
row = 4
 
direct = os.path.dirname(os.getcwd())
dirname = direct + '/4_data/'+ target

if not os.path.exists(dirname):
  os.mkdir(dirname)

result = lk.search_lightcurve(f'{target}', mission = 'TESS')
sector = result[row].mission[0] 
 
s = int(sector.replace('TESS Sector ', ''))
print(f'{target} Sector: {s}')
lc = result[row].download()
lc.to_fits(dirname + f'/{target}_{s}.fits')

t = lc.time.to_value('jd', 'long')
f = lc.flux

plt.plot(t-2457000.0 , f, '.k')
plt.xlabel(r'Time [BJD$_{\rm TDB} - 2457000$]')
plt.ylabel(r'Flux [e$^{-}$/s]')
plt.show()


