


# TESS Transit Timing of 278 Planets [![DOI](https://zenodo.org/badge/452065386.svg)](https://zenodo.org/badge/latestdoi/452065386)


This repository contains source code for the *"TESS Transit Timing of 278 Planets"* paper. The repository provides code (1) to process *TESS* light curves (with 2, 10, or 30-minute sampling) to extract transit times, and (2) to automatically download arXiv articles that contain a given planet's name and a number greater than 2,000,000. A detailed description of the code can be found in the paper. 


## Processing *TESS* light curves 
The code relies on the following Python packages:  numpy, pandas, scikit-learn, astropy, matplotlib, batman-package, scipy, lmfit, lxml, arxiv. 
You may install the packages via pip or, if you use anaconda, you may run the following commands:
```
module load anaconda
conda create -n tt python=3.7 numpy pandas scikit-learn astropy matplotlib batman-package scipy lightkurve lmfit lxml arxiv
conda activate tt
```
Then clone this repository to your machine:

```
git clone https://github.com/transit-timing/transit-timing.github.io.git
```
 To process *TESS* light curves stored in the *~/data/* folder and extract mid-transit times, first go to the directory containing the code:

```
cd 0_tt
```
Inside *main.py* script, you may specify the name of the planet whose *TESS* light curve you would like to process. Then, you can run the script:
```
python3 main.py
```
## Output
The results are stored in the *~/output* folder. A typical folder that this script produces for a given target contains the following:

1. Figure of the  *TESS* light curve that was processed (transits that were processed to extract mid-transit times are shown in red):

 
<img src="/2_figures/WASP_012/WASP-012_Sector_20_a_TimeSeries.png" alt="drawing" width="600"/>


2. Figure with the individual transits that were processed (best-fit transit model is shown in red):

 
<img src="/2_figures/WASP_012/WASP-012_Sector_20_b_IndividualTransitsWithFit.png" alt="drawing" width="600"/>


4. Folded light curve:

 
<img src="/2_figures/WASP_012/WASP-012_Sector_20_b_FoldedLightCurve.png" alt="drawing" width="600"/>

 

 

5. O-C diagram showing timing residuals. The timing residuals were calculated as follows: first, a linear model was fitted to the transit times extracted from *TESS* data. Then, the linear model was subtracted from the transit times at each epoch:

 
<img src="/2_figures/WASP_012/WASP-012_Sector_20_c_TimingResiduals.png" alt="drawing" width="600"/>

6. *planet_name_results.txt* contains extracted mid-transit times and their uncertainties.
7. *planet_name_log.txt* contains intermediate outputs of the code, such as the found best-fit transit model parameters and different statistics of the produced fits.


**Note:** in order to process *TESS* light curves, you may want to first download the relevant *.fits* files. To download a given light curve, you may do the following:

```
cd 0_download_data
```
In the  *download_single_target.py* script, specify the name of the target and the light curve that you would like to download. Then run the script:
```
python3 download_single_target.py
```


## Figures

O-C diagrams for each of the 382 systems are available [here](https://github.com/transit-timing/transit-timing/tree/main/2_figures/o_c).

## Tables

1. A single CSV file with all of the collected transit times for each target is available at [https://github.com/transit-timing/transit-timing/blob/main/3_database/table4.csv](https://github.com/transit-timing/transit-timing/blob/main/3_database/table3.csv)
2. A single CSV file with all of the target ephemerides is available at [https://github.com/transit-timing/transit-timing/blob/main/3_database/table3.csv](https://github.com/transit-timing/transit-timing/blob/main/3_database/table3.csv)
# Scraping literature
A semi-automated procedure was used to search the extensive literature on transiting planets to find all the previously reported times for all the systems in our sample. The details of our approach to scraping the data can be found in Section 4 of the paper. The basic steps are as follows.
1. Download bibcodes of all papers related to a given object from the NASA Astrophysics Data System database:
    - Go to https://ui.adsabs.harvard.edu/search;
    - Do cone search using RaDec coordinates with 3" radius. 
For example, you may execute the following search:    
```
object:"20h30m54.13s 6d25m46.4s:3"
```
Alternatively, you can type *object:"planet_name"* in the search bar of the NASA Astrophysics Data System database to retrieve articles for a planet with a given name *planet_name*. For example, you may type:
```
object:"WASP-12"
```
Then click "Export" -> "in BibTeX" and download the .bib file containing all of the papers referencing a given planet. Place the downloaded .bib file into *~/arxiv/article_database* folder and rename the file to whatever the name of the planet is (e.g. WASP-12.bib)

2. Run **bibcode2arxiv_id.py** to create a .txt file containing a list of arXiv IDs of papers related to a given object whose name is *planet_name*. The arXiv IDs are extracted from the .bib file downloaded in the previous step:
```
cd 1_scrape_literature
```
```  
python3 bibcode2arxiv_id.py -p planet_name
```
For example, if you downloaded a .bib file for WASP-12, then you would run
```  
python3 bibcode2arxiv_id.py -p WASP-12
```
3. Execute **run.py** to download *.tex* source files of arXiv papers retrieved in the previous step. The script untars the source files and transfer all .tex files to a different directory.
``` 
python3 -u run.py -p WASP-12
``` 
For a given planet, all tar.gz are stored in ~/arxiv_dir/{planet_name} directory.
For a given planet, all untarred files are stored in ~/untar_dir/{planet_name} directory.
4. Run **read_tex.py** to find numbers > 2,000,000; if such numbers and the name of the planet are present in a given paper, a PDF of the paper is transferred to ~/{planet_name} directory for a final manual review.

```
python3 read_tex.py -s wasp -id 12
```


# Citation
If you find the code useful, consider citing it:
[![DOI](https://zenodo.org/badge/452065386.svg)](https://zenodo.org/badge/latestdoi/452065386)

