# DEIMoS-Snakemake
I have made this to work with Sapelo2
to get this to work you must first 
ml Miniconda3/23.5.2-0
ml Anaconda3/2023.09-0
ml Mamba/23.1.0-4


When you want to get the deimos conda envirnment working you must use "mamba env create -f environment.yml". 
to activate envirnment you must use "source activate deimos"
once envirnment is active you can then run this command "pip install -e ." (might only have to do this once.


Explination on how the program works:

1 Rule mzml2hdf: This rule uses the deimos.load() function on all mzml files to turn them into h5 files. Specifically it seems to save MS levels. This is for both qc and sample files. The accession is the variable used here and it can be changed in the config file under "accession"

2.1 Rule factorize_qc: This rule will use the qc files output from mzml2hdf and use deimos.build_factors()

2.2 Rule threshold_qc: This rule uses the output of mzml2hdf and uses deimos.threshold() on it. Parameters can be changed in the config file under "threshold".

3 Rule smooth_qc: This rule uses the output of factorize_qc and threshold_qc to use the funvtion deimos.filters.smooth().

4 Rule peakpick_qc: This rule uses the output of factorize_qc and smooth_qc to run the function deimos.peakpick.presistant_homology().

5.Rule align_qc_qc: This rule uses the output of peakpick_qc and aligns the qc files to the other qc files (seems to also include itself). Creates plots. Check them to ensure things are corrcect. Can be found in output/qc/qc_transform/plots

6 Rule align_qc_sample: This rule uses the sample output from mzml2hdf and align_qc_qc to find the closest date_time_index out of the qc files. It then applies that transformation on the sample and saves it. 

7.1 Rule actorize_sample:

7.2 Rule threshold_sample:

8 Rule smooth_sample:

9 Rule peakpick_sample:

10 Rule Downselect_peaks:

11 Rule tba: This takes the downselected peak files. uses deimos.load to load them into a dictionary and saves the dictionary as a good csv file. Data can be found in output/samples/final_csvs
