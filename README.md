# DEIMoS-Snakemake
I have made this to work with Sapelo2
to get this to work you must first 
ml Miniconda3/23.5.2-0
ml Anaconda3/2023.09-0
ml Mamba/23.1.0-4

Then activate your snakemake conda envirnment as in the snakemake tutorial.

When you want to get the deimos conda envirnment working you must use "mamba env create -f environment.yml". The envirnment does not build with regualr conda. It MUST be mamba.
to activate envirnment you must use "source activate deimos"
once envirnment is active you can then run this command "pip install -e ." (might only have to do this once.)

code to create a DAG: 
snakemake -s auto_qc.smk --dag | dot -Tsvg > dag.svg

Code to run program (2 cores so you dont load too much into memory at once)
snakemake -s auto_qc.smk --cores 2

Directory breakdown:
auto_qc.smk
input
--- qc (where qc files go)
---samples (where sample files go)
output
--- qc (where all qc outputs go)
--- samples (where you can find final outputs and intermediates)



Explination on how the program works:

1 Rule mzml2hdf: This rule uses the deimos.load() function on all mzml files to turn them into h5 files. Specifically it seems to save MS levels. This is for both qc and sample files. The accession is the variable used here and it can be changed in the config file under "accession". The output of this command is saved in {sample_type}/parsed

2.1 Rule factorize_qc: This rule will use the qc files output from mzml2hdf and use deimos.build_factors()

2.2 Rule threshold_qc: This rule uses the output of mzml2hdf and uses deimos.threshold() on it. Parameters can be changed in the config file under "threshold".

3 Rule smooth_qc: This rule uses the output of factorize_qc and threshold_qc to use the funvtion deimos.filters.smooth().

4 Rule peakpick_qc: This rule uses the output of factorize_qc and smooth_qc to run the function deimos.peakpick.presistant_homology().

5.Rule align_qc_qc: This rule uses the output of peakpick_qc and aligns the qc files to all qc files. Creates plots. Check them to ensure things are corrcect. Can be found in output/qc/qc_transform/plots

6 Rule align_qc_sample: This rule uses the sample output from mzml2hdf and align_qc_qc to find the closest date_time_index out of the qc files. It then applies that transformation on the sample and saves it. Essentially I think it is moving peaks over to match the sample to the qc file.

7.1 Rule factorize_sample: This rule uses the output of align_qc_sample to run deimos.build_factors() on it.

7.2 Rule threshold_sample: This rule uses the output of align_qc_sample to run deimos.threshold() on it. This uses the 'threshold' form the config file.

8 Rule smooth_sample: This rule uses the outputs of factorize_sample and threshold_sample to use the function deimos.filters.smooth().

9 Rule peakpick_sample: This Rule uses the output of factorize_sample and smooth_sample to run the function deimos.peakpick.persistent_homology().

10 Rule Downselect_peaks: This Rule uses the output of peakpick_sample to run deimos.multi_sample_partition() and partitions.map(deimos.alignment.agglomerative_clustering()). It then does a few other transformations and finally drops all but the most intense peaks.

11 Rule save_to_csv: This takes the downselected peak files. uses deimos.load to load them into a dictionary and saves the dictionary as a good csv file. Data can be found in output/samples/final_csvs
