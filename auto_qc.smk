from datetime import datetime
import deimos
import glob
import h5py
import matplotlib.pyplot as plt
import numpy as np
from os.path import *
import pandas as pd
import pymzml



def get_ionization_mode(path):
    acc = deimos.get_accessions(path)
    if 'positive scan' in acc.keys():
        return 'POS'
    elif 'negative scan' in acc.keys():
        return 'NEG'
    else:
        raise ValueError('Ionization mode could not be inferred from metadata.')


def get_datetime(path):
    return datetime.strptime(pymzml.run.Reader(path).info['start_time'],
                             '%Y-%m-%dT%H:%M:%SZ')


def find_closest_datetime_index(query_datetime, datetime_list):
    deltas = np.array([abs(dt - query_datetime) for dt in datetime_list])
    idx = np.argmin(deltas) # Finding where the minimum lies on the array deltas

    return idx


#Load in the config file
configfile: "config.yaml"
# Change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# Filenames
sample_fns = [basename(x) for x in sorted(glob.glob(join('input', 'samples', '*.*')))]
qc_fns = [basename(x) for x in sorted(glob.glob(join('input', 'qc', '*.*')))]
tune_fns = [basename(x) for x in sorted(glob.glob(join('input', 'tune', '*.*')))]

# Identifiers
sample_ids = [splitext(splitext(x)[0])[0] for x in sample_fns]
qc_ids = [splitext(splitext(x)[0])[0] for x in qc_fns]
tune_ids = [splitext(splitext(x)[0])[0] for x in tune_fns]

# Subdirectories
sample_dirs = ['samples' for x in sample_fns]
qc_dirs = ['qc' for x in qc_fns]
tune_dirs = ['tune' for x in tune_fns]

# Ionization modes
sample_modes = [get_ionization_mode(join('input', 'samples', x)) for x in sample_fns]
qc_modes = [get_ionization_mode(join('input', 'qc', x)) for x in qc_fns]
tune_modes = [get_ionization_mode(join('input', 'tune', x)) for x in tune_fns]

# Column types
sample_columns = ['RP' if '_rp_' in x.lower() else 'HILIC' for x in sample_fns]
qc_columns = ['RP'  if '_rp_' in x.lower() else 'HILIC' for x in qc_fns]

# Datetimes
qc_dts = [get_datetime(join('input', 'qc', x)) for x in qc_fns]
tune_dts = [get_datetime(join('input', 'tune', x)) for x in tune_fns]

# Lookup dictionaries
fn_lookup = {k: v for k, v in zip(sample_ids + qc_ids + tune_ids,
                                  sample_fns + qc_fns + tune_fns)}
sample_type_lookup = {k: v for k, v in zip(sample_ids + qc_ids + tune_ids,
                                           sample_dirs + qc_dirs + tune_dirs)}
mode_lookup = {k: v for k, v in zip(sample_ids + qc_ids + tune_ids,
                                    sample_modes + qc_modes + tune_modes)}                                
column_lookup = {k: v for k, v in zip(sample_ids + qc_ids,
                                      sample_columns + qc_columns)}

# Collect all outputs
rule all:
    input:
        # OG ouptut
        # expand(join('output', '{sample_type}', 'downselected', '{id}.h5'),
        #        id=sample_ids, sample_type='samples')
        expand(join('output', '{sample_type}', 'final_csvs', '{id}.csv'),
                id=sample_ids, sample_type='samples')


rule mzml2hdf:
    input:
        lambda wildcards: join('input', sample_type_lookup[wildcards.id], fn_lookup[wildcards.id])
    output:
        join('output', '{sample_type}', 'parsed', '{id}.h5')
    run:
        # Read/parse mzml
        data = deimos.load(input[0], accession=config['accession'])

        # Enumerate MS levels
        for k, v in data.items():
            # Save as hdf5
            deimos.save(output[0], v, key=k, mode='a')


# Build factors
rule factorize_qc:
    input:
        rules.mzml2hdf.output
    output:
        join('output', '{sample_type}', 'factors', '{id}.npy')
    wildcard_constraints:
        sample_type='qc'
    run:
        # Get keys
        keys = list(h5py.File(input[0], 'r').keys())

        factors = {}
        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[0], key=k, columns=config['dims'] + ['intensity'])

            # Build factors
            factors[k] = deimos.build_factors(data, dims=config['dims'])

        # Save factors
        np.save(output[0], factors)


# Threshold data by intensity
rule threshold_qc:
    input:
        rules.mzml2hdf.output
    output:
        join('output', '{sample_type}', 'thresholded', '{id}.h5')
    wildcard_constraints:
        sample_type='qc'
    run:
        # Get keys
        keys = list(h5py.File(input[0], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[0], key=k, columns=config['dims'] + ['intensity'])

            # Threshold
            data = deimos.threshold(data, threshold=config['threshold'])

            # Save
            deimos.save(output[0], data, key=k, mode='a')   


# Smooth data
rule smooth_qc:
    input:
        rules.factorize_qc.output,
        rules.threshold_qc.output
    output:
        join('output', '{sample_type}', 'smoothed', '{id}.h5')
    wildcard_constraints:
        sample_type='qc'
    run:
        # Load factors
        factors = np.load(input[0], allow_pickle=True).item()

        # Get keys
        keys = list(h5py.File(input[1], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[1], key=k, columns=config['dims'] + ['intensity'])

            # Perform smoothing
            data = deimos.filters.smooth(data,
                                         factors=factors[k],
                                         dims=config['dims'],
                                         iterations=config['smooth']['iters'],
                                         radius=config['smooth']['radius'])

            # Save
            deimos.save(output[0], data, key=k, mode='a')


# Perform peak detection
rule peakpick_qc:
    input:
        rules.factorize_qc.output,
        rules.smooth_qc.output
    output:
        join('output', '{sample_type}', 'peakpicked', '{id}.h5')
    wildcard_constraints:
        sample_type='qc'
    run:
        # Load factors
        factors = np.load(input[0], allow_pickle=True).item()

        # Get keys
        keys = list(h5py.File(input[1], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[1], key=k, columns=config['dims'] + ['intensity'])

            # Perform peakpicking
            peaks = deimos.peakpick.persistent_homology(data,
                                                        factors=factors[k],
                                                        dims=config['dims'],
                                                        radius=config['peakpick']['radius'])

            # Save
            deimos.save(output[0], peaks, key=k, mode='a')


# Align quality control samples
rule align_qc_qc:
    input:
        source=rules.peakpick_qc.output,
        targets=expand(join('output', '{sample_type}', 'peakpicked', '{id}.h5'),
                       id=qc_ids, sample_type='qc')
    output:
        join('output', '{sample_type}', 'qc_transform', '{id}.npy'),
        join('output', '{sample_type}', 'qc_transform', 'plots', '{id}.png')
    wildcard_constraints:
        sample_type='qc'
    run:
        # Metadata
        # Remember wildcards.id refers to {id}
        column_type = column_lookup[wildcards.id]
        ionization_mode = mode_lookup[wildcards.id]

        # Possible targets
        qc_targets = [x for i, x in enumerate(input.targets)
                      if (qc_columns[i] == column_type)
                      & (qc_modes[i] == ionization_mode)]
        
        # Get keys
        keys = list(h5py.File(input.source[0], 'r').keys())

        # Figure init
        fig, ax = plt.subplots(len(config['dims']), len(keys),
                                   dpi=300, figsize=(3 * len(config['dims']), 3 * len(keys)),
                                   facecolor='w')

        # Enumerate MS levels
        transform = {}
        for i, k in enumerate(keys):
            # Plot title
            ax[0, i].set_title(k)

            # MS level
            transform[k] = {}

            # Load peak data
            source_peaks = deimos.load(input.source[0], key=k)
            target_peaks = deimos.load(qc_targets[0], key=k)

            # Threshold
            source_peaks = deimos.threshold(source_peaks,
                                            threshold=config['qc']['threshold'])
            target_peaks = deimos.threshold(target_peaks,
                                            threshold=config['qc']['threshold'])

            # Match
            a, b = deimos.alignment.tolerance(source_peaks, target_peaks,
                                              dims=config['dims'],
                                              tol=config['qc']['tol'],
                                              relative=config['qc']['relative'])
            print(f"Value of 'a': {a}")
            print(f"Value of 'b': {b}")
            print(f"i is: {keys}")
            # Fit alignment
            for j, dim in enumerate(config['dims']):
                transform[k][dim] = deimos.alignment.fit_spline(a, b, align=dim, kernel='linear')

                # Plot scatter
                ax[j, i].scatter(a[dim], b[dim], s=4, c='C0')

                # Plot fit
                newx = np.linspace(0, a[dim].max(), 1000)
                ax[j, i].plot(newx, transform[k][dim](newx),
                              c='black', linewidth=1, linestyle='--')

                # Axis labels
                ax[j, i].set_xlabel(dim, fontweight='bold')
                ax[j, i].set_ylabel(dim, fontweight='bold')

        plt.tight_layout()
        plt.savefig(output[1])
        plt.close()

        # Save transform
        np.save(output[0], transform)


rule align_qc_sample:
    input:
        rules.mzml2hdf.output,
        expand(join('output', '{sample_type}', 'qc_transform', '{id}.npy'),
               id=qc_ids, sample_type='qc')
    output:
        join('output', '{sample_type}', 'qc_aligned', '{id}.h5')
    wildcard_constraints:
        sample_type='samples'
    run:
        # Get datetime of sample
        dt = get_datetime(join('input', 'samples', fn_lookup[wildcards.id]))

        # Metadata
        column_type = column_lookup[wildcards.id]
        ionization_mode = mode_lookup[wildcards.id]

        # Find relevant QCs
        # TODO: Function to return relevant IDs and indices?
        relevant_qc_indices = [i for i, x in enumerate(qc_dts)
                               if (qc_columns[i] == column_type)
                               & (qc_modes[i] == ionization_mode)]

        # Find index of closest datetime
        idx = find_closest_datetime_index(dt, [qc_dts[i] for i in relevant_qc_indices])

        # Load transform for closest QC
        # TODO: Does expand on input preserve order? If so, can replace below.
        transform = np.load(join('output', 'qc', 'qc_transform', qc_ids[relevant_qc_indices[idx]] + '.npy'),
                            allow_pickle=True).item()
        
        # Get keys
        keys = list(h5py.File(input[0], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[0], key=k, columns=config['dims'] + ['intensity'])

            # Apply alignment
            for dim in config['dims']:
                data[dim] = transform[k][dim](data[dim].values)
            
            # Save
            deimos.save(output[0], data, key=k)


# Build factors
rule factorize_sample:
    input:
        rules.align_qc_sample.output
    output:
        join('output', '{sample_type}', 'factors', '{id}.npy')
    wildcard_constraints:
        sample_type='samples'
    run:
        # Get keys
        keys = list(h5py.File(input[0], 'r').keys())

        factors = {}
        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[0], key=k, columns=config['dims'] + ['intensity'])

            # Build factors
            factors[k] = deimos.build_factors(data, dims=config['dims'])

        # Save factors
        np.save(output[0], factors)


# Threshold data by intensity
rule threshold_sample:
    input:
        rules.align_qc_sample.output
    output:
        join('output', '{sample_type}', 'thresholded', '{id}.h5')
    wildcard_constraints:
        sample_type='samples'
    run:
        # Get keys
        keys = list(h5py.File(input[0], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[0], key=k, columns=config['dims'] + ['intensity'])

            # Threshold
            data = deimos.threshold(data, threshold=config['threshold'])

            # Save
            deimos.save(output[0], data, key=k, mode='a')   


# Smooth data
rule smooth_sample:
    input:
        rules.factorize_sample.output,
        rules.threshold_sample.output
    output:
        join('output', '{sample_type}', 'smoothed', '{id}.h5')
    wildcard_constraints:
        sample_type='samples'
    run:
        # Load factors
        factors = np.load(input[0], allow_pickle=True).item()

        # Get keys
        keys = list(h5py.File(input[1], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[1], key=k, columns=config['dims'] + ['intensity'])

            # Perform smoothing
            data = deimos.filters.smooth(data,
                                         factors=factors[k],
                                         dims=config['dims'],
                                         iterations=config['smooth']['iters'],
                                         radius=config['smooth']['radius'])

            # Save
            deimos.save(output[0], data, key=k, mode='a')

# Perform peak detection

rule peakpick_sample:
    input:
        rules.factorize_sample.output,
        rules.smooth_sample.output

    output:
        join('output', '{sample_type}', 'peakpicked', '{id}.h5')
    wildcard_constraints:
        sample_type='samples'
    run:
        # Load factors
        factors = np.load(input[0], allow_pickle=True).item()

        # Get keys
        keys = list(h5py.File(input[1], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[1], key=k, columns=config['dims'] + ['intensity'])

            # Perform peakpicking
            peaks = deimos.peakpick.persistent_homology(data,
                                                        factors=factors[k],
                                                        dims=config['dims'],
                                                        radius=config['peakpick']['radius'])

            # Save
            deimos.save(output[0], peaks, key=k, mode='a')

# This seems to keep peaks that are more intense
# Merge close peaks
rule downselect_peaks:
    input:
        rules.peakpick_sample.output
    output:
        join('output', '{sample_type}', 'downselected', '{id}.h5')
    wildcard_constraints:
        sample_type='samples'
    run:
        # Get keys
        keys = list(h5py.File(input[0], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            peaks = deimos.load(input[0], key=k)

            # Partition data
            partitions = deimos.multi_sample_partition(peaks,
                                                       split_on=config['downselect']['partition']['split_on'],
                                                       size=config['downselect']['partition']['size'],
                                                       tol=config['downselect']['partition']['tol'])

            # Map clustering
            peaks = partitions.map(deimos.alignment.agglomerative_clustering,
                                   dims=[x + '_weighted' for x in config['dims']],
                                   tol=config['downselect']['cluster']['tol'],
                                   relative=config['downselect']['cluster']['relative'])

            # Reindex over clusters, partitions
            peaks['cluster'] = peaks.groupby(by=['partition_idx', 'cluster']).ngroup().reset_index(drop=True)
    
            # Drop all but most intense
            peaks = peaks.sort_values(by='intensity', ascending=False).drop_duplicates(subset='cluster')
            peaks = peaks.drop(columns=['partition_idx', 'cluster'])
            
            # Save
            deimos.save(output[0], peaks, key=k, mode='a')


# Saves the output to an easy to read csv file
rule save_to_csv:
    input:
        rules.downselect_peaks.output
    output:
        join('output', '{sample_type}', 'final_csvs', '{id}.csv')
    wildcard_constraints:
        sample_type='samples'
    run:
        # Get keys
        keys = list(h5py.File(input[0], 'r').keys())

        # Enumerate MS levels
        for k in keys:
            # Load data
            data = deimos.load(input[0], key=k, columns=config['dims'] + ['intensity'])
        
        # Add column names
        column_names = config['dims'] + ['intensity']

        # Convert the list of column names to a comma-separated string
        column_names_str = ','.join(column_names)
        header = f"{column_names_str}"

        # save data to numpy data frame
        np.savetxt(output[0], data, delimiter=',', header=header, comments='')
