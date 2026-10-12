#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
import time
import multiprocessing as mp
import os
import warnings
import sys
import argparse


# -----------------------------
# Functions
# -----------------------------

def get_num_processes():
    # mp.cpu_count() returns the NODE's total core count, not what this
    # job was actually allocated -- on a shared HPC node under SLURM, that
    # means spawning far more worker processes than reserved cores,
    # oversubscribing CPU and slowing every worker down (or contending with
    # other jobs on the same node). os.sched_getaffinity(0) reports the
    # actual scheduling affinity mask (respects cgroup/SLURM cpuset
    # restrictions on Linux); SLURM_CPUS_PER_TASK is a fallback for
    # environments where affinity isn't cgroup-restricted directly.
    # sched_getaffinity is Linux-only (AttributeError on macOS/Windows),
    # hence the fallback chain down to cpu_count() for local/dev use.
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        pass
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        return int(slurm_cpus)
    return mp.cpu_count()


def calculate_f_for_row_chunk(row_chunk_data, groups):
    f_stats_chunk = []
    grouped_data = {}
    for key in groups:
        grouped_data[key] = row_chunk_data[groups[key]].values

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore') # This will silence warnings about rows that cannot generate f-stats
        f_stats_chunk = f_oneway(*grouped_data.values(), nan_policy="omit", axis=1)

    return f_stats_chunk.statistic


#-----------------------------
# Import args
#-----------------------------

parser = argparse.ArgumentParser(description="Calculate F-statistics on a mapped methylation file.")

parser.add_argument("--map", dest="map_file", type=str, required=True,
                    help="Path to the input mapped methylation file (TSV)")
parser.add_argument("--manifest", dest="manifest_file", type=str, required=True,
                    help="Path to the manifest file (TXT)")
parser.add_argument("--out_file", dest="output_file", type=str, required=True,
                    help="Path to the output file for F-statistics (TSV)")
args = parser.parse_args()

data_path = args.map_file
manifest_path = args.manifest_file
output_path = args.output_file

data = pd.read_csv(data_path, sep="\t", header=0)
manifest = pd.read_csv(manifest_path, sep="\t", header=0)

#------------------------------
# Main script
#------------------------------

if __name__ == "__main__":
    #print version of scipy (must be 1.15.0 or higher)
    print("scipy version: ", sys.modules['scipy'].__version__)

    start_time = time.time()

    # Get the list of all data columns once
    all_data_columns = data.columns.tolist()

    groups = {}
    for i in range(len(manifest)):
        if manifest.iloc[i, 0] not in groups:
            groups[manifest.iloc[i, 0]] = []
        groups[manifest.iloc[i, 0]].append(manifest.iloc[i, 1])

    # Fail fast, before spawning any worker process, if the manifest
    # references sample columns that don't actually exist in --map --
    # previously this surfaced as an opaque KeyError deep inside a
    # multiprocessing worker (mp.Pool.starmap), which is a much harder
    # failure to diagnose than a clear message up front.
    required_cols = sorted({col for cols in groups.values() for col in cols})
    missing_cols = [c for c in required_cols if c not in all_data_columns]
    if missing_cols:
        raise ValueError(
            f"--manifest references {len(missing_cols)} column(s) not present in --map: {missing_cols}\n"
            f"Columns present in --map: {all_data_columns}"
        )

    # --- Parallel Processing ---
    num_processes = get_num_processes()
    chunk_size = 100000 # Define a suitable chunk size (experiment with this)
    num_rows = len(data)
    chunks = []

    # Split the data into chunks
    for i in range(0, num_rows, chunk_size):
        chunk_indices = range(i, min(i + chunk_size, num_rows))
        chunk_data = data.iloc[chunk_indices] # Get numpy array for data chunk
        chunks.append((chunk_data, groups)) # Pass chunk with groups as a tuple

    print(f"Splitting {num_rows} rows into {len(chunks)} chunks, using {num_processes} process(es)...")

    # Use a multiprocessing Pool
    with mp.Pool(processes=num_processes) as pool:
        results = pool.starmap(calculate_f_for_row_chunk, chunks)

    # Flatten the list of results from chunks
    fstats = [item for sublist in results for item in sublist]
    fstats = np.array(fstats) # Convert back to numpy array

    end_time = time.time()
    print("Parallel F-statistic calculation time: ", end_time - start_time, " seconds")

    # --- Post-processing (Keep as in your original script) ---
    num_nans = np.sum(np.isnan(fstats))
    print("Number of NAs in F-statistic: ", num_nans)

    # Inf F-statistics arise from zero within-group variance with nonzero
    # between-group variance -- the theoretically CLEANEST possible marker,
    # not a degenerate case to discard. The original script never counted
    # these, so there was no visibility into whether they exist or whether
    # a downstream finite-value filter silently drops them.
    num_infs = np.sum(np.isinf(fstats))
    print("Number of Infs in F-statistic (zero within-group variance -- the strongest possible markers, verify these survive any downstream finite-value filter): ", num_infs)

    fstats_df = pd.DataFrame(fstats)
    fstats_df = fstats_df.round(4)

    #Merge the first three columns of data with fstats
    fstats_df = pd.concat([data.iloc[:, :3], fstats_df], axis=1)

    # Save the F-statistic to a file
    fstats_df.to_csv(output_path, sep="\t", index=False, header=False)
