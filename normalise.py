import pandas
import numpy as np


def normalise(infile, outfile, min_to_back=True):
    """In our dataset the fluorescence sometimes goes below the
    autofluorescence, for comparison we take the minimum
    value of fluorescence as the background. This function performs
    the normalisation taking the minimum flourescence as the
    background fluorescence, if and only if min_to_back is
    True. Otherwise, the normalisation is performed using only the
    autofluorescence plasmid. """
    df = pandas.read_csv(infile, index_col=0)
    print("starting normalisation")

    if min_to_back:
        df = pandas.merge(df,
                          df.groupby([
                              'backbone', 'strain'
                          ]).volume_decomposed_log_mean_gfp.min().to_frame(),
                          on=["strain", "backbone"],
                          how='outer',
                          suffixes=("", "_min"))
    else:
        df = pandas.merge(df,
                          df.query('plasmid=="1201"').groupby([
                              'backbone', 'strain'
                          ]).volume_decomposed_log_mean_gfp.mean().to_frame(),
                          on=["strain", "backbone"],
                          how='outer',
                          suffixes=("", "_min"))
    df = pandas.merge(df,
                      df.query('plasmid=="1717"').groupby([
                          'backbone', 'strain'
                      ]).volume_decomposed_log_mean_gfp.mean().to_frame(),
                      on=["strain", "backbone"],
                      how='outer',
                      suffixes=("", "_standard"))
    df["rrpu"] = ((np.exp(df.volume_decomposed_log_mean_gfp) -
                   np.exp(df.volume_decomposed_log_mean_gfp_min)) /
                  (np.exp(df.volume_decomposed_log_mean_gfp_standard) -
                   np.exp(df.volume_decomposed_log_mean_gfp_min)))
    compactdf = df.copy()
    compactdf = compactdf.drop(columns=[
        'filename', "date", "real_time", "log_mean_v_mean",
        "volume_decomposed_log_mean_gfp_min",
        "volume_decomposed_log_mean_gfp_standard"])
    compactdf.to_csv(outfile)
    return compactdf
