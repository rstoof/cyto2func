import numpy as np
import pandas
import fcsparser
from datetime import datetime
import scipy.stats as st
from scipy.optimize import curve_fit
import multiprocessing

# fitting to binormal##########################################################


# functions
def kdedata(data):
    """Makes a 2D histogram of log transformed data"""
    [fl, vl] = data
    fluormin = np.percentile(fl, 5)
    fluormax = np.percentile(fl, 95)
    volmin = np.percentile(vl, 5)
    volmax = np.percentile(vl, 95)
    xx, yy = np.mgrid[fluormin:fluormax:160j, volmin:volmax:160j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([fl, vl])
    kernel = st.gaussian_kde(values)

    cpus = multiprocessing.cpu_count()
    splits = np.array_split(positions, cpus, axis=1)
    pool = multiprocessing.Pool(processes=cpus)
    partial_results = pool.map(kernel.pdf, splits)
    results = np.concatenate(partial_results)
    f = np.reshape(results.T, xx.shape)
    pool.terminate()
    return [xx, yy, f]


def fit2binormal(kdedat):
    """Fits a 2D-normal distribution to 2D histogram of log transformed data"""
    [xx, yy, f] = kdedat
    xy_mesh = [xx, yy]
    guess_vals = [np.max(f), 3, 3, 2, 2, .8]

    # perform the fit, making sure to flatten the noisy data for the fit routine
    fit_params, cov_mat = curve_fit(normal_2d,
                                    xy_mesh,
                                    np.ravel(f),
                                    p0=guess_vals,
                                    maxfev=100000,
                                    bounds=([0, -3, -3, 0, 0,
                                             .0], [50, 15, 15, 5, 5, 1]))
    # calculate fit parameter errors from covariance matrix
    fit_errors = np.sqrt(np.diag(cov_mat))
    # manually calculate R-squared goodness of fit
    fit_residual = f - normal_2d(xy_mesh, *fit_params).reshape(
        np.outer(xx[:, 0], yy[0]).shape)
    fit_Rsquared = 1 - np.var(fit_residual) / np.var(f)
    return np.append(fit_params, fit_Rsquared)


def normal_2d(xy_mesh, amp, xc, yc, sigma_x, sigma_y, rho):
    """Defines a 2D normal with covariance rho"""
    # unpack 1D list into 2D x and y coords
    (x, y) = xy_mesh
    # make the 2D Gaussian matrix
    gauss = amp * np.exp(-(
        (x - xc)**2 / (sigma_x**2) - (2 * rho * (x - xc) * (y - yc)) /
        (sigma_x * sigma_y) + (y - yc)**2 /
        (sigma_y**2))) / (2 * np.pi * sigma_x * sigma_y * np.sqrt(1 - rho**2))
    # flatten the 2D Gaussian down to 1D
    return np.ravel(gauss)


def gate_data(data):
    return data[(data["SSC_H"] > np.exp(2.5)) &
                (data["SSC_A"] > 0) &
                (data["FSC_H"] > np.exp(1.5)) &
                (data["GFP_H"] > 0) &
                (data["SSC_A"] > data["SSC_H"]) &
                (data["SSC_H"] > data["FSC_H"])]


def decompose(out_file, desc_file, data_directory, channels, gate=True):
    """Perform the decomposition and fitting, writing the results to
    csv files"""
    df = pandas.read_csv(desc_file)
    datearr = []
    droparr = []
    fitarr = []

    print("Start fitting to bi log normal distribution")
    print("This may take some time...")

    for index, row in df.iterrows():
        try:
            path = data_directory + row.filename
            meta, data = fcsparser.parse(path, reformat_meta=True)
            data.columns = [x.strip().replace('-', '_') for x in data.columns]

            if gate:
                data = gate_data(data)

            date = datetime.strptime(f"{meta['$DATE']} {meta['$BTIM']}", '%Y-%b-%d %H:%M:%S')
            datearr.append(date)
            data = data.query(f"{channels[0]}>0 and {channels[1]}>0")
            data2 = kdedata([np.log(data[channels[1]]), np.log(data[channels[0]])])
            fits = fit2binormal(data2)
            fitarr.append(fits)
            print(index, fits)
        except FileNotFoundError:
            print(f"Failed to find file {data_directory + row.filename}")
            droparr.append(index)

    df = df.drop(droparr)
    df.insert(6, 'real_time', datearr)

    fitmufl = [fit[1] for fit in fitarr]
    fitmuv = [fit[2] for fit in fitarr]
    fitstdfl = [fit[3] for fit in fitarr]
    fitstdv = [fit[4] for fit in fitarr]
    fitrho = [fit[5] for fit in fitarr]
    fitgoodness = [fit[6] for fit in fitarr]

    df.insert(7, "log_mean_gfp", fitmufl)
    df.insert(8, "log_mean_v", fitmuv)
    df.insert(9, "log_std_gfp", fitstdfl)
    df.insert(10, "log_std_v", fitstdv)
    df.insert(11, "log_rho", fitrho)
    df.insert(12, "fit_goodness", fitgoodness)
    df.insert(13, "std_gfp_correct",
              np.sqrt(1 - np.power(np.array(df.log_rho), 2)) * np.array(df.log_std_gfp))
    df2 = pandas.merge(df,
                       df.groupby(['backbone',
                                   'strain']).log_mean_v.mean().to_frame(),
                       on=["strain", "backbone"],
                       how='outer',
                       suffixes=("", "_mean"))

    df2["volume_decomposed_log_mean_gfp"] = (
        df2["log_mean_gfp"] * df2["log_std_v"] - df2["log_std_gfp"] *
        df2["log_mean_v"] * df2["log_rho"] + df2["log_std_gfp"] *
        df2["log_mean_v_mean"] * df2["log_rho"]) / df2["log_std_v"]

    df2.to_csv(out_file)
    return df2
