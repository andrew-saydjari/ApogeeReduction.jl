import h5py
import matplotlib.pylab as plt
import numpy as np
from astropy.table import Table, join
import logging
from sklearn.model_selection import train_test_split
import numpy.lib.recfunctions as rfn
import xgboost as xgb
import json
import os
from matplotlib.colors import LogNorm
plt.style.use('/uufs/chpc.utah.edu/common/home/u6033276/mpl_styles/mystyle.mplstyle')


def data_loader(starLineCof_file: str,
                exp_file: str,
                training_data_file: str,
                logger):
    """
    Load the data for model training
    """
     # load data
    coeff_hf = h5py.File(starLineCof_file, 'r')  # star coefficient
    exp_hf = h5py.File(exp_file, 'r')  # exposure data
    save_data = np.load(training_data_file)['save_data']  # training data

    # select on matches and high snr
    ev_aspcap = np.isin(exp_hf['sdss_id'][:],
                        save_data['sdss_id'][save_data['sdss_id'] != -1])
    snr = exp_hf['snr'][:]
    snr_min = 5.
    snr_max = np.nanpercentile(snr, 99.)
    logger.info(f"Clipping SNR values to {snr_min:.2f} < SNR < {snr_max:.2f}")
    snr_good = (snr > snr_min) & (snr < snr_max)
    sdss_id_good = exp_hf['sdss_id'][:] != -1
    ev_labels = ev_aspcap & snr_good & sdss_id_good
    labelsi = np.where(ev_labels)[0]

    # grab the coefficients
    coeff_aspcap = np.zeros((len(labelsi), coeff_hf['x_starLineCof_v0'].shape[1]))

    logger.info(f"Reading in Coefficients")
    nchunk = 10000
    for i in range(0, len(labelsi), nchunk):
        coeff_aspcap[i: i+nchunk] = coeff_hf['x_starLineCof_v0'][labelsi[i: i+nchunk]]
    sdss_id_aspcap = exp_hf['sdss_id'][ev_labels]
    snr_aspcap = exp_hf['snr'][ev_labels]


    # clip extreme coefficients
    # lower, upper = np.nanpercentile(coeff_aspcap, [0.1, 99.9], axis=0)
    # outlier_mask = np.any((coeff_aspcap < lower) | (coeff_aspcap > upper), axis=1)
    # keep_idx = np.where(~outlier_mask)[0]

    # merge all of the data
    logger.info(f"Merging Data")
    label_tab = Table(save_data)
    coeff_tab = Table({'starLineCof': coeff_aspcap,
                       'sdss_id': sdss_id_aspcap,
                       'snr': snr_aspcap})
    coeff_label = join(coeff_tab, label_tab, keys='sdss_id')

    # make sure none have sdss_id == -1 !!!
    coeff_label = coeff_label[coeff_label['sdss_id'] != -1]
    logger.info(f"Have {len(coeff_label)} spectra to use for training")
    return coeff_label


def test_train_data(coeff_label: Table,
                    Ncoeff: int,
                    label_names: list,
                    stratify: bool = True):
    """
    create test train split
    """
    # stratify based on snr
    if stratify:
        snr_bins = np.nanquantile(coeff_label['snr'], np.linspace(0, 1, 21))
        snr_binned = np.digitize(coeff_label['snr'], snr_bins[1:-1])
    else:
        snr_binned = None
    X_train, X_test, y_train, y_test = train_test_split(
        np.column_stack((coeff_label['starLineCof'][:, :Ncoeff], coeff_label['snr'])),
        rfn.structured_to_unstructured(coeff_label[label_names].as_array()),
        test_size=0.2, random_state=42, stratify=snr_binned
    )

    # remove data with nans
    ev_train = np.all(np.isfinite(X_train), axis=1)
    ev_test = np.all(np.isfinite(X_test), axis=1)

    X_train = X_train[ev_train]
    y_train = y_train[ev_train]

    X_test = X_test[ev_test]
    y_test = y_test[ev_test]
    
    return X_train, y_train, X_test, y_test


def test_vs_prediction_plot(y_test: np.ndarray,
                            predictions: np.ndarray,
                            label_name: str,
                            save_dir: str,
                            snr: np.ndarray = None):
    """
    make a plot of test vs predictions
    """
    if snr is None:
        plt.figure(figsize=(12, 10))
        bins = np.linspace(y_test.min(), y_test.max(), 100)
        plt.hist2d(y_test, predictions, bins=bins, norm=LogNorm(), cmap='inferno')
        plt.colorbar(label='N')
        plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()],
                '--', c='r')
        plt.grid()
        rmse_per_target = np.sqrt(np.nanmean((y_test - predictions) ** 2))
        plt.title(f"{label_name}: RMSE = {rmse_per_target:.3f}")
        plt.xlabel('Test Data')
        plt.ylabel('Predictions')
        plt.savefig(f'{save_dir}/plots/test_predictions_{label_name}.png')
        plt.close()
    else:
        plt.figure(figsize=(12, 10))
        binsx = np.linspace(snr.min(), snr.max(), 100)
        maxx = np.nanpercentile(abs(y_test - predictions), 95)
        binsy = np.linspace(-maxx, maxx, 100)
        plt.hist2d(snr, y_test - predictions, bins=[binsx, binsy],
                   norm=LogNorm(), cmap='inferno')
        plt.colorbar(label='N')
        plt.axhline(0, linestyle='-', c='r')
        plt.grid()
        rmse_per_target = np.sqrt(np.nanmean((y_test - predictions) ** 2))
        plt.title(f"{label_name}: RMSE = {rmse_per_target:.3f}")
        plt.ylabel('Test Data - Predictions')
        plt.xlabel('SNR')
        plt.savefig(f'{save_dir}/plots/test_predictions_vs_snr_{label_name}.png')
        plt.close()


def kiel_diagram(y_test: np.ndarray,
                 predictions: np.ndarray,
                 label_names: list,
                 save_dir: str,
                 fe_h=False):
    """
    make a kiel diagram for test and predicted
    """
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 10))

    binsx = np.linspace(4000, 7500, 100)
    binsy = np.linspace(0, 5, 100)

    if fe_h:
        H_weights, xedges, yedges = np.histogram2d(y_test[:, label_names.index('teff')],
                                                y_test[:, label_names.index('logg')],
                                                bins=[binsx, binsy],
                                                weights=y_test[:, label_names.index('fe_h')])

        H_counts, _, _ = np.histogram2d(y_test[:, label_names.index('teff')],
                                        y_test[:, label_names.index('logg')],
                                        bins=(xedges, yedges))

        weighted_average = H_weights / H_counts


        res = ax1.imshow(weighted_average.T, origin='lower', aspect='auto',
                         extent=(binsx.min(), binsx.max(), binsy.min(), binsy.max()), cmap='inferno',
                        vmin=-1, vmax=0.3)
        plt.colorbar(res, label='[Fe/H]', ax=ax1)
    else:
        res = ax1.hist2d(y_test[:, label_names.index('teff')], y_test[:, label_names.index('logg')],
                        bins=[binsx, binsy], norm=LogNorm(), cmap='inferno')
        plt.colorbar(res[-1], label='N', ax=ax1)
    ax1.grid()
    ax1.set_title('Test Data')
    ax1.set_xlabel('Teff')
    ax1.set_ylabel('log(g)')
    ax1.invert_xaxis()
    ax1.invert_yaxis()

    if fe_h:
        H_weights, xedges, yedges = np.histogram2d(predictions[:, label_names.index('teff')],
                                                predictions[:, label_names.index('logg')],
                                                bins=[binsx, binsy],
                                                weights=predictions[:, label_names.index('fe_h')])

        H_counts, _, _ = np.histogram2d(predictions[:, label_names.index('teff')],
                                        predictions[:, label_names.index('logg')], bins=(xedges, yedges))

        weighted_average = H_weights / H_counts


        res = ax2.imshow(weighted_average.T, origin='lower', aspect='auto',
                         extent=(binsx.min(), binsx.max(), binsy.min(), binsy.max()), cmap='inferno',
                        vmin=-1, vmax=0.3)
        plt.colorbar(res, label='[Fe/H]', ax=ax2)
    else:
        res = ax2.hist2d(predictions[:, label_names.index('teff')], predictions[:, label_names.index('logg')],
                        bins=[binsx, binsy],
                        norm=LogNorm(vmin=res[3].norm.vmin, vmax=res[3].norm.vmax),
                        cmap='inferno')
        plt.colorbar(res[-1], label='N', ax=ax2)
    ax2.grid()
    ax2.set_title('Predictions')
    ax2.set_xlabel('Teff')
    ax2.set_ylabel('log(g)')
    ax2.invert_xaxis()
    ax2.invert_yaxis()
    if fe_h:
        plt.savefig(f'{save_dir}/plots/kiel_diagram_fe_h.png')
    else:
        plt.savefig(f'{save_dir}/plots/kiel_diagram.png')
    plt.close()


def alpha_fe_plot(y_test: np.ndarray,
                  predictions: np.ndarray,
                  label_names: list,
                  save_dir: str):
    """
    Make plot of alpha/M vs Fe/H for test and predicted
    """
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 10))

    binsx = np.linspace(-2, 0.5, 100)
    binsy = np.linspace(-0.4, 0.4, 100)

    res = ax1.hist2d(y_test[:, label_names.index('fe_h')],
                     y_test[:, label_names.index('raw_alpha_m_atm')],
                     bins=[binsx, binsy], norm=LogNorm(), cmap='inferno')
    plt.colorbar(res[-1], label='N', ax=ax1)
    ax1.grid()
    ax1.set_title('Test Data')
    ax1.set_xlabel('Fe/H')
    ax1.set_ylabel('alpha/M')

    res = ax2.hist2d(predictions[:, label_names.index('fe_h')],
                     predictions[:, label_names.index('raw_alpha_m_atm')],
                     bins=[binsx, binsy],
                     norm=LogNorm(vmin=res[3].norm.vmin, vmax=res[3].norm.vmax),
                     cmap='inferno')
    plt.colorbar(res[-1], label='N', ax=ax2)
    ax2.grid()
    ax2.set_title('Predictions')
    ax2.set_xlabel('Fe/H')
    ax2.set_ylabel('alpha/M')
    plt.savefig(f'{save_dir}/plots/alpha_m_vs_fe_h.png')
    plt.close()


def compute_label_statistics(true_labels, inferred_labels, label_names):
    """Compute bias, scatter, and MAD for each label."""
    stats = {}
    for i, name in enumerate(label_names):
        diff = inferred_labels[:, i] - true_labels[:, i]
        valid = np.isfinite(diff)
        bias = np.median(diff[valid])
        scatter = np.std(diff[valid])
        mad = np.median(np.abs(diff[valid] - bias))
        stats[name] = {'bias': bias, 'scatter': scatter, 'mad': mad}
    return stats
