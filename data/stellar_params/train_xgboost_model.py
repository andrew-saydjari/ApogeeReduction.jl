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
from model_utils import (data_loader, test_train_data, test_vs_prediction_plot,
                         kiel_diagram, alpha_fe_plot, compute_label_statistics)
plt.style.use('/uufs/chpc.utah.edu/common/home/u6033276/mpl_styles/mystyle.mplstyle')


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('xgboost_model.log'),
        logging.StreamHandler()  # Also print to console
    ]
)
logger = logging.getLogger(__name__)


def defaul_model_params():
    """
    default parameters used for the model
    """
    common_params = {
        "verbosity": 1,
        "n_jobs": 16,
        "learning_rate": 0.05,
        "n_estimators": 3000,
        "max_depth": 10,
        "subsample": 0.8,
        "colsample_bytree": 0.9,
        "colsample_bylevel": 0.8,
        "reg_alpha": 0.01,
        "reg_lambda": 1.0,
        "gamma": 0.5,
        "early_stopping_rounds": 150,
        "min_child_weight": 30,
    }
    return common_params


def train_model(X_train: np.ndarray,
                y_train: np.ndarray,
                X_test: np.ndarray,
                y_test: np.ndarray,
                label_names: list,
                save_dir: str,
                max_weight_perc: float = 80.,
                train_w_snr: bool = True,
                snr_weight: bool = True):
    """
    Model training step
    """
    models = {}
    predictions = np.zeros_like(y_test)

    if snr_weight:
        # create sample weights inversed by snr
        # hist, bin_edges = np.histogram(X_train[:, -1], bins=20)
        # bin_idx = np.digitize(X_train[:, -1], bin_edges[:-1])

        # weights = 1.0 / hist[bin_idx - 1]
        # weights /= np.mean(weights)  # normalize
        # weights = np.minimum(weights, max_weight)
        weights = X_train[:, -1] / np.nanmean(X_train[:, -1])
        max_weight = np.nanpercentile(weights, max_weight_perc)
        weights = np.minimum(weights, max_weight)
    else:
        weights = None

    # create directory to save
    os.makedirs(f'{save_dir}', exist_ok=True)
    os.makedirs(f'{save_dir}/plots', exist_ok=True)

    # save defaults used for run
    common_params = defaul_model_params()
    with open(f'{save_dir}/xbgoost_params.json', 'w') as f:
        json.dump(common_params, f, indent=4)
    logger.info(f"Running models with parameters: {str(common_params)}")

    for i, target_name in enumerate(label_names):
        logger.info(f"Training model for target {i} ({target_name})")
        y_train_i = y_train[:, i]
        y_val_i = y_test[:, i]

        model = xgb.XGBRegressor(**common_params)
        if train_w_snr:
            model.fit(
                X_train,
                y_train_i,
                eval_set=[(X_test, y_val_i)],
                sample_weight=weights
            )
            models[target_name] = model
            predictions[:, i] = model.predict(X_test)
        else:
            model.fit(
                X_train[:, :-1],
                y_train_i,
                eval_set=[(X_test[:, :-1], y_val_i)],
                sample_weight=weights
            )
            models[target_name] = model
            predictions[:, i] = model.predict(X_test[:, :-1])
        rmse_per_target = np.sqrt(np.nanmean((y_val_i - predictions[:, i]) ** 2))
        logger.info(f"RMSE for {target_name} = {rmse_per_target:.4f}")

        # save model
        model.save_model(f"{save_dir}/model_{target_name}.json")
    return models, predictions


if __name__ == '__main__':
    # load and prep the data
    starLineCof_file = 'arMADGICS_out_x_starLineCof_v0.h5'
    exp_file = 'exposures.h5'
    training_data_file = 'aspcap_training_data.npz'
    coeff_label = data_loader(starLineCof_file,
                              exp_file,
                              training_data_file,
                              logger)

    # get training and testing data
    label_names = ['teff', 'raw_alpha_m_atm', 'fe_h', 'logg']
    Ncoeff = 50
    train_w_snr = False
    snr_weight = True
    logger.info(f"Training with Ncoeff = {Ncoeff}, train_w_snr = {train_w_snr}, snr_weight = {snr_weight}")
    X_train, y_train, X_test, y_test = test_train_data(
        coeff_label, Ncoeff, label_names, stratify=False)

    # train the model
    save_dir = 'xgboost_models'
    models, predictions = train_model(X_train,
                                      y_train,
                                      X_test,
                                      y_test,
                                      label_names,
                                      save_dir,
                                      train_w_snr=train_w_snr,
                                      snr_weight=snr_weight)

    # print stats
    stats = compute_label_statistics(y_test, predictions, label_names)
    logger.info("Summary Statistics (Test Predictions)")
    for name in label_names:
        s = stats[name]
        logger.info(f"  {name:8s}: bias={s['bias']:+.4f}, scatter={s['scatter']:.4f}, MAD={s['mad']:.4f}")


    # save some nice plots to show results
    # 1-to-1 plots
    for i in range(len(label_names)):   
        test_vs_prediction_plot(y_test[:, i],
                                predictions[:, i],
                                label_names[i],
                                save_dir)
        test_vs_prediction_plot(y_test[:, i],
                                predictions[:, i],
                                label_names[i],
                                save_dir,
                                snr=X_test[:, -1])

    # kiel diagram
    kiel_diagram(y_test,
                 predictions,
                 label_names,
                 save_dir,
                 fe_h=False)

    # keil diagram with Fe/H
    kiel_diagram(y_test,
                 predictions,
                 label_names,
                 save_dir,
                 fe_h=True)

    # alpha/M vs Fe/H
    alpha_fe_plot(y_test,
                  predictions,
                  label_names,
                  save_dir)
