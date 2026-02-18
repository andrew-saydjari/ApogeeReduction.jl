import h5py
import matplotlib.pylab as plt
import numpy as np
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.cm as cm


if __name__ == '__main__':
    # load the RV valid data
    # based on HARPS data where RV stable == sigma_rv < 5 m/s
    res = np.load('rv_stable_data.npz')
    tri_ids = res['tri_ids']  # sdss_ids for targets
    tri_rv_data = res['tri_rv_data']  # mean RVs from Trifonov, et. al 2020

    # load the exposures file
    exp_hf = h5py.File('../data/exposures.h5', 'r')

    # grab columns and setup flags for good data to plot
    sdss_id = exp_hf['sdss_id'][:]
    v_rad_flags = exp_hf['v_rad_flags'][:]
    v_rad_good = v_rad_flags == 0  # cut on rv flag
    snr = exp_hf['snr'][:]
    snr_good = snr > 40.  # only high SNR
    MG = exp_hf['gaia_phot_g_mean_mag'][:] + 5 * np.log10(1e-3 * exp_hf['gaia_parallax'][:]) + 5
    MG_good = (MG > 3.26) & (MG < 6.85)  # choose sources roughly mid F to mid K
    mjd = exp_hf['mjd'][:]
    mjd_bad = (mjd > 59850) & (mjd < 60050)  # know bad rvs in this region

    ev_good = v_rad_good & snr_good & MG_good & ~mjd_bad

    obs = np.char.decode(exp_hf['observatory'][:], encoding='UTF-8')

    # create the plots
    for obsi in ['apo', 'lco']:
        ev_valid = np.isin(tri_ids, sdss_id[ev_good & (obs == obsi)])

        cmap = mpl.colormaps['inferno']
        vmin_val = exp_hf['adjfiberindx'][obs == obsi].min()
        vmax_val = exp_hf['adjfiberindx'][obs == obsi].max()
        norm = mcolors.Normalize(vmin=vmin_val, vmax=vmax_val)
        mappable = cm.ScalarMappable(norm=norm, cmap='inferno')

        f, ax1 = plt.subplots(1, 1, figsize=(16, 8))
        all_rvs = []
        for rv_sdss_id, rv in zip(tri_ids[ev_valid], tri_rv_data[ev_valid]):
            ev = (sdss_id == rv_sdss_id) & ev_good & (obs == obsi)
            mjd_ev = exp_hf['mjd'][ev]
            y_ev = exp_hf['v_rad'][ev] - rv
            yerr_ev = exp_hf['e_v_rad'][ev]
            fib_ev = exp_hf['adjfiberindx'][ev]
            for ax in [ax1]:
                ax.scatter(mjd_ev, y_ev, marker='.', c=fib_ev, cmap=cmap, norm=norm)
                for x_i, y_i, e_i, f_i in zip(mjd_ev, y_ev, yerr_ev, fib_ev):
                    ax.errorbar(x_i, y_i, yerr=e_i, fmt='none', capsize=2, ecolor=cmap(norm(f_i)))
            all_rvs += list(exp_hf['v_rad'][ev] - rv)
        for ax in [ax1]:
            ax.grid()
            ax.set_ylabel('RV (HARPS Mean Subtracted) [km/s]')
            ax.set_xlabel('MJD')
            f.colorbar(mappable, ax=ax, label='adjusted_fiber_index')
        all_rvs = np.array(all_rvs)

        ax1.set_title(r'%s, $\sigma_{RV}$ = %.3f m/s' % (obsi.upper(), np.std(all_rvs) * 1000))
        ax1.set_ylim(-1.1 * np.max(abs(all_rvs)), 1.1 * np.max(abs(all_rvs)))

        plt.savefig(f'manual_mwm_validation_rv_apogee_rv_stable_{obsi}.png',
                    dpi=200, bbox_inches='tight')
        plt.close()
