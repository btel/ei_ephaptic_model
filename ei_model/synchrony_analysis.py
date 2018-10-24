from scipy import signal
import numpy as np

def calc_correlation_matrix(solution, threshold):
    corr_mat = np.corrcoef((solution>threshold).T)
    ij = np.arange(corr_mat.shape[0])
    corr_mat[ij, ij] = 0
    return corr_mat
    
def plot_correlation_matrix(corr_mat, corr_mean):

    corr_max = np.abs(corr_mat).max()
    fig_corr_ring = plt.figure(facecolor='none')
    ax = plt.subplot(111)
    mappable = ax.matshow(corr_mat, cmap=plt.cm.BrBG_r, vmin=-corr_max, vmax=corr_max)
    cb = plt.colorbar(mappable)
    cb.set_label('correlation coefficient')
    ax.text(0.9, 0.8, '$\\bar{{r}}={:.3f}$'.format(corr_mean), size=15,
              transform=ax.transAxes, ha='right')
    return fig_corr_ring

def calc_pli(solution, dt, n_subsamp=100, filter_order=5, cutoff_frequency=1):
    """Calculate the phase locking index (PLI)"""
    # subsample signal
    y = solution[::n_subsamp, :].copy()
    
    # filter signal
    sampling_frequency = 1 / (dt * n_subsamp)
    b, a = signal.iirfilter(filter_order,  cutoff_frequency * 2 / sampling_frequency, btype='lowpass')
    y_filtered = signal.filtfilt(b, a, y[:], axis=0)
    phases = np.angle(signal.hilbert(y_filtered, axis=0))
    
    # calculate pairwise phase locking indice
    pli = np.abs(np.mean(np.exp((phases[:, :, None] - phases[:, None, :])*1j), axis=0))
    
    # zero out the diagonal for better visualisation
    ij = np.arange(pli.shape[0])
    pli[ij, ij] = np.NaN
    return pli

def plot_pli_matrix(pli, corr_mean):

    fig_corr_ring = plt.figure(facecolor='none')
    ax = plt.subplot(111)
    mappable = ax.matshow(pli)
    cb = plt.colorbar(mappable)
    cb.set_label('phase locking index')
    ax.text(0.9, 0.8, '$\\bar{{\Phi}}={:.3f}$'.format(corr_mean), size=15,
              transform=ax.transAxes, ha='right')
    return fig_corr_ring
