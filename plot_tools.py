"""
Functions used by LPNE for visualization.

For now this really just has one useful function: factor_gridplot, which
allows you to generate an 'upper triangular' plot of an electome factor.
"""
__date__ = "March 2021"


from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np
import os
from re import findall

from data_tools import load_data, feature_mat



def factor_gridplot(diags, offdiags1, offdiags2, areaList=None, freqs=1, \
    labels=('',)):
    """
    Plots a factor in an upper triangular grid format

    Parameters
    ----------
    diags : numpy.ndarray
        Contains features that should be plotted in subfigures along the
        diagonal (e.g. PSD). Shape: [AxF] where A = # of areas; F = # of
        frequencies
    offdiags1 : numpy.ndarray
        Contains first set of features to be plotted in off-diagonal subfigures.
        Should have same scale as diag features. Shape: [AxAxF]
    offdiags2 : numpy.ndarray
        Contains second set of features to be plotted in off-diagonal
        subfigures. Can have different scaling than diag features.
        Shape: [AxAxF]
    areaList : list of str
        List of areas that features are generated from.
    freqs : list or int
        List of frequencies or scalar frequency sampling spacing (Hz). ...
    labels: tuple
        Tuple of 2 y-axis text labels. First label is applied to diags and
        offdiags1. Second label applies to offdiags2.
    """
    YLIM_SCALE = 1.15
    R_LABEL_SIZE = 12
    R_LABEL_X_OFF = -27
    R_LABEL_Y_OFF = 0.5

    A,F = diags.shape

    # calculate ylim values for plots
    max_val1 = np.max([np.max(abs(diags)), np.max(abs(offdiags1))])
    max_val2 = np.max(abs(offdiags2))
    ylims1 = YLIM_SCALE*np.asarray([0, 1])
    ylims2 = YLIM_SCALE*np.asarray([-max_val2/max_val1, max_val2/max_val1])

    # generate 'upper triangular' grid of subplots
    diag_ax = np.full(A, None, dtype=object)
    offd_ax1 = np.full((A,A), None, dtype=object)
    offd_ax2 = np.full((A,A), None, dtype=object)
    fig1 = plt.figure(constrained_layout=False)
    spec1 = GridSpec(ncols=A, nrows=A, figure=fig1)
    for k in range(A):
        diag_ax[k] = fig1.add_subplot(spec1[k,k])
    for k in range(A):
        for l in range(k+1,A):
            offd_ax1[k,l] = fig1.add_subplot(spec1[k,l])
            offd_ax2[k,l] = offd_ax1[k,l].twinx()

    # convert frequencies to array if scalar frequency resolution given
    if np.isscalar(freqs):
        freqs = np.arange(1,F+1,freqs)

    # fill diagonal plots
    for k in range(A):
        diag_ax[k].plot(freqs, diags[k]/max_val1)
        diag_ax[k].set_ylim(ylims1)

    # fill off diagonal plots
    for k in range(A):
        for l in range(k+1,A):
            offd_ax1[k,l].plot(freqs, offdiags1[k,l]/max_val1, color='k')
            offd_ax2[k,l].plot(freqs, offdiags2[k,l]/max_val1, color='r')
            offd_ax1[k,l].set_ylim(ylims1)
            offd_ax2[k,l].set_ylim(ylims2)


    # handle legends, fonts, etc.. for diagonal plots
    for k in range(A):
        diag_ax[k].set_xlabel('Hz')
        diag_ax[k].set_ylabel(labels[0], color='k')
        freq_ticks = (freqs[0], int(freqs[-1]/2), freqs[-1])
        diag_ax[k].set_xticks(freq_ticks)
        if k == 0 and areaList is not None:
            diag_ax[0].set_title(areaList[0], fontsize=R_LABEL_SIZE, \
                    fontweight='bold')
        diag_ax[k].text(R_LABEL_X_OFF, R_LABEL_Y_OFF, areaList[k],
                        fontsize=R_LABEL_SIZE, fontweight='bold')

    # handle legends, fonts, etc.. for off diagonal plots
    for k in range(A):
        for l in range(k+1,A):
            offd_ax1[k,l].axes.yaxis.set_ticks([])
            offd_ax1[k,l].axes.xaxis.set_ticks([])
            offd_ax2[k,l].axes.xaxis.set_ticks([])
            if k == 0 and areaList is not None:
                offd_ax1[k,l].set_title(areaList[l], fontweight='bold')
            if l == (A-1):
                offd_ax2[k,l].set_ylabel(labels[1], color='r')
                offd_ax2[k,l].yaxis.set_label_position("right")
                offd_ax2[k,l].set_yticks((-1,0,1))
            else:
                offd_ax2[k,l].axes.yaxis.set_ticks([])

    plt.show()


def get_factor_features(feature_str, feature_labels, factor):
    """
    Return the indicated subset of factor features.

    Parameters
    ----------
    feature_str : str
        string indicating type of features to extract
    feature_labels : list of str
        list of labels indicating the features represented in 'factor'. These
        are generated by 'saveFeatures.m' and contained in the labels loaded in
        by 'load_data'
    factor : numpy.ndarray
        components vector, or loadings, from a single factor of a linear factor
        model. Shape = [Fx1] where F = # number of features

    Returns
    -------
    factor_features : numpy.ndarray
        The specified subset of the components vector.
    """
    f_list = [findall(feature_str, fl) for fl in feature_labels]
    f_idx = [bool(x) for x in f_list]
    return factor[f_idx]


def ld_plot_features(factor, labels):
    """
    Extract features to generate grid plot of factor composed of power
    and linear directionality. This was the format for the LD factor
    respresentation that Kaf has requested.

    Parameters
    ----------
    factor: numpy.ndarray
        components vector, or loadings, from a single factor of a linear
        factor model. Shape = [Fx1] where F = # number of features.
    labels: dictionary
        labels dictionary should match the dictionary returned by
        data_tools.load_data when the data used to generate the factor
        is loaded.

    Returns
    -------
    pow_mat: numpy.ndarray
        Array of PSD values extracted from factor. Shape = [AxF] where
        A = # of areas; F = # of frequencies.
    sync_mat: numpy.ndarray
        Array of 'synchrony' values representing the sum of linear
        directionality in both directions for a pair of regions.
        Shape = [AxAxF] where A = # of areas; F = # of frequencies.
    diff_mat: numpy.ndarray
        Array of 'difference' values representing the difference of
        linear directionality in both directions for a pair of regions.
        Shape = [AxAxF] where A = # of areas; F = # of frequencies.
    """
    A = len(labels['area'])
    F = len(labels['f'])

    pow_mat = np.zeros((A,F))
    sync_mat = np.zeros((A,A,F))
    diff_mat = np.zeros((A,A,F))

    feature_labels = \
            np.hstack((labels['powerFeatures'], labels['ldFeatures']))

    # compile power features
    for k, a in enumerate(labels['area']):
        pow_features = '^' + a + ' [0-9]{1,3}$'
        pow_mat[k] = get_factor_features(pow_features, feature_labels, factor)

        # compile directionality features
        for l, b in enumerate(labels['area'][k+1:]):
            ld1_features = '^' + a + '->' + b + ' [0-9]{1,3}$'
            ld1 = get_factor_features(ld1_features, feature_labels, factor)

            ld2_features = '^' + b + '->' + a + ' [0-9]{1,3}$'
            ld2 = get_factor_features(ld2_features, feature_labels, factor)

            this_sync = ld1 + ld2
            this_diff = ld2 - ld1

            b_idx = k+1+l
            sync_mat[k,b_idx] = this_sync
            diff_mat[k,b_idx] = this_diff

    return pow_mat, sync_mat, diff_mat



if __name__ == '__main__':
    # Make an example plot with the test data located on teams
    # You will need to download those TeST.mat and TeST.json files from
    # teams and save them in this local repository location for this
    # script to work.
    from sklearn.decomposition import NMF

    power, ld, labels = load_data(os.path.join('testData', 'TeST.mat'), \
            feature_list=['power', 'directionality'], f_bounds=(1,50))
    print("Label keys:", list(labels.keys()))

    X, _ = feature_mat(labels, power, ld)

    nmf_model = NMF(n_components=4, init='nndsvda').fit(X)

    comp = np.squeeze(nmf_model.components_[0])

    pow_mat, sync_mat, diff_mat = ld_plot_features(comp, labels)
    factor_gridplot(pow_mat, sync_mat, diff_mat, areaList=labels['area'], \
            freqs=np.arange(1,51), labels=('Power','Diff.'))



###
