"""
Functions used by LPNE for visualization.

For now this really just has two useful functions:
- factor_gridplot, which allows you to generate an 'upper triangular'
    plot of an electome factor.
- factor_full_gridplot, which generates a 'full matrix' representation
    of a factor, with upper and lower triangular sections representing
    different directions.
"""
__date__ = "April 2021"


from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np
import os
from re import findall

from data_tools import load_data, feature_mat

YLIM_SCALE = 1.1
R_LABEL_SIZE = 12
R_LABEL_X_OFF = -35
R_LABEL_Y_OFF = 0.5

def base_gridplot(A, full=True):
    """ Generates a figure made up of a grid of subplots. """
    diag_ax = np.full(A, None, dtype=object)
    offd_ax1 = np.full((A,A), None, dtype=object)
    fig1 = plt.figure(constrained_layout=False)
    spec1 = GridSpec(ncols=A, nrows=A, figure=fig1)
    for k in range(A):
        diag_ax[k] = fig1.add_subplot(spec1[k,k])

    # Handle off-diagonal subplots differently for full or upper
    # triangular plots
    if full:
        for k in range(A):
            for m in range(A):
                if m != k:
                    offd_ax1[k,m] = fig1.add_subplot(spec1[k,m])

        return (diag_ax, offd_ax1)
    else:
        offd_ax2 = np.full((A,A), None, dtype=object)
        for k in range(A):
            for m in range(k+1,A):
                offd_ax1[k,m] = fig1.add_subplot(spec1[k,m])
                offd_ax2[k,m] = offd_ax1[k,m].twinx()

        return (diag_ax, offd_ax1, offd_ax2)


def plot_diags(diag_ax, diags, freqs, max_val):
    """ Fill diagonal plots """
    ylims = YLIM_SCALE*np.asarray([0, 1])
    for k in range(diags.shape[0]):
        diag_ax[k].plot(freqs, diags[k]/max_val, color='k')
        diag_ax[k].set_ylim(ylims)
        diag_ax[k].axes.xaxis.set_ticks([])
        diag_ax[k].axes.yaxis.set_ticks([])


def set_xaxis(ax, freqs):
    """ Set x axis labels and ticks """
    ax.set_xlabel('Hz')
    freq_ticks = (freqs[0], int(freqs[-1]/2), freqs[-1])
    ax.set_xticks(freq_ticks)

def set_yaxis(ax, labels, area):
    """ Set y axis labels """
    ax.set_ylabel(labels, color='k')
    ax.text(R_LABEL_X_OFF, R_LABEL_Y_OFF, area, fontsize=R_LABEL_SIZE, fontweight='bold')


def factor_full_gridplot(diags, offdiags, areaList=None, freqs=1, ylabel=''):
    """
    Plots a factor in a full grid format.

    Parameters
    ----------
    diags : numpy.ndarray
        Contains features that should be plotted in subfigures along the
        diagonal (e.g. PSD). Shape: [AxF] where A = # of areas; F = # of
        frequencies
    offdiags : numpy.ndarray
        Contains the set of features to be plotted in off-diagonal subfigures.
        Should have same scale as diag features. Shape: [AxAxF]
    areaList : list of str
        List of areas that features are generated from.
    freqs : list or int
        List of frequencies or scalar frequency sampling spacing (Hz). ...
    ylabel: string
        y-axis text label, applied to diags and offdiags.
    """

    # convert frequencies to array if scalar frequency resolution given
    if np.isscalar(freqs):
        freqs = np.arange(1,F+1,freqs)

    A = diags.shape[0]
    diag_ax, offd_ax = base_gridplot(A, full=True)
    max_val = np.max([np.max(abs(diags)), np.max(abs(offdiags))])
    plot_diags(diag_ax, diags, freqs, max_val)

    # Fill off diagonal plots
    ylims = YLIM_SCALE*np.asarray([0, 1])
    for k in range(A):
        for m in range(A):
            if m != k:
                offd_ax[k,m].plot(freqs, offdiags[k,m]/max_val, color='k')
                offd_ax[k,m].set_ylim(ylims)
                offd_ax[k,m].axes.yaxis.set_ticks([])
                offd_ax[k,m].axes.xaxis.set_ticks([])

    # set x axis labels and ticks for bottom row
    set_xaxis(diag_ax[-1], freqs)
    for k in range(A-1):
        set_xaxis(offd_ax[-1,k], freqs)

    # set y axis labels for left column
    set_yaxis(diag_ax[0], ylabel, areaList[0])
    for k in range(1,A):
        set_yaxis(offd_ax[k,0], ylabel, areaList[k])

    # set area title for columns
    diag_ax[0].set_title(areaList[0], fontsize=R_LABEL_SIZE, fontweight='bold')
    for k in range(1,A):
        offd_ax[0,k].set_title(areaList[k], fontsize=R_LABEL_SIZE, fontweight='bold')

    plt.show()


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

    # convert frequencies to array if scalar frequency resolution given
    if np.isscalar(freqs):
        freqs = np.arange(1,F+1,freqs)

    A = diags.shape[0]
    diag_ax, offd_ax1, offd_ax2 = base_gridplot(A, full=False)
    max_val1 = np.max([np.max(abs(diags)), np.max(abs(offdiags1))])
    plot_diags(diag_ax, diags, freqs, max_val1)

    # calculate ylim values for plots
    max_val2 = np.max(abs(offdiags2))
    ylims1 = YLIM_SCALE*np.asarray([0, 1])
    ylims2 = YLIM_SCALE*np.asarray([-max_val2/max_val1, max_val2/max_val1])

    # fill off diagonal plots
    for k in range(A):
        for m in range(k+1,A):
            offd_ax1[k,m].plot(freqs, offdiags1[k,m]/max_val1, color='k')
            offd_ax2[k,m].plot(freqs, offdiags2[k,m]/max_val1, color='r')
            offd_ax1[k,m].set_ylim(ylims1)
            offd_ax2[k,m].set_ylim(ylims2)

    # handle legends, fonts, etc.. for diagonal plots
    for k in range(A):
        set_xaxis(diag_ax[k], freqs)
        set_yaxis(diag_ax[k], labels[0], areaList[k])
        if k == 0 and areaList is not None:
            diag_ax[0].set_title(areaList[0], fontsize=R_LABEL_SIZE, \
                    fontweight='bold')

    # handle legends, fonts, etc.. for off diagonal plots
    for k in range(A):
        for m in range(k+1,A):
            offd_ax1[k,m].axes.yaxis.set_ticks([])
            offd_ax1[k,m].axes.xaxis.set_ticks([])
            offd_ax2[k,m].axes.xaxis.set_ticks([])
            if k == 0 and areaList is not None:
                offd_ax1[k,m].set_title(areaList[m], fontsize=R_LABEL_SIZE,
                                        fontweight='bold')
            if m == (A-1):
                offd_ax2[k,m].set_ylabel(labels[1], color='r')
                offd_ax2[k,m].yaxis.set_label_position("right")
                offd_ax2[k,m].set_yticks((-1,0,1))
            else:
                offd_ax2[k,m].axes.yaxis.set_ticks([])

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


def ld_plot_features(factor, labels, sync_diff=True):
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
    sync_diff: boolean
        If true, return the sum (i.e. 'sync') and difference of
        directionality features. Otherwise just return a matrix
        containing the directionality features.

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
    if sync_diff:
        sync_mat = np.zeros((A,A,F))
        diff_mat = np.zeros((A,A,F))
    else:
        ld_mat = np.zeros((A,A,F))

    feature_labels = \
            np.hstack((labels['powerFeatures'], labels['ldFeatures']))

    # compile power features
    for k, a in enumerate(labels['area']):
        pow_features = '^' + a + ' [0-9]{1,3}$'
        pow_mat[k] = get_factor_features(pow_features, feature_labels, factor)

        # compile directionality features
        for m, b in enumerate(labels['area'][k+1:]):
            ld1_features = '^' + a + '->' + b + ' [0-9]{1,3}$'
            ld1 = get_factor_features(ld1_features, feature_labels, factor)

            ld2_features = '^' + b + '->' + a + ' [0-9]{1,3}$'
            ld2 = get_factor_features(ld2_features, feature_labels, factor)

            b_idx = k+1+m
            if sync_diff:
                sync_mat[k, b_idx] = ld1 + ld2
                diff_mat[k, b_idx] = ld2 - ld1
            else:
                ld_mat[k, b_idx] = ld1
                ld_mat[b_idx, k] = ld2
    if sync_diff:
        return pow_mat, sync_mat, diff_mat
    else:
        return pow_mat, ld_mat



if __name__ == '__main__':
    # Make an example plot with the test data located on teams
    # You will need to download those TeST.mat and TeST.json files from
    # teams and save them in this local repository location for this
    # script to work.
    from sklearn.decomposition import NMF

    power, ld, labels = load_data(os.path.join('testData', 'TeST.mat'), \
            feature_list=['power', 'directionality_pairwise'], f_bounds=(1,50))

    X, _ = feature_mat(labels, power, ld)

    nmf_model = NMF(n_components=5, init='nndsvdar').fit(X)

    comp = np.squeeze(nmf_model.components_[3])

    pow_mat, ld_mat = ld_plot_features(comp, labels, sync_diff=False)
    factor_full_gridplot(pow_mat, ld_mat, areaList=labels['area'], \
            freqs=np.arange(1,51), ylabel=('Pow. Dens.'))



###
