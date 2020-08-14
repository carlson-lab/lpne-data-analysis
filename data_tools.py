'''
LIST OF FUNCTIONS
load_data: Loads and extracts data
getLabels: extracts labels from dictionary
data_subset: extracts subset of data
combine_data: combines datasets
get_X: generate feature matrix from list of feature arrays
'''
import json
import h5py
import numpy as np
from sklearn.preprocessing import LabelEncoder
from copy import deepcopy

def load_data(filename, f_bounds=(1,56), feature_list=['power', 'coherence', 'granger']):
    """ Loads and extracts data from a JSON file with preprocessed data.
    
    Loads the data from the file, then extracts each field, takes only the data
    within the frequency bounds specified by the input f_bounds, and transforms each
    data matrix to be a 2-dimensional matrix where axis 0 is the time window and axis 1
    iterates first by brain area and then by frequency.
    
    INPUTS
    filename: name of the .mat file containing the data
        and labels variables. The .mat file should contain the fields you
        want to load (e.g. 'power', 'coherence', 'granger'), and the corresponding .json
        file contains labels.
    feature_list: list of strings indicating which variables to load from the .mat file
    f_bounds: frequency bounds to analyze; default is set to between 1 and 120 Hz, inclusive.
    
    OUTPUTS
    power: Transformed matrix of power values within the frequency range given by f_bounds
        in the form MxN, where M is time windows and N is a combination of brain area and
        frequency, iterating first frequency then by area.
    coherence: Transformed matrix of coherency values within the frequency range given by
        f_bounds; same dimensions as power. In this case, brain area refers to the pair of brain
        areas being compared.
    granger: Transformed matrix of exponentiated Granger causality values 
        within the frequency range given by f_bounds; same dimensions as power. 
        In this case, brain area refers to the pair of brain areas being compared.
    instant: Transformed matrix of exponentiated instantaneous causality values 
        within the frequency range given by f_bounds; same dimensions as power. 
        In this case, brain area refers to the pair of brain areas being compared.
    labels: Structure containing labeling information for the data. 
            FIELDS:
            'windows': Structure containing information on windows.
                FIELDS:
                'task': Task being done.
                'mouse': Name of mouse.
            'powerFeatures': list of string labels describing the
                    features represented in power.
            'cohFeatures': list of string labels describing the features 
                represented in coherence
            'gcFeatures': list of string labels describing the features
                represented in granger
            'instFeatures': list of string labels describing the features
                represented in instant
    """
    jfile = filename.replace('.mat','.json') 
    with open(jfile) as f:
        labels = json.load(f)
    if len(feature_list) == 0:
        return (labels, )

    features = list()
    with h5py.File(filename, 'r') as file:
        for ft in feature_list:
            features.append(list(file[ft]))

    # only take the data from frequency values within bounds given
    (fLow, fHigh) = f_bounds
    fIdx = [k for (k, f) in enumerate(labels['f']) if fLow <= f <= fHigh]
    labels['f'] = np.asarray(labels['f'])[fIdx]

    # convert to array, invert axes, take power, coherency, gc data at
    # indices specified from frequency
    for k,ft in enumerate(feature_list):
        if ft == 'power':
            features[k] = np.asarray(features[k])
            features[k] = np.swapaxes(features[k], 2,0)
            features[k] = features[k][fIdx,:,:]

            # reshape each nd array to matrix after reshaping, axis 0
            # is windows; axis 1 iterates through frequency first,
            # then channel
            a, b, c = features[k].shape
            features[k] = features[k].reshape(a*b, c, order='F').T

            if 'powerFeatures' in labels.keys():
                # reshape corresponding array of feature labels
                # MAKE SURE THESE OPERATIONS CORRESPOND TO OPERATIONS ON ACTUAL FEATURES ABOVE
                pf = np.asarray(labels['powerFeatures'])
                pf = pf[fIdx]
                labels['powerFeatures'] = pf.reshape(a*b, order='F')

            if 'powVersion' in labels.keys():
                p_version = labels['powVersion']
                print('version {0} used to calcuate power features'.format(p_version))
            else:
                print('Power features calculated using unknown version')

        if ft == 'coherence':
            features[k] = np.asarray(features[k])
            features[k] = np.swapaxes( features[k], 1,2)
            features[k] = np.swapaxes( features[k], 0,3)
            features[k] = np.swapaxes( features[k], 2,3)
            features[k] = features[k][fIdx,:,:,:].astype('float64')

            # collect indices of upper triangular portion of brain region x brain region area matrix
            r1, c1 = np.triu_indices( features[k].shape[-1], k=1)
            features[k] = features[k][..., r1,c1]

            features[k] = np.swapaxes(features[k],0,1)
            a, b, c = features[k].shape
            features[k] = features[k].reshape(a, b*c, order='F')

            if 'cohFeatures' in labels.keys():
                # reshape corresponding array of feature labels
                # MAKE SURE THESE OPERATIONS CORRESPOND TO OPERATIONS ON ACTUAL FEATURES ABOVE
                cf = np.asarray(labels['cohFeatures'])
                cf = np.swapaxes(cf, 1,2)
                cf = cf[fIdx,:,:]
                cf = cf[:,r1,c1]
                labels['cohFeatures'] = cf.reshape(b*c, order='F')

            if 'cohVersion' in labels.keys():
                c_version = labels['cohVersion']
                print('version {0} used to calcuate coherence features'.format(c_version))
            else:
                print('Coherence features calculated using unknown version')


        if ft == 'granger':
            gcFIdx = [k+1 for k in fIdx]

            gcArray = np.asarray(features[k])
            gcArray = gcArray[:, gcFIdx, :]
            gcArray = np.minimum(np.exp(gcArray), 10)
            features[k] = np.transpose(gcArray, (1,2,0))

            a,b,c = features[k].shape
            features[k] = features[k].reshape(a*b, c, order='F').T

            if 'gcFeatures' in labels.keys():
                # reshape corresponding array of feature labels
                # MAKE SURE THESE OPERATIONS CORRESPOND TO OPERATIONS ON ACTUAL FEATURES ABOVE
                gf = np.asarray(labels['gcFeatures'])
                gf = gf[:, gcFIdx].T
                labels['gcFeatures'] = gf.reshape(a*b, order='F')

            if 'gcVersion' in labels.keys():
                g_version = labels['gcVersion']
                print('version {0} used to calcuate granger features'.format(g_version))
            else:
                print('Granger features calculated using unknown version')


        if ft == 'instant':
            gcFIdx = [k+1 for k in fIdx]

            instArray = np.asarray(features[k])
            instArray = instArray[:, gcFIdx, :]
            instArray = np.minimum(np.exp(instArray), 10)
            features[k] = np.transpose(instArray, (1,2,0))
            a,b,c = features[k].shape
            features[k] = features[k].reshape(a*b, c, order='F').T

            if 'instFeatures' in labels.keys():
                # reshape corresponding array of feature labels
                # MAKE SURE THESE OPERATIONS CORRESPOND TO OPERATIONS ON ACTUAL FEATURES ABOVE
                ift = np.asarray(labels['instFeatures'])
                ift = ift[:, gcFIdx].T
                labels['instFeatures'] = ift.reshape(a*b, order='F')


        if ft == 'xFft':
            # account for double precision roundoff error
            tol = 1e-6
            sIdx = [k for (k, s) in enumerate(labels['s']) if (fLow-tol) <= s <= (fHigh+tol)]
            labels['s'] = np.asarray(labels['s'])[sIdx]

            xArray = np.asarray(features[k])
            W,C,S = xArray.shape

            features[k] = np.array([[[xArray[w,c,s][0] + 1j*xArray[w,c,s][1]
                                      for s in sIdx] for c in range(C)]
                                    for w in range(W)], dtype=np.complex64)

            if 'fftVersion' in labels.keys():
                f_version = labels['fftVersion']
                print('version {0} used to calcuate FFT'.format(f_version))

    if 'preprocessVersion' in labels.keys():
        pp_version = labels['preprocessVersion']
        print('Version {0} of preprocessing used'.format(pp_version))
        print('Make sure feature versions listed above match those used for any other '
              'dataset in the same project')
    else:
        print('Using data that was preprocessed with unknown preprocessing version. '
              'Please make sure all datasets in the same project were preprocessed '
              'the same way.')
        
    features.append(labels)

    return tuple(features)


def get_labels(labels, variable_name = 'task'):
    ''' Generates multivariate labels.

    INPUT
    labels: labels variable from preprocessed data
    variableName: name of the field you want for y variables

    OUTPUT
    y: task labels for each window '''

    task_strings = labels['windows'][variable_name]
    y = LabelEncoder().fit_transform(task_strings)    
    return y 


def data_subset(condition, *args):
    """ Returns a subset of the given data.

    INPUTS
    x: numpy array of data (WxF) where W is number of windows and F is number of features
    labels: labels variable from preprocessed data
    condition: boolean list/vector of length W
    args: optional arguments which may be numpy arrays with first dimension
        equal to W, or may be labels variable from preprocessed data
    """
    K = len(args)
    if K == 0 :
        return None

    # take subset for any additional arguments passed
    feats = list()
    for k in range(K):
        # special case, if tuple contains dictionaries, assume they are labels
        if type(args[k]) is dict:
            lab_copy = deepcopy(args[k])
            # iterate over each key in 'windows' dictionary
            for key, value in lab_copy['windows'].items():
                lab_copy['windows'][key] = np.asarray(value)[condition]
            feats.append(lab_copy)
        else:
            this_arr = np.asarray(args[k])
            feats.append(this_arr[condition])

    if K == 1:
        return feats[0]
    else:
        return tuple(feats)


def combine_data(*args):
    """ Combines multiple data subsets together, passed in as tuples.
    Assumes subsets are compatible (i.e. contain same feature space)
    Example: X, labels = combine_data((X1,X2,X3), (labels1,labels2,labels3))
    """
    K = len(args)
    if K == 0:
        return None
    
    features = list()
    for k in range(K):
        # special case, if tuple contains dictionaries, assume they are labels
        if type(args[k][0]) is dict:
            labels = deepcopy(args[k][0])
            for key in labels['windows']:
                values = [x for l in args[k] for x in l['windows'][key]]
                labels['windows'][key] = values
            features.append(labels)
        else:
            # combine any additional tuples passed in
            features.append(np.concatenate(args[k]))

    if K == 1:
        return features[0]
    else:
        return tuple(features)


def get_X(weights, feature_list, return_weights=False):
    rms = [np.sqrt(np.mean(f**2)) for f in feature_list]
    new_weights = weights / rms
    weighted_feat = [f*w for f,w in zip(feature_list, new_weights)]
    X = np.concatenate(weighted_feat, axis=1)
    if return_weights:
        return (X, new_weights)
    else:
        return X


def get_weights(group_list, subset_idx=None):
    """ For balancing sample weights. group_list should be a list of group
    labels. The first grouping listed will be balanced across the whole
    dataset. The next grouping will be relatively balanced within each
    individual group given by the first grouping. That pattern continues
    iteratively so that each grouping listed is relatively balanced within
    individual groups from the grouping that comes before it in the list.
    """
    if subset_idx is None:
        subset_idx = np.full(len(group_list[0]), True)
    else:
        subset_idx = np.array(subset_idx, dtype=np.bool_)

    group_labels = np.asarray(group_list[0])
    these_group_labels = group_labels[subset_idx]
    group_names = set(these_group_labels)

    weights = np.zeros(sum(subset_idx))
    for g in group_names:
        g_idx = np.asarray([x == g for x in group_labels])
        gsub_idx = g_idx[subset_idx]

        if len(group_list) > 1:
            subgroup_weights = get_weights(group_list[1:], subset_idx & g_idx)
            weights[gsub_idx] = subgroup_weights / sum(gsub_idx)
        else:
            weights[gsub_idx] = 1/sum(gsub_idx)

    weights = weights/np.mean(weights)
    return weights
