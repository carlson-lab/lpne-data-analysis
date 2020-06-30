'''
LIST OF FUNCTIONS
run_nested_cv
run_cv
run_split
get_balanced_groupings
gen_splits
get_performance
train_factor_model
make_folds
make_param_dict
split_features
plot_factors
'''

import numpy as np
from sklearn.decomposition import PCA, NMF
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import LabelEncoder, LabelBinarizer
from sklearn.svm import SVC
import collections
from sklearn.ensemble import RandomForestClassifier
from data_tools import get_X
import re
import matplotlib.pyplot as plt

RAND_STATE = 42

def run_nested_cv(featureList, y, group, parameters, folds=0, metric='auc',
                  modelOpts={'reduction_method':'NMF','classifier':'logistic'},
                  standardizeFeatures=True, split=None):
    ''' Function to run nested cross-validation. Will use averaged performance
        across cross-validation splits to find the best hyperparameters.
    
    INPUT:
    X: feature data (for example power, coherency, Grainger
        causality). Array with first dimension equal to the number of
        windows. Remaining dimension(s) iterates over freatures, most
        likely sorted by frequency and brain region, iterating first
        by brain region and then frequency.
    y: target data. WxL one-hot array, where W is # of windows and L
        is # of classes.
    group: pre-split folds, dividing up data into sections for cross
        validation. List or array of labels with length equal to number of
        windows. Each window should have a label corresponding a fold.
    metric: performance metric used, either 'accuracy', 'precision', or 'auc'.
    folds: number of folds to use for cross-validation. If 0, uses 
        hold-one-out scheme.
    modelOpts: dictionary which may contain one of the two following
            pairs of fields:
                'classifier': either 'logistic','forest', or 'svm'
                'reduction_method': either 'PCA' or 'NMF'
            OR
                'train func': function for training a custom
                    model. Takes parameters 'data, options'.  data
                    should be a tuple of values (X,y) continaing the
                    training data, which are the same form as X and y
                    described above. options should be a dictionary of
                    options or parameters for training the model.
                'eval func': function for evaluating the performance
                    of a custom model. Takes parameters 'model, data,
                    metric'. model should contain a trained
                    model. Data contains the data to be evaluated and
                    takes the same form as in 'train func'. metric can
                    be any of the options given in the
                    metric field described above.
    parameters: dictionary which may contain any fields required
            for training the desired model. If not using a custom
            model, the following fields are expected:
            'dimensions': number of features to use for dimensionality
                reduction
            'reg_class': regularization values to validate for
                classifier. See train_factor_model for details.
            'reg_factor': regularization values to validate for
                factor model. See train_factor_model for details.
    OUTPUTS:
    myResults: dictionary containing results. Fields vary depending on
        layer number
            'performance': vector of test set performances
            'fold results': list of dictionaries containing
            'performance', 'model', 'performance all models', and
            'parameters' for all test sets in the cross validation
        '''
    if split is None:
        split = gen_splits(group, folds)
    nGroups = len(split)

    if standardizeFeatures:
        for k,f in enumerate(featureList):
            featureList[k] = f / np.std(f)
        
    custom = 'train func' in modelOpts
    modelOpts['custom'] = custom

    splitResults = [None]*nGroups
    performance = np.zeros(nGroups)
    perf_data= [None]*nGroups
        
    for k,s in enumerate(split):
        trIdx = s['train']
        tsIdx = s['test']
        featuresTrain = split_features(featureList, trIdx)
        featuresTest = split_features(featureList, tsIdx)
        yTrain = y[trIdx]
        yTest = y[tsIdx]
        groupTrain = group[trIdx]

        # for each test mouse, run 1 layer and return results
        if 'validation split' in s.keys():
            vs = s['validation split']
        else:
            vs = None
        splitResults[k] = run_cv(featuresTrain, yTrain, groupTrain, parameters=parameters,
                                 folds=folds, metric=metric, modelOpts=modelOpts,
                                 split=vs, standardizeFeatures=False)

        s['validation split'] = splitResults[k]['split']
        thisModel = splitResults[k]['model']
        theseWeights = splitResults[k]['parameters']['feature_weights']

        XTest = get_X(theseWeights, featuresTest)

        if custom:
            performance[k], perf_data[k] = modelOpts['eval func'](thisModel, (XTest, yTest),
                                                                  metric)
        else:
            performance[k], perf_data[k] = get_performance(thisModel, XTest, yTest, metric)

    myResults = dict()
    myResults['performance'] = performance
    myResults['fold results'] = splitResults
    myResults['performance data'] = perf_data
    myResults['split'] = split
    
    return myResults


def run_cv(featureList, y, group, parameters, folds=0, split=None, metric='auc',
           modelOpts={'reduction_method':'NMF','classifier':'logistic'},
           standardizeFeatures=True):
    ''' Function to run cross-validation. Data is split into K sets; each set 
    is then rotated as a holdout set and the average performance is evaluated 
    to determine the best combination of hyperparameters. Will use averaged 
    performance across cross-validation splits to find the best hyperparameters.
    
    INPUT:
    X: feature data (for example power, coherency, Grainger
        causality). Array with first dimension equal to the number of
        windows. Remaining dimension(s) iterates over freatures, most
        likely sorted by frequency and brain region, iterating first
        by brain region and then frequency.
    y: target data. WxL one-hot array, where W is # of windows and L
        is # of classes.
    group: pre-split folds, dividing up data into sections for cross
        validation. List or array of labels with length equal to number of
        windows. Each window should have a label corresponding a fold.
    metric: performance metric used, either 'accuracy', 'precision', or 'auc'.
    folds: number of folds to use for cross-validation. If 0, uses 
        hold-one-out scheme.
    split: information about which data to use based upon prior
        generation of test/train/validation splits.  Consists of a
        list with length = number of splits; each index contains a
        dictionary with entries 'train', 'val', and 'test'. Each entry
        links to the windows corresponding to that category.
    modelOpts: dictionary which may contain one of the two following
            pairs of fields:
                'classifier': either 'logistic','forest', or 'svm'
                'reduction_method': either 'PCA' or 'NMF'
            OR
                'train func': function for training a custom
                    model. Takes parameters 'data, options'.  data
                    should be a tuple of values (X,y) continaing the
                    training data, which are the same form as X and y
                    described above. options should be a dictionary of
                    options or parameters for training the model.
                'eval func': function for evaluating the performance
                    of a custom model. Takes parameters 'model, data,
                    metric'. model should contain a trained
                    model. Data contains the data to be evaluated and
                    takes the same form as in 'train func'. metric can
                    be any of the options given in the
                    metric field described above.
    parameters: dictionary which may contain any fields required
            for training the desired model. If not using a custom
            model, the following fields are expected:
            'dimensions': number of features to use for dimensionality
                reduction
            'reg_class': regularization values to validate for
                classifier. See train_factor_model for details.
            'reg_factor': regularization values to validate for
                factor model. See train_factor_model for details.
    OUTPUTS:
    myResults: dictionary containing results.
            'performance': best validation set performance 
            'model': model with best validation performance
            'performance all models': array of all validation set performances
            'parameters': parameters associated with best validation performance
        '''
    if split is None:
        split = gen_splits(group, folds)
    nGroups = len(split)

    if standardizeFeatures:
        for k,f in enumerate(featureList):
            featureList[k] = f / np.std(f)
    
    # check if custom model being used
    custom = 'train func' in modelOpts
    modelOpts['custom'] = custom

    # initialize storage for results of each split
    ps = [len(x) for x in parameters.values() if isinstance(x,list)]
    pSize = [len(split)] + ps
    performance = np.zeros(pSize)
    perfData = [np.full(ps, None)]*len(split)

    for k,s in enumerate(split):
        trIdx = s['train']
        tsIdx = s['test']

        featuresTrain = split_features(featureList, trIdx)
        featuresTest = split_features(featureList, tsIdx)
        yTrain = y[trIdx]
        yTest = y[tsIdx]
        
        oneTrainClass = len(np.unique(yTrain)) < 2
        oneTestClass = len(np.unique(yTest)) < 2
        if oneTestClass:
            print('Warning: single class present in test data!')
            performance[k,...] = np.nan
        elif oneTrainClass:
            print('Warning: single class present in train data!')
            performance[k,...] = np.nan
        else:
            performance[k,...], perfData[k] = run_split(featuresTrain, yTrain,
                                                        featuresTest, yTest,
                                                        parameters, modelOpts, metric)
        
    # collect information describing best performing model 
    meanPerf = np.nanmean(performance, axis=0)
    bestIdx = np.unravel_index(meanPerf.argmax(), meanPerf.shape)
    bestParams = make_param_dict(parameters, bestIdx)

    X = get_X(bestParams['feature_weights'], featureList)
    if custom:
        bestModel = modelOpts['train func']((X, y), bestParams)
    else:
        bestModel = train_factor_model(X, y, modelOpts, bestParams)
    
    for k,pd in enumerate(perfData):
        perfData[k] = pd[bestIdx]
        
    bestPerf = performance
    for k in bestIdx:
        bestPerf = bestPerf[:,k]

    myResults = dict()
    myResults['parameters'] = bestParams
    myResults['performance all models'] = performance
    myResults['model'] = bestModel
    myResults['performance'] = bestPerf
    myResults['performance data'] = perfData
    myResults['split'] = split
    return myResults


def run_split(featuresTr, yTr, featuresTs, yTs, parameters, modelOpts, metric):
    # assumes features in trainset have already been normalized wrt each other

    # iterate through all combinations of parameters and train/evaluate
    pSize = [len(x) for x in parameters.values() if isinstance(x,list)]
    performance = np.zeros(pSize)
    perfData = np.full(pSize, dict())
    for idx, _ in np.ndenumerate(performance):
        theseParams = make_param_dict(parameters, idx)

        # weight features to generate training data
        if 'feature_weights' in theseParams:
            weights = theseParams['feature_weights']
        else:
            weights = np.ones(featuresTr.shape)

        XTr = get_X(weights, featuresTr)
        XTs = get_X(weights, featuresTs)

        if modelOpts['custom']:
            model = modelOpts['train func']((XTr, yTr), theseParams)
            performance[idx], perfData[idx] = modelOpts['eval func'](model, (XTs, yTs),
                                                                     metric)
        else:
            model = train_factor_model(XTr, yTr, modelOpts, theseParams)
            performance[idx], perfData[idx] = get_performance(model, XTs, yTs, metric)
    return (performance, perfData)


def get_balanced_groupings(labels, groupKey = 'genotype'):
    
    # Convert mouse names to a value and get a set of all unique names
    mousenames = labels['windows']['mouse']
    mouseId = LabelEncoder().fit_transform(mousenames)
    namesSet = list(set(mouseId))
    # Get a set of all unique group names 
    groupVals = LabelEncoder().fit_transform(labels['windows'][groupKey])
    groupSet = set(groupVals)
    
    # Initialize
    groupCounts = np.zeros(len(groupSet))
    groupIds = [None]*len(groupSet)
    
    # Separate out into unique groups and count how many in each group
    for i,g in enumerate(groupSet):
        groupIds[i] = np.unique(mouseId[groupVals==g])
        np.random.shuffle(groupIds[i])
        groupCounts[i] = len(groupIds[i])
    
    # Find the minimum number of groups to make
    minimum = min(groupCounts)
    idx = np.argwhere(groupCounts==minimum)[0]
    
    # Find how many of each group per new mouse grouping
    groupCountsNorm = np.floor(groupCounts/minimum)
    # Find how many need to be left over
    remainder = np.mod(groupCounts,minimum)
    
    # Initialize your final groupings
    finalGroups = np.asarray([-1]*len(mousenames))
    usedGroups = [0]*len(groupSet)
    
    
    # Make as few groups as possible
    # Iterate over number of final groups, then over number of groups
    # you're balancing
    for k in range(int(minimum)):
        # for each balancing category
        for j,g in enumerate(groupSet):
            # check if you need to add an extra mouse
            if remainder[j]>0:
                micePerGroup = groupCountsNorm[j]+1
                remainder[j]-=1
            else:
                micePerGroup = groupCountsNorm[j]
            # Add the right number of mice, determined earlier, to this group
            for m in range(int(micePerGroup)):
                thisMouseIdx = mouseId==groupIds[j][usedGroups[j]]
                usedGroups[j]+=1
                finalGroups[thisMouseIdx] = k
                
    return finalGroups


def gen_splits(group, folds=0):
    '''Generates nested splits using test mouse holdout and cross
        validation mouse holdout.

    INPUT
    group: list-like with length equal to the number of windows that
        need to be split up. Each element should have a label
        assigning the corresponding window to a group, where groups of
        windows will remain together across all splits. This may be
        desirable when all the windows belonging to the same subject
        should be grouped together, among other cases.
    folds: number of folds to use for cross-validation. If 0, uses 
        hold-one-out scheme.

    OUTPUT
    splits: information about which data to use based upon prior
        generation of test/train/validation splits.  Consists of a
        list with length = number of splits; each index contains a
        dictionary with entries 'train', 'val', and 'test'. Each entry
        links to the windows corresponding to that category.
    nGroups: total number of groups; can be used if creating
        double-layer splits
    '''
    # change to numerical format, convert to a list
    groupId = LabelEncoder().fit_transform(group)
    groupList = np.unique(groupId)
    nGroups = len(groupList)

    if folds:
        # count windows in each group, then sort indicies by size
        groupSize = np.zeros(nGroups)
        for g in groupList:
            groupSize[g] = np.sum(groupId==g)
        sizeIdx = np.argsort(groupSize)

        # count how many groups go in each fold, and how many folds
        # have 1 'extra' group
        groupsPerFold = nGroups//folds
        nBigFolds = nGroups%folds

        # get indicies associated with each fold
        nBigFGroups = (groupsPerFold+1)*nBigFolds
        bigFolds = make_folds(sizeIdx[:nBigFGroups], groupsPerFold+1)
        littleFolds = make_folds(sizeIdx[nBigFGroups:], groupsPerFold)

        # create foldId array like groupId, but for folds
        Folds = bigFolds + littleFolds
        foldId = np.asarray([k for g in groupId for k,f in enumerate(Folds) if g in f])
        foldList = np.unique(foldId)
        nFolds = len(foldList)
    else:
        foldId = groupId
        foldList = groupList
        nFolds = nGroups

    # create list of dictionaries corresponding to each split
    splits=[];
    for testFold in foldList:
        testMice = foldId==testFold
        trainMice = foldId!=testFold
        splits.append(dict())
        splits[-1]['test'] = testMice
        splits[-1]['train'] = trainMice
            
    return splits


def get_performance(model, X, y, metric):
    factorModel, classifier = model
    perfData = dict()
    
    try:
        scores = factorModel.transform(X)

        if metric == 'accuracy':
            performance = classifier.score(scores, y)

        else:
            perfData['d_func'] = classifier.decision_function(scores)
            lb = LabelBinarizer()
            y_1hot = lb.fit_transform(y)

            if metric =='auc':
                perfData['ovr_auc'] = roc_auc_score(y_1hot, perfData['d_func'], average=None)
                performance = np.mean(perfData['ovr_auc'])
            elif metric == 'precision':
                perfData['ovr_avg_prec'] = average_precision_score(y_1hot, perfData['d_func'],
                                                                   average=None)
                performance = np.mean(perfData['ovr_avg_prec'])

            
    except:
        print("Warning: performance could not be evaluated!")
        performance = np.nan
        
    return (performance, perfData)


def train_factor_model(X, y, modelOptions, parameters):
    N = parameters['dimensions']
    C = parameters['reg_class']
    alpha = parameters['reg_factor']
    reductionMethod = modelOptions['reduction_method']
    classifierMethod = modelOptions['classifier']


    
    # Run dimensionality reduction on the given data
    if reductionMethod=='NMF':
        factorModel = NMF(N, alpha=alpha, init='nndsvd', shuffle=True, random_state=RAND_STATE)
    elif reductionMethod=='PCA':
        factorModel = PCA(N, random_state=RAND_STATE)
    scores = factorModel.fit_transform(X)
    
    # just fits a classifier to the data given the regularization
    if classifierMethod=='logistic':
        classifier=LogisticRegression(penalty='l1', C=C, solver='saga', random_state=RAND_STATE,
                                      class_weight='balanced')
    elif classifierMethod == 'svm':
        classifier=SVC(C=C, kernel='rbf', class_weight='balanced', random_state=RAND_STATE)     
    elif classifierMethod =='forest':
        classifier=RandomForestClassifier(max_depth=C, random_state=RAND_STATE,
                                          class_weight='balanced')
    else:
        raise ValueError('Classifier type must be "logistic", "svm", or "forest".')
        
    classifier.fit(scores,y) 
    print('trained model:', parameters)

    return (factorModel, classifier)

def make_folds(idx, foldSize):
    '''divide group indicies (sorted by size) into folds
    
    INPUT 
    idx: list of indicies. Sorted by size, such that the first
        index corresponds to the smallest group.
    foldSize: number of groups to sort into each fold
   
    OUTPUT
    folds: list of tuples. Each tuple contains indicies from idx that
        are to be grouped into the same fold.
    '''
    idx = idx.tolist()
    folds = []
    
    takeMid = foldSize%2
    nOuter = foldSize//2

    while idx:
        thisFold = []
        if takeMid:
            midIdx = idx.pop(len(idx)//2)
            thisFold.append(midIdx)
        for k in range(nOuter):
            thisFold.append(idx.pop(0))
            thisFold.append(idx.pop(-1))
        folds.append(tuple(thisFold))
            
    return folds


def make_param_dict(parameters, idx):
    theseParams = {key:parameters[key] for key in parameters.keys()
                   if not isinstance(parameters[key],list)}
    listKeys = [key for key in parameters.keys() if isinstance(parameters[key],list)]
    for k,key in enumerate(listKeys):
        theseParams[key] = parameters[key][idx[k]]

    return theseParams


def split_features(featureList, idx):
    newList = list()

    for f in featureList:
        newList.append(f[idx])

    return newList


def plot_factors(factors, plot_feature, factor_no, feature_labels):
    plot_feature = '^' + plot_feature + ' [0-9]{1,3}$'
    plot_list = [re.findall(plot_feature, fl) for fl in feature_labels]
    plot_idx = [bool(x) for x in plot_list]
    freq = [float(re.findall('[0-9]{1,3}$', pl[0])[0]) for pl in plot_list if pl]
    
    plt.plot(freq, factors[factor_no, plot_idx])
    plt.xlabel('Frequency (Hz)')
    plt.show()
