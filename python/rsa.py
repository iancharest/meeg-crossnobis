"""Module containing functions for temporal RSA"""
import mne
import numpy as np
from mne.cov import compute_whitener
from joblib.parallel import Parallel, delayed
from scipy.spatial.distance import cdist
from sklearn.model_selection import StratifiedShuffleSplit


def correlation(x, y, *args):
    x_ = x.mean(axis=0)
    y_ = y.mean(axis=0)
    return np.arctanh(1. - cdist(x_, y_, metric='correlation'))


def euclidean(x, y, *args):
    x_ = x.mean(axis=0)
    y_ = y.mean(axis=0)
    return cdist(x_, y_, metric='euclidean')


def _compute_fold(epoch, targets, train, test, metric_fx=correlation,
                  cv_multi_normalize=None):
    """Computes pairwise metric across time for one fold

    Arguments
    ---------
    epoch : instance of mne.Epoch
    targets : array (n_trials,)
        target (condition) for each trials
    train : array-like of int
        indices of the training data
    test : array-like of int
        indices of the testing data
    metric_fx : function(x, y, targets_train, targets_test)
        any function that returns a scalar given two arrays.
        This condition must hold: metric_fx(x, y) == metric_fx(y, x)
    cv_multi_normalize : str | None (default None)
        Multivariate normalize the trials before computing the distance between
        pairwise conditions.
        Valid values are 'epoch' | 'baseline' | None:
            - 'epoch' computes the covariance matrix on the entire epoch;
            - 'baseline' uses only the baseline condition; requires to pass an
              array times

    Returns
    -------
    rdms: array (n_pairwise_targets, n_times, n_times)
    """

    targets = np.asarray(targets)

    if cv_multi_normalize not in ('epoch', 'baseline', None):
        raise ValueError(
            "cv_multi_normalize must be one of {0}".format(
                ('epoch', 'baseline', None)))

    # make sure we don't have tuple as labels but strings or integers
    # to avoid possible bug in scikit-learn
    for t in targets:
        assert (isinstance(t, (str, int)))

    # preallocate array
    unique_targets = np.unique(targets)
    n_unique_targets = len(unique_targets)
    n_times = len(epoch.times)
    n_triu = n_unique_targets * (n_unique_targets - 1) / 2 + n_unique_targets
    rdms = np.zeros((n_triu, n_times, n_times))
    # get training and testing data
    epoch_train = epoch.copy()[train]
    targets_train = targets[train]
    assert(len(epoch_train) == len(targets_train))
    epoch_test = epoch.copy()[test]
    targets_test = targets[test]
    assert(len(epoch_test) == len(targets_test))

    if cv_multi_normalize:
        tmax = 0. if cv_multi_normalize == 'baseline' else None
        cov_train = mne.compute_covariance(epoch_train,
                                           tmax=tmax, method='shrunk')
        W_train, ch_names = compute_whitener(cov_train, epoch_train.info)
        # whiten both training and testing set
        epo_data_train = np.array([np.dot(W_train, e)
                                   for e in epoch_train.get_data()])
        epo_data_test = np.array([np.dot(W_train, e)
                                  for e in epoch_test.get_data()])
    else:
        epo_data_train = epoch_train.get_data()
        epo_data_test = epoch_test.get_data()

    # compute pairwise metric
    idx = 0
    for i, target1 in enumerate(unique_targets):
        mask_target1 = np.where([t == target1 for t in targets_train])[0]
        assert (len(mask_target1) > 0)
        epo_data_target1 = epo_data_train[mask_target1]
        for j, target2 in enumerate(unique_targets[i:]):
            mask_target2 = np.where([t == target2 for t in targets_test])[0]
            assert (len(mask_target2) > 0)
            epo_data_target2 = epo_data_test[mask_target2]
            # now loop through time
            for t1 in range(n_times):
                for t2 in range(t1, n_times):
                    # XXX: in this way we can also pass a classifier
                    rdms[idx, t1, t2] = \
                        metric_fx(epo_data_target1[..., t1][None, :],
                                  epo_data_target2[..., t2][None, :],
                                  targets_train, targets_test)
                    rdms[idx, t2, t1] = rdms[idx, t1, t2]
            idx += 1
    return rdms


def compute_temporal_rdm(epoch, targets, metric_fx=correlation,
                         cv=StratifiedShuffleSplit(n_splits=10, test_size=0.5),
                         cv_multi_normalize=None,
                         n_jobs=1):
    """Computes pairwise metric across time

    Arguments
    ---------
    epoch : instance of mne.Epoch
    targets : array (n_trials,)
        target (condition) for each trials
    metric_fx : function(x, y, train_targets, test_targets)
        any function that returns a scalar given two arrays.
        This condition must hold: metric_fx(x, y) == metric_fx(y, x)
    cv : sklearn cross-validator
    cv_multi_normalize : str | None (default None)
        Multivariate normalize the trials before computing the distance between
        pairwise conditions.
        Valid values are 'epoch' | 'baseline' | None:
            - 'epoch' computes the covariance matrix on the entire epoch;
            - 'baseline' uses only the baseline condition; requires to pass an
              array times
    n_jobs : int (default 1)

    Returns
    -------
    rdm: array (n_pairwise_targets, n_times, n_times)
        the cross-validated RDM over time
    """
    splits = cv.split(targets, targets)
    rdm = Parallel(n_jobs=n_jobs)(
        delayed(_compute_fold)
        (epoch, targets, train, test,
         metric_fx=metric_fx, cv_multi_normalize=cv_multi_normalize)
        for train, test in splits)

    rdm = np.stack(rdm, axis=-1).mean(axis=-1)
    return rdm
