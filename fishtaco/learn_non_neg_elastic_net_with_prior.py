__author__ = 'ManorLab'

################################################################################################################
# LearnNonNegElasticWithPrior.py: Learns a cross-validation Non Negative Elastic Net model with a Prior
# from given features
#
# notes:
# - does NOT normalize/centers X, so X needs to be normalized/centered before sending to function
# - does NOT center Y, so Y needs to be centered before sending to function
#
# LearnNonNegElasticWithPrior.learn(cov_train, res_train, params):
#   Learns the model
#
# cov_train: the covariate matrix of size NxP, with N samples and P covariates
#
# res_train: the response vector of size N, with N samples
#
# params: a dictionary containing optional params
#   params.['covariates_prior']: a vector of priors for the different features in the range [0,1]
#   params.['class_subfeatures']: a vector of binary case/control labels for the samples, dividing them to two classes.
#                                 If this parameter is given when learning, 2*P covariates will be created.
#                                 The P covariates are zero for controls and copied from cov_train for
#                                 the cases, and the second P covariates are zero for cases and copied from cov_train
#                                 for the controls. This matrix will replace the original, so the returned model
#                                 will have 2*P weights
#   params.['num_cv']: the number of folds in the internal cross-validation (default: 5)
#   params.['l1_ratio']: the ratio of L1 penalty in the regularization in the range [0,1] (default: 0.5)
#
# returns the elastic-net model and the vector of validation R^2
#
#
################################################################################################################

from sklearn import cross_validation, linear_model
import numpy as np


def learn(cov_train, res_train, params={}):
    np.set_printoptions(precision=2, suppress=False, linewidth=200)  # nicer output

    # update train covariates by the given params
    if 'covariates_prior' in params.keys() and params['covariates_prior'] is not None:
        cov_train = cov_train * params['covariates_prior']

    if 'class_subfeatures' in params.keys() and params['class_subfeatures'] is not None:
        # cases subfeatures
        cases = params['class_subfeatures'] > 0
        case_subfeatures = np.zeros_like(cov_train)
        case_subfeatures[cases, :] = cov_train[cases, :]
        # control subfeatures
        controls = params['class_subfeatures'] == 0
        control_subfeatures = np.zeros_like(cov_train)
        control_subfeatures[controls, :] = cov_train[controls, :]
        # update cov_train
        cov_train = np.hstack((case_subfeatures, control_subfeatures))

    # more params
    if 'num_cv' in params.keys():
        num_cv = params['num_cv']
    else:
        num_cv = 5

    if 'l1_ratio' in params.keys():
        l1_ratio = params['l1_ratio']
    else:
        l1_ratio = 0.5

    k_fold = cross_validation.KFold(len(res_train), n_folds=num_cv, shuffle=True)

    best_validation_rsqr = np.zeros(num_cv)
    best_validation_alpha = np.zeros(num_cv)

    for inner_k, (inner_train, inner_validation) in enumerate(k_fold):
        #print(inner_k)
        cov_inner_train = cov_train[inner_train, :]
        cov_inner_validation = cov_train[inner_validation, :]
        response_inner_train = res_train[inner_train]
        response_inner_validation = res_train[inner_validation]

        alphas_positive_enet, coefs_positive_enet, _ = linear_model.enet_path(cov_inner_train, response_inner_train, l1_ratio=l1_ratio,
                                                                              fit_intercept=False, normalize=False, positive=True, return_models=False)

        num_alphas = len(alphas_positive_enet)

        prediction_validation = np.dot(coefs_positive_enet.transpose(), cov_inner_validation.transpose())

        rep_res_val = np.repeat(response_inner_validation, num_alphas).reshape(len(response_inner_validation), num_alphas).transpose()
        rep_mean_val = np.repeat(np.mean(response_inner_validation), len(response_inner_validation)*num_alphas).reshape(len(response_inner_validation), num_alphas).transpose()

        sos_residual = np.sum((prediction_validation - rep_res_val) ** 2, axis=1)
        sos_original = np.sum((rep_res_val - rep_mean_val) ** 2, axis=1)

        rep_validation_rsqr = np.array(1 - (sos_residual / sos_original))

        sorted_ind = np.argsort(rep_validation_rsqr)[::-1]
        best_validation_rsqr[inner_k] = rep_validation_rsqr[sorted_ind[0]]
        best_validation_alpha[inner_k] = alphas_positive_enet[sorted_ind[0]]

    mean_best_alpha = np.mean(best_validation_alpha)
    #print(mean_best_alpha)
    #print(best_validation_rsqr)

    # now learn one unified model on the given data using the mean_best_alpha
    enet = linear_model.ElasticNet(l1_ratio=l1_ratio, alpha=mean_best_alpha, fit_intercept=False,
                                   normalize=False, positive=True)

    enet.fit(cov_train, res_train)

    return enet, best_validation_rsqr


def test(enet_model, cov_test, res_test, params={}):

    # update test covariates by the given params
    if 'covariates_prior' in params.keys() and params['covariates_prior'] is not None:
        cov_test = cov_test * params['covariates_prior']

    if 'class_subfeatures' in params.keys() and params['class_subfeatures'] is not None:
        # cases subfeatures
        cases = params['class_subfeatures'] > 0
        case_subfeatures = np.zeros_like(cov_test)
        case_subfeatures[cases, :] = cov_test[cases, :]
        # control subfeatures
        controls = params['class_subfeatures'] == 0
        control_subfeatures = np.zeros_like(cov_test)
        control_subfeatures[controls, :] = cov_test[controls, :]
        # update cov_train
        cov_test = np.hstack((case_subfeatures, control_subfeatures))

    prediction = enet_model.predict(cov_test)

    sos_residual = np.sum((prediction - res_test) ** 2)
    sos_original = np.sum((res_test - np.mean(res_test)) ** 2)

    test_rsqr = 1 - (sos_residual / sos_original)

    return prediction, test_rsqr
