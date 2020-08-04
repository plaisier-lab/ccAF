##########################################################
## ccAF:  classifier.py                                 ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

# Pandas gives lots of future warnings, suppressing them for now
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#########################################
## Load Python packages for classifier ##
#########################################

# General
import pandas as pd
import scanpy as sc
from scipy.sparse import isspmatrix

# ACTINN
import ccAF.actinn
import ccAF.actinn as actinn

#################
## Classifier ##
#################

class Classifier_ACTINN:
    """A class designed to facilitate ACTINN classifier construction
    and use. Can also be used for """
    def __init__(self, train, label, learning_rate = 0.0001, num_epochs = 200, minibatch_size = 128, print_cost = True, output_probability = False):
        self.train = train
        self.label = label
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.minibatch_size = minibatch_size
        self.print_cost = print_cost
        self.output_probability = output_probability
        self.label = label
        self.classifier, self.label_to_type_dict, self.genes = self.__build_classifier()

    # Prepare data
    def __prep_data(self, data):
        # Make indicies unique for
        data.var_names_make_unique()
        # Remove all genes with zero counts
        sc.pp.filter_genes(data, min_cells=1)
        if isspmatrix(data.X):
            return pd.DataFrame(data.X.todense(), index = data.obs_names, columns = data.var_names).T
        else:
            return pd.DataFrame(data.X, index = data.obs_names, columns = data.var_names).T
    
    # Prepare test data for predicting
    def __prep_predict_data(self, data):
        missing = set(self.genes).difference(data.index)
        if len(missing)>0:
            data = pd.concat([data, pd.DataFrame(0,index=missing, columns = data.columns)])
        return data.loc[list(self.genes)]

    # Build classifier
    def __build_classifier(self):
        train = self.train
        # Convert into pandas DataFrame
        train_data = self.__prep_data(train)
        labels = self.train.obs[self.label]
        clf, label_to_type_dict, genes = actinn.train_model(train_data, labels, learning_rate = self.learning_rate, num_epochs = self.num_epochs, minibatch_size = self.minibatch_size, print_cost = self.print_cost)
        return clf, label_to_type_dict, genes

    # Predict labels with rejection
    def predict_labels(self, newData):
        test_data = self.__prep_data(newData)
        pred_data = self.__prep_predict_data(test_data)
        labels = actinn.predict_new_data(pred_data, self.classifier, self.label_to_type_dict, self.genes)
        return list(labels['celltype'])

    # Predict labels with rejection
    def predict_probs(self, newData, axis = -1):
        test_data = self.__prep_data(newData)
        pred_data = self.__prep_predict_data(test_data)
        probs = actinn.predict_probabilities(pred_data, self.classifier, self.genes, axis=axis)
        probs.index = pd.Index([self.label_to_type_dict[i] for i in list(probs.index)])
        return probs

