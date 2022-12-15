#!/usr/bin/env python3

## Fewer nodes in inter layers. Less layers in u models.

import subprocess
import os
import itertools
import numpy as np
from sklearn.naive_bayes import GaussianNB
import tensorflow as tf
import random
import time
import math
# import autokeras as ak
import sys
import argparse

# import keras_tuner as kt

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, train_test_split
from tensorboard.plugins.hparams import api as hp
# from keras.models import load_model
# import pandas as pd

import lazypredict

# from lazypredict import Supervised
# Importing Regression class from lazypredcit
# from lazypredict.Supervised import LazyRegressor
# Importing Classification class from lazypredict
from lazypredict.Supervised import LazyClassifier

# Argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument('-m', type=str, required=True)
parser.add_argument('-na', type=str, required=True)
# parser.add_argument('-ml',type=str,required=True)
parser.add_argument('--h')
args = parser.parse_args()

if args.h:
    print('Herzlich Willkommen to ML_reconstruction program!')
    print('A program to reconstruct phylogenetic trees.')
    print('USAGE: python3 Simulations_and_ML_DNA_and_AA_final.py [-na <string>] [-m <string>]')
    print('Parameters: ')
    print('-na :\tNucleotide or Amino acids data. Type DNA for nucleotide. Type AA for Amino acids')
    print('-m:\tSelection of the evolution models.')
    print('\tList of the available nucleotide substitution models: JC,K2P,F81,F84,HKY,GTR')
    print('\tList of the available amino acids subtitution models: JTT, LG, WAG_OLD, WAG, WAG_STAR, DA')
    print('-h:\tDisplays usage information and exit')
    sys.exit(0)

mode_grid_search = 1

python_random_seed = 13

os.environ['MKL_NUM_THREADS'] = '10'
os.environ['GOTO_NUM_THREADS'] = '10'
os.environ['OMP_NUM_THREADS'] = '10'
os.environ['openmp'] = 'True'

# tf.config.threading.set_inter_op_parallelism_threads(10)
# tf.config.threading.set_intra_op_parallelism_threads(10)

# tf.compat.v1.Session(config
# config = tf.compat.v1.ConfigProto(device_count={"CPU": 120},
#                                  inter_op_parallelism_threads=120,
#                                  intra_op_parallelism_threads=120)
# sess = tf.compat.v1.Session(config=config)


dir = "/home/nkulikov/PycharmProjects/test-ml-1/"
tmpdir = "/dev/shm/Evolution_models/"
runID = str(os.getpid() * random.randint(1, 100000))

## Switches:
predict_mode = 0  ## 0: normal, 1: short internal branches,
random.seed(python_random_seed)

# global variables:
global_replicates = 1000
global_replicates_short_internalbranch = 1000
total_replicates = 3 * global_replicates + 3 * global_replicates_short_internalbranch
global_seqlen = 100000

prediction_replicates = 10
global_predict_seqlen = 1000

global_min_branchlength = 0.1
global_max_branchlength = 0.5
global_min_short_internal_branchlength = 0.001
global_max_short_internal_branchlength = 0.02

base_models_choice_list = ['JC', 'K2P', 'F81', 'F84', 'HKY', 'GTR']
base_models_choice_aa = ['JTT', 'LG', 'WAG_OLD', 'WAG', 'WAG_STAR', 'DAY']
# Parameters for nucleotide substitution models
# All models
global_min_pinv = 0
global_max_pinv = 0.50
global_min_shape = 0.01
global_max_shape = 4
# K2P, F84, HKY, GTR
global_min_tstv = 1.0
global_max_tstv = 3.0
# F81, F84, HKY, GTR
global_min_basefreq = 0.2
global_max_basefreq = 0.3
# GTR
global_min_rrates_GTR = 0.1
global_max_rrates_GTR = 1.0

iqtree_predict_mode = 0

print("Parameters used:")
if iqtree_predict_mode == 1:
    print("Testing reconstruction success with ML, MP, NJ, not with ANNs.")

print("prediction_replicates_DNA:                  ", prediction_replicates)
print("prediction_replicates_AA:                  ", prediction_replicates)
print("global_predict_seqlen_DNA:                  ", global_predict_seqlen)
print("global_predict_seqlen_AA:                  ", global_predict_seqlen)
print("global_min_branchlength:                ", global_min_branchlength)
print("global_max_branchlength:                ", global_max_branchlength)
print("global_min_short_internal_branchlength: ", global_min_short_internal_branchlength)
print("global_max_short_internal_branchlength: ", global_max_short_internal_branchlength)

print("global_min_pinv:                        ", global_min_pinv)
print("global_max_pinv:                        ", global_max_pinv)
print("global_min_shape:                       ", global_min_shape)
print("global_max_shape:                       ", global_max_shape)
print("base_models_choice_list:                ", base_models_choice_list)
print("base_models_choice_aa:                ", base_models_choice_aa)


##### Grid search parameters:

# HP_NUM_UNITS = hp.HParam('num_units', hp.Discrete([16, 32]))
# HP_DROPOUT = hp.HParam('dropout', hp.RealInterval(0.1, 0.2))

# lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
#    initial_learning_rate=1e-2,
#    decay_steps=10000,
#    decay_rate=0.9)


# opt_SGD_1   = tf.keras.optimizers.SGD()
# opt_SGD_2   = tf.keras.optimizers.SGD(lr=0.001)
# opt_SGD_3   = tf.keras.optimizers.SGD(lr=0.005)
# opt_SGD_4   = tf.keras.optimizers.SGD(learning_rate=lr_schedule)
# opt_SGD_5   = tf.keras.optimizers.SGD(learning_rate=0.0005)
# opt_SGD_6   = tf.keras.optimizers.SGD(learning_rate=0.0003)
# opt_SGD_7   = tf.keras.optimizers.SGD(learning_rate=0.0001)
# opt_SGD_8   = tf.keras.optimizers.SGD(learning_rate=0.0005, momentum=0.9)
# opt_SGD_9   = tf.keras.optimizers.SGD(learning_rate=0.0005, momentum=0.9, nesterov=True)
# ad   = tf.keras.optimizers.Adam()

def get_optimizer(name):
    if name == "opt_SGD_1":
        return opt_SGD_1
    if name == "opt_SGD_2":
        return opt_SGD_2
    if name == "opt_SGD_3":
        return opt_SGD_3
    if name == "opt_SGD_4":
        return opt_SGD_4

    if name == "opt_SGD_5":
        return opt_SGD_5
    if name == "opt_SGD_6":
        return opt_SGD_6
    if name == "opt_SGD_7":
        return opt_SGD_7
    if name == "opt_SGD_8":
        return opt_SGD_8
    if name == "opt_SGD_9":
        return opt_SGD_9
    if name == "ad":
        return ad


### HPARAM
# HP_ACTIVATION = hp.HParam('activation',   hp.Discrete(['tanh']))
# HP_DROPOUT    = hp.HParam('dropout',      hp.Discrete([0.3]))
# HP_OPTIMIZER  = hp.HParam('optimizer',    hp.Discrete(['ad']))
# METRIC_ACCURACY = 'accuracy'

# HPARAMS = [
#    HP_LAYERS,
#    HP_KERNEL_SIZE,
#    HP_DENSE_LAYERS,
#    HP_ACTIVATION,
#    HP_DROPOUT,
#    HP_OPTIMIZER,
# ]

##########

# param_grid = {
#                'batch_size':           np.asarray([256]),
#                'epochs' :              np.asarray([50]),
# #               'hidden_layers' :      np.asarray([(136,68,34), # best
# #                                        (68,34,17),
# #                                        (272,136,68),
# #                                        (136,68,34,17),
# #                                        (136,68,34,17,8),]),
#                'optimizer' :           np.asarray([tf.keras.optimizers.Adam(lr=0.001, decay=0.001, epsilon=1e-4),
#                                                    tf.keras.optimizers.Nadam(lr=0.001, decay=0.001, epsilon=1e-4),
#                                                    tf.keras.optimizers.RMSprop(),  ## Important: strange error message if parentheses are left out.
#                                                    tf.keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.99, epsilon=1e-4, amsgrad=False),
#                                                    tf.keras.optimizers.Adamax(learning_rate=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-04)  ]),
#                'dropout_rate' :        np.asarray([0.3]),       #, 0.4
#                'activation' :          np.asarray([tf.keras.layers.ReLU()])     #   ('relu', 'elu', 'softplus')
#              }

EPOCHS = 1000
LEN = 160000


### Evaluation routines: (moved here in Test2):
def highest_likelihood(predictions):
    topo = []
    for prediction in predictions:
        if prediction[0] > prediction[1]:
            if prediction[0] > prediction[2]:
                topo.append(0)
            else:  # predict[2] >= prediction[0] > prediction[1]:
                topo.append(2)
        else:  # predict[1] >= prediction[0]
            if prediction[1] > prediction[2]:
                topo.append(1)
            else:  # prediction[2] >= predict[1] >= prediction[0]
                topo.append(2)
    return topo


def which_topology(predictions, threshold):
    topo = []
    for prediction in predictions:
        if prediction[0] > threshold:
            topo.append(0)
        elif prediction[1] > threshold:
            topo.append(1)
        elif prediction[2] > threshold:
            topo.append(2)
        else:
            topo.append("?")
    return topo


def count_topologies(predictions, threshold):
    topo = which_topology(predictions, threshold)
    return [topo.count(0) / len(topo), topo.count(1) / len(topo), topo.count(2) / len(topo)]


def count_tolologies_by_likelihood(predictions):
    topo = highest_likelihood(predictions)
    return [topo.count(0) / len(topo), topo.count(1) / len(topo), topo.count(2) / len(topo)]


def compute_mean_accuracy_of_predictions(predictions, topology_type):
    prediction_list = []
    for prediction in predictions:
        prediction_list.append(prediction[topology_type])
    L = [float(prediction) for prediction in prediction_list if prediction]
    acc = (sum(L) / len(L) if L else '-') * 100
    return acc


def print_highest_likelihood_evaluation(model, d0, d1, d2):
    pred0 = model.predict(d0)
    pred1 = model.predict(d1)
    pred2 = model.predict(d2)
    print(count_tolologies_by_likelihood(pred0))
    print(count_tolologies_by_likelihood(pred1))
    print(count_tolologies_by_likelihood(pred2))


def find_frequency(freq_file):
    with open(freq_file, "r") as frequency_file:
        frequency_list = []
        #        frequency_file.readline()
        for line in frequency_file:
            column = line.split("\t")
            if len(column[0]) == 4:
                frequency_list.append(np.float32(column[1]))
            else:
                f_e = open('file_error.txt', 'a')
                f_e.write(line)
                f_e.write('\n')
                f_e.close()
    return frequency_list


def make_tree_file(file_name, topology, B1, B2, B3, B4, B5, B6, seqlen):
    with open(file_name, "w") as tree_file:
        if topology == 0:
            string = "1.0 {0} Model_JC ((A:{1},B:{2}):{3},(C:{4},D:{5}):{6})\n".format(seqlen, B1, B2, B3, B4, B5, B6)
        if topology == 1:
            string = "1.0 {0} Model_JC ((A:{1},C:{2}):{3},(B:{4},D:{5}):{6})\n".format(seqlen, B1, B2, B3, B4, B5, B6)
        if topology == 2:
            string = "1.0 {0} Model_JC ((A:{1},D:{2}):{3},(C:{4},B:{5}):{6})\n".format(seqlen, B1, B2, B3, B4, B5, B6)
        tree_file.write(string)


def make_tree_file_model(file_name, topology, B1, B2, B3, B4, B5, B6, seqlen, modelname):
    with open(file_name, "w") as tree_file:
        if topology == 0:
            string = "1.0  {0} {1} ((A:{2},B:{3}):{4},(C:{5},D:{6}):{7})\n".format(seqlen, modelname, B1, B2, B3, B4,
                                                                                   B5, B6)
        if topology == 1:
            string = "1.0  {0} {1} ((A:{2},C:{3}):{4},(B:{5},D:{6}):{7})\n".format(seqlen, modelname, B1, B2, B3, B4,
                                                                                   B5, B6)
        if topology == 2:
            string = "1.0  {0} {1} ((A:{2},D:{3}):{4},(C:{5},B:{6}):{7})\n".format(seqlen, modelname, B1, B2, B3, B4,
                                                                                   B5, B6)
        print(string)
        tree_file.write(string)


def make_model_file(file_name, model, alpha, pinv, tstv, basefreq_1, basefreq_2, GTR_A_C, GTR_A_G, GTR_A_T, GTR_C_G,
                    GTR_C_T, GTR_G_T, model_name):
    with open(file_name, "w") as model_file:
        if model == 'JC':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'K2P':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'tstv: {0}\n'.format(tstv)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'F81':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1, basefreq_2, basefreq_1, basefreq_2)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'HKY':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'tstv: {0}\n'.format(tstv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1, basefreq_2, basefreq_1, basefreq_2)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'F84':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'tstv: {0}\n'.format(tstv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1, basefreq_2, basefreq_1, basefreq_2)
            model_file.write(string)
            model_file.write("end model\n")
        elif model == 'GTR':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'tstv: {0}\n'.format(tstv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1, basefreq_2, basefreq_1, basefreq_2)
            model_file.write(string)
            string = 'rrates: {} {} {} {} {} {}\n'.format(GTR_A_C, GTR_A_G, GTR_A_T, GTR_C_G, GTR_C_T, GTR_G_T)
            model_file.write(string)
            model_file.write("end model\n")
        elif model in base_models_choice_aa:
            model_file.write("begin aa-model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n '.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            model_file.write("end model\n")
            print("Model file written:")


def simulate_random_lengths_random_free_model(topology, seed, base_model_list, min_brlen, max_brlen, min_internal_brlen,
                                              max_internal_brlen, min_pinv, max_pinv,
                                              min_shape, max_shape, seqlen, base_model_list_aa, global_min_tstv,
                                              global_max_tstv,
                                              global_min_rrates_GTR, global_max_rrates_GTR):
    model_filename = tmpdir + "Modelfile_free" + runID + "_" + str(seed) + ".txt"
    seq_filename = tmpdir + "sim-sequences" + runID + "_" + str(seed) + ".txt"
    sim_filename = tmpdir + "sim_tree" + runID + "_" + str(seed) + ".txt"

    if args.na == 'DNA':
        model = args.m
        modelname = model + '_model'
    elif args.na == 'AA':
        model = args.m
        modelname = model + '_model'
        # if model != "JC":
        #    print("So far, only the JC model is implemented in this function.\n")
        #   sys.exit(1)

    shape = random.uniform(min_shape, max_shape)
    pinv = random.uniform(min_pinv, max_pinv)

    r1 = random.uniform(min_brlen, max_brlen)
    r2 = random.uniform(min_brlen, max_brlen)
    r3 = random.uniform(min_internal_brlen, max_internal_brlen)
    r4 = random.uniform(min_brlen, max_brlen)
    r5 = random.uniform(min_brlen, max_brlen)
    #    r6 = random.uniform(min_brlen,max_brlen)
    r3 = r3 / 2
    r6 = r3
    tstv = random.uniform(global_min_tstv, global_max_tstv)
    basefreq_1 = np.random.normal(0.5, 0.03, 1) / 2
    while basefreq_1 > 0.3 or basefreq_1 < 0.2:
        basefreq_1 = np.random.normal(0.5, 0.03, 1) / 2
    basefreq_2 = (1 - 2 * basefreq_1) / 2
    basefreq_1 = float(basefreq_1)
    basefreq_2 = float(basefreq_2)
    GTR_A_C = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_A_G = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_A_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_C_G = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_C_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_G_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)

    print("Calling make_model_file: ", model_filename, " ", model, " ", shape, " ", pinv, " ", modelname)
    make_model_file(model_filename, model, shape, pinv, tstv, basefreq_1, basefreq_2, GTR_A_C, GTR_A_G, GTR_A_T,
                    GTR_C_G, GTR_C_T, GTR_G_T, modelname)
    make_tree_file_model(sim_filename, topology, r1, r2, r3, r4, r5, r6, seqlen, modelname)

    all_sim_filename = "all_sim_tree" + runID + ".txt"
    os.system("echo " + str(seed) + " >> " + all_sim_filename)
    os.system("cat " + sim_filename + " >> " + all_sim_filename)
    command = "PolyMoSim-v1.1.3b -s " + str(
        seed) + " -m " + model_filename + " -t " + sim_filename + " -f site_pattern_freq_relative_fill -n 1 1> " + seq_filename
    subprocess.run([command], shell=True)
    individual_frequencies = find_frequency(seq_filename)
#    if len(individual_frequencies) != 256:
#        print('Sequence file_name: ', seq_filename)
#        print('Shape of the temporary array: ', len(individual_frequencies))
#        sys.exit(0)
    os.remove(seq_filename)
    os.remove(sim_filename)
    os.remove(model_filename)
    return individual_frequencies


####################################
# my main program starts here:
####################################

topology_list = []

topology_array = np.empty((total_replicates, 1))
frequency_array = np.empty((total_replicates, 256))

print(topology_array.shape)
print(frequency_array.shape)

### Frequency files:
if (args.m == 'JC'):
    global_seqlen = 100000
    global_replicates = 100000
    global_replicates_short_internalbranch = 100000
    frequency_simulations_filename = 'trainFrequency' + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.m in base_models_choice_aa:
    global_seqlen = 1000000
    global_replicates = 1000
    global_replicates_short_internalbranch = 1000
    frequency_simulations_filename = 'trainFrequency' + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.m == 'K2P':
    global_seqlen = 100000
    global_replicates = 100000
    global_replicates_short_internalbranch = 100000
    frequency_simulations_filename = 'trainFrequency' + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
        global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.m == 'F81':
    global_seqlen = 100000
    global_replicates = 100000
    global_replicates_short_internalbranch = 100000
    frequency_simulations_filename = 'trainFrequency' + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.m == 'F84' or args.m == 'HKY':
    global_seqlen = 100000
    global_replicates = 100000
    global_replicates_short_internalbranch = 100000
    frequency_simulations_filename = 'trainFrequency' + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.m == 'GTR':
    global_seqlen = 100000
    global_replicates = 100000
    global_replicates_short_internalbranch = 100000
    frequency_simulations_filename = 'trainFrequency' + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) + "_pseed_" + str(
        python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.m + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(
        global_max_rrates_GTR) + "_pseed_" + str(python_random_seed) + ".npy"
print("")
print("Filenames training:")
print(frequency_simulations_filename)
print(topology_simulations_filename)
print("")
sys.stdout.flush()

# time before simulate
timeBeforeSimulate = time.time()
if os.path.isfile(frequency_simulations_filename):
    frequency_array = np.load(frequency_simulations_filename)
    topology_array = np.load(topology_simulations_filename)

else:
    if args.na == 'DNA':
        global_seqlen = 100000
        global_replicates = 100000
        global_replicates_short_internalbranch = 100000
        total_replicates = 3 * global_replicates + 3 * global_replicates_short_internalbranch
        frequency_array = np.empty((total_replicates, 256))
        for i in range(0, global_replicates):
            #        temporary_array = simulate_random_lengths_randomJK_model(0,i,global_min_branchlength,global_max_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(0)
        for i in range(global_replicates, 2 * global_replicates):
            #        temporary_array = simulate_random_lengths_randomJK_model(1,i,global_min_branchlength,global_max_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa,
                                                                        global_min_tstv, global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(1)
        for i in range(2 * global_replicates, 3 * global_replicates):
            #        temporary_array = simulate_random_lengths_randomJK_model(2,i,global_min_branchlength,global_max_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, \
                                                                        global_min_tstv, global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)
            frequency_array[i] = temporary_array
            topologies = topology_list.append(2)

        for i in range(3 * global_replicates, 3 * global_replicates + global_replicates_short_internalbranch):
            #        temporary_array = simulate_random_short_internal_branch_randomJK_model(0,i,global_min_branchlength,global_max_branchlength,global_min_short_internal_branchlength,global_max_short_internal_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_short_internal_branchlength,
                                                                        global_max_short_internal_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa,
                                                                        global_min_tstv, global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(0)
        for i in range(3 * global_replicates + global_replicates_short_internalbranch,
                       3 * global_replicates + 2 * global_replicates_short_internalbranch):
            #        temporary_array = simulate_random_short_internal_branch_randomJK_model(1,i,global_min_branchlength,global_max_branchlength,global_min_short_internal_branchlength,global_max_short_internal_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_short_internal_branchlength,
                                                                        global_max_short_internal_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(1)
        for i in range(3 * global_replicates + 2 * global_replicates_short_internalbranch,
                       3 * global_replicates + 3 * global_replicates_short_internalbranch):
            #        temporary_array = simulate_random_short_internal_branch_randomJK_model(2,i,global_min_branchlength,global_max_branchlength,global_min_short_internal_branchlength,global_max_short_internal_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_short_internal_branchlength,
                                                                        global_max_short_internal_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(2)

        topology_array = topology_list

        np.save(frequency_simulations_filename, frequency_array)
        np.save(topology_simulations_filename, topology_array)



    elif args.na == 'AA':
        global_seqlen = 1000000
        global_replicates = 10
        global_replicates_short_internalbranch = 10
        total_replicates = 3 * global_replicates + 3 * global_replicates_short_internalbranch
        frequency_array = np.empty((total_replicates, 160000))
        for i in range(0, global_replicates):
            #        temporary_array = simulate_random_lengths_randomJK_model(0,i,global_min_branchlength,global_max_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(0)
        for i in range(global_replicates, 2 * global_replicates):
            #        temporary_array = simulate_random_lengths_randomJK_model(1,i,global_min_branchlength,global_max_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(1)
        for i in range(2 * global_replicates, 3 * global_replicates):
            #        temporary_array = simulate_random_lengths_randomJK_model(2,i,global_min_branchlength,global_max_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(2)
        for i in range(3 * global_replicates, 3 * global_replicates + global_replicates_short_internalbranch):
            #        temporary_array = simulate_random_short_internal_branch_randomJK_model(0,i,global_min_branchlength,global_max_branchlength,global_min_short_internal_branchlength,global_max_short_internal_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_short_internal_branchlength,
                                                                        global_max_short_internal_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(0)

        for i in range(3 * global_replicates + global_replicates_short_internalbranch,
                       3 * global_replicates + 2 * global_replicates_short_internalbranch):
            #        temporary_array = simulate_random_short_internal_branch_randomJK_model(1,i,global_min_branchlength,global_max_branchlength,global_min_short_internal_branchlength,global_max_short_internal_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_short_internal_branchlength,
                                                                        global_max_short_internal_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(1)
        for i in range(3 * global_replicates + 2 * global_replicates_short_internalbranch,
                       3 * global_replicates + 3 * global_replicates_short_internalbranch):
            #        temporary_array = simulate_random_short_internal_branch_randomJK_model(2,i,global_min_branchlength,global_max_branchlength,global_min_short_internal_branchlength,global_max_short_internal_branchlength)
            temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                        global_min_branchlength,
                                                                        global_max_branchlength,
                                                                        global_min_short_internal_branchlength,
                                                                        global_max_short_internal_branchlength,
                                                                        global_min_pinv, global_max_pinv,
                                                                        global_min_shape,
                                                                        global_max_shape, global_seqlen,
                                                                        base_models_choice_aa, global_min_tstv,
                                                                        global_max_tstv,
                                                                        global_min_rrates_GTR, global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(2)
        topology_array = topology_list
        np.save(frequency_simulations_filename, frequency_array)
        np.save(topology_simulations_filename, topology_array)

print('Shape of the frequency array: ', frequency_array.shape)
print(frequency_array)

print('Shape of the topology array: ', frequency_array.shape)
print(topology_array)
# time after
timeAfterSimulate = time.time()

temporary_0 = np.empty((prediction_replicates, 256))
temporary_1 = np.empty((prediction_replicates, 256))
temporary_2 = np.empty((prediction_replicates, 256))

# time before simulate data
timeBeforeSimulateData = time.time()

python_random_seed += 1
random.seed(python_random_seed)

if predict_mode == 0:
    if (args.m == 'JC') or (args.m in base_models_choice_aa):
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"
    elif args.m == 'K2P':
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"
    elif args.m == 'F81':
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"
    elif args.m == 'F84' or args.m == 'HKY':
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"
    elif args.m == 'GTR':
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(
            global_max_rrates_GTR) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(
            global_max_rrates_GTR) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(
            global_max_rrates_GTR) + "_pseed_" + str(python_random_seed) + ".npy"
if predict_mode == 1:
    if (args.m == 'JC') or (args.m in base_models_choice_aa):
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(
            python_random_seed) + ".npy"
    elif args.m == 'K2P':
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(
            python_random_seed) + ".npy"
    elif args.m == 'F81':
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(
            global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(
            global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(
            global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"
    elif args.m == 'F84' or args.m == 'HKY':
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"
    elif args.m == 'GTR':
        temporary0_filename = "predict0_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(
            global_max_rrates_GTR) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(
            global_max_rrates_GTR) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.m + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
            global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(
            global_max_rrates_GTR) + "_pseed_" + str(
            python_random_seed) + ".npy"

print("Filenames for prediction files:")
print("predict0: ", temporary0_filename)
print("predict1: ", temporary1_filename)
print("predict2: ", temporary2_filename)
sys.stdout.flush()

#    temporary0_filename = "predict0_mode_randMod" + str(predict_mode) + "_" + str(prediction_replicates) + "_" + str(global_min_branchlength) + "_" + str(global_max_branchlength) + "__" + str(global_min_short_internal_branchlength) + "_" + str(global_max_short_internal_branchlength) + "_pseed_" + str(python_random_seed)+ ".npy"
#    temporary1_filename = "predict1_mode_randMod" + str(predict_mode) + "_" + str(prediction_replicates) + "_" + str(global_min_branchlength) + "_" + str(global_max_branchlength) + "__" + str(global_min_short_internal_branchlength) + "_" + str(global_max_short_internal_branchlength) + "_pseed_" + str(python_random_seed)+ ".npy"
#    temporary2_filename = "predict2_mode_randMod" + str(predict_mode) + "_" + str(prediction_replicates) + "_" + str(global_min_branchlength) + "_" + str(global_max_branchlength) + "__" + str(global_min_short_internal_branchlength) + "_" + str(global_max_short_internal_branchlength) + "_pseed_" + str(python_random_seed)+ ".npy"


if os.path.isfile(temporary0_filename):
    temporary_0 = np.load(temporary0_filename, allow_pickle=True)
    temporary_1 = np.load(temporary1_filename, allow_pickle=True)
    temporary_2 = np.load(temporary2_filename, allow_pickle=True)
else:
    if args.na == 'DNA':
        temporary_0 = np.empty((prediction_replicates, 256))
        temporary_1 = np.empty((prediction_replicates, 256))
        temporary_2 = np.empty((prediction_replicates, 256))
        if predict_mode == 0:
            for j in range(0, prediction_replicates):
                temporary_0[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)
                temporary_1[j] = simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)
                temporary_2[j] = simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)

        if predict_mode == 1:
            for j in range(0, prediction_replicates):
                temporary_0[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_short_internal_branchlength,
                                                                           global_max_short_internal_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)

                temporary_1[j] = simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_short_internal_branchlength,
                                                                           global_max_short_internal_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)
                temporary_2[j] = simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_short_internal_branchlength,
                                                                           global_max_short_internal_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)

        np.save(temporary0_filename, temporary_0)
        np.save(temporary1_filename, temporary_1)
        np.save(temporary2_filename, temporary_2)

    if args.na == 'AA':
        temporary_0 = np.empty((prediction_replicates, 160000))
        temporary_1 = np.empty((prediction_replicates, 160000))
        temporary_2 = np.empty((prediction_replicates, 160000))
        if predict_mode == 0:
            for j in range(0, prediction_replicates):
                temporary_0[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)

                temporary_1[j] = simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)
                temporary_2[j] = simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)

        if predict_mode == 1:

            for j in range(0, prediction_replicates):
                temporary_0[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_short_internal_branchlength,
                                                                           global_max_short_internal_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)

                temporary_1[j] = simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_short_internal_branchlength,
                                                                           global_max_short_internal_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)
                temporary_2[j] = simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                           global_min_branchlength,
                                                                           global_max_branchlength,
                                                                           global_min_short_internal_branchlength,
                                                                           global_max_short_internal_branchlength,
                                                                           global_min_pinv, global_max_pinv,
                                                                           global_min_shape, global_max_shape,
                                                                           global_predict_seqlen, base_models_choice_aa,
                                                                           global_min_tstv, global_max_tstv,
                                                                           global_min_rrates_GTR, global_max_rrates_GTR)

        np.save(temporary0_filename, temporary_0)
        np.save(temporary1_filename, temporary_1)
        np.save(temporary2_filename, temporary_2)

print('Shape of the temporary_0: ', temporary_0)
print(temporary_0)

print('Shape of the temporary_1: ', temporary_1)
print(temporary_1)

print('Shape of the temporary_2: ', temporary_2)
print(temporary_2)

# ANN training starts here

sc = StandardScaler()
frequency_array = sc.fit_transform(frequency_array)
# frequency_train = sc.fit_transform(frequency_train)
# frequency_test  = sc.transform(frequency_test)

frequency_train, frequency_test, topology_train, topology_test = train_test_split(frequency_array, topology_array,
                                                                                  test_size=0.04, random_state=42)

frequency_train = np.asarray(frequency_train)
topology_train = np.asarray(topology_train)

frequency_test = np.asarray(frequency_test)
topology_test = np.asarray(topology_test)

# topology_train = tf.keras.utils.to_categorical(topology_train, 3)
# topology_test = tf.keras.utils.to_categorical(topology_test, 3)

print('Shape of the frequency_train: ', frequency_train.shape)
print('Shape of the frequency_test: ', frequency_test.shape)
print('Shape of the topology_train: ', topology_train.shape)
print('Shape of the topology_test: ', topology_test.shape)

# topology_train = topology_train.flatten()
# topology_test = topology_test.flatten()

# print('Shape of the topology_train: ', topology_train.shape)
# print('Shape of the topology_test: ', topology_test.shape)

# lazypredict

# multiple_ML_model = LazyClassifier(verbose=0,ignore_warnings=True,predictions=True)
multiple_ML_model = LazyClassifier(verbose=0, ignore_warnings=False, custom_metric=None, predictions=True)

models, predictions = multiple_ML_model.fit(frequency_train, frequency_test, topology_train, topology_test)

print(models)