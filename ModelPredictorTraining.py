#!/usr/bin/env python3

## Fewer nodes in inter layers. Less layers in u models.

import subprocess
import os
import itertools
import numpy as np
from   sklearn.naive_bayes import GaussianNB
import tensorflow as tf
import random
import time
import math
#import autokeras as ak
import sys
import argparse

#import keras_tuner as kt

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, train_test_split
from tensorboard.plugins.hparams import api as hp



# Argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument('-sequence_type', type=str, default='DNA', required=False)
parser.add_argument('-neural_network', type=str, default='b1', required=False)
parser.add_argument('-seqlen', type=int, default=100000, required=False)
parser.add_argument('-seqnum', type=int, default=100000, required=False)
parser.add_argument('-epochs', type=int, default=500, required=False)
parser.add_argument('-predict_mode', type=int, default=0, required=False)
parser.add_argument('-predict_seqlen', type=int, default=30000, required=False)
parser.add_argument('-predict_seqnum', type=int, default=1000, required=False)
parser.add_argument('-info', type=int, required=False)

args = parser.parse_args()

if args.info == 1:
    print('Welcome to program that predicts substitution models!')
    print('USAGE: python3 ModelPredictorTraining.py [-sequence-type <string>] [-neural-network <string>] \n [-sequences-length <integer>], [-reps <integer>] [-EPOCHS <integer>]')
    print('Optional parameters: ')
    print('-sequence_type :\tNucleotide or Amino acids data. Enter DNA for nucleotide sequences. Enter AA for amino acid sequences. Default: DNA ')
    print('-neural_network:\tDefault: b1 ')
    print('List of the available neural networks: b1, b2, b3, b10, u, cu, cud')
    print('-seqlen :\tThe length of the training sequences. Default: 100000')
    print('-seqnum :\tThe number of sequences used for training. Actual number is multiplied by 4 for amino acids and by 5 for nucleotide sequences. Default: 100000')
    print('-epochs :\tThe number of epochs is a hyperparameter that defines the number times that the learning algorithm will work through the entire training dataset. Default: 500 ')
    print('-predict_seqlen :\tThe length of the prediction sequences. Default: 30000 ')
    print('-predict_seqnum :\tThe number of the prediction sequences. Default: 1000 ')
    print('-predict_mode :\tType of predicted sequences either normal internal branches(0) or short internal branches(1). Default: 0')
    print('List of the available nucleotide substitution models: JC,K2P,F81,HKY,GTR')
    print('List of the available amino acids substitution models: JTT, LG, WAG, DAY')

    print('-info:\tDisplays usage information and exit')
    sys.exit(0)

python_random_seed = 13

os.environ['MKL_NUM_THREADS'] = '10'
os.environ['GOTO_NUM_THREADS'] = '10'
os.environ['OMP_NUM_THREADS'] = '10'
os.environ['openmp'] = 'True'

tf.config.threading.set_inter_op_parallelism_threads(10)
tf.config.threading.set_intra_op_parallelism_threads(10)

tmpdir = os.getcwd() + "/"
runID = str(os.getpid()*random.randint(1,100000))

## Switches:
predict_mode = args.predict_mode ## 0: normal, 1: short internal branches,
random.seed(python_random_seed)

#global variables:

global_seqlen = args.seqlen
global_replicates = args.seqnum
global_replicates_short_internalbranch = args.seqnum
prediction_replicates = args.predict_seqnum
global_predict_seqlen = args.predict_seqlen

total_replicates = 3*global_replicates + 3*global_replicates_short_internalbranch

if args.sequence_type == 'DNA' :
    frequency_array = np.empty((5 * total_replicates, 256))
    evo_mod_array = np.empty((5 * total_replicates, 1))
if args.sequence_type == 'AA' :
    frequency_array = np.empty((4 * total_replicates, 160000))
    evo_mod_array = np.empty((4 * total_replicates, 1))

global_min_branchlength = 0.1
global_max_branchlength = 0.5
global_min_short_internal_branchlength = 0.001
global_max_short_internal_branchlength = 0.02


base_models_choice_list = ['JC','K2P','F81','HKY','GTR']
base_models_choice_aa = ['JTT', 'LG', 'WAG', 'DAY']
# Parameters for nucleotide substitution models
# All models
global_min_pinv   = 0
global_max_pinv   = 0.50
global_min_shape  = 0.01
global_max_shape  = 4
# K2P, F84, HKY, GTR
global_min_tstv = 1.0
global_max_tstv = 3.0
# F81, F84, HKY, GTR
global_min_basefreq = 0.2
global_max_basefreq = 0.3
# GTR
global_min_rrates_GTR = 0.1
global_max_rrates_GTR = 1.0

print("Parameters used:")


print("prediction_replicates_DNA:                  ",prediction_replicates)
print("prediction_replicates_AA:                  ",prediction_replicates)
print("global_predict_seqlen_DNA:                  ",global_predict_seqlen)
print("global_predict_seqlen_AA:                  ",global_predict_seqlen)
print("global_min_branchlength:                ",global_min_branchlength)
print("global_max_branchlength:                ",global_max_branchlength)
print("global_min_short_internal_branchlength: ",global_min_short_internal_branchlength)
print("global_max_short_internal_branchlength: ",global_max_short_internal_branchlength )

print("global_min_pinv:                        ", global_min_pinv)
print("global_max_pinv:                        ", global_max_pinv)
print("global_min_shape:                       ", global_min_shape)
print("global_max_shape:                       ", global_max_shape)
print("base_models_choice_list:                ", base_models_choice_list)
print("base_models_choice_aa:                ", base_models_choice_aa)



##### Grid search parameters:


lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
    initial_learning_rate=1e-2,
    decay_steps=10000,
    decay_rate=0.9)

ad   = tf.keras.optimizers.Adam()

def get_optimizer(name):

    if name == "ad":
        return ad

### HPARAM
HP_ACTIVATION = hp.HParam('activation',   hp.Discrete(['tanh']))
HP_DROPOUT    = hp.HParam('dropout',      hp.Discrete([0.3]))
HP_OPTIMIZER  = hp.HParam('optimizer',    hp.Discrete(['ad']))
METRIC_ACCURACY = 'accuracy'

HPARAMS = [
    HP_ACTIVATION,
    HP_DROPOUT,
    HP_OPTIMIZER,
]


### Evaluation routines: (moved here in Test2):
def highest_likelihood (predictions):
    topo = []
    if args.sequence_type == 'DNA':
        for prediction in predictions:
            if prediction[0] > prediction[1] and prediction[0] > prediction[2] and prediction[0] > prediction[3] and prediction[0] > prediction[4]:
                topo.append('JC')
            if prediction[1] > prediction[0] and prediction[1] > prediction[2] and prediction[1] > prediction[3] and prediction[1] > prediction[4]:
                topo.append('K2P')
            if prediction[2] > prediction[0] and prediction[2] > prediction[1] and prediction[2] > prediction[3] and prediction[2] > prediction[4]:
                topo.append('F81')
            if prediction[3] > prediction[0] and prediction[3] > prediction[1] and prediction[3] > prediction[2] and prediction[3] > prediction[4]:
                topo.append('HKY')
            if prediction[4] > prediction[0] and prediction[4] > prediction[1] and prediction[4] > prediction[2] and prediction[4] > prediction[3]:
                topo.append('GTR')
    if args.sequence_type == 'AA':
        for prediction in predictions:
            if prediction[0] > prediction[1] and prediction[0] > prediction[2] and prediction[0] > prediction[3] and \
                    prediction[0] > prediction[4]:
                topo.append('JTT')
            if prediction[1] > prediction[0] and prediction[1] > prediction[2] and prediction[1] > prediction[3] and \
                    prediction[1] > prediction[4]:
                topo.append('LG')
            if prediction[2] > prediction[0] and prediction[2] > prediction[1] and prediction[2] > prediction[3] and \
                    prediction[2] > prediction[4]:
                topo.append('WAG')
            if prediction[4] > prediction[0] and prediction[4] > prediction[1] and prediction[4] > prediction[2] and \
                    prediction[4] > prediction[3]:
                topo.append('DAY')

    return topo

def count_evo_mod_by_likelihood(predictions):
    topo = highest_likelihood (predictions)
    if args.sequence_type == 'DNA' :
        return [topo.count('JC')/len(topo), topo.count('K2P')/len(topo), topo.count('F81')/len(topo), topo.count('HKY')/len(topo), topo.count('GTR')/len(topo)]
    if args.sequence_type == 'AA' :
        return [topo.count('JTT')/len(topo), topo.count('LG')/len(topo), topo.count('WAG')/len(topo), topo.count('DAY')/len(topo)]


def print_highest_likelihood_evaluation(model, d0, d1, d2,d3,d4):
    pred0 = model.predict(d0)
    pred1 = model.predict(d1)
    pred2 = model.predict(d2)
    pred3 = model.predict(d3)
    print(count_evo_mod_by_likelihood(pred0))
    print(count_evo_mod_by_likelihood(pred1))
    print(count_evo_mod_by_likelihood(pred2))
    print(count_evo_mod_by_likelihood(pred3))
    if args.sequence_type == 'DNA':
        pred4 = model.predict(d4)
        print(count_evo_mod_by_likelihood(pred4))


def find_frequency(freq_file):
    with open(freq_file,"r") as frequency_file:
        frequency_list = []
        for line in frequency_file:
            column = line.split("\t")
            if len(column[0]) == 4:
                frequency_list.append(column[1])
            else :
                f_e = open('file_error.txt', 'a')
                f_e.write(line)
                f_e.write('\n')
                f_e.close()
    return frequency_list


def make_tree_file_model(file_name,topology,B1,B2,B3,B4,B5,B6, seqlen, modelname):
    with open(file_name,"w") as tree_file:
        if topology == 0:
            string = "1.0  {0} {1} ((A:{2},B:{3}):{4},(C:{5},D:{6}):{7})\n".format(seqlen, modelname, B1,B2,B3,B4,B5,B6)
        if topology ==1 :
            string = "1.0  {0} {1} ((A:{2},C:{3}):{4},(B:{5},D:{6}):{7})\n".format(seqlen, modelname, B1,B2,B3,B4,B5,B6)
        if topology ==2:
            string = "1.0  {0} {1} ((A:{2},D:{3}):{4},(C:{5},B:{6}):{7})\n".format(seqlen, modelname, B1,B2,B3,B4,B5,B6)
        print(string)
        tree_file.write(string)


def make_model_file(file_name, model, alpha, pinv, tstv,basefreq_1,basefreq_2,GTR_A_C,GTR_A_G,GTR_A_T,GTR_C_G,GTR_C_T,GTR_G_T, model_name):
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
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1,basefreq_2,basefreq_1,basefreq_2)
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
        elif model == 'GTR':
            model_file.write("begin model\n")
            model_file.write("name: {0}\n".format(model_name))
            model_file.write('modeltype: {0}\n'.format(model))
            string = "shape: {0}\n".format(alpha)
            model_file.write(string)
            string = "pinv: {0}\n".format(pinv)
            model_file.write(string)
            string = 'basefreq: {0}  {1}  {2}  {3}\n'.format(basefreq_1, basefreq_2, basefreq_1, basefreq_2)
            model_file.write(string)
            string = 'rrates: {} {} {} {} {} {}\n'.format(GTR_A_C,GTR_A_G,GTR_A_T,GTR_C_G,GTR_C_T,GTR_G_T)
            model_file.write(string)
            model_file.write("end model\n")
        elif model in base_models_choice_aa :
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
                                              min_shape, max_shape, seqlen, base_model_list_aa,global_min_tstv,global_max_tstv,
                                              global_min_rrates_GTR,global_max_rrates_GTR, md):
    model_filename = tmpdir + "Modelfile_free" + runID + "_" + str(seed) + ".txt"
    seq_filename = tmpdir + "sim-sequences" + runID + "_" + str(seed) + ".txt"
    sim_filename = tmpdir + "sim_tree" + runID + "_" + str(seed) + ".txt"

    if args.sequence_type == 'DNA':
        model = md
        modelname = model + '_model'
    elif args.sequence_type == 'AA':
        model = md
        modelname = model + '_model'

    shape = random.uniform(min_shape, max_shape)
    pinv = random.uniform(min_pinv, max_pinv)

    r1 = random.uniform(min_brlen, max_brlen)
    r2 = random.uniform(min_brlen, max_brlen)
    r3 = random.uniform(min_internal_brlen, max_internal_brlen)
    r4 = random.uniform(min_brlen, max_brlen)
    r5 = random.uniform(min_brlen, max_brlen)
    r3 = r3 / 2
    r6 = r3
    tstv = random.uniform(global_min_tstv, global_max_tstv)
    basefreq_1 = np.random.normal(0.5,0.03,1)/2
    while basefreq_1 > 0.3 or basefreq_1 < 0.2:
        basefreq_1 = np.random.normal(0.5,0.03,1)/2
    basefreq_2 = (1 - 2*basefreq_1)/2
    basefreq_1 = float(basefreq_1)
    basefreq_2 = float(basefreq_2)
    GTR_A_C = random.uniform(global_min_rrates_GTR,global_max_rrates_GTR)
    GTR_A_G = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_A_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_C_G = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_C_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)
    GTR_G_T = random.uniform(global_min_rrates_GTR, global_max_rrates_GTR)

    print("Calling make_model_file: ", model_filename, " ", model, " ", shape, " ", pinv, " ", modelname)
    make_model_file(model_filename, model, shape, pinv, tstv,basefreq_1,basefreq_2,GTR_A_C,GTR_A_G,GTR_A_T,GTR_C_G,GTR_C_T,GTR_G_T, modelname)
    make_tree_file_model(sim_filename, topology, r1, r2, r3, r4, r5, r6,  seqlen, modelname)

    all_sim_filename = "all_sim_tree" + runID + ".txt"
    os.system("echo " + str(seed) + " >> " + all_sim_filename)
    os.system("cat " + sim_filename + " >> " + all_sim_filename)
    command = "PolyMoSim-v1.1.4 -s " + str(
        seed) + " -m " + model_filename + " -t " + sim_filename + " -f site_pattern_freq_relative_fill -n 1 1> " + seq_filename
    subprocess.run([command], shell=True)
    individual_frequencies = find_frequency(seq_filename)
    os.remove(seq_filename)
    os.remove(sim_filename)
    os.remove(model_filename)
    return individual_frequencies



####################################
# my main program starts here:
####################################

topology_list = []
evo_mod_list = []

frequency_simulations_filename = 'trainFrequency' + '_all_evo_mod' + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
        global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"

evo_mod_simulations_filename = "trainEvoModClass" + '_all_evo_mod'+ '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) +  "_pseed_" + str(python_random_seed) + ".npy"
### Frequency files:

print("Filenames training:")
print(frequency_simulations_filename)
print(evo_mod_simulations_filename)

sys.stdout.flush()

if os.path.isfile(frequency_simulations_filename):
    frequency_array = np.load(frequency_simulations_filename)
    evo_mod_array= np.load(evo_mod_simulations_filename)

else:
    if args.sequence_type == 'DNA':

        for j in base_models_choice_list:
            if j == 'JC':
                counter = 0
                x = 0
            elif j == 'K2P':
                counter = total_replicates
                x = 1
            elif j == 'F81':
                counter = 2*total_replicates
                x = 2
            elif j == 'HKY':
                counter = 3*total_replicates
                x = 3
            elif j == 'GTR':
                counter = 4*total_replicates
                x = 4
            for i in range(0, global_replicates):
                temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                              global_min_rrates_GTR,global_max_rrates_GTR, j)

                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)
            for i in range(global_replicates, 2 * global_replicates):
                temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,
                                                                    global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR, j)

                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)
            for i in range(2 * global_replicates, 3 * global_replicates):
                temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,\
                                                                    global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR, j)
                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)

            for i in range(3 * global_replicates, 3 * global_replicates + global_replicates_short_internalbranch):
                temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,
                                                                    global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR, j)

                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)
            for i in range(3 * global_replicates + global_replicates_short_internalbranch,
                   3 * global_replicates + 2 * global_replicates_short_internalbranch):
                 temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR, j)

                 frequency_array[i+counter] = temporary_array
                 evo_mod_list.append(x)
            for i in range(3 * global_replicates + 2 * global_replicates_short_internalbranch,
                   3 * global_replicates + 3 * global_replicates_short_internalbranch):
                 temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR, j)

                 frequency_array[i+counter] = temporary_array
                 evo_mod_list.append(x)

        print('Freq_shape: ', frequency_array.shape)
        evo_mod_array = np.array(evo_mod_list)
        evo_mod_array = np.reshape(evo_mod_array,(5*total_replicates,1))
        print('Evo_mod_shape', evo_mod_array.shape)

        np.save(frequency_simulations_filename, frequency_array)
        np.save(evo_mod_simulations_filename, evo_mod_array)


    elif args.sequence_type == 'AA':

        for j in base_models_choice_aa:
            if j == 'JTT':
                counter = 0
                x = 0
            elif j == 'LG':
                counter = total_replicates
                x = 1
            elif j == 'WAG':
                counter = 2 * total_replicates
                x = 2
            elif j == 'DAY':
                counter = 3 * total_replicates
                x = 3
            for i in range(0, global_replicates):

                temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR,j)

                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)
            for i in range(global_replicates, 2 * global_replicates):

                temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR,j)

                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)

            for i in range(2 * global_replicates, 3 * global_replicates):

                temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR,j)
                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)


            for i in range(3 * global_replicates, 3 * global_replicates + global_replicates_short_internalbranch):

                temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR,j)
                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)



            for i in range(3 * global_replicates + global_replicates_short_internalbranch,
                   3 * global_replicates + 2 * global_replicates_short_internalbranch):

                temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR,j)

                frequency_array[i+counter] = temporary_array
                evo_mod_list.append(x)
            for i in range(3 * global_replicates + 2 * global_replicates_short_internalbranch,
                   3 * global_replicates + 3 * global_replicates_short_internalbranch):

                 temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR,j)
                 frequency_array[i+counter] = temporary_array
                 evo_mod_list.append(x)

        print('Freq_shape: ', frequency_array.shape)
        evo_mod_array = np.array(evo_mod_list)
        evo_mod_array = np.reshape(evo_mod_array,(5*total_replicates,1))
        print('Evo_mod_shape', evo_mod_array.shape)

        np.save(frequency_simulations_filename, frequency_array)
        np.save(evo_mod_simulations_filename, evo_mod_array)


print('Shape of the frequency array: ', frequency_array.shape)
print(frequency_array)


temporary_JC = np.empty((prediction_replicates, 256))
temporary_K2P = np.empty((prediction_replicates, 256))
temporary_F81 = np.empty((prediction_replicates, 256))
temporary_HKY = np.empty((prediction_replicates, 256))
temporary_GTR = np.empty((prediction_replicates, 256))

temporary_JTT = np.empty((prediction_replicates, 160000))
temporary_LG = np.empty((prediction_replicates, 160000))
temporary_WAG = np.empty((prediction_replicates, 160000))
temporary_DAY = np.empty((prediction_replicates, 160000))


# time before simulate data
timeBeforeSimulateData = time.time()

python_random_seed += 1
random.seed(python_random_seed)

if predict_mode == 0:
    if args.sequence_type == 'DNA' :

        temporary_JC_filename = "predict_JC" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_K2P_filename = "predict_K2P" + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_F81_filename = "predict_F81"  + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_HKY_filename = "predict_HKY" + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_GTR_filename = "predict_GTR" + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"
    elif args.sequence_type == 'AA' :
        temporary_JTT_filename = "predict_JTT" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_LG_filename = "predict_LG" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_WAG_filename = "predict_WAG" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"


        temporary_DAY_filename = "predict_DAY" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

if predict_mode == 1:
  if args.sequence_type == 'DNA' :
    temporary_JC_filename = "predict_JC" + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

    temporary_K2P_filename = "predict_K2P" + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) +  '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) +"_pseed_" + str(python_random_seed) + ".npy"

    temporary_F81_filename = "predict_F81" + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

    temporary_HKY_filename = "predict_HKY" + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

    temporary_GTR_filename = "predict_GTR" + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) + "_pseed_" + str(python_random_seed) + ".npy"
  elif args.sequence_type == 'AA':
        temporary_JTT_filename = "predict_JTT" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_LG_filename = "predict_LG" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_WAG_filename = "predict_WAG" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary_DAY_filename = "predict_DAY" +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

print("Filenames for prediction files:")
if  args.sequence_type == 'DNA' :
    print("predict_JC: ", temporary_JC_filename)
    print("predict_K2P: ", temporary_K2P_filename)
    print("predict_F81: ", temporary_F81_filename)
    print("predict_HKY: ", temporary_HKY_filename)
    print("predict_GTR: ", temporary_GTR_filename)
elif args.sequence_type == 'AA' :
    print("predict_JTT: ", temporary_JTT_filename)
    print("predict_LG: ", temporary_LG_filename)
    print("predict_WAG: ", temporary_WAG_filename)
    print("predict_DAY: ", temporary_DAY_filename)
sys.stdout.flush()


if os.path.isfile(temporary_JC_filename):
    if args.sequence_type == 'DNA' :
        temporary_JC = np.load(temporary_JC_filename,allow_pickle=True)
        temporary_K2P = np.load(temporary_K2P_filename,allow_pickle=True)
        temporary_F81 = np.load(temporary_F81_filename,allow_pickle=True)
        temporary_HKY = np.load(temporary_HKY_filename, allow_pickle=True)
        temporary_GTR = np.load(temporary_GTR_filename, allow_pickle=True)

    elif args.sequence_type == 'AA' :
        temporary_JTT = np.load(temporary_JTT_filename,allow_pickle=True)
        temporary_LG = np.load(temporary_LG_filename,allow_pickle=True)
        temporary_WAG = np.load(temporary_WAG_filename,allow_pickle=True)
        temporary_DAY = np.load(temporary_DAY_filename, allow_pickle=True)
else:
     if args.sequence_type == 'DNA':
        if predict_mode == 0:
            base_models_choice_list = ['JC','K2P','F81','HKY','GTR']
            for i in base_models_choice_list:
                for j in range(0, prediction_replicates):
                    if i == 'JC':
                        temporary_JC[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR,i)
                    if i == 'K2P':
                        temporary_K2P[j] = simulate_random_lengths_random_free_model(0, 20000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)
                    if i == 'F81':
                        temporary_F81[j] = simulate_random_lengths_random_free_model(0, 30000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)

                    if i == 'HKY':
                        temporary_HKY[j] = simulate_random_lengths_random_free_model(0, 50000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)
                    if i == 'GTR':
                        temporary_GTR[j] = simulate_random_lengths_random_free_model(0, 60000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)


        if predict_mode == 1:

            for i in base_models_choice_list:
                for j in range(0, prediction_replicates):
                    if i == 'JC':
                        temporary_JC[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength, global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR,i)
                    if i == 'K2P':
                        temporary_K2P[j] = simulate_random_lengths_random_free_model(0, 20000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_short_internal_branchlength,
                                                                            global_max_short_internal_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)
                    if i == 'F81':
                        temporary_F81[j] = simulate_random_lengths_random_free_model(0, 30000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_short_internal_branchlength,
                                                                            global_max_short_internal_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)

                    if i == 'HKY':
                        temporary_HKY[j] = simulate_random_lengths_random_free_model(0, 50000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_short_internal_branchlength,
                                                                            global_max_short_internal_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)
                    if i == 'GTR':
                        temporary_GTR[j] = simulate_random_lengths_random_free_model(0, 60000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_short_internal_branchlength,
                                                                            global_max_short_internal_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)


        np.save(temporary_JC_filename, temporary_JC)
        np.save(temporary_K2P_filename, temporary_K2P)
        np.save(temporary_F81_filename, temporary_F81)
        np.save(temporary_HKY_filename, temporary_HKY)
        np.save(temporary_GTR_filename, temporary_GTR)

     if args.sequence_type == 'AA':
        if predict_mode == 0:
            for i in base_models_choice_aa:
                for j in range(0, prediction_replicates):
                    if i == 'JTT':
                        temporary_JTT[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR,i)
                    if i == 'LG':
                        temporary_LG[j] = simulate_random_lengths_random_free_model(0, 20000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)
                    if i == 'WAG':
                        temporary_WAG[j] = simulate_random_lengths_random_free_model(0, 30000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)

                    if i == 'DAY':
                        temporary_DAY[j] = simulate_random_lengths_random_free_model(0, 60000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)


        if predict_mode == 1:
            for i in base_models_choice_aa:
                for j in range(0, prediction_replicates):
                    if i == 'JTT':
                        temporary_JTT[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength, global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR,i)
                    if i == 'LG':
                        temporary_LG[j] = simulate_random_lengths_random_free_model(0, 20000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_short_internal_branchlength,
                                                                            global_max_short_internal_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)
                    if i == 'WAG':
                        temporary_WAG[j] = simulate_random_lengths_random_free_model(0, 30000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_short_internal_branchlength,
                                                                            global_max_short_internal_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)

                    if i == 'DAY':
                        temporary_DAY[j] = simulate_random_lengths_random_free_model(0, 60000000 + j, base_models_choice_list,
                                                                            global_min_branchlength,
                                                                            global_max_branchlength,
                                                                            global_min_short_internal_branchlength,
                                                                            global_max_short_internal_branchlength,
                                                                            global_min_pinv, global_max_pinv,
                                                                            global_min_shape, global_max_shape,
                                                                            global_predict_seqlen,
                                                                            base_models_choice_aa, global_min_tstv,
                                                                            global_max_tstv,
                                                                            global_min_rrates_GTR,
                                                                            global_max_rrates_GTR, i)

        np.save(temporary_JTT_filename, temporary_JTT)
        np.save(temporary_LG_filename, temporary_LG)
        np.save(temporary_WAG_filename, temporary_WAG)
        np.save(temporary_DAY_filename, temporary_DAY)

if  args.sequence_type == 'DNA':
    print('Shape of the temporary_JC: ', temporary_JC.shape)
    print(temporary_JC)
    print('Shape of the temporary_K2P: ', temporary_K2P.shape)
    print(temporary_K2P)
    print('Shape of the temporary_F81: ', temporary_F81.shape)
    print(temporary_F81)
    print('Shape of the temporary_HKY: ', temporary_HKY.shape)
    print(temporary_HKY)
    print('Shape of the temporary_GTR: ', temporary_GTR.shape)
    print(temporary_GTR)
elif args.sequence_type == 'AA':
    print('Shape of the temporary_JTT: ', temporary_JTT.shape)
    print(temporary_JTT)
    print('Shape of the temporary_LG: ', temporary_LG.shape)
    print(temporary_LG)
    print('Shape of the temporary_WAG: ', temporary_WAG.shape)
    print(temporary_WAG)
    print('Shape of the temporary_DAY: ', temporary_DAY.shape)
    print(temporary_DAY)
# ANN training starts here

sc = StandardScaler()
frequency_array = sc.fit_transform(frequency_array)

frequency_train, frequency_test, evo_mod_train, evo_mod_test = train_test_split(frequency_array, evo_mod_array,
                                                                                  test_size=0.04, random_state=42)

frequency_train = np.asarray(frequency_train)
evo_mod_train = np.asarray(evo_mod_train)

frequency_test = np.asarray(frequency_test)
evo_mod_test = np.asarray(evo_mod_test)

### User provided arguments:
EPOCHS = args.epochs

if args.sequence_type == 'DNA':
    LEN = 256
    evo_mod_train = tf.keras.utils.to_categorical(evo_mod_train, 5)
    evo_mod_test = tf.keras.utils.to_categorical(evo_mod_test, 5)
elif args.sequence_type == 'AA' :
    LEN = 160000
    evo_mod_train = tf.keras.utils.to_categorical(evo_mod_train, 4)
    evo_mod_test = tf.keras.utils.to_categorical(evo_mod_test, 4)

### Neural Networks

def branched_model_1(hparam):
    input = tf.keras.Input(shape=(LEN,))
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x)

    y = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x)
    y = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(y)
    y = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(y)
    if args.sequence_type == 'DNA':
        y = tf.keras.layers.Dense(units=5, activation='softmax')(y)
    elif args.sequence_type == 'AA':
        y = tf.keras.layers.Dense(units=4, activation='softmax')(y)
    ann = tf.keras.Model(inputs=[input], outputs=y)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_2(hparam):
    input = tf.keras.Input(shape=(LEN,))
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x)

    x_ann = tf.keras.Model(inputs=input, outputs=x)

    y = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    y = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(y)
    y = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(y)
    y = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(y)
    y = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(y)

    y_ann = tf.keras.Model(inputs=input, outputs=y)

    combined = tf.keras.layers.concatenate([x_ann.output, y_ann.output])

    z = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(combined)
    z = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(z)
    z = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(combined)
    z = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(z)
    z = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(z)
    if args.sequence_type == 'DNA':
        z = tf.keras.layers.Dense(units=5, activation='softmax')(z)
    elif args.sequence_type == 'AA':
        z = tf.keras.layers.Dense(units=4, activation='softmax')(z)
    ann = tf.keras.Model(inputs=[input], outputs=z)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_3(hparam):
    input = tf.keras.Input(shape=(LEN,))
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x)

    x_ann = tf.keras.Model(inputs=input, outputs=x)

    y = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    y = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(y)
    y = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(y)
    y = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(y)
    y = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(y)

    y_ann = tf.keras.Model(inputs=input, outputs=y)

    w = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    w = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(w)
    w = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(w)
    w = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(w)
    w = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(w)
    w_ann = tf.keras.Model(inputs=input, outputs=w)

    combined = tf.keras.layers.concatenate([x_ann.output, y_ann.output, w_ann.output])

    z = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(combined)
    z = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(z)
    z = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(z)
    z = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(z)
    z = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(z)
    if args.sequence_type == 'DNA':
        z = tf.keras.layers.Dense(units=5, activation='softmax')(z)
    elif args.sequence_type == 'AA':
        z = tf.keras.layers.Dense(units=4, activation='softmax')(z)
    ann = tf.keras.Model(inputs=[input], outputs=z)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_10(hparam):
    input = tf.keras.Input(shape=(LEN,))
    x1 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x1 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x1)
    x1 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x1)
    x1 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x1)
    x1 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x1)
    x1 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x1)
    x1 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x1)
    x1_ann = tf.keras.Model(inputs=input, outputs=x1)

    x2 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x2 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x2)
    x2 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x2)
    x2 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x2)
    x2 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x2)
    x2_ann = tf.keras.Model(inputs=input, outputs=x2)

    x3 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x3 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x3)
    x3 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x3)
    x3_ann = tf.keras.Model(inputs=input, outputs=x3)

    x4 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x4 = tf.keras.layers.Dropout(0.5)(x4)
    x4 = tf.keras.layers.Dense(units=16, activation=hparam[HP_ACTIVATION])(x4)
    x4 = tf.keras.layers.Dropout(0.5)(x4)
    x4 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x4)
    x4_ann = tf.keras.Model(inputs=input, outputs=x4)

    x5 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(input)
    x5 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x5)
    x5 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x5)
    x5_ann = tf.keras.Model(inputs=input, outputs=x5)

    x6 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(input)
    x6 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x6)
    x6 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x6)
    x6_ann = tf.keras.Model(inputs=input, outputs=x6)

    x7 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(input)
    x7 = tf.keras.layers.Dropout(0.5)(x7)
    x7 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x7)
    x7_ann = tf.keras.Model(inputs=input, outputs=x7)

    x8 = tf.keras.layers.Dense(units=111, activation=hparam[HP_ACTIVATION])(input)
    x8 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x8)
    x8 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x8)
    x8_ann = tf.keras.Model(inputs=input, outputs=x8)

    x9 = tf.keras.layers.Dense(units=77, activation="relu")(input)
    x9 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x9)
    x9 = tf.keras.layers.Dense(units=32, activation="relu")(x9)
    x9_ann = tf.keras.Model(inputs=input, outputs=x9)

    x10 = tf.keras.layers.Dense(units=155, activation="relu")(input)
    x10 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x10)
    x10 = tf.keras.layers.Dense(units=64, activation="relu")(x10)
    x10 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x10)
    x10 = tf.keras.layers.Dense(units=32, activation="relu")(x10)
    x10_ann = tf.keras.Model(inputs=input, outputs=x10)

    combined = tf.keras.layers.concatenate(
        [x1_ann.output, x2_ann.output, x3_ann.output, x4_ann.output, x5_ann.output, x6_ann.output, x7_ann.output,
         x8_ann.output, x9_ann.output, x10_ann.output])

    z = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(combined)
    z = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(z)
    z = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(z)
    z = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(z)
    z = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(z)
    if args.sequence_type == 'DNA':
        z = tf.keras.layers.Dense(units=5, activation='softmax')(z)
    elif args.sequence_type == 'AA':
        z = tf.keras.layers.Dense(units=4, activation='softmax')(z)
    ann = tf.keras.Model(inputs=[input], outputs=z)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_1u(hparam):
    input = tf.keras.Input(shape=(LEN,))
    x3 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(
        input)
    x4 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x3)
    x5 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x4)
    x6d = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x5)
    x7 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x6d)
    x7d = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x7)
    x8 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x7d)
    x9 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x8)
    x10 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x9)
    x11 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x10)
    x14 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x11)
    x15 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x14)
    x16 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x15)
    x17 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x16)
    x18 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x17)
    x19 = tf.keras.layers.Dense(units=16, activation=hparam[HP_ACTIVATION])(x18)
    if args.sequence_type == 'DNA':
        x20 = tf.keras.layers.Dense(units=5, activation='softmax')(x19)
    elif args.sequence_type == 'AA':
        x20 = tf.keras.layers.Dense(units=4, activation='softmax')(x19)
    ann = tf.keras.Model(inputs=[input], outputs=x20)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])

    return ann


def branched_model_cu(hparam):
    input = tf.keras.Input(shape=(LEN,))
    x3 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x4 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x3)
    x5 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x4)
    x6d = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x5)
    x7 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x6d)
    x7d = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x7)
    x8 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x7d)

    y0 = tf.keras.layers.concatenate([x7, x8])

    x9 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(y0)

    y1 = tf.keras.layers.concatenate([x5, x9])

    x10 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(y1)

    y2 = tf.keras.layers.concatenate([x10, x3])

    x11 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(y2)
    x14 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x11)
    x15 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x14)
    x16 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x15)
    x17 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x16)
    x19 = tf.keras.layers.Dense(units=16, activation=hparam[HP_ACTIVATION])(x17)
    if args.sequence_type == 'DNA':
        x20 = tf.keras.layers.Dense(units=5, activation='softmax')(x19)
    elif args.sequence_type == 'AA':
        x20 = tf.keras.layers.Dense(units=4, activation='softmax')(x19)
    ann = tf.keras.Model(inputs=[input], outputs=x20)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_cud(hparam):
    input = tf.keras.Input(shape=(LEN,))
    x3 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(input)
    x4 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x3)
    x5 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x4)
    x6d = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x5)
    x7 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x6d)
    x7d = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x7)
    x8 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x7d)

    y0 = tf.keras.layers.concatenate([x7, x8])

    x9 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(y0)

    y1 = tf.keras.layers.concatenate([x5, x9])

    x10 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(y1)

    y2 = tf.keras.layers.concatenate([x10, x3])

    x11 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(y2)
    x14 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x11)
    x15 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x14)
    x16 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x15)
    x17 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x16)
    x19 = tf.keras.layers.Dense(units=16, activation=hparam[HP_ACTIVATION])(x17)
    if args.sequence_type == 'DNA':
        x20 = tf.keras.layers.Dense(units=5, activation='softmax')(x19)
    elif args.sequence_type == 'AA':
        x20 = tf.keras.layers.Dense(units=4, activation='softmax')(x19)

    ann = tf.keras.Model(inputs=[input], outputs=x20)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def run_hparam_one_model_one_hset(model_func, X_train, Y_train, X_test, Y_test, hparams, do_probabilities=False):
    model = model_func(hparams)

    logdir = "model_b10"
    callback = tf.keras.callbacks.TensorBoard(logdir, profile_batch=0)
    hparams_callback = hp.KerasCallback(logdir, hparams)

    result = model.fit(x=X_train, y=Y_train, batch_size=256, epochs=EPOCHS, validation_data=(X_test, Y_test),
                       callbacks=[callback, hparams_callback], verbose=2)


    return model  ### result #, pred


def run_hparam_on_grid(model_func, X_train, Y_train, X_test, Y_test):
    with tf.summary.create_file_writer("tf_hparams_log.log").as_default():
        hp.hparams_config(hparams=HPARAMS, metrics=[hp.Metric(METRIC_ACCURACY, display_name='Accuracy')])

    session_num = 0

    for activation in HP_ACTIVATION.domain.values:
        for dropout in HP_DROPOUT.domain.values:
            for optimizer in HP_OPTIMIZER.domain.values:
                hparams = {HP_ACTIVATION: activation,
                           HP_DROPOUT: dropout,
                           HP_OPTIMIZER: optimizer, }
                run_name = "run-%d" % session_num
                print('--- Starting trial: %s' % run_name)
                print({h.name: hparams[h] for h in hparams})
                fitted_model = run_hparam_one_model_one_hset(model_func, X_train, Y_train, X_test, Y_test, hparams)
                print("Hparams:", hparams)
                if args.sequence_type == 'DNA':
                    print_highest_likelihood_evaluation(fitted_model, sc.transform(temporary_JC), sc.transform(temporary_K2P),
                                                    sc.transform(temporary_F81),
                                                    sc.transform(temporary_HKY), sc.transform(temporary_GTR))
                    fitted_model.save('DNAmodelpred' + '.' + args.neural_network)
                elif args.sequence_type == 'AA':
                        print_highest_likelihood_evaluation(fitted_model, sc.transform(temporary_JTT),
                                                            sc.transform(temporary_LG),
                                                            sc.transform(temporary_WAG), sc.transform(temporary_DAY))
                        fitted_model.save('AAmodelpred' + '.' + args.neural_network)
                session_num += 1


############
### Main:
############


if args.neural_network == 'b1':
    run_hparam_on_grid(branched_model_1,
                       frequency_train, evo_mod_train,
                       frequency_test, evo_mod_test
              )
if args.neural_network == 'b2':
    run_hparam_on_grid(branched_model_2,
                       frequency_train, evo_mod_train,
                       frequency_test, evo_mod_test
              )
if args.neural_network == 'b3':
    run_hparam_on_grid(branched_model_3,
                       frequency_train, evo_mod_train,
                       frequency_test, evo_mod_test
             )
if args.neural_network == 'b10':
    run_hparam_on_grid(branched_model_10,
                       frequency_train, evo_mod_train,
                       frequency_test, evo_mod_test
                       )
if args.neural_network == 'u':
    run_hparam_on_grid(branched_model_1u,
                       frequency_train, evo_mod_train,
                       frequency_test, evo_mod_test
             )
if args.neural_network == 'cu':
    run_hparam_on_grid(branched_model_cu,
                       frequency_train, evo_mod_train,
                       frequency_test, evo_mod_test
             )
if args.neural_network == 'cud':
    run_hparam_on_grid(branched_model_cud,
                       frequency_train, evo_mod_train,
                       frequency_test, evo_mod_test
        )

sys.exit(0)


