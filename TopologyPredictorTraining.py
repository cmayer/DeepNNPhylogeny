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
import sys
import argparse


from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, train_test_split
from tensorboard.plugins.hparams import api as hp



# Argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument('-sequence_type', type=str, choices=['DNA','AA'],
                    help="Compulsory argument. Nucleotide or Amino acids data. Enter DNA for nucleotide sequences. Enter AA for amino acid sequences.", required=True)
parser.add_argument('-substitution_model',type=str, choices=['JC','K2P','F81','F84','HKY','GTR','JTT','LG','WAG_OLD','WAG','WAG_STAR','DAY'],
                    help="Compulsory argument. JC,K2P,F81,F84,HKY,GTR are nucleotide models. JTT,LG,WAG_OLD,WAG,WAG_STAR,DAY are amino acid models.", required=True)
parser.add_argument('-seqlen', type=int, default=None, help = "Optional argument. The length of the training sequences. Default for DNA: 100000\tDefault for amino acids: 1000000", required=False)
parser.add_argument('-seqnum', type=int, default=None,
                    help = "Optional argument. The number of sequences used for training. Actual number is multiplied by 6. Default for DNA: 100000\tDefault for amino acids: 1000", required=False)
parser.add_argument('-neural_network', type=str, default='b1', choices=['b1','b2','b3','b10','u','cu','cud'],
                    help="Optional argument. Which neural network should be trained. Default: b1", required=False)
parser.add_argument('-epochs', type=int, default=500,
                    help="Optional parameter. The number of epochs is a hyperparameter that defines the number times that the learning algorithm will work through the entire training dataset. Default: 500", required=False)
parser.add_argument('-predict_mode', type=int, default=0, choices=[0,1],
                    help="Optional argument. Type of predicted sequences either normal internal branches(0) or short internal branches(1). Default: 0", required=False)
parser.add_argument('-predict_seqlen', type=int, default=1000,help="Optional argument. The length of the prediction sequences. Default: 1000", required=False)
parser.add_argument('-predict_seqnum', type=int, default=1000,help="Optional argument. The number of the prediction sequences. Default: 1000", required=False)

args = parser.parse_args()

python_random_seed = 13

os.environ['MKL_NUM_THREADS'] = '56'
os.environ['GOTO_NUM_THREADS'] = '56'
os.environ['OMP_NUM_THREADS'] = '56'
os.environ['openmp'] = 'True'

tf.config.threading.set_inter_op_parallelism_threads(1)
tf.config.threading.set_intra_op_parallelism_threads(1)

tmpdir = os.getcwd() + "/"
runID = str(os.getpid()*random.randint(1,100000))

## Switches:
predict_mode = args.predict_mode ## 0: normal, 1: short internal branches,
random.seed(python_random_seed)

#global variables:
if (args.sequence_type == 'DNA' and args.seqlen == None and args.seqnum == None):
    global_seqlen = 100000
    global_replicates = 100000
    global_replicates_short_internalbranch = 100000
elif (args.sequence_type == 'AA' and args.seqlen == None and args.seqnum == None):
    global_seqlen = 1000000
    global_replicates = 1000
    global_replicates_short_internalbranch = 1000
elif (args.sequence_type == 'DNA' and args.seqlen == None):
    global_seqlen = 100000
    global_replicates = args.seqnum
    global_replicates_short_internalbranch = args.seqnum
elif (args.sequence_type == 'AA' and args.seqlen == None):
    global_seqlen = 1000000
    global_replicates = args.seqnum
    global_replicates_short_internalbranch = args.seqnum
elif (args.sequence_type == 'DNA' and args.seqnum == None):
    global_seqlen = args.seqlen
    global_replicates = 100000
    global_replicates_short_internalbranch = 100000
elif (args.sequence_type == 'AA' and args.seqnum == None):
    global_seqlen = args.seqlen
    global_replicates = 1000
    global_replicates_short_internalbranch = 1000
else:
    global_seqlen = args.seqlen
    global_replicates = args.seqnum
    global_replicates_short_internalbranch = args.seqnum

prediction_replicates = args.predict_seqnum
global_predict_seqlen = args.predict_seqlen

total_replicates = 3*global_replicates + 3*global_replicates_short_internalbranch

if args.sequence_type == 'DNA' :
    frequency_array = np.empty((total_replicates, 256))
    evo_mod_array = np.empty((total_replicates, 1))
if args.sequence_type == 'AA' :
    frequency_array = np.empty((total_replicates, 160000))
    evo_mod_array = np.empty((total_replicates, 1))


global_min_branchlength = 0.1
global_max_branchlength = 0.5
global_min_short_internal_branchlength = 0.001
global_max_short_internal_branchlength = 0.02


base_models_choice_list = ['JC','K2P','F81','F84','HKY','GTR']
base_models_choice_aa = ['JTT', 'LG', 'WAG_OLD', 'WAG', 'WAG_STAR', 'DAY']
# Parameters for nucleotide substitution models
# All models
global_min_pinv   = 0
global_max_pinv   = 0.50
global_min_shape  = 0.01
global_max_shape  = 4
# K2P, F84, HKY
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

##########

### Evaluation routines: 
def highest_likelihood (predictions):
    topo = []
    for prediction in predictions:
        if prediction[0] > prediction[1]:
            if prediction[0] > prediction[2]:
                topo.append(0)
            else:  # predict[2] >= prediction[0] > prediction[1]:
                topo.append(2)
        else:      # predict[1] >= prediction[0]
            if prediction[1] > prediction[2]:
                topo.append(1)
            else:  # prediction[2] >= predict[1] >= prediction[0]
                topo.append(2)
    return topo

def which_topology (predictions,threshold):
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


def count_tolologies_by_likelihood(predictions):
    topo = highest_likelihood (predictions)
    return [topo.count(0)/len(topo), topo.count(1)/len(topo), topo.count(2)/len(topo)]


def print_highest_likelihood_evaluation(model, d0, d1, d2):
    pred0 = model.predict(d0)
    pred1 = model.predict(d1)
    pred2 = model.predict(d2)
    print(count_tolologies_by_likelihood(pred0))
    print(count_tolologies_by_likelihood(pred1))
    print(count_tolologies_by_likelihood(pred2))



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
                                              global_min_rrates_GTR,global_max_rrates_GTR):
    model_filename = tmpdir + "Modelfile_free" + runID + "_" + str(seed) + ".txt"
    seq_filename = tmpdir + "sim-sequences" + runID + "_" + str(seed) + ".txt"
    sim_filename = tmpdir + "sim_tree" + runID + "_" + str(seed) + ".txt"


    model = args.substitution_model
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

### Frequency files:
if (args.substitution_model == 'JC') :

    frequency_simulations_filename = 'trainFrequency' + args.substitution_model +'_R_' + str(global_replicates) + "_B_" + str(
    global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
    global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
    global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(global_max_pinv) + "_a_" + str(
    global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.substitution_model +'_R_' + str(global_replicates) + "_B_" + str(
    global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
    global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
    global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(global_max_pinv) + "_a_" + str(
    global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.substitution_model in base_models_choice_aa :

    frequency_simulations_filename = 'trainFrequency' + args.substitution_model +'_R_' + str(global_replicates) + "_B_" + str(
    global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
    global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
    global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(global_max_pinv) + "_a_" + str(
    global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.substitution_model +'_R_' + str(global_replicates) + "_B_" + str(
    global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
    global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
    global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(global_max_pinv) + "_a_" + str(
    global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.substitution_model == 'K2P' :

    frequency_simulations_filename = 'trainFrequency' + args.substitution_model + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
        global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.substitution_model + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) +  "_pseed_" + str(python_random_seed) + ".npy"
elif args.substitution_model == 'F81' :

    frequency_simulations_filename = 'trainFrequency' + args.substitution_model + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.substitution_model + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.substitution_model == 'F84' or args.substitution_model == 'HKY' :

    frequency_simulations_filename = 'trainFrequency' + args.substitution_model + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.substitution_model + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"
elif args.substitution_model == 'GTR' :

    frequency_simulations_filename = 'trainFrequency' + args.substitution_model + '_R_' + str(global_replicates) + "_B_" + str(
        global_min_branchlength) + "_" + str(global_max_branchlength) + "_r_" + str(
        global_replicates_short_internalbranch) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(
        global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(
        global_min_tstv) + '_' + str(
        global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
        global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) + "_pseed_" + str(python_random_seed) + ".npy"

    topology_simulations_filename = "trainTopoClass" + args.substitution_model + '_R_' + str(global_replicates) + "_B_" + str(
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

if os.path.isfile(frequency_simulations_filename):
    frequency_array = np.load(frequency_simulations_filename)
    topology_array = np.load(topology_simulations_filename)

else:
    if args.sequence_type == 'DNA':

        for i in range(0, global_replicates):
            temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                              global_min_rrates_GTR,global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(0)
        for i in range(global_replicates, 2 * global_replicates):
            temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,
                                                                    global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(1)
        for i in range(2 * global_replicates, 3 * global_replicates):
            temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,\
                                                                    global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)
            frequency_array[i] = temporary_array
            topologies = topology_list.append(2)

        for i in range(3 * global_replicates, 3 * global_replicates + global_replicates_short_internalbranch):
            temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,
                                                                    global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(0)
        for i in range(3 * global_replicates + global_replicates_short_internalbranch,
                   3 * global_replicates + 2 * global_replicates_short_internalbranch):
             temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

             frequency_array[i] = temporary_array
             topologies = topology_list.append(1)
        for i in range(3 * global_replicates + 2 * global_replicates_short_internalbranch,
                   3 * global_replicates + 3 * global_replicates_short_internalbranch):
             temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

             frequency_array[i] = temporary_array
             topologies = topology_list.append(2)

        topology_array = topology_list

        np.save(frequency_simulations_filename, frequency_array)
        np.save(topology_simulations_filename, topology_array)



    elif args.sequence_type == 'AA':

        for i in range(0, global_replicates):
            temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(0)
        for i in range(global_replicates, 2 * global_replicates):
            temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(1)
        for i in range(2 * global_replicates, 3 * global_replicates):
            temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(2)
        for i in range(3 * global_replicates, 3 * global_replicates + global_replicates_short_internalbranch):
            temporary_array = simulate_random_lengths_random_free_model(0, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(0)
            
        for i in range(3 * global_replicates + global_replicates_short_internalbranch,
                   3 * global_replicates + 2 * global_replicates_short_internalbranch):
            temporary_array = simulate_random_lengths_random_free_model(1, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

            frequency_array[i] = temporary_array
            topologies = topology_list.append(1)
        for i in range(3 * global_replicates + 2 * global_replicates_short_internalbranch,
                   3 * global_replicates + 3 * global_replicates_short_internalbranch):
             temporary_array = simulate_random_lengths_random_free_model(2, i, base_models_choice_list,
                                                                    global_min_branchlength, global_max_branchlength,
                                                                    global_min_short_internal_branchlength,
                                                                    global_max_short_internal_branchlength,
                                                                    global_min_pinv, global_max_pinv, global_min_shape,
                                                                    global_max_shape, global_seqlen,
                                                                    base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                    global_min_rrates_GTR,global_max_rrates_GTR)

             frequency_array[i] = temporary_array
             topologies = topology_list.append(2)
        topology_array = topology_list
        np.save(frequency_simulations_filename, frequency_array)
        np.save(topology_simulations_filename, topology_array)

print('Shape of the frequency array: ', frequency_array.shape)
print(frequency_array)

print('Shape of the topology array: ', frequency_array.shape)
print(topology_array)

if args.sequence_type == 'DNA':
    temporary_0 = np.empty((prediction_replicates, 256))
    temporary_1 = np.empty((prediction_replicates, 256))
    temporary_2 = np.empty((prediction_replicates, 256))
elif args.sequence_type == 'AA':
    temporary_0 = np.empty((prediction_replicates, 160000))
    temporary_1 = np.empty((prediction_replicates, 160000))
    temporary_2 = np.empty((prediction_replicates, 160000))



python_random_seed += 1
random.seed(python_random_seed)

if predict_mode == 0:
    if (args.substitution_model == 'JC') or (args.substitution_model in base_models_choice_aa) :
        temporary0_filename = "predict0_" + args.substitution_model +"_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
        global_min_pinv) + "_" + str(global_max_pinv) + "_pseed_" + str(python_random_seed) + ".npy"
    elif args.substitution_model == 'K2P':
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) +"_pseed_" + str(python_random_seed) + ".npy"
    elif args.substitution_model == 'F81':
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) +"_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) +"_pseed_" + str(python_random_seed) + ".npy"
    elif args.substitution_model == 'F84' or args.substitution_model == 'HKY' :
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) +"_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) +"_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) +"_pseed_" + str(python_random_seed) + ".npy"
    elif args.substitution_model == 'GTR':
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) +"_pseed_" + str(python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) +'_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) +"_pseed_" + str(python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_i_" + str(
            global_min_pinv) + "_" + str(global_max_pinv) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) +"_pseed_" + str(python_random_seed) + ".npy"
if predict_mode == 1:
    if (args.substitution_model == 'JC') or (args.substitution_model in base_models_choice_aa) :
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(
        python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(
        python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
        prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
        global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
        global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
        global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + "_pseed_" + str(
        python_random_seed) + ".npy"
    elif args.substitution_model == 'K2P':
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) +  '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + "_pseed_" + str(
            python_random_seed) + ".npy"
    elif args.substitution_model == 'F81' :
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"
    elif args.substitution_model == 'F84' or args.substitution_model == 'HKY':
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + "_pseed_" + str(
            python_random_seed) + ".npy"
    elif args.substitution_model == 'GTR':
        temporary0_filename = "predict0_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary1_filename = "predict1_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) + "_pseed_" + str(
            python_random_seed) + ".npy"

        temporary2_filename = "predict2_" + args.substitution_model + "_L_" + str(global_predict_seqlen) + "_r_" + str(
            prediction_replicates) + "_B_" + str(global_min_branchlength) + "_" + str(
            global_max_branchlength) + "_b_" + str(global_min_short_internal_branchlength) + "_" + str(
            global_max_short_internal_branchlength) + "_i_" + str(global_min_pinv) + "_" + str(
            global_max_pinv) + "_a_" + str(global_min_shape) + "_" + str(global_max_shape) + '_tstv_' + str(global_min_tstv) + '_' + str(
            global_max_tstv) + '_bf_' + str(global_min_basefreq) + '_' + str(
            global_max_basefreq) + '_rt_' + str(global_min_rrates_GTR) + '_' + str(global_max_rrates_GTR) + "_pseed_" + str(
            python_random_seed) + ".npy"

print("Filenames for prediction files:")
print("predict0: ", temporary0_filename)
print("predict1: ", temporary1_filename)
print("predict2: ", temporary2_filename)
sys.stdout.flush()


if os.path.isfile(temporary0_filename):
    temporary_0 = np.load(temporary0_filename,allow_pickle=True)
    temporary_1 = np.load(temporary1_filename,allow_pickle=True)
    temporary_2 = np.load(temporary2_filename,allow_pickle=True)
else:
     if args.sequence_type == 'DNA':
        temporary_0 = np.empty((prediction_replicates, 256))
        temporary_1 = np.empty((prediction_replicates, 256))
        temporary_2 = np.empty((prediction_replicates, 256))
        if predict_mode == 0:
            for j in range(0, prediction_replicates):
                temporary_0[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                temporary_1[j] = simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                temporary_2[j] = simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

        if predict_mode == 1:
            for j in range(0, prediction_replicates):
                temporary_0[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

                temporary_1[j] = simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                temporary_2[j] = simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

        np.save(temporary0_filename, temporary_0)
        np.save(temporary1_filename, temporary_1)
        np.save(temporary2_filename, temporary_2)

     if args.sequence_type == 'AA':
        temporary_0 = np.empty((prediction_replicates, 160000))
        temporary_1 = np.empty((prediction_replicates, 160000))
        temporary_2 = np.empty((prediction_replicates, 160000))
        if predict_mode == 0:
            for j in range(0, prediction_replicates):
                temporary_0[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

                temporary_1[j] = simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                temporary_2[j] = simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

        if predict_mode == 1:

            for j in range(0, prediction_replicates):
                temporary_0[j] = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

                temporary_1[j] = simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                temporary_2[j] = simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)


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

frequency_train, frequency_test, topology_train, topology_test = train_test_split(frequency_array, topology_array,
                                                                                  test_size=0.04, random_state=42)

frequency_train = np.asarray(frequency_train)
topology_train = np.asarray(topology_train)

frequency_test = np.asarray(frequency_test)
topology_test = np.asarray(topology_test)

topology_train = tf.keras.utils.to_categorical(topology_train, 3)
topology_test = tf.keras.utils.to_categorical(topology_test, 3)

EPOCHS = args.epochs
### Put ak code here if used:
if args.sequence_type == 'DNA':
    LEN = 256
elif args.sequence_type == 'AA' :
    LEN = 160000



def branched_model_1(hparam, normlayer):
    input = tf.keras.Input(shape=(LEN,))
    inputn = normlayer(input)

    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x)

    y = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x)
    y = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(y)
    y = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(y)
    y = tf.keras.layers.Dense(units=3, activation='softmax')(y)
    ann = tf.keras.Model(inputs=[input], outputs=y)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_2(hparam, normlayer):
    input = tf.keras.Input(shape=(LEN,))
    inputn = normlayer(input)

    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x)

    x_ann = tf.keras.Model(inputs=input, outputs=x)

    y = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
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
    z = tf.keras.layers.Dense(units=3, activation='softmax')(z)
    ann = tf.keras.Model(inputs=[input], outputs=z)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_3(hparam, normlayer):
    input = tf.keras.Input(shape=(LEN,))
    inputn = normlayer(input)

    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(x)
    x = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x)
    x = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x)

    x_ann = tf.keras.Model(inputs=input, outputs=x)

    y = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
    y = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(y)
    y = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(y)
    y = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(y)
    y = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(y)

    y_ann = tf.keras.Model(inputs=input, outputs=y)

    w = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
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
    z = tf.keras.layers.Dense(units=3, activation='softmax')(z)
    ann = tf.keras.Model(inputs=[input], outputs=z)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_10(hparam, normlayer):
    input = tf.keras.Input(shape=(LEN,))
    inputn = normlayer(input)

    x1 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
    x1 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x1)
    x1 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x1)
    x1 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x1)
    x1 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(x1)
    x1 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x1)
    x1 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x1)
    x1_ann = tf.keras.Model(inputs=input, outputs=x1)

    x2 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
    x2 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x2)
    x2 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(x2)
    x2 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x2)
    x2 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x2)
    x2_ann = tf.keras.Model(inputs=input, outputs=x2)

    x3 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
    x3 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x3)
    x3 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x3)
    x3_ann = tf.keras.Model(inputs=input, outputs=x3)

    x4 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
    x4 = tf.keras.layers.Dropout(0.5)(x4)
    x4 = tf.keras.layers.Dense(units=16, activation=hparam[HP_ACTIVATION])(x4)
    x4 = tf.keras.layers.Dropout(0.5)(x4)
    x4 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x4)
    x4_ann = tf.keras.Model(inputs=input, outputs=x4)

    x5 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(inputn)
    x5 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x5)
    x5 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x5)
    x5_ann = tf.keras.Model(inputs=input, outputs=x5)

    x6 = tf.keras.layers.Dense(units=64, activation=hparam[HP_ACTIVATION])(inputn)
    x6 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x6)
    x6 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x6)
    x6_ann = tf.keras.Model(inputs=input, outputs=x6)

    x7 = tf.keras.layers.Dense(units=128, activation=hparam[HP_ACTIVATION])(inputn)
    x7 = tf.keras.layers.Dropout(0.5)(x7)
    x7 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x7)
    x7_ann = tf.keras.Model(inputs=input, outputs=x7)

    x8 = tf.keras.layers.Dense(units=111, activation=hparam[HP_ACTIVATION])(inputn)
    x8 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x8)
    x8 = tf.keras.layers.Dense(units=32, activation=hparam[HP_ACTIVATION])(x8)
    x8_ann = tf.keras.Model(inputs=input, outputs=x8)

    x9 = tf.keras.layers.Dense(units=77, activation="relu")(inputn)
    x9 = tf.keras.layers.Dropout(hparam[HP_DROPOUT])(x9)
    x9 = tf.keras.layers.Dense(units=32, activation="relu")(x9)
    x9_ann = tf.keras.Model(inputs=input, outputs=x9)

    x10 = tf.keras.layers.Dense(units=155, activation="relu")(inputn)
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
    z = tf.keras.layers.Dense(units=3, activation='softmax')(z)
    ann = tf.keras.Model(inputs=[input], outputs=z)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_1u(hparam, normlayer):
    input = tf.keras.Input(shape=(LEN,))
    inputn = normlayer(input)

    x3 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(
        inputn) 
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
    x20 = tf.keras.layers.Dense(units=3, activation="softmax")(x19)

    ann = tf.keras.Model(inputs=[input], outputs=x20)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])

    return ann


def branched_model_cu(hparam, normlayer):
    input = tf.keras.Input(shape=(LEN,))
    inputn = normlayer(input)

    x3 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
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
    x20 = tf.keras.layers.Dense(units=3, activation="softmax")(x19)

    ann = tf.keras.Model(inputs=[input], outputs=x20)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def branched_model_cud(hparam, normlayer):
    input = tf.keras.Input(shape=(LEN,))
    inputn = normlayer(input)

    x3 = tf.keras.layers.Dense(units=256, activation=hparam[HP_ACTIVATION])(inputn)
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
    x20 = tf.keras.layers.Dense(units=3, activation="softmax")(x19)

    ann = tf.keras.Model(inputs=[input], outputs=x20)
    optimizer = get_optimizer(hparam[HP_OPTIMIZER])
    ann.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
    return ann


def run_hparam_one_model_one_hset(model_func, X_train, Y_train, X_test, Y_test, hparams, normlayer, do_probabilities=False):
    model = model_func(hparams, normlayer)

    logdir = 'TopPred' + args.substitution_model + '_' + args.neural_network + '_' + str(args.epochs) + '.log' 
    callback = tf.keras.callbacks.TensorBoard(logdir, profile_batch=0)
    hparams_callback = hp.KerasCallback(logdir, hparams)

    result = model.fit(x=X_train, y=Y_train, batch_size=256, epochs=EPOCHS, validation_data=(X_test, Y_test),
                       callbacks=[callback, hparams_callback], verbose=2)

    return model  ### result 


def run_hparam_on_grid(model_func, X_train, Y_train, X_test, Y_test, normlayer):

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
                fitted_model = run_hparam_one_model_one_hset(model_func, X_train, Y_train, X_test, Y_test, hparams, normlayer)
                print("Hparams:", hparams)
                print_highest_likelihood_evaluation(fitted_model, temporary_0, temporary_1,
                                                    temporary_2)

                session_num += 1
                if args.substitution_model == 'DAY':
                    fitted_model.save('TopPred' + 'Dayhoff' + '_' + str(args.epochs) + '.' + args.neural_network)
                    print('The name of the saved NN: TopPredDayhoff_' + str(args.epochs) + '.' + args.neural_network)
                else:
                    fitted_model.save('TopPred' + args.substitution_model  + '_' + str(args.epochs) + '.' + args.neural_network)
                    print('The name of the saved NN: TopPred' + args.substitution_model  + '_' + str(args.epochs) + '.' + args.neural_network)

############
### Main:
############

normlayer = tf.keras.layers.experimental.preprocessing.Normalization()
normlayer.adapt(frequency_train)

if args.neural_network == 'b1':
    run_hparam_on_grid(branched_model_1,
                       frequency_train, topology_train,
                       frequency_test, topology_test, normlayer
            )
elif args.neural_network == 'b2':
    run_hparam_on_grid(branched_model_2,
                       frequency_train, topology_train,
                       frequency_test, topology_test, normlayer
            )
elif args.neural_network == 'b3':
    run_hparam_on_grid(branched_model_3,
                       frequency_train, topology_train,
                       frequency_test, topology_test, normlayer
            )
elif args.neural_network == 'b10':
    run_hparam_on_grid(branched_model_10,
                       frequency_train, topology_train,
                       frequency_test, topology_test, normlayer
                       )
elif args.neural_network == 'u':
    run_hparam_on_grid(branched_model_1u,
                       frequency_train, topology_train,
                       frequency_test, topology_test, normlayer
            )
elif args.neural_network == 'cu':
    run_hparam_on_grid(branched_model_cu,
                       frequency_train, topology_train,
                       frequency_test, topology_test, normlayer
            )
elif args.neural_network == 'cud':
    run_hparam_on_grid(branched_model_cud,
                       frequency_train, topology_train,
                       frequency_test, topology_test, normlayer
            )

sys.exit(0)
