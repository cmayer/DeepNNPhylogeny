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
from keras.models import load_model


# Argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument('-m',type=str, required=True)
parser.add_argument('-na', type=str, required=True)
parser.add_argument('-bl', type=int, required=True)
parser.add_argument('-top', type=int, required=True)
args = parser.parse_args()



python_random_seed = 13

os.environ['MKL_NUM_THREADS'] = '10'
os.environ['GOTO_NUM_THREADS'] = '10'
os.environ['OMP_NUM_THREADS'] = '10'
os.environ['openmp'] = 'True'

tf.config.threading.set_inter_op_parallelism_threads(10)
tf.config.threading.set_intra_op_parallelism_threads(10)


dir   = "/home/nkulikov/PycharmProjects/test-ml-1/"
tmpdir = "/home/nkulikov/PycharmProjects/time_count/"
runID = str(os.getpid()*random.randint(1,100000))

## Switches:
predict_mode = args.bl ## 0: normal, 1: short internal branches,
random.seed(python_random_seed)

#global variables:
global_replicates = 1
global_replicates_short_internalbranch = 1
total_replicates = 3*global_replicates + 3*global_replicates_short_internalbranch
global_seqlen = 100000

prediction_replicates = 1
global_predict_seqlen = 100000

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
# K2P, F84, HKY, GTR
global_min_tstv = 1.0
global_max_tstv = 3.0
# F81, F84, HKY, GTR
global_min_basefreq = 0.2
global_max_basefreq = 0.3
# GTR
global_min_rrates_GTR = 0.1
global_max_rrates_GTR = 1.0


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
    seq_filename = tmpdir + "100000_sim-sequences" + '_' + args.m + '_' + 'bl' + '_' + str(args.bl) + '_' + 'top' + '_' + str(args.top) + ".fas"
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
    command = "PolyMoSim-v1.1.3b -s " + str(
        seed) + " -m " + model_filename + " -t " + sim_filename + " -f fasta -n 1 1> " + seq_filename
    subprocess.run([command], shell=True)
    os.remove(all_sim_filename)
    os.remove(sim_filename)
    os.remove(model_filename)
    print('Done!')



####################################
# my main program starts here:
####################################

if args.na == 'DNA':
        if predict_mode == 0:
            for j in range(0, prediction_replicates):
                if args.top == 0:
                    simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                elif args.top == 1:
                    simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                elif args.top == 2:
                    simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

        if predict_mode == 1:
            for j in range(0, prediction_replicates):
                if args.top == 0:
                    simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

                elif args.top == 1:
                    simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                elif args.top == 2:
                    simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)



if args.na == 'AA':
        if predict_mode == 0:
            for j in range(0, prediction_replicates):
                if args.top == 0:
                    simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

                elif args.top == 1:
                    simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                elif args.top == 2:
                    simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)

        if predict_mode == 1:

            for j in range(0, prediction_replicates):
                if args.top == 0:
                    simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                elif args.top == 1:

                    simulate_random_lengths_random_free_model(1, 20000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)
                elif args.top == 2:
                    simulate_random_lengths_random_free_model(2, 30000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength,
                                                                       global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR)


