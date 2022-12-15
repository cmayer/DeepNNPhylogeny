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
import re

#import keras_tuner as kt

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, train_test_split
from tensorboard.plugins.hparams import api as hp



# Argparse arguments

parser = argparse.ArgumentParser()
#parser.add_argument('-m',type=str, required=True)
parser.add_argument('-na', type=str, required=True)
parser.add_argument('-bl', type=int, required=True)
parser.add_argument('--h')
args = parser.parse_args()

if args.h :
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

tf.config.threading.set_inter_op_parallelism_threads(10)
tf.config.threading.set_intra_op_parallelism_threads(10)

#tf.compat.v1.Session(config
#config = tf.compat.v1.ConfigProto(device_count={"CPU": 120},
#                                  inter_op_parallelism_threads=120,
#                                  intra_op_parallelism_threads=120)
#sess = tf.compat.v1.Session(config=config)


dir   = "/home/nkulikov/PycharmProjects/test-ml-1/"
tmpdir = "/home/nkulikov/PycharmProjects/ML_NJ_iqtree/evo_mod/"
runID = str(os.getpid()*random.randint(1,100000))

## Switches:
predict_mode = args.bl ## 0: normal, 1: short internal branches,
random.seed(python_random_seed)

#global variables:
global_replicates = 1000
global_replicates_short_internalbranch = 1000
total_replicates = 3*global_replicates + 3*global_replicates_short_internalbranch
global_seqlen = 100000

prediction_replicates = 1000
global_predict_seqlen = 30000

global_min_branchlength = 0.1
global_max_branchlength = 0.5
global_min_short_internal_branchlength = 0.001
global_max_short_internal_branchlength = 0.02


base_models_choice_list = ['JC','K2P','F81','HKY','GTR']
base_models_choice_aa = ['JTT', 'LG', 'WAG', 'WAG_STAR', 'DAY']
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

iqtree_predict_mode = 1

print("Parameters used:")
if iqtree_predict_mode == 1:
    print("Testing reconstruction success with ML, MP, NJ, not with ANNs.")



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


### Evaluation routines: (moved here in Test2):


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



def iqtree(sim_file):
    command =  'iqtree2.1.3 -mset JC,K2P,F81,HKY,GTR -mrate I,G,I+G -s ' + sim_file # + '>' +  sim_file + '.treefile'
    subprocess.run([command], shell=True)
    f = open(sim_file + '.iqtree', "r")
    letters = ''
    flag = 0
    if args.na == 'DNA':
      for i in f:
        if i.startswith('Best-fit'):
            print('We found')
            for j in i:
                if j == ':':
                    flag = 1
                if flag == 1 :
                    if j != ' ' and j != ':' and j != '\t':
                        letters +=j
                        if len(letters) >= 3:
                            break

    print('letters:', letters)
    if letters[0:2] == 'JC':
        print('my JC:', letters[0:2])
        topolog = 'JC'

    elif letters[0:3] == 'K2P':
        topolog = 'K2P'

    elif letters[0:3] == 'F81':
        topolog = 'F81'

    elif letters[0:3] == 'HKY':
        topolog = 'HKY'

    elif letters[0:3] == 'GTR':
        print('my GTR:', letters[0:3])
        topolog = 'GTR'


    else :
        topolog = '?'




    os.remove(sim_file + '.treefile')
    os.remove(sim_file + '.mldist')
    os.remove(sim_file + '.log')
    os.remove(sim_file + '.iqtree')
    os.remove(sim_file + '.model.gz')
    os.remove(sim_file + '.ckp.gz')
    os.remove(sim_file + '.bionj')

    return topolog


def count_correct_trees (topo) :
    if args.na == 'DNA' :
        return [topo.count('JC')/len(topo), topo.count('K2P')/len(topo), topo.count('F81')/len(topo), topo.count('HKY')/len(topo), topo.count('GTR')/len(topo)]
    if args.na == 'AA' :
        return [topo.count('JTT')/len(topo), topo.count('LG')/len(topo), topo.count('WAG')/len(topo), topo.count('WAG_STAR')/len(topo), topo.count('DAY')/len(topo)]


def simulate_random_lengths_random_free_model(topology, seed, base_model_list, min_brlen, max_brlen, min_internal_brlen,
                                              max_internal_brlen, min_pinv, max_pinv,
                                              min_shape, max_shape, seqlen, base_model_list_aa,global_min_tstv,global_max_tstv,
                                              global_min_rrates_GTR,global_max_rrates_GTR, md):
    model_filename = tmpdir + "Modelfile_free" + runID + "_" + str(seed) + ".txt"
    seq_filename = tmpdir + "sim-sequences" + runID + "_" + str(seed) + ".fa"
    sim_filename = tmpdir + "sim_tree" + runID + "_" + str(seed) + ".txt"

    if args.na == 'DNA':
        model = md
        modelname = model + '_model'
    elif args.na == 'AA':
        model = md
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

    command =  'PolyMoSim-v1.1.3b -s ' + str(seed) + ' -m ' + model_filename + ' -t ' + sim_filename + ' -n 1 -f fasta 1> ' + seq_filename
    subprocess.run([command], shell=True)
    iqtree_reconstruction = iqtree(seq_filename)
    os.remove(seq_filename)
    os.remove(sim_filename)
    os.remove(model_filename)

    return iqtree_reconstruction



####################################
# my main program starts here:
####################################


# time before simulate data
timeBeforeSimulateData = time.time()

python_random_seed += 1
random.seed(python_random_seed)

if args.na == 'DNA':
    temporary_JC = []
    temporary_K2P = []
    temporary_F81 = []
    temporary_HKY = []
    temporary_GTR = []
    if predict_mode == 0:
            base_models_choice_list = ['JC','K2P','F81','HKY','GTR']
            for i in base_models_choice_list:
                for j in range(0, prediction_replicates):
                    if i == 'JC':
                        temporary = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR,i)
                        temporary_JC.append(temporary)
                    if i == 'K2P':
                        temporary = simulate_random_lengths_random_free_model(0, 20000000 + j, base_models_choice_list,
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
                        temporary_K2P.append(temporary)
                    if i == 'F81':
                        temporary = simulate_random_lengths_random_free_model(0, 30000000 + j, base_models_choice_list,
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
                        temporary_F81.append(temporary)

                    if i == 'HKY':
                        temporary = simulate_random_lengths_random_free_model(0, 50000000 + j, base_models_choice_list,
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
                        temporary_HKY.append(temporary)
                    if i == 'GTR':
                        temporary = simulate_random_lengths_random_free_model(0, 60000000 + j, base_models_choice_list,
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
                        temporary_GTR.append(temporary)



    if predict_mode == 1:
            base_models_choice_list = ['JC','K2P','F81','HKY','GTR']
            for i in base_models_choice_list:
                for j in range(0, prediction_replicates):
                    if i == 'JC':
                        temporary = simulate_random_lengths_random_free_model(0, 10000000 + j, base_models_choice_list,
                                                                       global_min_branchlength, global_max_branchlength,
                                                                       global_min_short_internal_branchlength, global_max_short_internal_branchlength,
                                                                       global_min_pinv, global_max_pinv,
                                                                       global_min_shape, global_max_shape,
                                                                       global_predict_seqlen,base_models_choice_aa,global_min_tstv,global_max_tstv,
                                                                       global_min_rrates_GTR,global_max_rrates_GTR,i)
                        temporary_JC.append(temporary)

                    if i == 'K2P':
                        temporary = simulate_random_lengths_random_free_model(0, 20000000 + j, base_models_choice_list,
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
                        temporary_K2P.append(temporary)


                    if i == 'F81':
                        temporary = simulate_random_lengths_random_free_model(0, 30000000 + j, base_models_choice_list,
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
                        temporary_F81.append(temporary)


                    if i == 'HKY':
                        temporary = simulate_random_lengths_random_free_model(0, 50000000 + j, base_models_choice_list,
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
                        temporary_HKY.append(temporary)

                    if i == 'GTR':
                        temporary = simulate_random_lengths_random_free_model(0, 60000000 + j, base_models_choice_list,
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

                        temporary_GTR.append(temporary)

print('temporary_JC: ', temporary_JC)

print('temporary_K2P: ', temporary_K2P)

print('temporary_F81: ', temporary_F81)

print('temporary_HKY: ', temporary_HKY)

print('temporary_GTR: ', temporary_GTR)



print('Maximum likelihood reconstructions')
if args.na == 'DNA':
    print(count_correct_trees(temporary_JC))
    print(count_correct_trees(temporary_K2P))
    print(count_correct_trees(temporary_F81))
    print(count_correct_trees(temporary_HKY))
    print(count_correct_trees(temporary_GTR))

#if args.na == 'AA':
#    print(count_correct_trees(temporary_JTT))
#    print(count_correct_trees(temporary_LG))
#    print(count_correct_trees(temporary_WAG))
#    print(count_correct_trees(temporary_WAG_STAR))
#    print(count_correct_trees(temporary_DAY))