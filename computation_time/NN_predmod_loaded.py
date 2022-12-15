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
from numpy import argmax
#import keras_tuner as kt

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, train_test_split
from tensorboard.plugins.hparams import api as hp
from keras.models import load_model


# Argparse arguments

parser = argparse.ArgumentParser()
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

####################################
# my main program starts here:
####################################
sc = StandardScaler()



# load model
f = open("Time_NN_predmod_2.txt", "w")

path_of_the_directory = '/home/nkulikov/PycharmProjects/time_count/'
ext = ('.fas')
for files in os.listdir(path_of_the_directory):
    if files.endswith(ext):

        if 'F81_bl_0_top_0' or 'JC_bl_0_top_0' or 'K2P_bl_0_top_0' or 'F84_bl_0_top_0' or 'HKY_bl_0_top_0' or 'GTR_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + '.' + j)
                start = time.time()
                for i in range(100):
                    command = './quartet-pattern-counter-v1.1 ' + files + ' /dev/shm/out.npy'
                    subprocess.run([command], shell=True)
                    frequency_array = np.load('/dev/shm/out.npy')
                    frequency_array = np.reshape(frequency_array,(1,-1))
                    freq_ar = sc.fit_transform(frequency_array)
                    prediction = model.predict(freq_ar)
                    x = argmax(prediction)
                    y = x.item()
                    if y == 0:
                        print('Substitution model: JC')
                    elif y == 1:
                        print('Substitution model: K2P')
                    elif y == 2:
                        print('Substitution model: F81')
                    elif y == 3:
                        print('Substitution model: HKY')
                    else:
                        print('Substitution model: GTR')
            end = time.time()
            total_time = end - start
            f.write('File name: ')
            f.write(files)
            f.write('\n')
            f.write('Total time: ')
            f.write(str(total_time))
            f.write('\n')
            f.write('Average time: ')
            f.write(str(total_time / 100))
            f.write('\n')
            f.write('NN model: ')
            f.write(j)
            f.write('\n')
            f.write('\n')

        if 'JTT_bl_0_top_0' or 'LG_bl_0_top_0' or 'DAY_bl_0_top_0' or 'WAG_bl_0_top_0' in files:
            for j in ['b1','b2', 'b3', 'b10', 'u', 'cu', 'cud']:
                model = load_model('evomodel' + '.' + j)
                start = time.time()
                for i in range(100):
                    command = './quartet-pattern-counter-v1.1 -p ' + files + ' /dev/shm/out.npy'
                    subprocess.run([command], shell=True)
                    frequency_array = np.load('/dev/shm/out.npy')
                    frequency_array = np.reshape(frequency_array, (1, -1))
                    freq_ar = sc.fit_transform(frequency_array)
                    prediction = model.predict(freq_ar)
                    x = argmax(prediction)
                    y = x.item()
                    if y == 0:
                        print('Substitution model: JTT')
                    elif y == 1:
                        print('Substitution model: LG')
                    elif y == 2:
                        print('Substitution model: WAG')
                    else:
                        print('Substitution model: Dayhoff')

                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Total time: ')
                f.write(str(total_time))
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

f.close()

