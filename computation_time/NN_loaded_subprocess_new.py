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
from Bio import Phylo
from io import StringIO


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

def tree_topology(topology,m):
    if topology == 0:
#        string = "1.0 Model_" + m + " ((A,B),(C,D))"
        return "((A,B),(C,D))"
    elif topology == 1:
#        string = "1.0 Model_" + m + " ((A,C),(B,D))"
        return "((A,C),(B,D))"
    elif topology == 2:
#        string = "1.0 Model_" + m + " ((A,D),(C,B))"
        return "((A,D),(C,B))"




# load model
f = open("Time_NN_top_new.txt", "w")
k = open('Newick_trees.txt', 'w')

path_of_the_directory = '/home/nkulikov/PycharmProjects/time_count/'
ext = ('.fas')
for files in os.listdir(path_of_the_directory):
    if files.endswith(ext):

        if 'F81_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + 'F81' + '.' + j)
                m = 'F81'
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
                    tree = tree_topology(y,m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'HKY_bl_0_top_0' in files:
            for j in ['b1','b2', 'b3', 'b10', 'u', 'cu', 'cud']:
                model = load_model('model' + 'HKY' + '.' + j)
                m = 'HKY'
                start = time.time()
                for i in range(100):
                    command = './quartet-pattern-counter-v1.1 ' + files + ' /dev/shm/out.npy'
                    subprocess.run([command], shell=True)
                    frequency_array = np.load('/dev/shm/out.npy')
                    frequency_array = np.reshape(frequency_array, (1, -1))
                    freq_ar = sc.fit_transform(frequency_array)
                    prediction = model.predict(freq_ar)
                    x = argmax(prediction)
                    y = x.item()
                    tree = tree_topology(y, m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'K2P_bl_0_top_0' in files :
            for j in ['b1','b2', 'b3', 'b10', 'u', 'cu', 'cud']:
                model = load_model('model' + 'K2P' + '.' + j)
                m = 'K2P'
                start = time.time()
                for i in range(100):
                    command = './quartet-pattern-counter-v1.1 ' + files + ' /dev/shm/out.npy'
                    subprocess.run([command], shell=True)
                    frequency_array = np.load('/dev/shm/out.npy')
                    frequency_array = np.reshape(frequency_array, (1, -1))
                    freq_ar = sc.fit_transform(frequency_array)
                    prediction = model.predict(freq_ar)
                    x = argmax(prediction)
                    y = x.item()
                    tree = tree_topology(y, m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Total time: ')
                f.write(str(total_time))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'JC_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + 'JC' + '.' + j)
                m = 'JC'
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
                    tree = tree_topology(y,m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'F84_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + 'F84' + '.' + j)
                m = 'F84'
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
                    tree = tree_topology(y,m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'GTR_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + 'GTR' + '.' + j)
                m = 'GTR'
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
                    tree = tree_topology(y,m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'JTT_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + 'JTT' + '.' + j)
                m = 'JTT'
                start = time.time()
                for i in range(100):
                    command = './quartet-pattern-counter-v1.1 -p ' + files + ' /dev/shm/out.npy'
                    subprocess.run([command], shell=True)
                    frequency_array = np.load('/dev/shm/out.npy')
                    frequency_array = np.reshape(frequency_array,(1,-1))
                    freq_ar = sc.fit_transform(frequency_array)
                    prediction = model.predict(freq_ar)
                    x = argmax(prediction)
                    y = x.item()
                    tree = tree_topology(y,m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'LG_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + 'LG' + '.' + j)
                m = 'LG'
                start = time.time()
                for i in range(100):
                    command = './quartet-pattern-counter-v1.1 -p ' + files + ' /dev/shm/out.npy'
                    subprocess.run([command], shell=True)
                    frequency_array = np.load('/dev/shm/out.npy')
                    frequency_array = np.reshape(frequency_array,(1,-1))
                    freq_ar = sc.fit_transform(frequency_array)
                    prediction = model.predict(freq_ar)
                    x = argmax(prediction)
                    y = x.item()
                    tree = tree_topology(y,m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'DAY_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + 'DAY' + '.' + j)
                m = 'DAY'
                start = time.time()
                for i in range(100):
                    command = './quartet-pattern-counter-v1.1 -p ' + files + ' /dev/shm/out.npy'
                    subprocess.run([command], shell=True)
                    frequency_array = np.load('/dev/shm/out.npy')
                    frequency_array = np.reshape(frequency_array,(1,-1))
                    freq_ar = sc.fit_transform(frequency_array)
                    prediction = model.predict(freq_ar)
                    x = argmax(prediction)
                    y = x.item()
                    tree = tree_topology(y,m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')

        if 'WAG_bl_0_top_0' in files:
            for j in ['b1','b2','b3','b10','u','cu','cud']:
                model = load_model('model' + 'WAG' + '.' + j)
                m = 'WAG'
                start = time.time()
                for i in range(100):
                    command = './quartet-pattern-counter-v1.1 -p ' + files + ' /dev/shm/out.npy'
                    subprocess.run([command], shell=True)
                    frequency_array = np.load('/dev/shm/out.npy')
                    frequency_array = np.reshape(frequency_array,(1,-1))
                    freq_ar = sc.fit_transform(frequency_array)
                    prediction = model.predict(freq_ar)
                    x = argmax(prediction)
                    y = x.item()
                    tree = tree_topology(y,m)
                    k.write(tree)
                    k.write('\n')
                end = time.time()
                total_time = end - start
                f.write('File name: ')
                f.write(files)
                f.write('\n')
                f.write('Average time: ')
                f.write(str(total_time / 100))
                f.write('\n')
                f.write('NN model: ')
                f.write(j)
                f.write('\n')
                f.write('\n')
#    print('Prediction', prediction)
f.close()
k.close()
