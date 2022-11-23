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
parser.add_argument('-sequence_type', type=str, choices=["DNA","AA"],
                    help="Compulsory argument. Nucleotide or Amino acids data. Enter DNA for nucleotide sequences. Enter AA for amino acid sequences.", required=True)
parser.add_argument('-NN_name',type=str, help="Compulsory argument. Enter the name of the NN folder.", required=True)
parser.add_argument('-alignment_file',type=str, help="Compulsory argument. Enter the name of the multiplealignment file.",required=True)

args = parser.parse_args()

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
str = args.alignment_file + '_' + args.NN_name + '_' + 'substitution_model.txt'
str = str.replace("/", "")
alignment_file = args.alignment_file

f = open(str, "w")
# load model
model = load_model(args.NN_name)
if args.sequence_type == 'DNA':
    if os.path.isfile(alignment_file):
        command = './quartet-pattern-counter-v1.1 ' + alignment_file + ' /dev/shm/out.npy'
        subprocess.run([command], shell=True)
        frequency_array = np.load('/dev/shm/out.npy')
        print(frequency_array)
        frequency_array = np.reshape(frequency_array,(1,-1))
        freq_ar = sc.fit_transform(frequency_array)
        prediction = model.predict(freq_ar)
        print(prediction)
        x = argmax(prediction)
        print(x)
        y = x.item()
        print(y)
        if y == 0:
            print('JC')
            f.write('JC')
        elif y == 1:
            print('K2P')
            f.write('K2P')
        elif y == 2:
            print('F81')
            f.write('F81')
        elif y == 3:
            print('HKY')
            f.write('HKY')
        else:
            print('GTR')
            f.write('GTR')
    else:
        print("The multiplealignment file does not exist!")
        print("Please try again.")
        sys.exit()
elif args.sequence_type == 'AA':
    if os.path.isfile(alignment_file):
        command = './quartet-pattern-counter-v1.1 -p ' + alignment_file + ' /dev/shm/out.npy'
        subprocess.run([command], shell=True)
        frequency_array = np.load('/dev/shm/out.npy')
        frequency_array = np.reshape(frequency_array, (1, -1))
        freq_ar = sc.fit_transform(frequency_array)
        prediction = model.predict(freq_ar)
        x = argmax(prediction)
        y = x.item()
        if y == 0:
            print('JTT')
            f.write('JTT')
        elif y == 1:
            print('LG')
            f.write('LG')
        elif y == 2:
            print('WAG')
            f.write('WAG')
        else:
            print('DAY')
            f.write('DAY')
    else:
        print("The multiplealignment file does not exist!")
        print("Please try again.")
        sys.exit()
f.close()