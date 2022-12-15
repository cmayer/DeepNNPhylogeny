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
f = open("Time_BioNJ.txt", "w")
path_of_the_directory = '/home/nkulikov/PycharmProjects/time_count/'
ext = ('.fas')
for files in os.listdir(path_of_the_directory):
    if files.endswith(ext):
       if 'F81' in files:
           start = time.time()
           for i in range(1000):
               command = "iqtree2.1.3 -fast -mset JC,K2P,F81,F84,HKY,GTR -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'F84' in files:
           start = time.time()
           for i in range(1000):
               command = "iqtree2.1.3 -fast -mset JC,K2P,F81,F84,HKY,GTR -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'JC' in files:
           start = time.time()
           for i in range(1000):
               command = "iqtree2.1.3 -fast -mset JC,K2P,F81,F84,HKY,GTR -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'K2P' in files:
           start = time.time()
           for i in range(1000):
               command = "iqtree2.1.3 -fast -mset JC,K2P,F81,F84,HKY,GTR -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'GTR' in files:
           start = time.time()
           for i in range(1000):
               command = "iqtree2.1.3 -fast -mset JC,K2P,F81,F84,HKY,GTR -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'HKY' in files:
           start = time.time()
           for i in range(1000):
               command = "iqtree2.1.3 -fast -mset JC,K2P,F81,F84,HKY,GTR -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'JTT' in files:
           start = time.time()
           for i in range(100):
               command = "iqtree2.1.3 -fast -mset JTT,LG,WAG,Dayhoff -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'LG' in files:
           start = time.time()
           for i in range(100):
               command = "iqtree2.1.3 -fast -mset JTT,LG,WAG,Dayhoff -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'WAG_bl' in files:
           start = time.time()
           for i in range(100):
               command = "iqtree2.1.3 -fast -mset JTT,LG,WAG,Dayhoff -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')
       elif 'DAY' in files:
           start = time.time()
           for i in range(100):
               command = "iqtree2.1.3 -fast -mset JTT,LG,WAG,Dayhoff -redo -s " + files
               subprocess.run([command], shell=True)
           end = time.time()
           total_time = end - start
           f.write('File name: ')
           f.write(files)
           f.write('\n')
           f.write('Total time: ')
           f.write(str(total_time))
           f.write('\n')
           f.write('Average time: ')
           f.write(str(total_time / 1000))
           f.write('\n')

f.close()