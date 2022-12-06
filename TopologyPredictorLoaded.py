#!/usr/bin/env python3

## Fewer nodes in inter layers. Less layers in u models.

import subprocess
import os
import itertools
import numpy as np
from   sklearn.naive_bayes import GaussianNB
import tensorflow as tf
import random
import math

import sys
import argparse
from numpy import argmax

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, train_test_split
from tensorboard.plugins.hparams import api as hp
from keras.models import load_model
from Bio import Phylo
from io import StringIO


# Argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument('-sequence_type', type=str, choices=["DNA","AA"],
                    help="Compulsory argument. Nucleotide or Amino acids data. Enter DNA for nucleotide sequences. Enter AA for amino acid sequences.", required=True)
parser.add_argument('-substitution_model',type=str, choices=['JC','K2P','F81','F84','HKY','GTR','JTT','LG','WAG_OLD','WAG','WAG_STAR','DAY'],
                    help="Compulsory argument. JC,K2P,F81,F84,HKY,GTR are nucleotide models. JTT,LG,WAG_OLD,WAG,WAG_STAR,DAY are amino acid models.", required=True)
parser.add_argument('-NN_name',type=str, help="Compulsory argument. Enter the name of the NN folder.", required=True)
parser.add_argument('-alignment_file',type=str, help="Compulsory argument. Enter the name of the multiplealignment file.",required=True)

args = parser.parse_args()

python_random_seed = 13

os.environ['MKL_NUM_THREADS'] = '10'
os.environ['GOTO_NUM_THREADS'] = '10'
os.environ['OMP_NUM_THREADS'] = '10'
os.environ['openmp'] = 'True'

tf.config.threading.set_inter_op_parallelism_threads(10)
tf.config.threading.set_intra_op_parallelism_threads(10)


def tree_topology(topology,m):
    if topology == 0:
        string = "1.0 Model_" + m + " ((A,B),(C,D))"
        print(string)
        return "((A,B),(C,D))"
    elif topology == 1:
        string = "1.0 Model_" + m + " ((A,C),(B,D))"
        print(string)
        return "((A,C),(B,D))"
    elif topology == 2:
        string = "1.0 Model_" + m + " ((A,D),(C,B))"
        print(string)
        return "((A,D),(C,B))"


model = load_model(args.NN_name)
alignment_file = args.alignment_file
sc = StandardScaler()

if args.sequence_type == 'DNA':
    command = './quartet-pattern-counter-v1.1 ' + alignment_file + ' /dev/shm/out.npy'
elif args.sequence_type == 'AA':
    command = './quartet-pattern-counter-v1.1 -p ' + alignment_file + ' /dev/shm/out.npy'
subprocess.run([command], shell=True)
frequency_array = np.load('/dev/shm/out.npy')
frequency_array = np.reshape(frequency_array,(1,-1))
#freq_ar = sc.fit_transform(frequency_array)
prediction = model.predict(frequency_array)
x = argmax(prediction)
print(prediction)
y = x.item()
tree = tree_topology(y,args.substitution_model)
tree = Phylo.read(StringIO(tree), 'newick')
Phylo.write(tree, "tree_topology_" + args.alignment_file + '_' + args.NN_name +".nwk", "newick")





