#!/usr/bin/env python3

## Import libraries

import subprocess
import os
import itertools
import numpy as np
import tensorflow as tf
import random
import math
import sys
import argparse
from numpy import argmax
from keras.models import load_model
from configparser import ConfigParser



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

config = ConfigParser()
config_file = "DeepNNPhylogeny.config"

def quartet_pattern(q):
    command = q
    output = subprocess.check_output([command], shell=True)
    output = str(output)
    if output.startswith('Usage:'): 
        pass
    else:
        print('quartet-pattern-counter-v1.1 was not found!')
        sys.exit()

def check_read_config():
    first_check = os.getcwd() + "/" + config_file
    home_dir = os.path.expanduser("~") + "/" 
    second_check = home_dir + config_file
    third_check = home_dir + "DeepNNPhylogeny-main/" + config_file 
    if os.path.isfile(first_check):
        print("The config file: ", first_check, " was found!")
        config.read(first_check)
        return config
    elif os.path.isfile(second_check):
        print("The config file: ", second_check, " was found!")
        config.read(second_check)
        return config
    elif os.path.isfile(third_check):
        print("The config file: ", third_check, " was found!")
        config.read(third_check)
        return config
    else:     
        print("Configuration file was not found!")
        sys.exit()


def search_for_the_NN(conf):
    if os.path.isdir(args.NN_name):
        model = load_model(args.NN_name)
        return model
    elif (os.path.isdir(args.NN_name) == False):
        NN_name = args.NN_name.replace("/", "")
        config_NN = conf["NN-Search-Path"]
        for i in range(0,13):
            pathway_to_NN = config_NN[i]
            pathway_to_NN = pathway_to_NN + NN_name
            if os.path.isdir(pathway_to_NN):
                model = load_model(pathway_to_NN)
                break
        return model
    else: 
        print("Neural network model was not found!")
        sys.exit()
    
####################################
# my main program starts here:
####################################

str = args.alignment_file + '_' + args.NN_name + '_' + 'substitution_model.txt'
str = str.replace("/", "")

alignment_file = args.alignment_file

f = open(str, "w")

# load model

# Check whether the config file exists 
config = check_read_config()

# Check whether the NN exist and load model 
model = search_for_the_NN(config)

if args.sequence_type == 'DNA':
    quartet_pattern("quartet-pattern-counter-v1.1")
    if os.path.isfile(alignment_file):
        command = 'quartet-pattern-counter-v1.1 ' + alignment_file + " " + os.getcwd() + "/out.npy"
        path =  os.getcwd() + "/out.npy"
        subprocess.run([command], shell=True)
        frequency_array = np.load(path)
        frequency_array = np.reshape(frequency_array,(1,-1))
        prediction = model.predict(frequency_array)
        print('softmax values: ', prediction)
        x = argmax(prediction)
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
        command = 'quartet-pattern-counter-v1.1 -p ' + alignment_file + " " + os.getcwd() + "/out.npy"
        path =  os.getcwd() + "/out.npy"
        subprocess.run([command], shell=True)
        frequency_array = np.load(path)
        frequency_array = np.reshape(frequency_array, (1, -1))
        prediction = model.predict(frequency_array)
        print('softmax values: ', prediction)
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