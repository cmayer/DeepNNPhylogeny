#!/usr/bin/env python3

## Import libraries

import subprocess
import os
import numpy as np
import tensorflow as tf
import sys
import argparse
from numpy import argmax
from keras.models import load_model
from configparser import ConfigParser
# These libraries will be needed when tree length prediction is introduced 
#from Bio import Phylo
#from io import StringIO 


# Argparse arguments

parser = argparse.ArgumentParser()
parser.add_argument('-sequence_type', type=str, choices=["DNA","AA"],
                    help="Compulsory argument. Nucleotide or Amino acids data. Enter DNA for nucleotide sequences. Enter AA for amino acid sequences.", required=True)
parser.add_argument('-substitution_model',type=str, choices=['JC','K2P','F81','F84','HKY','GTR','JTT','LG','WAG_OLD','WAG','WAG_STAR','DAY'],
                    help="Compulsory argument. JC,K2P,F81,F84,HKY,GTR are nucleotide models. JTT,LG,WAG_OLD,WAG,WAG_STAR,DAY are amino acid models.", required=True)
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


def is_quartet_counter_available():
    try:
        subprocess.check_output(['which', 'quartet-pattern-counter-v1.1'])
        return True
    except subprocess.CalledProcessError:
        return False

def read_alignment_file(alignment_file):
    f = open(alignment_file, "r")
    return_list = [] 
    for line in f :
        line = line.rstrip()
        if ">" in line : 
            line = line.replace(">","")
            line = line.split(" ")
            return_list.append(line[0])
    return return_list          

def tree_topology(topology,m):
    if topology == 0:
        header_list = read_alignment_file(args.alignment_file)
        string = "((" + str(header_list[0]) + ","  + str(header_list[1]) + "),(" + str(header_list[2]) + "," + str(header_list[3]) + "))"
        print("Model:",m," Tree topology: ",string)
        return string
    elif topology == 1:
        header_list = read_alignment_file(args.alignment_file)
        string = "((" + str(header_list[0]) + "," + str(header_list[2]) + "),(" + str(header_list[1]) + "," + str(header_list[3]) + "))"
        print("Model:",m," Tree topology: ",string)
        return string
    elif topology == 2:
        header_list = read_alignment_file(args.alignment_file)
        string = "((" + str(header_list[0]) + "," + str(header_list[3]) + "),(" + str(header_list[2]) + "," + str(header_list[1]) + "))"
        print("Model:",m," Tree topology: ",string)
        return string

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

# Check whether the config file exists 
config = check_read_config()
# Check whether the NN exist and load model 
model = search_for_the_NN(config)

model = load_model(args.NN_name)
alignment_file = args.alignment_file

# Check the quartet-pattern-counter-v1.1 availability
if is_quartet_counter_available():
    print("quartet-pattern-counter-v1.1 is available")
else:
    print("quartet-pattern-counter-v1.1 is not available")


if args.sequence_type == 'DNA':
    command = 'quartet-pattern-counter-v1.1 ' + alignment_file + " " + os.getcwd() + "/out.npy"
elif args.sequence_type == 'AA':
    command = 'quartet-pattern-counter-v1.1 -p ' + alignment_file + " " + os.getcwd() + "/out.npy"
subprocess.run([command], shell=True)
path = os.getcwd() + "/out.npy"
frequency_array = np.load(path)
frequency_array = np.reshape(frequency_array,(1,-1))
prediction = model.predict(frequency_array)
x = argmax(prediction)
print('The softmax values of topologies: ', prediction)
y = x.item()
tree = tree_topology(y,args.substitution_model)
temporary_str = args.alignment_file.replace(".fas","")
neural_name = args.NN_name.replace("/","")
os.remove("out.npy")
out_file =  "tree_topology_" + temporary_str + '_' + neural_name +".nwk"
f = open(out_file, "w")
f.write(tree)
f.write("\n")
f.close()

# These lines will be needed when tree length prediction is introduced 
# tree = Phylo.read(StringIO(tree), 'newick')
# Phylo.write(tree, "tree_topology_" + temporary_str + '_' + neural_name +".nwk", "newick")
