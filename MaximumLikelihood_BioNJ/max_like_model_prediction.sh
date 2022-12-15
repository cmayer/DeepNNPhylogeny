#!/bin/bash


python3 ML_iqtree_evo_mod.py -na DNA -bl 0 1> result_ML_evo_mod.txt 2> log_ML_evo_mod.log
python3 ML_iqtree_evo_mod.py -na DNA -bl 1 1> result_ML_evo_mod_sb.txt 2> log_ML_evo_mod_sb.log

