#!/bin/bash


python3 ML_NJ_iqtree.py -na DNA -m JC -bl 0 -alg bionj 1> result_JC_NJ.txt 2> log_JC_NJ.log
python3 ML_NJ_iqtree.py -na DNA -m K2P -bl 0 -alg bionj 1> result_K2P_NJ.txt 2> log_K2P_NJ.log
python3 ML_NJ_iqtree.py -na DNA -m F81 -bl 0 -alg bionj 1> result_F81_NJ.txt 2> log_F81_NJ.log
python3 ML_NJ_iqtree.py -na DNA -m F84 -bl 0 -alg bionj 1> result_F84.txt 2> log_F84_NJ.log
python3 ML_NJ_iqtree.py -na DNA -m HKY -bl 0 -alg bionj 1> result_HKY_NJ.txt 2> log_HKY_NJ.log
python3 ML_NJ_iqtree.py -na DNA -m GTR -bl 0 -alg bionj 1> result_GTR_NJ.txt 2> log_GTR_NJ.log

