#!/bin/bash


python3 ML_NJ_iqtree.py -na DNA -m JC -bl 0 -alg ml  1> result_JC_max_like.txt 2> log_JC_max_like.log
python3 ML_NJ_iqtree.py -na DNA -m K2P -bl 0 -alg ml 1> result_K2P_max_like.txt 2> log_K2P_max_like.log
python3 ML_NJ_iqtree.py -na DNA -m F81 -bl 0 -alg ml 1> result_F81_max_like.txt 2> log_F81_max_like.log
python3 ML_NJ_iqtree.py -na DNA -m F84 -bl 0 -alg ml 1> result_F84_max_like.txt 2> log_F84_max_like.log
python3 ML_NJ_iqtree.py -na DNA -m HKY -bl 0 -alg ml 1> result_HKY_max_like.txt 2> log_HKY_max_like.log
python3 ML_NJ_iqtree.py -na DNA -m GTR -bl 0 -alg ml 1> result_GTR_max_like.txt 2> log_GTR_max_like.log
