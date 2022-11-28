#!/bin/bash

while getopts s:n:a: flag
do
    case $flag in
        s) sequence_type=${OPTARG};;
        n) NN_name=${OPTARG};;
        a) alignment_file=${OPTARG};;
    esac
done
echo "${sequence_type}"
python3 ModelPredictorLoaded.py -sequence_type $sequence_type -NN_name $NN_name -alignment_file $alignment_file 
var1="_substitution_model.txt"
underscore="_"
substitution_model_file_name="$alignment_file$underscore$NN_name$var1"
echo "$substitution_model_file_name"
tag=$(tail -n 1 ${substitution_model_file_name}) 

echo "The substitution model is: $tag"
echo "Please enter the name of the topology predictor Neural network: "
read NeuralNetworkName  

python3 TopologyPredictorLoaded.py -sequence_type $sequence_type -NN_name $NeuralNetworkName -alignment_file $alignment_file -substitution_model $tag 
