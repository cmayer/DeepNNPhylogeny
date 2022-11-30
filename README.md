# DeepNNPhylogeny
A Deep Neural Networks for phylogenetic tree reonstructions.

## Setup
Conda environment installation 
```
conda create --name name_of_the_conda_environment python=3.9
conda activate name_of_the_conda_environment
conda install tensorflow
conda install sklearn
conda install -c anaconda scikit-learn 
conda install matplotlib pandas (optional)
conda install -c anaconda tensorflow-datasets 
conda install -c conda-forge kerastuner 
conda install -c conda-forge tensorflow-hub
```

## Make the programs executable 
```
unzip DeepNNPhylogeny-main.zip
cd DeepNNPhylogeny-main/
chmod u+x ModelPred_TopPred.sh
cd quartet-pattern-counter-v1.1/
./compile
cp quartet-pattern-counter-v1.1 /pathway/DeepNNPhylogeny-main/
```

## Neural networks training 
Quick start with the default parameters. 
For model prediction training run:
```
python3 ModelPredictorTraining.py -sequence_type DNA
```
for nucleotide sequences
or 
```
python3 ModelPredictorTraining.py -sequence_type AA
```
for amino acid sequences 

For model topology training run:
```
python3 TopologyPredictorTraining.py  -sequence_type * -substitution_model **
```
where * is DNA or AA, and <br />
** is 'JC','K2P','F81','F84','HKY','GTR' - nucleotide substitution models <br />
** is 'JTT','LG','WAG_OLD','WAG','WAG_STAR','DAY' - amino acid substitution models

To see all available parameters, their description and usage, run: 
```
python3 ModelPredictorTraining.py --help
python3 TopologyPredictorTraining.py --help
```
## Trained neural networks usage 

To use the trained substitution neural network run: 
```
python3 ModelPredictorLoaded.py -sequence_type * -NN_name ** -alignment_file ***
```
where * DNA or AA, and <br />
** is a name of the substitution model neural network predictor folder  <br />
*** is a name of the multiplealignments file <br />
For trained topology neural network run: 
```
python3 TopologyPredictorLoaded.py -sequence_type * -NN_name ** -alignment_file *** -substitution_model ****
```
where * DNA or AA, and <br />
** is a name of the substitution model neural network predictor folder  <br />
*** is a name of the multiplealignments file <br />
**** is 'JC','K2P','F81','F84','HKY','GTR' - nucleotide substitution models <br />
**** is 'JTT','LG','WAG_OLD','WAG','WAG_STAR','DAY' - amino acid substitution models <br />
<br />
It is possible to run both neural networks sequentially in the same program. <br />
First predict the substitution model, and then 
predict the tree topology based on the predicted substitution model. 
```
./ModelPred_TopPred.sh  -s * -n ** -a *** 
```
Where * is DNA or AA <br />
** is a name of the substitution model neural network predictor folder  <br />
*** is a name of the multiplealignments file <br />
After substituion model prediction it will ask you for an input. You should enter the name of the topology prediction. 
