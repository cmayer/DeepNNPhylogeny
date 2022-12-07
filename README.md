# DeepNNPhylogeny
A Deep Neural Networks for phylogenetic tree reonstructions. Currently it is limited to quartet trees. On quartet trees it performs in most cases as good as the maximum likelihood method.

## Setup
Conda environment installation 
```
conda create --name name_of_the_conda_environment python tensorflow -c anaconda scikit-learn
conda activate name_of_the_conda_environment
```

## On Linux and Mac OS X the programs can be installed as follows:
```
- Download the DeepNNPhylogeny archive or clone the github repository locally.
- If you downloaded the archive: type on the command line: 
unzip DeepNNPhylogeny-main.zip
cd DeepNNPhylogeny-main/
chmod u+x ModelPred_TopPred.sh
cd quartet-pattern-counter-v1.1/
./compile
cp quartet-pattern-counter-v1.1 /pathway/DeepNNPhylogeny-main/
```

Topology and model predictions/classifications require the quartet-pattern-counter to be in your system path or in the directory you run the python programs in.
(Needs to be verified.)

## Pretrained neural networks can be downloaded on Dryad. You can train your own neural networks as follows: 
Quick start with the default parameters.

Training neural networks is done by simulating a large number of data sets and using the pattern frequency vectors together with the known topology
to train a neural network to classify the topology or the model of sequence evolution. Simulations are conducted with the software [PolyMoSim avaiable on github](https://github.com/cmayer/PolyMoSim). PolyMoSim has to be installed and must be avaiable in the system path in order to run training tasks. The PolyMoSim software is not required if you only want to predict/classify models or topologies.

In order to train a neural network for model prediction run:
```
python3 ModelPredictorTraining.py -sequence_type DNA
```
for nucleotide sequences
or 
```
python3 ModelPredictorTraining.py -sequence_type AA
```
for amino acid sequences 

In order to train a neural network for topology prediction run:
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
## Predicting models or topologies for user specified alignemtns using pretrained models: 

To use the trained substitution neural network run: 

Download pre-trained neural networks from [xxx](https://www.dryad.com) and

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
