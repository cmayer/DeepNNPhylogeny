# DeepNNPhylogeny
In this project we developed Deep Neural Networks for phylogenetic tree reconstructions. Currently it is limited to quartet trees. On quartet trees it performs in most cases as good as the maximum likelihood method, which means that it is not significantly better or worse than maximum likelyhood.
In very few scenarios the neural network is significantly inferior than the maximum likelihood method, but the differences are so small that we think that neural networks and training can be improved such that both methods perform equally good.<br>

This repository contains the software for training und using neural networks for the following tasks:
- model prediction for nucleotide alignments
- topology prediction for nucleotide alignments
- model prediction for amino acid alignments
- topology prediction for amino acid alignments

## Installing the required machine learning libraries using the Anaconda package manager:
It is recommended to create a conda environment for your tensorflow module: 
```
conda create --name name_of_the_conda_environment python tensorflow scikit-learn
conda activate name_of_the_conda_environment
```

## On Linux and Mac OS X the DeepNNPhylogeny software can be installed as follows:
```
- Download the DeepNNPhylogeny archive or clone the github repository locally.
- If you downloaded the archive: type on the command line:
unzip DeepNNPhylogeny-main.zip
cd DeepNNPhylogeny-main/
chmod u+x ModelPred_TopPred.sh
cd quartet-pattern-counter-v1.1_src/
make
# If make is not installed, type: chmod u+x compile.sh; ./compile.sh
cp quartet-pattern-counter-v1.1 /pathway/DeepNNPhylogeny-main/
```

## Training your own neural networks:

Pre-trained neural networks can be downloaded from Dryad. If you prefer to use pre-trained models, you can skip this section.<br>

**You can train your own neural networks as follows:<br>**

Training neural networks is done by simulating a large number of data sets and using the pattern frequency vectors together with the known topology
to train a neural network to classify the topology or the model of sequence evolution. Simulations are conducted with the software [PolyMoSim available on github](https://github.com/cmayer/PolyMoSim). PolyMoSim has to be installed and must be available in the system path in order to run training tasks. The PolyMoSim software is not required if you only want to predict/classify models or topologies.

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

## Topology and evolutionary model predictions/classifications:

Topology and model predictions require the quartet-pattern-counter program to be in your system path or in the directory you run the python programs in.
(Needs to be verified.)

You can use pre-trained neural networks that can be downloaded from Dryad or you can train your own neural networks.

## Predicting models of sequence evolution for user specified alignments using pre-trained models: 

Download pre-trained neural networks from [xxx](https://www.dryadcom) or use a model you have trained yourself.

```
python3 ModelPredictorLoaded.py -sequence_type (*) -NN_name (**) -alignment_file (***)
```
where<br>
(*)   is DNA or AA, and <br>
(**)  is a name of the substitution model neural network predictor folder  <br>
(***) is a name of the multiple-sequence-alignment file <br>
For trained topology neural network run: 
```
python3 TopologyPredictorLoaded.py -sequence_type * -NN_name ** -alignment_file *** -substitution_model ****
```
where * DNA or AA, and <br />
** is a name of the substitution model neural network predictor folder  <br />
*** is a name of the multiple-sequence-alignment file <br />
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
*** is a name of the multiple-sequence-alignment file <br />
After substitution model prediction it will ask you for an input. You should enter the name of the topology prediction.


