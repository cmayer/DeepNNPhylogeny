# DeepNNPhylogeny
In this project we developed Deep Neural Networks for phylogenetic tree reconstruction and evolutionary model selection. Currently the networks are limited to quartet trees. On quartet trees the neural network classifiers perform in most cases as good as the maximum likelihood method, which means that it is not significantly better or worse than maximum likelihood.
In very few scenarios the neural network is marginally but significantly inferior to the maximum likelihood method, but we think that neural networks and training can be improved such that both methods perform equally good.<br>

This repository contains the software for training und using neural networks for the following tasks:
- model prediction for nucleotide alignments
- model prediction for amino acid alignments
- topology prediction for nucleotide alignments
- topology prediction for amino acid alignments

## Installing the required machine learning libraries using the Anaconda package manager:
It is recommended to create a conda environment for your tensorflow module: 
```
conda create --name name_of_the_conda_environment 
conda activate name_of_the_conda_environment
conda install tensorflow
```

## On Linux and Mac OS X the DeepNNPhylogeny software can be installed as follows:

- Download the DeepNNPhylogeny archive or clone the github repository locally.
- If you downloaded the archive: type on the command line:
```
unzip DeepNNPhylogeny-main.zip
```
Next compile the quartet-pattern-counter-v1.1 program:
```
cd DeepNNPhylogeny-main/
chmod u+x ModelPred_TopPred.sh
cd quartet-pattern-counter-v1.1_src/
make
# If "make" is not installed, type: chmod u+x compile.sh; ./compile.sh
```
- Make sure that you copy the compiled quartet-pattern-counter-v1.1 program to a folder that is listed in your $PATH variable so that your system can always find it, or copy it to the folder you want to use it from, or specify the full path to the program.

## Training your own neural networks:

Pre-trained neural networks can be downloaded from [DryAd](https://doi.org/10.5061/dryad.ksn02v783). If you prefer to use pre-trained models, you can skip next section.<br>

**You can train your own neural networks as follows:<br>**

Training neural networks is done by simulating a large number of data sets and using the pattern frequency vectors together with the known topology
to train a neural network to classify the topology or the model of sequence evolution. Simulations are conducted with the software [PolyMoSim available on github](https://github.com/cmayer/PolyMoSim). PolyMoSim has to be installed and must be available in the system path in order to run training tasks. The PolyMoSim software is not required if you only want to predict/classify models or topologies.

In order to train a neural network for model prediction run:
```
python3 ModelPredictorTraining.py -sequence_type DNA
```
for nucleotide sequences,
or 
```
python3 ModelPredictorTraining.py -sequence_type AA
```
for amino acid sequences. 

In order to train a neural network with default parameters for topology prediction run:
```
python3 TopologyPredictorTraining.py  -sequence_type (*) -substitution_model (**)
```
where (\*) is DNA or AA, and <br />
(\*\*) is 'JC', 'K2P', 'F81', 'F84', 'HKY', 'GTR' - nucleotide substitution models <br />
(\*\*) is 'JTT', 'LG', 'WAG_OLD', 'WAG','WAG_STAR', 'DAY' - amino acid substitution models


To see all available parameters, their description and usage, run: 
```
python3 ModelPredictorTraining.py --help
python3 TopologyPredictorTraining.py --help
```

Simulating amino acid data sets takes much longer than simulating nucleotide data sets.
For a large number of amino acid replicates, we recommend to use multiprocessing library (will be added to the program soon) and to conduct the training on a computer with a large number of core. 


## Topology and evolutionary model predictions/classifications:

Topology and model predictions require the quartet-pattern-counter program to be in your system path or in the directory you run the python programs in.

You can use pre-trained neural networks that can be downloaded from [DryAd](https://doi.org/10.5061/dryad.ksn02v783) or you can train your own neural networks.

## Predicting models of sequence evolution for user specified alignments using pre-trained models: 

Download pre-trained neural networks from [DryAd](https://doi.org/10.5061/dryad.ksn02v783) or use a model you have trained yourself.
The programs search for the neural network in the working directory, specified path to the neural network, or DeepNNPhylogeny.config file.
Make sure that you have placed a DeepNNPhylogeny.config file in the working directory, the home directory, or the DeepNNPhylogeny-main folder placed in the home directory. 
The default content of the DeepNNPhylogeny.config is the DryAD folders contained pre-trained NNs. 

```
python3 ModelPredictorLoaded.py -sequence_type (*) -NN_name (**) -alignment_file (***)
```
where (\*) is DNA or AA <br>
(\*\*)  is a name of the substitution model neural network predictor folder  <br>
(\*\*\*) is a name of the multiple-sequence-alignment file <br>
For trained topology neural network run: 
```
python3 TopologyPredictorLoaded.py -sequence_type (*) -NN_name (**) -alignment_file (***) -substitution_model (****)
```
where (\*) DNA or AA <br />
(\*\*) is a name of the substitution model neural network predictor folder  <br />
(\*\*\*) is a name of the multiple-sequence-alignment file <br />
(\*\*\*\*) is 'JC','K2P','F81','F84','HKY','GTR' - nucleotide substitution models <br />
(\*\*\*\*) is 'JTT','LG','WAG_OLD','WAG','WAG_STAR','DAY' - amino acid substitution models <br />
<br />
It is possible to run both neural networks sequentially in the same program. <br />
First predict the substitution model, and then 
predict the tree topology based on the predicted substitution model. 
```
./ModelPred_TopPred.sh  -s (*) -n (**) -a (***) 
```
where (\*) is DNA or AA <br />
(\*\*) is a name of the substitution model neural network predictor folder  <br />
(\*\*\*) is a name of the multiple-sequence-alignment file <br />
After substitution model prediction it will ask you for an input. You should enter the name of the topology prediction NN.

## How to cite DeepNNPhylogeny:
Machine learning can be as good as maximum likelihood when reconstructing phylogenetic trees and determining the best evolutionary model on four taxon alignments

https://www.biorxiv.org/content/10.1101/2023.07.12.548770v1.article-metrics



## Frequently asked questions <a id="Frequently-aksed-questions"></a>
No questions have been asked so far.

