# DeepNNPhylogeny
A Deep Neural Networks for phylogenetic tree reonstructions.

+## Setup
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

+## Make the programs executable 

unzip DeepNNPhylogeny-main.zip
cd DeepNNPhylogeny-main/
chmod u+x ModelPred_TopPred.sh
cd quartet-pattern-counter-v1.1/
./compile
cp quartet-pattern-counter-v1.1 /pathway/DeepNNPhylogeny-main/

## Neural networks training 

