# DNA_Data_Pipeline
A Modular Data Pipeline for DNA Data Storage (MSc Project)

## Installation and Setup

### Pipeline
To create pipeline environment:
```
git clone --recursive https://github.com/ryan-tsk/DNA_Data_Pipeline
conda create --name DNApipeline python=3.8
conda activate DNApipeline
pip install -r requirements.txt
git clone https://github.com/Omer-Sella/turboDNA
git clone https://github.com/Omer-Sella/Euclid
```

### Bonito
To create Bonito environment:
```
conda create --name Bonito python=3.8
conda activate bonito
cd bonito
pip install -r requirements.txt
bonito download --models
```
To use the latest Bonito version instead of the test pipeline version:
```
pip install ont-bonito
```

### Chiron 
To create Chiron environment:
```
conda create --name Bonito python=3.7
pip install Chiron
pip install tensorflow-gpu==1.15
pip install protobuf==3.20.*
conda install cudatoolkit=10.0
conda install -c anaconda cudnn=7
```
The above set up is for GPU-based Chiron (recommended).
For CPU base, replace the tensorflow line as follow:
```
pip install tensorflow-gpu==1.15
```
Note that cudatoolkit and cudnn is not needed for CPU base
### Alignment Programs

The alignment programs supported by this pipeline are ClustalW2, ClustalO, MUSCLE and PRANK  

1. Make a directory for each of the programs
2. Change directory - download and install each binary
3. Export relevant binary path

```
mkdir alignments/bin/clustalo
cd alignments/bin/clustalo
wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64
export PATH=$PATH:/alignments/bin/clustalo
```

## User Guide

### Configuration and Modules
Place all function files in one modules folder  
Refer to template.yaml for more information on creating Configuration files

### Running the Pipeline
Use main.py to run the pipeline:
```
python3 main.py -c=config/config_A1.yaml -s=True
```
Refer to main.py --help for supported arguments





