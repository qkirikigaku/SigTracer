# SigTracer

SigTracer provides a clone decomposition based on mutation signatures for one tumor sample like below:

![sample-1_clone](https://user-images.githubusercontent.com/26032572/111058359-73c12080-84d1-11eb-9b54-2bdeb1af945f.png)

## Quick start

We have prepared Docker image file that completed the installation of SigTracer.
You can easily download from https://hub.docker.com/r/taromss/sigtracer , and run the container with:

```
docker run -it taromss/sigtracer
cd /root/SigTracer
```

If you use this Docker image, following manual installation is not necessary at all.

## Manual installation

You can install SigTracer in the local environment with the following procedure.

### Requirements
SigTracer is implemented in C++ and Python 3.5 with some packages.

For C++, you need the **Boost** package to calculate digamma and gamma functions.
  
According to your environment, please install the package following the official site ( https://www.boost.org/ ) and add the path to the compiler.

For Python, you need some packages listed below:
* numpy
* scipy
* matplotlib
* seaborn

### Compile
Clone the repository and execute the command:
```
git clone https://github.com/qkirikigaku/SigTracer
cd SigTracer
make compile
```

## Running example with simulation data
We have prepared a simulated tumor sample following the generative process of SigTracer, and you can try to run SigTracer to reproduce the result with a simple command:

```
python SigTracer.py example
```

If your environment is capable to run with parallel core, you can easily parallelized experiments by changing `SigTracer.py` as below:

```python:SigTracer.py
num_processes=8 # FIXME according to your environment.
```

The execution result will be output to the `result/example` directory.

## Usage with an actual data

You should prepare the following input directry.

### Input directory (`SigTracer/data`) structure

```
data
├-- ref
│    ├-- MT_SBS.txt
│    └-- signature_probability_v3.1.csv
└-- ${exp_name}
       ├-- purity.csv
       ├-- ref_sig.csv
       ├-- ${sample_name1}.csv
       ├-- ${sample_name2}.csv
       ....
       └-- ${sample_nameS}.csv
```

The directory, `data/${exp_name}`, have three types of files, `purity.csv`, `ref_sig.csv`, and `${sample_name}.csv` for each sample.
Please refer to `data/example_real` directory as an example of acutal sequenced data, and the directory includes:

* `purity.csv`

It shows the purity, the proportion of the reads from cancer cells in all the sequenced reads, and it contains a header line and purity values for each sample:
```
Sample_name,purity
${sample_name1},0.9929135909899457
${sample_name2},0.9934101600184388
.....
${sample_nameS},0.5705353718426889
```

* `ref_sig.csv`

It shows the signature active in each sample:
```
${sample_name1},SBS5,SBS7a,SBS7c,SBS8,SBS9,SBS10a,SBS10b
${sample_name2},SBS3,SBS6,SBS7a
.....
${sample_nameS},SBS1,SBS3,SBS4,SBS5,SBS6,SBS7a,SBS9,SBS10a
```

* `${sample_name}.csv`

For one sample, namely ${sample_name}, all the mutations are listed in this file:
```
mutation_id,chromosome,position,ref_counts,var_counts,normal_cn,minor_cn,mut_cn,major_cn,total_cn,trinucleotide,annotation
mut_1,1,5,257,21,2,1,1,2,3,61,SBS59-target
mut_2,1,10,224,87,2,0,1,4,4,44,SBS13-target
.....
mut_2000,1,10000,64,13,2,1,1,1,2,51,SBS17b-target
```

With the actual sequenced data, since we do not know the values of `mut_cn`, please fill in the zeros as needed (these values are only used for evaluation with the simulation data).
`normal_cn` and `major_cn` are necessary for prediction, and they should be estimated by other tools and filled with estimators.
`annotation` column shows the mutated region for each mutation, and we listed non-existent gene names as examples here.

### Running SigTracer
Before running SigTracer, please fix the corresponding line in `SigTracer.py` to change to the mode for acutal data as following:
```
Simulation_flag = False # FIXME to change "False" when applying an actual data.
```

Using inputs in `data/example_real`, execute the following command resulting in parameter estimation:
```
python SigTracer.py example_real
```
Note that `data/example_real` includes 100 artificial samples, so it may take a lot of time.
For this reason, we recommend running the program in an environment where massively parallel computing is possible.

As is the above simulation, the execution result will be output to the `result/example_real` directory.

### Statistical test for measuring the relationship between mutations and signatures
Detailed information is described in our original paper:, you can extract the mutations related with certain signatures.

After parameter estimation with the above procedure, run the following command:

```
python TestClone.py example_real
```

The result will be output to the `result/example_real/tables/summary/test_mutation.tsv`.
