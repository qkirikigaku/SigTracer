# SigTracer

SigTracer provides a clone decomposition based on mutation signatures for one tumor sample like below:

![sample-1_clone](https://user-images.githubusercontent.com/26032572/111058359-73c12080-84d1-11eb-9b54-2bdeb1af945f.png)

## Quick start

## Manual insatll

### Requirements
SigTracer is implemented in C++ and Python 3.5 with some packages.

For C++, you need the `Boost` package to calculate digamma and gamma functions.
  
According to your environment, please install the package from the official site ( https://www.boost.org/ ) and add the path to the compiler.

For Python, you need some packages listed below:
* numpy
* scipy
* subprocess
* pickle
* itertools

### Compile
clone the repository and execute the command:
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

The execution result will be output to the `result` directory.
