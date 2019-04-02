## work flow
### (1) compile main_proj_2.cpp

the parameters are inside the source code file

however, only seed is the outside parameter, this example takes 33 as seed

```console
linuxServer~$: main_proj_2 33
```

### (2) copy main_proj_2 to folder data1 and submit job array

the output file is in format seed_sequence.dat
```console
...
002_00015.dat
002_00016.dat
002_00017.dat
...
```

### (3) analysis the data

In the folder dataAnalysis, there are two files.

data2.ipynb will read the step (2) folder and convert to best and error two data files.

analyticContinuation.ipynb is using Python to do Pade Regession, still not finished writing.

```
analyticContinuation.ipynb
data2.ipynb
```
