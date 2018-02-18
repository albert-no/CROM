# Running Example

## Parameter
[main.cpp](main.cpp) is the main code.
In `main()`, user can set parameters
```
fname : name of qscore file
xdim : block length. this will be a number of rows in subqscore files. It should be a power of 2.
num_x : a number of columns.
rd_param : tuning parameter of r-d function (1.4 is reasonable for xdim=65536).
R_enc : encoding rate
mode : `e` for encoding and `d` for decoding
verbose : whether printing intermediate steps
file_idx : number of subqscore files
num_dec_pts : number of target decoding rates
R_dec : target decoding rates
```

## Execution
```
make
./crom_example.o
```

## Clean up
To clean up the intermediate files, run the python code.
```
python cleanup.py
```
