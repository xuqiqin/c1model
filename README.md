# c1model
modeling of Pol II recycling after CTDP1 degradation
## R dependencies
```
R=3.6.1
cowplot
data.table
default
deSolve
dplyr
ggplot2
groHMM
pbmcapply
rjson
```
## Description
### Input
- a matrix computed by deepTools for pSer2 ChIP-Seq after CTDP1 degradation ```data/C1D.matrix.gz```
- a configuration file in ```data/config.txt```
### Usage
To reproduce the figure in the paper, follow the instructions in ```scripts/recycle.html```
