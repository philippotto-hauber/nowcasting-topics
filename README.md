# nowcasting-topics
This repo implements nowcasts of GDP based on topics extracted from daily newspapers

## Structure of the repo

- **model**: code for the EM-algorithm to estimate a mixed frequency factor model with daily, weekly, monthly and quarterly data. References are [Modugno (2011)](https://ideas.repec.org/p/ecb/ecbwps/20111324.html) and [Banbura et al. (2011)](https://ideas.repec.org/p/red/sed012/555.html). 

- **data**: notebook that illustrates how the raw topics are transformed prior to estimation; various plots. The subfolder **topics** contains the raw topics time series and their description. 
