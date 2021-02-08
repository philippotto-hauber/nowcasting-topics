# nowcasting-topics
This repo implements nowcasts of GDP based on topics extracted from daily newspapers

## Structure of the repo

- **model**: code for the EM-algorithm to estimate a mixed frequency factor model with daily, weekly, monthly and quarterly data. References are [Modugno (2011)](https://ideas.repec.org/p/ecb/ecbwps/20111324.html) and [Banbura et al. (2011)](https://ideas.repec.org/p/red/sed012/555.html). For the time being, the code can mix daily, monthly and quarterly series. For the latter two it distinguishes between stock and flow variables. Weekly data can both be lower frequency series or the highest frequency. In case of the former, it is straight-forward to adapt the code in the same for monthly and quarterly data to allow for weekly flows and stocks. Note, however, that this has not been implemented yet and the current code implicitly assumes that all weekly series are stocks. Regarding weekly data as the highest frequency, the weights for monthly and quarterly flow variables would have to be adjusted. Similarly, the
<img src="https://render.githubusercontent.com/render/math?math=\Xi_m"> 
and 
<img src="https://render.githubusercontent.com/render/math?math=\Xi_q"> 
cumulator variables that indicate the start of a new period, would also have to be adjusted as they are currently formulated in days. A useful guide on how to deal with mixing quarterly and weekly data is [Aruoba et al. (2011)](https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/ads/real-time-measurement-of-business-conditions14.pdf?la=en&hash=8CB33CA37D2F88A57F3622F4060D2CD5).

- **data**: notebook that illustrates how the raw topics are transformed prior to estimation; various plots. The subfolder **topics** contains the raw topics time series and their description. 
