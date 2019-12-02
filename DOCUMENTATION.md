
## Overview

The repository `adapta-fit-sads` contains the R code used in to fit Species Abundance Distributions (SADs) using an adaptive cutoff on the abundances. When using this code, please acknowledge the authors by citing  [Ser-Giacomi et al. (2018)](#references).



## Index
This documentation is organized as follows:

- [Input file format](#input-file-format)
- [Parameters setting](#parameters-setting)
	- [Command line arguments](#command-line-arguments)
	- [Parameters inside the code](#parameters-inside-the-code)	
- [ Running the code](#running-the-code)
- [Outputs format](#outputs-format)
- [References](#references)



## Input file format
The input file `tara_data` is a `.csv` table. To each row is associated a precise species (or OTU) and to each column a single sampling station. The entry *i-j* of the table corresponds to the abundance of the species *i* in the station *j*.



## Parameters setting
Here we list all the parameters that should be set to run the code. The place in the code where parameters should be inserted is denoted by `#INSERT`.


#### Command line arguments
- Two command line arguments should be given while executing the code. They are useful when launching several jobs on the same dataset:
```R
## how many samples are taken in each job 
job_number_sample = as.integer(args[1])

## counter for the initial sample number for a job
start_sample = as.integer(args[2])
``` 



#### Parameters inside the code

- The path of the input file:
```R
## reading tara data
tara_data = read.csv(#INSERT);
```

- The number of boots used by the Monte Carlo to calculate the p-value:
```R
## number of boots for p-value calculation
boot_reps = #INSERT
```

- The minimum value allowed for the cutoff on the abundances i.e. the value from which the loop over increasing `x_max` starts:
```R
    ## setting starting value for x_max
    start_x_max = as.integer(#INSERT)
```

- The maximum value allowed for the cutoff on the abundances i.e. the value at which the loop over increasing `x_max` ends:
```R
    ## setting end value for x_max
    end_x_max = as.integer(#INSERT)
```

- The size of the step among two consecutive values of `x_max` in the loop. Note that it could be also a nonlinear function of `x_max`:
```R
        ## stepping
        x_max = start_x_max + #INSERT
        counter = counter + 1
        vec_x_max[counter] = x_max
```

- The parameters needed by the generalized simulated annealing algorithm, please refer to [Xiang et al. (2013)](#references) for the details:
```R
        estimated_parameters = GenSA(c(#INSERT), L_func)$par
```


## Running the code
The code can be executed in the command line as follows:
``` bash
Rscript adapta-fit-sads.R <job_number_sample> <start_sample>
```
where `<job_number_sample` and `<start_sample>` are two integers.


## Outputs format



## References

[[Ser-Giacomi et al. 2018]](https://www.nature.com/articles/s41559-018-0587-2) Ser-Giacomi, E., Zinger, L., Malviya, S., De Vargas, C., Karsenti, E., Bowler, C., & De Monte, S. (2018). Ubiquitous abundance distribution of non-dominant plankton across the global ocean. *Nature Ecology & Evolution*, 2(8), 1243.

 [[Xiang et al. (2013)]](https://www.researchgate.net/profile/Sylvain_Gubian/publication/265058751_Generalized_Simulated_Annealing_for_Global_Optimization_The_GenSA_Package/links/53fdca890cf22f21c2f8470e/Generalized-Simulated-Annealing-for-Global-Optimization-The-GenSA-Package.pdf) Xiang, Y., Gubian, S., Suomela, B., & Hoeng, J. (2013). Generalized Simulated Annealing for Global Optimization: The GenSA Package. *R Journal*, 5(1).







