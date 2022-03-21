## Pogrom files for models in *"Investment slumps during financial crises: The real effects of credit supply"* from Alexandros Fakos, Plutarchos Sakellaris, and Tiago Tavares

This repository contains program codes and auxiliary files necessary to replicate the results presented in the paper *"Investment slumps during financial crises: The real effects of credit supply"* from Alexandros Fakos, Plutarchos Sakellaris, and Tiago Tavares. Most programs are written in Fortran and were tested for GFortran (GCC version 11.2), nvfortran (NVIDIA Compilers and Tools version 20.7), and MPICH (version 3.3.2). Additionally, some figures are generated with gnuplot (version 5.4). Different makefiles automatize compilation of the different programs. Additionally, this repository also includes some a few shell scripts used to send run MPI programs on an AWS server.

#### Replication of results in Table 1 and 2 - Model estimates
Set up first auxiliary functions and a fixed random sequence of random shocks by running `main.f90` (`make main`) and then `stats.f90` (`make stats_show`). The *simulated minimum distance* estimates summarized in tables 1 and 2 of the paper can be found by running either `smm_estimation.f90` (`make smm_estimation` for a serial version) or `smm_estimation_par.f90` (`make smm_estimation_par` with `Cluster_script.sh`, `Sync_in_server.sh`, and `run.sh` for a parallel version using MPI). Once the program finds the estimates, you can update the parameters used in `main.f90` and `main_abc.f90` to replicate figure 3 by running the programs generated by `make main` and `make main_abc` followed by `make gplot` and `make gplot_abc`. The programs generate the text files `Policyfnz_baseline.txt` and `Policyfnz_abc.txt` that generate figure 3 with the gnuplot file `figure_scripts/Leverage_z_baseline_abc.gp`.

#### Replication of results in Table 4 - Estimating the credit-supply series
The estimated series of credit-supply can be generated using the program `greek_crisis_exp.f90` (`make greek_crisis_exp`). The program contains two target leverage series (raw and corrected for provisions). The relevant part of the code for the two targets has to be uncommented in separate runs of the program in order to generate the two time series presented the table 4.

#### Replication of results in Table 5 - Decomposing the investment contraction
Once the program generated from `greek_crisis_exp.f90` is run, one should have output files `plot_data/Distributions_cal_estimation_PUTebkacc.txt` and `plot_data/Distributions_cal_estimation_PUTebkprov.txt`. To generate the results in table 5, update the code option inside `greek_crisis_exp.f90` (set `iflag_creditshock = 1` in line 96 to `iflag_creditshock = 0`) and run again the program. Now, an additional output file `plot_data/Distributions_cal_estimation_PUT.txt` is generated. The decompositions are based in the three output files contained in the folder `plot_data/` by comparing average and aggregate investment rates.

#### Replication of Figure 5 - Exporters and Non-exporters subsamples


[description later - plenty of options]

#### Program `smm_estimation.f90` and `smm_estimation_par.f90`

[description later - a few options and aux files need to edit both simultaneously]

#### Gnuplot scripts: inside folder `./figure_scripts`






