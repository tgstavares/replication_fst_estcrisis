CCMPI  = mpif90.mpich

#CC =gfortran
#CFLAGS = -O3 -fno-range-check -ffree-line-length-none -w -fopenmp -cpp

CC = nvfortran
CFLAGS =-acc -Minfo=accel -O3 -fast -Minline -ta=tesla:managed -cpp

LFLAGS = #/usr/local/lib/liblapack.a /usr/local/lib/libblas.a /usr/local/lib/libnlopt.a
INCLUDES = -I/usr/local/include/

SRS1 = tools/prec.f90 \
       tools/normal.f90 \
       tools/linear_op.f90 \
       tools/kind_module.f90 \
       tools/bobyqa.f90 \
       tools/brent.f90 \
       tools/asa047.f90 \
       tools/ACRS.f90 \
       tools/sort_tiago.f90

SRS2 = #$(wildcard tools/interp/*.f)
SRS4 = $(wildcard tools/stats/*.f)

SRS5 = tools/interpf90/bspline_kinds_module.f90 \
       tools/interpf90/bspline_sub_module.f90 \
       tools/interpf90/bspline_oo_module.f90 \
       tools/interpf90/bspline_module.f90

SRS6 = $(wildcard tools/pspline/*.f90)

SRS = $(SRS1) $(SRS4) $(SRS2) $(SRS5) $(SRS6)

OBJ= $(patsubst tools/%.f90,%.o, \
     $(patsubst tools/%.f,%.o, \
     $(patsubst tools/interpf90/%.f90,%.o, \
     $(patsubst tools/pspline/%.f90,%.o, \
     $(patsubst tools/stats/%.f,%.o,$(SRS))))))

SA1 = aux_routines/parameters.f90 \
      aux_routines/globals.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values.f90 \
      aux_routines/equilibrium.f90 \
      aux_routines/simulation.f90

SA1_nonexporters = aux_routines/parameters_nonexporters.f90 \
      aux_routines/globals.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values.f90 \
      aux_routines/equilibrium.f90 \
      aux_routines/simulation.f90

SA1_exporters = aux_routines/parameters_exporters.f90 \
      aux_routines/globals.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values.f90 \
      aux_routines/equilibrium.f90 \
      aux_routines/simulation.f90

SA_ebc = aux_routines/parameters_ebc.f90 \
      aux_routines/globals.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values.f90 \
      aux_routines/equilibrium.f90 \
      aux_routines/simulation.f90

SA_defval = aux_routines/parameters_defval.f90 \
      aux_routines/globals_defval.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_defval.f90 \
      aux_routines/equilibrium_defval.f90 \
      aux_routines/simulation_defval.f90

SA_defval_nonexporters = aux_routines/parameters_defval_nonexporters.f90 \
      aux_routines/globals_defval.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_defval.f90 \
      aux_routines/equilibrium_defval.f90 \
      aux_routines/simulation_defval.f90

SA_defcf_nonexporters = aux_routines/parameters_defcf_nonexporters.f90 \
      aux_routines/globals_defcf.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_defcf.f90 \
      aux_routines/equilibrium_defcf.f90 \
      aux_routines/simulation_defcf.f90

SA_defcf_exporters = aux_routines/parameters_defcf_exporters.f90 \
      aux_routines/globals_defcf.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_defcf.f90 \
      aux_routines/equilibrium_defcf.f90 \
      aux_routines/simulation_defcf.f90

SA_defval_exporters = aux_routines/parameters_defval_exporters.f90 \
      aux_routines/globals_defval.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_defval.f90 \
      aux_routines/equilibrium_defval.f90 \
      aux_routines/simulation_defval.f90

SA_defcf = aux_routines/parameters_defcf.f90 \
      aux_routines/globals_defcf.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_defcf.f90 \
      aux_routines/equilibrium_defcf.f90 \
      aux_routines/simulation_defcf.f90

SA_kutz = aux_routines/parameters_kutz.f90 \
      aux_routines/globals_kutz.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_kutz.f90 \
      aux_routines/equilibrium_kutz.f90 \
      aux_routines/simulation_kutz.f90

SA_entry = aux_routines/parameters_entry.f90 \
      aux_routines/globals_entry.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_entry.f90 \
      aux_routines/equilibrium_entry.f90 \
      aux_routines/simulation_entry.f90

SA_dis = aux_routines/parameters_dis.f90 \
      aux_routines/globals_dis.f90 \
      aux_routines/globals_montecarlo.f90 \
      aux_routines/process.f90 \
      aux_routines/process.f90 \
      aux_routines/values_dis.f90 \
      aux_routines/equilibrium_dis.f90 \
      aux_routines/simulation_dis.f90

$(OBJ): $(SRS)
	$(CC) $(CFLAGS) $(SRS) -c

main: $(OBJ) $(SA1) main.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_nonexporters: $(OBJ) $(SA1_nonexporters) main_nonexporters.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_exporters: $(OBJ) $(SA1_exporters) main_exporters.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_abc: $(OBJ) $(SA1) main_abc.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_ebc: $(OBJ) $(SA_ebc) main_ebc.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_defval: $(OBJ) $(SA_defval) main_defval.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_defval_nonexporters: $(OBJ) $(SA_defval_nonexporters) main_defval_nonexporters.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_defcf_nonexporters: $(OBJ) $(SA_defcf_nonexporters) main_defcf_nonexporters.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_defcf_exporters: $(OBJ) $(SA_defcf_exporters) main_defcf_exporters.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_defval_exporters: $(OBJ) $(SA_defval_exporters) main_defval_exporters.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_defcf: $(OBJ) $(SA_defcf) main_defcf.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_kutz: $(OBJ) $(SA_kutz) main_kutz.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_entry: $(OBJ) $(SA_entry) main_entry.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

main_dis: $(OBJ) $(SA_dis) main_dis.f90
	$(CC) $^ -o main $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show: $(OBJ) $(SA1) stats.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_nonexporters: $(OBJ) $(SA1_nonexporters) stats_nonexporters.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_exporters: $(OBJ) $(SA1_exporters) stats_exporters.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_abc: $(OBJ) $(SA1) stats.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_ebc: $(OBJ) $(SA_ebc) stats.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_defval: $(OBJ) $(SA_defval) stats_defval.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_defval_nonexporters: $(OBJ) $(SA_defval_nonexporters) stats_defval_nonexporters.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_defcf_nonexporters: $(OBJ) $(SA_defcf_nonexporters) stats_defcf_nonexporters.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_defcf_exporters: $(OBJ) $(SA_defcf_exporters) stats_defcf_exporters.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_defval_exporters: $(OBJ) $(SA_defval_exporters) stats_defval_exporters.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_defcf: $(OBJ) $(SA_defcf) stats_defcf.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_kutz: $(OBJ) $(SA_kutz) stats_kutz.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_entry: $(OBJ) $(SA_entry) stats_entry.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

stats_show_dis: $(OBJ) $(SA_dis) stats_dis.f90
	$(CC) $^ -o stats_show $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp: $(OBJ) $(SA1) greek_crisis_exp.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_utilization: $(OBJ) $(SA1) greek_crisis_exp_utilization.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_policies: $(OBJ) $(SA1) greek_crisis_exp_policies.f90
	$(CC) $^ -o greek_crisis_policies $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_nonexporters: $(OBJ) $(SA1_nonexporters) greek_crisis_exp_nonexporters.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_exporters: $(OBJ) $(SA1_exporters) greek_crisis_exp_exporters.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_abc: $(OBJ) $(SA1) greek_crisis_exp.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_ebc: $(OBJ) $(SA_ebc) greek_crisis_exp.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_defval: $(OBJ) $(SA_defval) greek_crisis_exp_defval.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_defval_nonexporters: $(OBJ) $(SA_defval_nonexporters) greek_crisis_exp_defval_nonexporters.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_defcf_nonexporters: $(OBJ) $(SA_defcf_nonexporters) greek_crisis_exp_defcf_nonexporters.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_defcf_exporters: $(OBJ) $(SA_defcf_exporters) greek_crisis_exp_defcf_exporters.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_defval_exporters: $(OBJ) $(SA_defval_exporters) greek_crisis_exp_defval_exporters.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_defcf: $(OBJ) $(SA_defcf) greek_crisis_exp_defcf.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_kutz: $(OBJ) $(SA_kutz) greek_crisis_exp_kutz.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_entry: $(OBJ) $(SA_entry) greek_crisis_exp_entry.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

greek_crisis_exp_dis: $(OBJ) $(SA_dis) greek_crisis_exp_dis.f90
	$(CC) $^ -o greek_crisis_exp $(CFLAGS) $(LFLAGS) $(INCLUDES)

gplot: $(OBJ) $(SA1) gplot.f90
	$(CC) $^ -o gplot $(CFLAGS) $(LFLAGS) $(INCLUDES)

gplot_abc: $(OBJ) $(SA1) gplot_abc.f90
	$(CC) $^ -o gplot $(CFLAGS) $(LFLAGS) $(INCLUDES)

gplot_defval: $(OBJ) $(SA_defval) gplot_defval.f90
	$(CC) $^ -o gplot $(CFLAGS) $(LFLAGS) $(INCLUDES)

gplot_defcf: $(OBJ) $(SA_defcf) gplot_defcf.f90
	$(CC) $^ -o gplot $(CFLAGS) $(LFLAGS) $(INCLUDES)

smm_estimation: $(OBJ) $(SA1) smm_estimation.f90
	$(CC) $^ -o smm_estimation $(CFLAGS) $(LFLAGS) $(INCLUDES)

smm_estimation_nonexporters: $(OBJ) $(SA1_nonexporters) smm_estimation_nonexporters.f90
	$(CC) $^ -o smm_estimation $(CFLAGS) $(LFLAGS) $(INCLUDES)

smm_estimation_exporters: $(OBJ) $(SA1_exporters) smm_estimation_exporters.f90
	$(CC) $^ -o smm_estimation $(CFLAGS) $(LFLAGS) $(INCLUDES)

smm_estimation_par: $(OBJ) $(SA1) smm_estimation_par.f90
	$(CCMPI) $^ -o smm_estimation_par $(CFLAGS) $(LFLAGS) $(INCLUDES)

smm_estimation_par_nonexporters: $(OBJ) $(SA1_nonexporters) smm_estimation_par_nonexporters.f90
	$(CCMPI) $^ -o smm_estimation_par $(CFLAGS) $(LFLAGS) $(INCLUDES)

smm_estimation_par_exporters: $(OBJ) $(SA1_exporters) smm_estimation_par_exporters.f90
	$(CCMPI) $^ -o smm_estimation_par $(CFLAGS) $(LFLAGS) $(INCLUDES)

clean:
	rm -f *.o *.mod *~ tools/*~ scripts/*~ aux_routines/*~ figure_scripts/*~ \
	fort.* *# main main_abc main_ebc main_defval \
	stats_show stats_show_defval \
	gplot greek_crisis_exp smm_estimation smm_estimation_par greek_crisis_policies




