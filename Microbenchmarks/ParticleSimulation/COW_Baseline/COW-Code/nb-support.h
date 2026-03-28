#ifndef NB_SUPPORT_H
#define NB_SUPPORT_H

#include "nb.h"

int 
allocate_memory(struct nb_simulation_opts *simopts, 
                struct nb_simulation_data *sim);

int 
prepare_initial_conditions_circular(struct nb_simulation_opts *simopts, 
                                    struct nb_simulation_data *sim);
        
int 
prepare_initial_conditions_random(struct nb_simulation_opts *simopts, 
                                  struct nb_simulation_data *sim);

int
run_simulation(struct nb_simulation_opts *simopts,
               struct nb_simulation_data *sim,
               struct nb_timing *timing);

int
run_perf(struct nb_simulation_opts *simopts,
         struct nb_simulation_data *sim,
         struct nb_perf *perf);

void 
perform_requested_output(struct nb_simulation_opts *simopts,
                         struct nb_simulation_data *sim,
                         int istep);

void output_frame(struct nb_simulation_opts *simopts,
                  struct nb_simulation_data *sim,
                  int istep);

#endif
