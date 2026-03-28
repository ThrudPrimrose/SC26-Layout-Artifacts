#ifndef NB_MAIN_H
#define NB_MAIN_H

#include <argp.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <linux/limits.h>
#include <stdlib.h>

// ---------------------------------------------------------------------------
//
//      Define compiler type
//
// ---------------------------------------------------------------------------

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#define HAVE_INTEL
#else
#define HAVE_GNU
#endif

// ---------------------------------------------------------------------------
//
//      Select double vs float
//
// ---------------------------------------------------------------------------
#if !defined(USE_FLOAT) && !defined(USE_DOUBLE)
#define USE_DOUBLE
#endif

#if defined(USE_FLOAT) && defined(USE_DOUBLE)
#error "Must choose one of USE_DOUBLE or USE_FLOAT"
#endif

#ifdef USE_FLOAT
typedef float real;
#else
typedef double real;
#endif

#ifdef USE_FLOAT
enum { NB_ALIGN = 128 };
#else
enum { NB_ALIGN = 256 };
#endif

#ifdef USE_FLOAT
#   define   INV(a)    (1.0f / (a))
#   define  SQRT(a)    sqrtf(a)
#   define RSQRT(a)    INV(SQRT(a))
#else
#   define   INV(a)    (1.0 / (a))
#   define  SQRT(a)    sqrt(a)
#   define RSQRT(a)    INV(SQRT(a))
#endif

// ---------------------------------------------------------------------------
//
//      Compiler specific functions or keywords
//
// ---------------------------------------------------------------------------

#if defined(HAVE_INTEL)

#define RESTRICT                restrict
#define NB_ASSUME_ALIGNED(var)  __assume_aligned((var), NB_ALIGN)
#define NB_MALLOC(var,sz)       (var) = _mm_malloc((sz),NB_ALIGN)

#elif defined(HAVE_GNU)

#define RESTRICT                __restrict__
#define NB_ASSUME_ALIGNED(var)  (var) = __builtin_assume_aligned((var), NB_ALIGN)
#define NB_MALLOC(var,sz)       posix_memalign((void**)&(var),NB_ALIGN,(sz))

#endif

// ---------------------------------------------------------------------------
//
//      Undefine things that weren't selected by the user
//
// ---------------------------------------------------------------------------

#ifndef USE_RESTRICT
#undef  RESTRICT
#define RESTRICT
#endif

#ifndef USE_ASSUME_ALIGNED
#undef  NB_ASSUME_ALIGNED
#define NB_ASSUME_ALIGNED(x) 
#endif

#ifndef USE_MM_MALLOC
#undef  NB_MALLOC
#define NB_MALLOC(var,sz) (var)=malloc(sz)
#endif

// ---------------------------------------------------------------------------
//
//      Cache tile size
//
// ---------------------------------------------------------------------------

#if !defined(TILE_SIZE) || TILE_SIZE <= 0 
#   undef  USE_TILES
#   undef  TILE_SIZE
#   define  TILE_SIZE 1
#else
#   define USE_TILES
#endif

// ---------------------------------------------------------------------------
//
//      Logging macro with verbosity level.
//
// ---------------------------------------------------------------------------

#define log(lvl,simopts,...) do { if ((simopts)->verbose >= lvl) fprintf(simopts->log_fp, __VA_ARGS__); } while (0)

// ---------------------------------------------------------------------------
//
//      Flags so the tests can fail and not abort, where a normal run would.
//
// ---------------------------------------------------------------------------

enum flag_type
{
    NO_FLAGS        =   0,
    EXIT_ON_ERROR   =   1,
    RETURN_ON_ERROR =   2,
};

// ---------------------------------------------------------------------------
//
//      Initial conditions
//
// ---------------------------------------------------------------------------

enum ic_type 
{ 
    IC_NONE,
    IC_CIRCULAR,
    IC_RANDOM,
};
enum { IC_COUNT = 2 };

// ---------------------------------------------------------------------------
//
//      Timing measurements
//
// ---------------------------------------------------------------------------

#define TIME_ACC(t) do { (t).dt += deltat(&(t).t0,&(t).t1); (t).n++; } while (0)
#define TIME(fn,t)  do { wallclock(&(t).t0); (fn); wallclock(&(t).t1); TIME_ACC((t)); } while (0)

struct sim_measure
{
    struct timeval t0, t1;
    double dt;
    int n;
};

struct nb_timing
{
    struct sim_measure sim;
    struct sim_measure drift, kick, accel;
    struct sim_measure io;
};

struct perf_measure
{ 
    struct timeval t0, t1;
    double dt;

    long nflops;

    long real_load;
    long real_store;
    long ptr_load;

    double arithmetic_intensity;    // [FLOP count / bytes moved]
    double FLOPS;
    double data_xfer;               // [GB]
    double data_rate;               // [GB/s]

}; 

// ---------------------------------------------------------------------------
//
//      Function specific measurements
//
// ---------------------------------------------------------------------------

struct nb_perf
{
    int n;
    struct perf_measure kick, drift, accel;
};

// ---------------------------------------------------------------------------
//
//      Raw command line arguments
//
// ---------------------------------------------------------------------------

struct nb_commandline_opts
{
    int run_sim;
    char *run_example;
    int run_tests;
    int run_perf;
    int verbose;
    int N;
    char *ic;
    int steps;
    real dt;
    int compute_energy;
    int status_freq_nsteps;
    int output_freq_nsteps;
    char *output_tag;
    char *status_output_file;
    char *log_output_file;
    char *perf_output_file;
    int overwrite;
    real rsoft;
    char *rng_seed;
};

// ---------------------------------------------------------------------------
//
//      Parsed and corrected options that are used.
//
// ---------------------------------------------------------------------------

struct nb_simulation_opts
{
    int run_sim;
    int run_tests;
    int run_perf;

    int verbose;
    int npart;
    int nranks;
    int nthreads;
    int nsteps;
    enum ic_type ic;
    real dt;
    int compute_E;
    int status_freq_nsteps;
    int output_freq_nsteps;
    char output_tag[PATH_MAX+1];
    char status_path[PATH_MAX+1];     FILE *status_fp;
    char log_path[PATH_MAX+1];        FILE *log_fp;
    char perf_path[PATH_MAX+1];       FILE *perf_fp;
    int overwrite;
    int first_output;
    real rsoft;
    long int rng_seed;
};

#if defined(USE_SOA)
#   include "nb-soa-data-layout.h"
#elif defined(USE_AOS)
#   include "nb-aos-data-layout.h"
#endif

// ---------------------------------------------------------------------------
//
//      Simulation data
//
// ---------------------------------------------------------------------------

struct energy
{
    double K,P;
};

struct nb_simulation_data
{
    PARTICLE_DATA;

    struct energy E, E0;
};

#include "nb-support.h"
#if defined(USE_SOA)
#   include "nb-soa.h"
#elif defined(USE_AOS)
#   include "nb-aos.h"
#endif

// ---------------------------------------------------------------------------
//
//      Function prototypes.
//
// ---------------------------------------------------------------------------

error_t parse_commandline_options(int argc, char *argv[], int flags, struct nb_commandline_opts *cmdlopts);
void prepare_commandline_argp(struct argp *argp);
error_t parse_opt (int key, char *arg, struct argp_state *state);

int 
check_option_validity(struct nb_commandline_opts *cmdlopts, int flags);

int 
assure_option_coherency(struct nb_commandline_opts *cmdlopts, struct nb_simulation_opts *simopts, int flags);

int 
prepare_initial_conditions(struct nb_simulation_opts *simopts,
                           struct nb_simulation_data *sim);
 

void
inject_example_configurations(struct nb_commandline_opts *cmdlopts);

void
first_touch(void *s, size_t size);

int 
imax(const int a, const int b);

int 
imin(const int a, const int b);

int 
time_to_report_status(struct nb_simulation_opts *simopts, int istep);

int 
time_to_output_frame(struct nb_simulation_opts *simopts, int istep);

void 
display_status(FILE *fp, struct nb_simulation_data *sim, int istep);

void 
display_startup_report(FILE *fp, struct nb_simulation_opts *simopts);

void 
wallclock(struct timeval *t);

double 
deltat(struct timeval *t1, struct timeval *t2);

void report_perf(FILE *fp, struct nb_perf *perf);

void 
report_sim_stats
(FILE *fp, struct nb_simulation_opts *simopts, struct nb_timing *t);

void compute_perf_stats(struct nb_perf *perf);

void compute_perf_measurements(int n, struct perf_measure *pm);

void 
write_perf_stats
(FILE *fp, struct nb_simulation_opts *simopts, struct nb_perf *perf);

int
nranks_from_environment();

int
nthreads_from_environment();

#endif
