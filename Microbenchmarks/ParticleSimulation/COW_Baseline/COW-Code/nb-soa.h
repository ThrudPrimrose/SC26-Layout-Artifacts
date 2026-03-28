#ifndef NB_SOA_SUPPORT_H
#define NB_SOA_SUPPORT_H

#include "nb.h"

void 
accel(int N, real *m, struct vec3 *r, struct vec3 *a, real rsoft);

void 
accel_tiled(int N, real *m, struct vec3 *r, struct vec3 *a, real rsoft);

void 
kick(real dt, int N, struct vec3 *v, struct vec3 *a);

void 
drift(real dt, int N, struct vec3 *r, struct vec3 *v);


void 
perf_drift(real dt, int N, struct vec3 *r, struct vec3 *v, struct perf_measure *pm);

void 
perf_kick(real dt, int N, struct vec3 *v, struct vec3 *a, struct perf_measure *pm);

void 
perf_accel(int N, real *m, struct vec3 *r, struct vec3 *a, real rsoft, struct perf_measure *pm);

void
compute_energies(int N, real *m, struct vec3 *r, struct vec3 *v, real rsoft, double *KE, double *PE);

double 
kinetic_energy(int N, real *m, struct vec3 *v);

double 
potential_energy(int N, real *m, struct vec3 *r, real rsoft);


#endif
