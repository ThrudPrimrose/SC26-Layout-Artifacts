#ifndef NB_AOS_KDA_H
#define NB_AOS_KDA_H

void accel(int N, struct particle *p, real rsoft);

void kick(real dt, int N, struct particle *p);

void drift(real dt, int N, struct particle *p);

void 
perf_drift
(real dt, int N, struct particle *p, struct perf_measure *pm);

void 
perf_kick
(real dt, int N, struct particle *p, struct perf_measure *pm);

void 
perf_accel
(int N, struct particle *p, real rsoft, struct perf_measure *pm);

void
compute_energies(int N, struct particle *p, real rsoft, double *KE, double *PE);

double 
kinetic_energy
(int N, struct particle *p);

double 
potential_energy
(int N, struct particle *p, real rsoft);

#endif
