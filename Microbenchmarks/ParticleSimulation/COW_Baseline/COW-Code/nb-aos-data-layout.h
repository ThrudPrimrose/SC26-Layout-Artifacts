#ifndef NB_AOS_DATA_LAYOUT_H
#define NB_AOS_DATA_LAYOUT_H

struct vec3
{
    real x,y,z;
};

struct particle
{
    real m;
    struct vec3 r,v,a;
};

#define PARTICLE_DATA struct particle *p

#endif
