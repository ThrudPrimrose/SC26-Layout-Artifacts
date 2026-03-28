#ifndef NB_SOA_DATA_LAYOUT_H
#define NB_SOA_DATA_LAYOUT_H

struct vec3
{
    union 
    {
        struct { real *x,*y,*z; };
        struct { real *w[3];  };
    };
};

#define PARTICLE_DATA      \
struct                     \
{                          \
    real *m;               \
    struct vec3 r,v,a;     \
}

#endif
