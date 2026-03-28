COMPILER?=gnu

ifeq ($(COMPILER),gnu)

    #OPTIMIZE=-DUSE_RESTRICT -DUSE_ASSUME_ALIGNED -DUSE_MM_MALLOC -DUSE_FLOAT -DUSE_TILES -DTILE_SIZE=8

    CC=gcc
    override CFLAGS += -Wall -Wpedantic -O3 -march=native 
    override CFLAGS += -fopenmp -DUSE_FIRST_TOUCH -DUSE_ASSUME_ALIGNED -DUSE_SIMD
    #-vec-optimized 
    override LDFLAGS=-lm

else ifeq ($(COMPILER),intel)

    CC=icc
    override CFLAGS  += -I.. -fno-inline-functions
    override CFLAGS  += -qopt-report=5
    override CFLAGS  += -qopt-report-phase=vec -qopt-report-phase=openmp
    override CFLAGS  += -qopt-report-routine=kick,drift,accel
    override CFLAGS  += -qopt-report-file=stdout
    override LDFLAGS +=

endif
UNITYFLAGS=-DUNITY_OUTPUT_COLOR -DUNITY_INCLUDE_DOUBLE

vpath %.c ..

SOA_SRC=unity.c nb-soa-tests.c nb-main.c nb-soa-support.c nb-soa-alloc.c nb-soa-kda.c
AOS_SRC=unity.c nb-aos-tests.c nb-main.c nb-aos-support.c nb-aos-alloc.c nb-aos-kda.c

nbody-%: $(SOA_SRC)
	$(CC) -I.. -DUSE_SOA $(CFLAGS) $(LDFLAGS) $(UNITYFLAGS) $? -o $@ > $@.optrpt -lm

nbody-soa-%: $(SOA_SRC)
	$(CC) -I.. -DUSE_SOA $(CFLAGS) $(LDFLAGS) $(UNITYFLAGS) $? -o $@ > $@.optrpt -lm

nbody-aos-%: $(AOS_SRC)
	$(CC) -I.. -DUSE_AOS $(CFLAGS) $(LDFLAGS) $(UNITYFLAGS) $? -o $@ > $@.optrpt -lm

default:
	@echo "Usage: make nbody-NAME"
	@echo "Usage: make nbody-soa-NAME"
	@echo "Usage: make nbody-aos-NAME"

clean: .PHONY
	rm -f *.optrpt
.PHONY:

.FORCE:
