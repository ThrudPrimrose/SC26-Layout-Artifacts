NVCC = nvcc
NVCCFLAGS = -O3 -std=c++17 -Xcompiler -fopenmp

TARGET = run_addusxx_g

SRCS = main.cu usxx_kernels.cu usxx_kernels_soa.cu

all: $(TARGET)

$(TARGET): $(SRCS) types.h data_loading.h usxx_kernels.h
	$(NVCC) $(NVCCFLAGS) -o $@ $(SRCS)

clean:
	rm -f $(TARGET)
