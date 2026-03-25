print("hipcc  -O3 -ffast-math -fPIC -Wall -Wextra  -DHIP_PLATFORM_AMD=1 -D__HIP_PLATFORM_AMD__=1  "
             "-Wno-unused-parameter -munsafe-fp-atomics --offload-arch=gfx942 "
             "-ffp-contract=fast -Wno-ignored-attributes -Wno-unused-result -std=c++17 transpose_hiptensor.cpp -o transpose_lib -lhiptensor")