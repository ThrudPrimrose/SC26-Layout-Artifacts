#pragma once
#include "types.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>

inline bool file_exists(const char* filename) {
    struct stat st;
    return stat(filename, &st) == 0;
}

inline long file_size(const char* filename) {
    struct stat st;
    if (stat(filename, &st) != 0) return -1;
    return st.st_size;
}

inline bool data_load_int(int& val, const char* filename) {
    if (!file_exists(filename)) {
        printf("File does not exist: %s\n", filename);
        return false;
    }
    FILE* f = fopen(filename, "rb");
    if (!f) { printf("Error opening file %s\n", filename); return false; }
    fread(&val, sizeof(int), 1, f);
    fclose(f);
    return true;
}

inline bool data_load_real(DP& val, const char* filename) {
    if (!file_exists(filename)) {
        printf("File does not exist: %s\n", filename);
        return false;
    }
    FILE* f = fopen(filename, "rb");
    if (!f) { printf("Error opening file %s\n", filename); return false; }
    fread(&val, sizeof(DP), 1, f);
    fclose(f);
    return true;
}

inline bool data_load_int_array(int* data, int count, const char* filename) {
    if (!file_exists(filename)) {
        printf("File does not exist: %s\n", filename);
        return false;
    }
    long fsize = file_size(filename);
    if (count != fsize / 4) {
        printf("Error of different file_size for %s. Array size: %d, in-file: %ld\n", filename, count, fsize);
        return false;
    }
    FILE* f = fopen(filename, "rb");
    if (!f) { printf("Error opening file %s\n", filename); return false; }
    fread(data, sizeof(int), count, f);
    fclose(f);
    return true;
}

inline bool data_load_real_array(DP* data, int count, const char* filename) {
    if (!file_exists(filename)) {
        printf("File does not exist: %s\n", filename);
        return false;
    }
    long fsize = file_size(filename);
    if (count != fsize / 8) {
        printf("Error of different file_size for %s. Array size: %d, in-file: %ld\n", filename, count, fsize);
        return false;
    }
    FILE* f = fopen(filename, "rb");
    if (!f) { printf("Error opening file %s\n", filename); return false; }
    fread(data, sizeof(DP), count, f);
    fclose(f);
    return true;
}

inline bool data_load_cmplx_array(Complex_DP* data, int count, const char* filename) {
    if (!file_exists(filename)) {
        printf("File does not exist: %s\n", filename);
        return false;
    }
    long fsize = file_size(filename);
    if (count != fsize / 16) {
        printf("Error of different file_size for %s. Array size: %d, in-file: %ld\n", filename, count, fsize);
        return false;
    }
    FILE* f = fopen(filename, "rb");
    if (!f) { printf("Error opening file %s\n", filename); return false; }
    fread(data, sizeof(Complex_DP), count, f);
    fclose(f);
    return true;
}
