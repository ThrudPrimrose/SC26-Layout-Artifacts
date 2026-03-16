
#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <cstring>
#include <filesystem>
#include "hdf5.h"
#include "codegen/cloudsc_py.h"
#include "sensitivity.h"

double* load_h5_array_double(hid_t file_id, std::string name, size_t max_size) {
    double* ptr = new double[max_size];
    for(size_t i=0; i<max_size; ++i) ptr[i] = 0.0;

    if (file_id >= 0 && H5Lexists(file_id, name.c_str(), H5P_DEFAULT) > 0) {
        hid_t dataset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
        hid_t space_id = H5Dget_space(dataset_id);
        hssize_t num_elements = H5Sget_simple_extent_npoints(space_id);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr);
        H5Sclose(space_id);
        H5Dclose(dataset_id);
    } else {
        for(size_t i=0; i<max_size; ++i) ptr[i] = 0.5;
    }
    return ptr;
}

int* load_h5_array_int(hid_t file_id, std::string name, size_t max_size) {
    int* ptr = new int[max_size];
    for(size_t i=0; i<max_size; ++i) ptr[i] = 0;

    if (file_id >= 0 && H5Lexists(file_id, name.c_str(), H5P_DEFAULT) > 0) {
        hid_t dataset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
        hid_t space_id = H5Dget_space(dataset_id);
        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr);
        H5Sclose(space_id);
        H5Dclose(dataset_id);
    } else {
        for(size_t i=0; i<max_size; ++i) ptr[i] = 1;
    }
    return ptr;
}

void save_h5_array_double(hid_t file_id, std::string name, double* ptr, size_t size, int klev, int klon, int nclv) {
    hsize_t dims[3];
    int rank = 0;
    if (size == (size_t)klon * (klev + 1)) { dims[0] = klev + 1; dims[1] = klon; rank = 2; }
    else if (size == (size_t)klon * klev) { dims[0] = klev; dims[1] = klon; rank = 2; }
    else if (size == (size_t)klon * klev * nclv) { dims[0] = nclv; dims[1] = klev; dims[2] = klon; rank = 3; }
    else if (size == (size_t)klon) { dims[0] = klon; rank = 1; }
    else { dims[0] = size; rank = 1; }

    hid_t space_id = H5Screate_simple(rank, dims, NULL);
    hid_t dataset_id = H5Dcreate2(file_id, name.c_str(), H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr);
    H5Dclose(dataset_id);
    H5Sclose(space_id);
}

void save_h5_array_int(hid_t file_id, std::string name, int* ptr, size_t size, int klon) {
    hsize_t dims[1] = { (hsize_t)size };
    hid_t space_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate2(file_id, name.c_str(), H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr);
    H5Dclose(dataset_id);
    H5Sclose(space_id);
}

int main(int argc, char** argv) {
    int num_steps = 1;
    if (argc > 1 && std::string(argv[1]) != "--save" && std::string(argv[1]) != "--sensitivity") num_steps = std::stoi(argv[1]);
    bool save_output = false;
    bool sensitivity_mode = false;
    double sens_eps = 1e-5;
    for(int i=1; i<argc; ++i) {
        if(std::string(argv[i]) == "--save") save_output = true;
        if(std::string(argv[i]) == "--sensitivity") sensitivity_mode = true;
        if(std::string(argv[i]) == "--sens-eps" && i+1 < argc) sens_eps = std::stod(argv[++i]);
    }

    int klon = 100, klev = 137, nclv = 5;
    int ncldqi = 2, ncldql = 1, ncldqr = 3, ncldqs = 4, ncldqv = 5;
    int kidia = 1, kfdia = 100;
    double ptsphy = 3600;

    std::string input_file = "input.h5";
    hid_t file_id = -1;
    if (std::filesystem::exists(input_file)) {
        std::cout << "Reading input data from: " << input_file << std::endl;
        file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    } else {
        std::cerr << "Warning: " << input_file << " not found, using dummy data." << std::endl;
    }

    std::cout << "Running CloudSC for " << num_steps << " steps." << std::endl;

    // Constants
    double ydcst_rcpd = 1004.70885783306744;
    double ydcst_rd = 287.059673666590697;
    double ydcst_retv = 0.607766731611463662;
    double ydcst_rg = 9.80664999999999942;
    double ydcst_rlmlt = 333700.0;
    double ydcst_rlstt = 2834500.0;
    double ydcst_rlvtt = 2500800.0;
    double ydcst_rtt = 273.160000000000025;
    double ydcst_rv = 461.524993308387877;
    double ydthf_r2es = 380.160870344284717;
    double ydthf_r3ies = 22.5869999999999997;
    double ydthf_r3les = 17.5019999999999989;
    double ydthf_r4ies = -0.699999999999999956;
    double ydthf_r4les = 32.1899999999999977;
    double ydthf_r5alscp = 17451123.2533625774;
    double ydthf_r5alvcp = 10497584.6816953104;
    double ydthf_r5ies = 6185.67582000000039;
    double ydthf_r5les = 4217.45694000000003;
    double ydthf_ralfdcp = 332.13601870666929;
    double ydthf_ralsdcp = 2821.21529824409345;
    double ydthf_ralvdcp = 2489.07927953742455;
    double ydthf_rkoop1 = 2.58300000000000018;
    double ydthf_rkoop2 = 0.00481160000000000008;
    double ydthf_rtice = 250.160000000000025;
    double ydthf_rticecu = 250.160000000000025;
    double ydthf_rtwat = 273.160000000000025;
    double ydthf_rtwat_rtice_r = 0.0434782608695652162;
    double ydthf_rtwat_rticecu_r = 0.0434782608695652162;
    int yrecldp_laericeauto = false;
    int yrecldp_laericesed = false;
    int yrecldp_laerliqautolsp = false;
    int yrecldp_laerliqcoll = false;
    int yrecldp_ncldtop = 15;
    int yrecldp_nssopt = 1;
    double yrecldp_ramid = 0.800000000000000044;
    double yrecldp_ramin = 1.00000000000000002e-08;
    double yrecldp_rccn = 125.0;
    double yrecldp_rcl_apb1 = 714000000000.0;
    double yrecldp_rcl_apb2 = 116000000.0;
    double yrecldp_rcl_apb3 = 241.599999999999994;
    double yrecldp_rcl_cdenom1 = 557000000000.0;
    double yrecldp_rcl_cdenom2 = 103000000.0;
    double yrecldp_rcl_cdenom3 = 204.0;
    double yrecldp_rcl_const1i = 3.62318801151369977e-06;
    double yrecldp_rcl_const1r = 1.38230076757950893;
    double yrecldp_rcl_const1s = 3.62318801151369977e-06;
    double yrecldp_rcl_const2i = 6283185.30717958603;
    double yrecldp_rcl_const2r = 2143.2299120517614;
    double yrecldp_rcl_const2s = 6283185.30717958603;
    double yrecldp_rcl_const3i = 596.999847583599831;
    double yrecldp_rcl_const3r = 0.634999999999999787;
    double yrecldp_rcl_const3s = 596.999847583599831;
    double yrecldp_rcl_const4i = 0.66666666666666663;
    double yrecldp_rcl_const4r = -0.200000000000000178;
    double yrecldp_rcl_const4s = 0.66666666666666663;
    double yrecldp_rcl_const5i = 0.921166666666666689;
    double yrecldp_rcl_const5r = 8685252.96508213319;
    double yrecldp_rcl_const5s = 0.921166666666666689;
    double yrecldp_rcl_const6i = 1.00000009489611852;
    double yrecldp_rcl_const6r = -4.79999999999999982;
    double yrecldp_rcl_const6s = 1.00000009489611852;
    double yrecldp_rcl_const7s = 90363515.7635107338;
    double yrecldp_rcl_const8s = 1.17566666666666664;
    double yrecldp_rcl_fac1 = 4146.9027898470631;
    double yrecldp_rcl_fac2 = 0.55555555555555558;
    double yrecldp_rcl_fzrab = -0.660000000000000031;
    double yrecldp_rcl_ka273 = 0.0240000000000000005;
    double yrecldp_rcl_kk_cloud_num_land = 300.0;
    double yrecldp_rcl_kk_cloud_num_sea = 50.0;
    double yrecldp_rcl_kkaac = 67.0;
    double yrecldp_rcl_kkaau = 1350.0;
    double yrecldp_rcl_kkbac = 1.14999999999999991;
    double yrecldp_rcl_kkbaun = -1.79000000000000004;
    double yrecldp_rcl_kkbauq = 2.4700000000000002;
    double yrecldp_rclcrit_land = 0.000550000000000000033;
    double yrecldp_rclcrit_sea = 0.000250000000000000005;
    double yrecldp_rcldiff = 3.00000000000000008e-06;
    double yrecldp_rcldiff_convi = 7.0;
    double yrecldp_rcldtopcf = 0.0100000000000000002;
    double yrecldp_rcovpmin = 0.100000000000000006;
    double yrecldp_rdensref = 1.0;
    double yrecldp_rdepliqrefdepth = 500.0;
    double yrecldp_rdepliqrefrate = 0.100000000000000006;
    double yrecldp_riceinit = 9.9999999999999998e-13;
    double yrecldp_rkconv = 0.000166666666666666661;
    double yrecldp_rkooptau = 10800.0;
    double yrecldp_rlcritsnow = 3.00000000000000008e-05;
    double yrecldp_rlmin = 1.00000000000000002e-08;
    double yrecldp_rnice = 0.0269999999999999997;
    double yrecldp_rpecons = 5.54725619859993018e-05;
    double yrecldp_rprc1 = 100.0;
    double yrecldp_rprecrhmax = 0.699999999999999956;
    double yrecldp_rsnowlin1 = 0.00100000000000000002;
    double yrecldp_rsnowlin2 = 0.0299999999999999989;
    double yrecldp_rtaumel = 7200.0;
    double yrecldp_rthomo = 235.160000000000025;
    double yrecldp_rvice = 0.130000000000000004;
    double yrecldp_rvrain = 4.0;
    double yrecldp_rvrfactor = 0.00508999999999999991;
    double yrecldp_rvsnow = 1.0;

    // Array allocations & loading
    int* ktype = load_h5_array_int(file_id, "KTYPE", klon * (klev + 1) * nclv);
    int* ldcum = load_h5_array_int(file_id, "LDCUM", klon * (klev + 1) * nclv);
    double* pa = load_h5_array_double(file_id, "PA", klon * (klev + 1) * nclv);
    double* pap = load_h5_array_double(file_id, "PAP", klon * (klev + 1) * nclv);
    double* paph = load_h5_array_double(file_id, "PAPH", klon * (klev + 1) * nclv);
    double* pccn = load_h5_array_double(file_id, "PCCN", klon * (klev + 1) * nclv);
    double* pclv = load_h5_array_double(file_id, "PCLV", klon * (klev + 1) * nclv);
    double* pcovptot = load_h5_array_double(file_id, "PCOVPTOT", klon * (klev + 1) * nclv);
    double* pdyna = load_h5_array_double(file_id, "PDYNA", klon * (klev + 1) * nclv);
    double* pdyni = load_h5_array_double(file_id, "PDYNI", klon * (klev + 1) * nclv);
    double* pdynl = load_h5_array_double(file_id, "PDYNL", klon * (klev + 1) * nclv);
    double* pfcqlng = load_h5_array_double(file_id, "PFCQLNG", klon * (klev + 1) * nclv);
    double* pfcqnng = load_h5_array_double(file_id, "PFCQNNG", klon * (klev + 1) * nclv);
    double* pfcqrng = load_h5_array_double(file_id, "PFCQRNG", klon * (klev + 1) * nclv);
    double* pfcqsng = load_h5_array_double(file_id, "PFCQSNG", klon * (klev + 1) * nclv);
    double* pfhpsl = load_h5_array_double(file_id, "PFHPSL", klon * (klev + 1) * nclv);
    double* pfhpsn = load_h5_array_double(file_id, "PFHPSN", klon * (klev + 1) * nclv);
    double* pfplsl = load_h5_array_double(file_id, "PFPLSL", klon * (klev + 1) * nclv);
    double* pfplsn = load_h5_array_double(file_id, "PFPLSN", klon * (klev + 1) * nclv);
    double* pfsqif = load_h5_array_double(file_id, "PFSQIF", klon * (klev + 1) * nclv);
    double* pfsqitur = load_h5_array_double(file_id, "PFSQITUR", klon * (klev + 1) * nclv);
    double* pfsqlf = load_h5_array_double(file_id, "PFSQLF", klon * (klev + 1) * nclv);
    double* pfsqltur = load_h5_array_double(file_id, "PFSQLTUR", klon * (klev + 1) * nclv);
    double* pfsqrf = load_h5_array_double(file_id, "PFSQRF", klon * (klev + 1) * nclv);
    double* pfsqsf = load_h5_array_double(file_id, "PFSQSF", klon * (klev + 1) * nclv);
    double* phrlw = load_h5_array_double(file_id, "PHRLW", klon * (klev + 1) * nclv);
    double* phrsw = load_h5_array_double(file_id, "PHRSW", klon * (klev + 1) * nclv);
    double* picrit_aer = load_h5_array_double(file_id, "PICRIT_AER", klon * (klev + 1) * nclv);
    double* plcrit_aer = load_h5_array_double(file_id, "PLCRIT_AER", klon * (klev + 1) * nclv);
    double* plsm = load_h5_array_double(file_id, "PLSM", klon * (klev + 1) * nclv);
    double* plu = load_h5_array_double(file_id, "PLU", klon * (klev + 1) * nclv);
    double* plude = load_h5_array_double(file_id, "PLUDE", klon * (klev + 1) * nclv);
    double* pmfd = load_h5_array_double(file_id, "PMFD", klon * (klev + 1) * nclv);
    double* pmfu = load_h5_array_double(file_id, "PMFU", klon * (klev + 1) * nclv);
    double* pnice = load_h5_array_double(file_id, "PNICE", klon * (klev + 1) * nclv);
    double* pq = load_h5_array_double(file_id, "PQ", klon * (klev + 1) * nclv);
    double* prainfrac_toprfz = load_h5_array_double(file_id, "PRAINFRAC_TOPRFZ", klon * (klev + 1) * nclv);
    double* pre_ice = load_h5_array_double(file_id, "PRE_ICE", klon * (klev + 1) * nclv);
    double* psnde = load_h5_array_double(file_id, "PSNDE", klon * (klev + 1) * nclv);
    double* psupsat = load_h5_array_double(file_id, "PSUPSAT", klon * (klev + 1) * nclv);
    double* pt = load_h5_array_double(file_id, "PT", klon * (klev + 1) * nclv);
    double* pvervel = load_h5_array_double(file_id, "PVERVEL", klon * (klev + 1) * nclv);
    double* pvfa = load_h5_array_double(file_id, "PVFA", klon * (klev + 1) * nclv);
    double* pvfi = load_h5_array_double(file_id, "PVFI", klon * (klev + 1) * nclv);
    double* pvfl = load_h5_array_double(file_id, "PVFL", klon * (klev + 1) * nclv);
    double* tendency_loc_a = load_h5_array_double(file_id, "TENDENCY_LOC_A", klon * (klev + 1) * nclv);
    double* tendency_loc_cld = load_h5_array_double(file_id, "TENDENCY_LOC_CLD", klon * (klev + 1) * nclv);
    double* tendency_loc_q = load_h5_array_double(file_id, "TENDENCY_LOC_Q", klon * (klev + 1) * nclv);
    double* tendency_loc_t = load_h5_array_double(file_id, "TENDENCY_LOC_T", klon * (klev + 1) * nclv);
    double* tendency_tmp_a = load_h5_array_double(file_id, "TENDENCY_TMP_A", klon * (klev + 1) * nclv);
    double* tendency_tmp_cld = load_h5_array_double(file_id, "TENDENCY_TMP_CLD", klon * (klev + 1) * nclv);
    double* tendency_tmp_q = load_h5_array_double(file_id, "TENDENCY_TMP_Q", klon * (klev + 1) * nclv);
    double* tendency_tmp_t = load_h5_array_double(file_id, "TENDENCY_TMP_T", klon * (klev + 1) * nclv);

    if (file_id >= 0) H5Fclose(file_id);

    std::cout << "Initializing DaCe CloudSC CPU program (klon=" << klon << ", klev=" << klev << ")..." << std::endl;
    auto *handle = __dace_init_cloudsc_py(klev, klon, ncldqi, ncldql, ncldqr, ncldqs, ncldqv, nclv);

    std::cout << "Executing CloudSC..." << std::endl;
    double total_time = 0;
    for (int s=0; s<num_steps; ++s) {
        auto start = std::chrono::high_resolution_clock::now();
        __program_cloudsc_py(handle, ktype, ldcum, pa, pap, paph, pccn, pclv, pcovptot, pdyna, pdyni, pdynl, pfcqlng, pfcqnng, pfcqrng, pfcqsng, pfhpsl, pfhpsn, pfplsl, pfplsn, pfsqif, pfsqitur, pfsqlf, pfsqltur, pfsqrf, pfsqsf, phrlw, phrsw, picrit_aer, plcrit_aer, plsm, plu, plude, pmfd, pmfu, pnice, pq, prainfrac_toprfz, pre_ice, psnde, psupsat, pt, pvervel, pvfa, pvfi, pvfl, tendency_loc_a, tendency_loc_cld, tendency_loc_q, tendency_loc_t, tendency_tmp_a, tendency_tmp_cld, tendency_tmp_q, tendency_tmp_t, kfdia, kidia, klev, klon, ncldqi, ncldql, ncldqr, ncldqs, ncldqv, nclv, ptsphy, ydcst_rcpd, ydcst_rd, ydcst_retv, ydcst_rg, ydcst_rlmlt, ydcst_rlstt, ydcst_rlvtt, ydcst_rtt, ydcst_rv, ydthf_r2es, ydthf_r3ies, ydthf_r3les, ydthf_r4ies, ydthf_r4les, ydthf_r5alscp, ydthf_r5alvcp, ydthf_r5ies, ydthf_r5les, ydthf_ralfdcp, ydthf_ralsdcp, ydthf_ralvdcp, ydthf_rkoop1, ydthf_rkoop2, ydthf_rtice, ydthf_rticecu, ydthf_rtwat, ydthf_rtwat_rtice_r, ydthf_rtwat_rticecu_r, yrecldp_laericeauto, yrecldp_laericesed, yrecldp_laerliqautolsp, yrecldp_laerliqcoll, yrecldp_ncldtop, yrecldp_nssopt, yrecldp_ramid, yrecldp_ramin, yrecldp_rccn, yrecldp_rcl_apb1, yrecldp_rcl_apb2, yrecldp_rcl_apb3, yrecldp_rcl_cdenom1, yrecldp_rcl_cdenom2, yrecldp_rcl_cdenom3, yrecldp_rcl_const1i, yrecldp_rcl_const1r, yrecldp_rcl_const1s, yrecldp_rcl_const2i, yrecldp_rcl_const2r, yrecldp_rcl_const2s, yrecldp_rcl_const3i, yrecldp_rcl_const3r, yrecldp_rcl_const3s, yrecldp_rcl_const4i, yrecldp_rcl_const4r, yrecldp_rcl_const4s, yrecldp_rcl_const5i, yrecldp_rcl_const5r, yrecldp_rcl_const5s, yrecldp_rcl_const6i, yrecldp_rcl_const6r, yrecldp_rcl_const6s, yrecldp_rcl_const7s, yrecldp_rcl_const8s, yrecldp_rcl_fac1, yrecldp_rcl_fac2, yrecldp_rcl_fzrab, yrecldp_rcl_ka273, yrecldp_rcl_kk_cloud_num_land, yrecldp_rcl_kk_cloud_num_sea, yrecldp_rcl_kkaac, yrecldp_rcl_kkaau, yrecldp_rcl_kkbac, yrecldp_rcl_kkbaun, yrecldp_rcl_kkbauq, yrecldp_rclcrit_land, yrecldp_rclcrit_sea, yrecldp_rcldiff, yrecldp_rcldiff_convi, yrecldp_rcldtopcf, yrecldp_rcovpmin, yrecldp_rdensref, yrecldp_rdepliqrefdepth, yrecldp_rdepliqrefrate, yrecldp_riceinit, yrecldp_rkconv, yrecldp_rkooptau, yrecldp_rlcritsnow, yrecldp_rlmin, yrecldp_rnice, yrecldp_rpecons, yrecldp_rprc1, yrecldp_rprecrhmax, yrecldp_rsnowlin1, yrecldp_rsnowlin2, yrecldp_rtaumel, yrecldp_rthomo, yrecldp_rvice, yrecldp_rvrain, yrecldp_rvrfactor, yrecldp_rvsnow);
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double>(end - start).count();

        if (save_output) {
            std::filesystem::create_directories("outputs_cpp");
            std::string out_name = "outputs_cpp/cpp_output_step_" + std::to_string(s) + ".h5";
            std::filesystem::path full_path = std::filesystem::absolute(out_name);
            hid_t out_id = H5Fcreate(out_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            save_h5_array_int(out_id, "KTYPE", ktype, klon, klon);
            save_h5_array_int(out_id, "LDCUM", ldcum, klon, klon);
            save_h5_array_double(out_id, "PA", pa, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PAP", pap, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PAPH", paph, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PCCN", pccn, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PCLV", pclv, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PCOVPTOT", pcovptot, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PDYNA", pdyna, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PDYNI", pdyni, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PDYNL", pdynl, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PFCQLNG", pfcqlng, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFCQNNG", pfcqnng, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFCQRNG", pfcqrng, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFCQSNG", pfcqsng, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFHPSL", pfhpsl, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFHPSN", pfhpsn, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFPLSL", pfplsl, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFPLSN", pfplsn, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFSQIF", pfsqif, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFSQITUR", pfsqitur, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFSQLF", pfsqlf, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFSQLTUR", pfsqltur, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFSQRF", pfsqrf, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PFSQSF", pfsqsf, klon * (klev + 1), klev, klon, nclv);
            save_h5_array_double(out_id, "PHRLW", phrlw, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PHRSW", phrsw, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PICRIT_AER", picrit_aer, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PLCRIT_AER", plcrit_aer, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PLSM", plsm, klon, klev, klon, nclv);
            save_h5_array_double(out_id, "PLU", plu, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PLUDE", plude, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PMFD", pmfd, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PMFU", pmfu, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PNICE", pnice, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PQ", pq, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PRAINFRAC_TOPRFZ", prainfrac_toprfz, klon, klev, klon, nclv);
            save_h5_array_double(out_id, "PRE_ICE", pre_ice, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PSNDE", psnde, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PSUPSAT", psupsat, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PT", pt, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PVERVEL", pvervel, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PVFA", pvfa, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PVFI", pvfi, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "PVFL", pvfl, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "TENDENCY_LOC_A", tendency_loc_a, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "TENDENCY_LOC_CLD", tendency_loc_cld, klon * klev * nclv, klev, klon, nclv);
            save_h5_array_double(out_id, "TENDENCY_LOC_Q", tendency_loc_q, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "TENDENCY_LOC_T", tendency_loc_t, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "TENDENCY_TMP_A", tendency_tmp_a, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "TENDENCY_TMP_CLD", tendency_tmp_cld, klon * klev * nclv, klev, klon, nclv);
            save_h5_array_double(out_id, "TENDENCY_TMP_Q", tendency_tmp_q, klon * klev, klev, klon, nclv);
            save_h5_array_double(out_id, "TENDENCY_TMP_T", tendency_tmp_t, klon * klev, klev, klon, nclv);
            H5Fclose(out_id);
            std::cout << "Saved results to: " << full_path << std::endl;
        }
    }

    std::cout << "Total execution time: " << total_time << "s" << std::endl;
    std::cout << "Average time per step: " << total_time / num_steps << "s" << std::endl;

    if (sensitivity_mode) {
        // Re-load fresh inputs (kernel may have mutated them)
        if (std::filesystem::exists(input_file))
            file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        auto reload = [&](double* dst, const char* name, size_t n) {
            if (file_id >= 0 && H5Lexists(file_id, name, H5P_DEFAULT) > 0) {
                hid_t ds = H5Dopen2(file_id, name, H5P_DEFAULT);
                H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst);
                H5Dclose(ds);
            }
        };
        size_t n2 = klon * klev, n3 = klon * (klev + 1), n4 = klon * klev * nclv;
        reload(pt, "PT", n2); reload(pq, "PQ", n2);
        reload(tendency_tmp_t, "TENDENCY_TMP_T", n2); reload(tendency_tmp_q, "TENDENCY_TMP_Q", n2);
        reload(tendency_tmp_a, "TENDENCY_TMP_A", n2); reload(tendency_tmp_cld, "TENDENCY_TMP_CLD", n4);
        reload(tendency_loc_t, "TENDENCY_LOC_T", n2); reload(tendency_loc_q, "TENDENCY_LOC_Q", n2);
        reload(tendency_loc_a, "TENDENCY_LOC_A", n2); reload(tendency_loc_cld, "TENDENCY_LOC_CLD", n4);
        reload(pvfa, "PVFA", n2); reload(pvfl, "PVFL", n2); reload(pvfi, "PVFI", n2);
        reload(pdyna, "PDYNA", n2); reload(pdynl, "PDYNL", n2); reload(pdyni, "PDYNI", n2);
        reload(phrsw, "PHRSW", n2); reload(phrlw, "PHRLW", n2);
        reload(pvervel, "PVERVEL", n2); reload(pap, "PAP", n2); reload(paph, "PAPH", n2);
        reload(plsm, "PLSM", n2); reload(plu, "PLU", n2); reload(plude, "PLUDE", n2);
        reload(psnde, "PSNDE", n2); reload(pmfu, "PMFU", n2); reload(pmfd, "PMFD", n2);
        reload(pa, "PA", n2); reload(pclv, "PCLV", n2); reload(psupsat, "PSUPSAT", n2);
        reload(plcrit_aer, "PLCRIT_AER", n2); reload(picrit_aer, "PICRIT_AER", n2);
        reload(pre_ice, "PRE_ICE", n2); reload(pccn, "PCCN", n2); reload(pnice, "PNICE", n2);
        reload(pcovptot, "PCOVPTOT", n2); reload(prainfrac_toprfz, "PRAINFRAC_TOPRFZ", n2);
        reload(pfcqlng, "PFCQLNG", n3); reload(pfcqnng, "PFCQNNG", n3);
        reload(pfcqrng, "PFCQRNG", n3); reload(pfcqsng, "PFCQSNG", n3);
        reload(pfsqlf, "PFSQLF", n3); reload(pfsqrf, "PFSQRF", n3);
        reload(pfsqif, "PFSQIF", n3); reload(pfsqsf, "PFSQSF", n3);
        reload(pfsqltur, "PFSQLTUR", n3); reload(pfsqitur, "PFSQITUR", n3);
        reload(pfplsl, "PFPLSL", n3); reload(pfplsn, "PFPLSN", n3);
        reload(pfhpsl, "PFHPSL", n3); reload(pfhpsn, "PFHPSN", n3);
        if (file_id >= 0) H5Fclose(file_id);

        auto save_copy = [](double* src, size_t n) {
            double* c = new double[n]; std::memcpy(c, src, n * sizeof(double)); return c;
        };
        // Backup ALL arrays (kernel may mutate any of them in-place)
        double* pt_bk = save_copy(pt, n2);
        double* pq_bk = save_copy(pq, n2);
        double* tendency_tmp_t_bk = save_copy(tendency_tmp_t, n2);
        double* tendency_tmp_q_bk = save_copy(tendency_tmp_q, n2);
        double* tendency_tmp_a_bk = save_copy(tendency_tmp_a, n2);
        double* tendency_tmp_cld_bk = save_copy(tendency_tmp_cld, n4);
        double* tendency_loc_t_bk = save_copy(tendency_loc_t, n2);
        double* tendency_loc_q_bk = save_copy(tendency_loc_q, n2);
        double* tendency_loc_a_bk = save_copy(tendency_loc_a, n2);
        double* tendency_loc_cld_bk = save_copy(tendency_loc_cld, n4);
        double* pvfa_bk = save_copy(pvfa, n2);
        double* pvfl_bk = save_copy(pvfl, n2);
        double* pvfi_bk = save_copy(pvfi, n2);
        double* pdyna_bk = save_copy(pdyna, n2);
        double* pdynl_bk = save_copy(pdynl, n2);
        double* pdyni_bk = save_copy(pdyni, n2);
        double* phrsw_bk = save_copy(phrsw, n2);
        double* phrlw_bk = save_copy(phrlw, n2);
        double* pvervel_bk = save_copy(pvervel, n2);
        double* pap_bk = save_copy(pap, n2);
        double* paph_bk = save_copy(paph, n2);
        double* plu_bk = save_copy(plu, n2);
        double* plude_bk = save_copy(plude, n2);
        double* psnde_bk = save_copy(psnde, n2);
        double* pmfu_bk = save_copy(pmfu, n2);
        double* pmfd_bk = save_copy(pmfd, n2);
        double* pa_bk = save_copy(pa, n2);
        double* pclv_bk = save_copy(pclv, n4);
        double* psupsat_bk = save_copy(psupsat, n2);
        double* plcrit_aer_bk = save_copy(plcrit_aer, n2);
        double* picrit_aer_bk = save_copy(picrit_aer, n2);
        double* pre_ice_bk = save_copy(pre_ice, n2);
        double* pccn_bk = save_copy(pccn, n2);
        double* pnice_bk = save_copy(pnice, n2);
        double* pcovptot_bk = save_copy(pcovptot, n2);
        double* pfsqlf_bk = save_copy(pfsqlf, n3);
        double* pfsqrf_bk = save_copy(pfsqrf, n3);
        double* pfsqif_bk = save_copy(pfsqif, n3);
        double* pfsqsf_bk = save_copy(pfsqsf, n3);
        double* pfcqlng_bk = save_copy(pfcqlng, n3);
        double* pfcqnng_bk = save_copy(pfcqnng, n3);
        double* pfcqrng_bk = save_copy(pfcqrng, n3);
        double* pfcqsng_bk = save_copy(pfcqsng, n3);
        double* pfhpsn_bk = save_copy(pfhpsn, n3);
        double* pfplsn_bk = save_copy(pfplsn, n3);
        double* pfsqltur_bk = save_copy(pfsqltur, n3);
        double* pfsqitur_bk = save_copy(pfsqitur, n3);
        double* pfhpsl_bk = save_copy(pfhpsl, n3);
        double* pfplsl_bk = save_copy(pfplsl, n3);
        double* prainfrac_toprfz_bk = save_copy(prainfrac_toprfz, n2);
        double* plsm_bk = save_copy(plsm, klon);

        auto restore = [&]() {
            std::memcpy(pt, pt_bk, n2 * sizeof(double));
            std::memcpy(pq, pq_bk, n2 * sizeof(double));
            std::memcpy(tendency_tmp_t, tendency_tmp_t_bk, n2 * sizeof(double));
            std::memcpy(tendency_tmp_q, tendency_tmp_q_bk, n2 * sizeof(double));
            std::memcpy(tendency_tmp_a, tendency_tmp_a_bk, n2 * sizeof(double));
            std::memcpy(tendency_tmp_cld, tendency_tmp_cld_bk, n4 * sizeof(double));
            std::memcpy(tendency_loc_t, tendency_loc_t_bk, n2 * sizeof(double));
            std::memcpy(tendency_loc_q, tendency_loc_q_bk, n2 * sizeof(double));
            std::memcpy(tendency_loc_a, tendency_loc_a_bk, n2 * sizeof(double));
            std::memcpy(tendency_loc_cld, tendency_loc_cld_bk, n4 * sizeof(double));
            std::memcpy(pvfa, pvfa_bk, n2 * sizeof(double));
            std::memcpy(pvfl, pvfl_bk, n2 * sizeof(double));
            std::memcpy(pvfi, pvfi_bk, n2 * sizeof(double));
            std::memcpy(pdyna, pdyna_bk, n2 * sizeof(double));
            std::memcpy(pdynl, pdynl_bk, n2 * sizeof(double));
            std::memcpy(pdyni, pdyni_bk, n2 * sizeof(double));
            std::memcpy(phrsw, phrsw_bk, n2 * sizeof(double));
            std::memcpy(phrlw, phrlw_bk, n2 * sizeof(double));
            std::memcpy(pvervel, pvervel_bk, n2 * sizeof(double));
            std::memcpy(pap, pap_bk, n2 * sizeof(double));
            std::memcpy(paph, paph_bk, n2 * sizeof(double));
            std::memcpy(plu, plu_bk, n2 * sizeof(double));
            std::memcpy(plude, plude_bk, n2 * sizeof(double));
            std::memcpy(psnde, psnde_bk, n2 * sizeof(double));
            std::memcpy(pmfu, pmfu_bk, n2 * sizeof(double));
            std::memcpy(pmfd, pmfd_bk, n2 * sizeof(double));
            std::memcpy(pa, pa_bk, n2 * sizeof(double));
            std::memcpy(pclv, pclv_bk, n4 * sizeof(double));
            std::memcpy(psupsat, psupsat_bk, n2 * sizeof(double));
            std::memcpy(plcrit_aer, plcrit_aer_bk, n2 * sizeof(double));
            std::memcpy(picrit_aer, picrit_aer_bk, n2 * sizeof(double));
            std::memcpy(pre_ice, pre_ice_bk, n2 * sizeof(double));
            std::memcpy(pccn, pccn_bk, n2 * sizeof(double));
            std::memcpy(pnice, pnice_bk, n2 * sizeof(double));
            std::memcpy(pcovptot, pcovptot_bk, n2 * sizeof(double));
            std::memcpy(pfsqlf, pfsqlf_bk, n3 * sizeof(double));
            std::memcpy(pfsqrf, pfsqrf_bk, n3 * sizeof(double));
            std::memcpy(pfsqif, pfsqif_bk, n3 * sizeof(double));
            std::memcpy(pfsqsf, pfsqsf_bk, n3 * sizeof(double));
            std::memcpy(pfcqlng, pfcqlng_bk, n3 * sizeof(double));
            std::memcpy(pfcqnng, pfcqnng_bk, n3 * sizeof(double));
            std::memcpy(pfcqrng, pfcqrng_bk, n3 * sizeof(double));
            std::memcpy(pfcqsng, pfcqsng_bk, n3 * sizeof(double));
            std::memcpy(pfhpsn, pfhpsn_bk, n3 * sizeof(double));
            std::memcpy(pfplsn, pfplsn_bk, n3 * sizeof(double));
            std::memcpy(pfsqltur, pfsqltur_bk, n3 * sizeof(double));
            std::memcpy(pfsqitur, pfsqitur_bk, n3 * sizeof(double));
            std::memcpy(pfhpsl, pfhpsl_bk, n3 * sizeof(double));
            std::memcpy(pfplsl, pfplsl_bk, n3 * sizeof(double));
            std::memcpy(prainfrac_toprfz, prainfrac_toprfz_bk, n2 * sizeof(double));
            std::memcpy(plsm, plsm_bk, klon * sizeof(double));
        };

        auto execute = [&]() {
            __program_cloudsc_py(handle, ktype, ldcum, pa, pap, paph, pccn, pclv, pcovptot, pdyna, pdyni, pdynl, pfcqlng, pfcqnng, pfcqrng, pfcqsng, pfhpsl, pfhpsn, pfplsl, pfplsn, pfsqif, pfsqitur, pfsqlf, pfsqltur, pfsqrf, pfsqsf, phrlw, phrsw, picrit_aer, plcrit_aer, plsm, plu, plude, pmfd, pmfu, pnice, pq, prainfrac_toprfz, pre_ice, psnde, psupsat, pt, pvervel, pvfa, pvfi, pvfl, tendency_loc_a, tendency_loc_cld, tendency_loc_q, tendency_loc_t, tendency_tmp_a, tendency_tmp_cld, tendency_tmp_q, tendency_tmp_t, kfdia, kidia, klev, klon, ncldqi, ncldql, ncldqr, ncldqs, ncldqv, nclv, ptsphy, ydcst_rcpd, ydcst_rd, ydcst_retv, ydcst_rg, ydcst_rlmlt, ydcst_rlstt, ydcst_rlvtt, ydcst_rtt, ydcst_rv, ydthf_r2es, ydthf_r3ies, ydthf_r3les, ydthf_r4ies, ydthf_r4les, ydthf_r5alscp, ydthf_r5alvcp, ydthf_r5ies, ydthf_r5les, ydthf_ralfdcp, ydthf_ralsdcp, ydthf_ralvdcp, ydthf_rkoop1, ydthf_rkoop2, ydthf_rtice, ydthf_rticecu, ydthf_rtwat, ydthf_rtwat_rtice_r, ydthf_rtwat_rticecu_r, yrecldp_laericeauto, yrecldp_laericesed, yrecldp_laerliqautolsp, yrecldp_laerliqcoll, yrecldp_ncldtop, yrecldp_nssopt, yrecldp_ramid, yrecldp_ramin, yrecldp_rccn, yrecldp_rcl_apb1, yrecldp_rcl_apb2, yrecldp_rcl_apb3, yrecldp_rcl_cdenom1, yrecldp_rcl_cdenom2, yrecldp_rcl_cdenom3, yrecldp_rcl_const1i, yrecldp_rcl_const1r, yrecldp_rcl_const1s, yrecldp_rcl_const2i, yrecldp_rcl_const2r, yrecldp_rcl_const2s, yrecldp_rcl_const3i, yrecldp_rcl_const3r, yrecldp_rcl_const3s, yrecldp_rcl_const4i, yrecldp_rcl_const4r, yrecldp_rcl_const4s, yrecldp_rcl_const5i, yrecldp_rcl_const5r, yrecldp_rcl_const5s, yrecldp_rcl_const6i, yrecldp_rcl_const6r, yrecldp_rcl_const6s, yrecldp_rcl_const7s, yrecldp_rcl_const8s, yrecldp_rcl_fac1, yrecldp_rcl_fac2, yrecldp_rcl_fzrab, yrecldp_rcl_ka273, yrecldp_rcl_kk_cloud_num_land, yrecldp_rcl_kk_cloud_num_sea, yrecldp_rcl_kkaac, yrecldp_rcl_kkaau, yrecldp_rcl_kkbac, yrecldp_rcl_kkbaun, yrecldp_rcl_kkbauq, yrecldp_rclcrit_land, yrecldp_rclcrit_sea, yrecldp_rcldiff, yrecldp_rcldiff_convi, yrecldp_rcldtopcf, yrecldp_rcovpmin, yrecldp_rdensref, yrecldp_rdepliqrefdepth, yrecldp_rdepliqrefrate, yrecldp_riceinit, yrecldp_rkconv, yrecldp_rkooptau, yrecldp_rlcritsnow, yrecldp_rlmin, yrecldp_rnice, yrecldp_rpecons, yrecldp_rprc1, yrecldp_rprecrhmax, yrecldp_rsnowlin1, yrecldp_rsnowlin2, yrecldp_rtaumel, yrecldp_rthomo, yrecldp_rvice, yrecldp_rvrain, yrecldp_rvrfactor, yrecldp_rvsnow);
        };

        std::vector<sensitivity::Field> sens_inputs = {
            {pt, (int)n2, "pt"},
            {pq, (int)n2, "pq"},
            {tendency_tmp_t, (int)n2, "tendency_tmp_t"},
            {tendency_tmp_q, (int)n2, "tendency_tmp_q"},
            {tendency_tmp_a, (int)n2, "tendency_tmp_a"},
            {tendency_tmp_cld, (int)n4, "tendency_tmp_cld"},
            {pvfa, (int)n2, "pvfa"},
            {pvfl, (int)n2, "pvfl"},
            {pvfi, (int)n2, "pvfi"},
            {pdyna, (int)n2, "pdyna"},
            {pdynl, (int)n2, "pdynl"},
            {pdyni, (int)n2, "pdyni"},
            {phrsw, (int)n2, "phrsw"},
            {phrlw, (int)n2, "phrlw"},
            {pvervel, (int)n2, "pvervel"},
            {pap, (int)n2, "pap"},
            {paph, (int)n2, "paph"},
            {plsm, klon, "plsm"},
            {plu, (int)n2, "plu"},
            {plude, (int)n2, "plude"},
            {psnde, (int)n2, "psnde"},
            {pmfu, (int)n2, "pmfu"},
            {pmfd, (int)n2, "pmfd"},
            {pa, (int)n2, "pa"},
            {pclv, (int)n4, "pclv"},
            {psupsat, (int)n2, "psupsat"},
            {plcrit_aer, (int)n2, "plcrit_aer"},
            {picrit_aer, (int)n2, "picrit_aer"},
            {pre_ice, (int)n2, "pre_ice"},
            {pccn, (int)n2, "pccn"},
            {pnice, (int)n2, "pnice"},
        };

        std::vector<sensitivity::Field> sens_outputs = {
            {tendency_loc_t, (int)n2, "tendency_loc_t"},
            {tendency_loc_q, (int)n2, "tendency_loc_q"},
            {tendency_loc_a, (int)n2, "tendency_loc_a"},
            {tendency_loc_cld, (int)n4, "tendency_loc_cld"},
            {pcovptot, (int)n2, "pcovptot"},
            {pfsqlf, (int)n3, "pfsqlf"},
            {pfsqrf, (int)n3, "pfsqrf"},
            {pfsqif, (int)n3, "pfsqif"},
            {pfsqsf, (int)n3, "pfsqsf"},
            {pfcqlng, (int)n3, "pfcqlng"},
            {pfcqnng, (int)n3, "pfcqnng"},
            {pfcqrng, (int)n3, "pfcqrng"},
            {pfcqsng, (int)n3, "pfcqsng"},
            {pfhpsn, (int)n3, "pfhpsn"},
            {pfplsn, (int)n3, "pfplsn"},
            {pfsqltur, (int)n3, "pfsqltur"},
            {pfsqitur, (int)n3, "pfsqitur"},
            {pfhpsl, (int)n3, "pfhpsl"},
            {pfplsl, (int)n3, "pfplsl"},
            {prainfrac_toprfz, klon, "prainfrac_toprfz"},
        };

        std::cerr << "Running sensitivity analysis (eps=" << sens_eps << ")...\n";
        sensitivity::run(restore, execute, sens_inputs, sens_outputs, sens_eps);

        delete[] pt_bk; delete[] pq_bk;
        delete[] tendency_tmp_t_bk; delete[] tendency_tmp_q_bk;
        delete[] tendency_tmp_a_bk; delete[] tendency_tmp_cld_bk;
        delete[] tendency_loc_t_bk; delete[] tendency_loc_q_bk;
        delete[] tendency_loc_a_bk; delete[] tendency_loc_cld_bk;
        delete[] pvfa_bk; delete[] pvfl_bk; delete[] pvfi_bk;
        delete[] pdyna_bk; delete[] pdynl_bk; delete[] pdyni_bk;
        delete[] phrsw_bk; delete[] phrlw_bk; delete[] pvervel_bk;
        delete[] pap_bk; delete[] paph_bk;
        delete[] plu_bk; delete[] plude_bk; delete[] psnde_bk;
        delete[] pmfu_bk; delete[] pmfd_bk;
        delete[] pa_bk; delete[] pclv_bk; delete[] psupsat_bk;
        delete[] plcrit_aer_bk; delete[] picrit_aer_bk;
        delete[] pre_ice_bk; delete[] pccn_bk; delete[] pnice_bk;
        delete[] pcovptot_bk;
        delete[] pfsqlf_bk; delete[] pfsqrf_bk; delete[] pfsqif_bk; delete[] pfsqsf_bk;
        delete[] pfcqlng_bk; delete[] pfcqnng_bk; delete[] pfcqrng_bk; delete[] pfcqsng_bk;
        delete[] pfhpsn_bk; delete[] pfplsn_bk;
        delete[] pfsqltur_bk; delete[] pfsqitur_bk;
        delete[] pfhpsl_bk; delete[] pfplsl_bk;
        delete[] prainfrac_toprfz_bk; delete[] plsm_bk;
    }

    __dace_exit_cloudsc_py(handle);
    return 0;
}
