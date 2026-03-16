/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "hash.h"

struct cloudsc_py_state_t {

};

void __program_cloudsc_py_internal(cloudsc_py_state_t*__state, int * __restrict__ ktype, int * __restrict__ ldcum, double * __restrict__ pa, double * __restrict__ pap, double * __restrict__ paph, double * __restrict__ pccn, double * __restrict__ pclv, double * __restrict__ pcovptot, double * __restrict__ pdyna, double * __restrict__ pdyni, double * __restrict__ pdynl, double * __restrict__ pfcqlng, double * __restrict__ pfcqnng, double * __restrict__ pfcqrng, double * __restrict__ pfcqsng, double * __restrict__ pfhpsl, double * __restrict__ pfhpsn, double * __restrict__ pfplsl, double * __restrict__ pfplsn, double * __restrict__ pfsqif, double * __restrict__ pfsqitur, double * __restrict__ pfsqlf, double * __restrict__ pfsqltur, double * __restrict__ pfsqrf, double * __restrict__ pfsqsf, double * __restrict__ phrlw, double * __restrict__ phrsw, double * __restrict__ picrit_aer, double * __restrict__ plcrit_aer, double * __restrict__ plsm, double * __restrict__ plu, double * __restrict__ plude, double * __restrict__ pmfd, double * __restrict__ pmfu, double * __restrict__ pnice, double * __restrict__ pq, double * __restrict__ prainfrac_toprfz, double * __restrict__ pre_ice, double * __restrict__ psnde, double * __restrict__ psupsat, double * __restrict__ pt, double * __restrict__ pvervel, double * __restrict__ pvfa, double * __restrict__ pvfi, double * __restrict__ pvfl, double * __restrict__ tendency_loc_a, double * __restrict__ tendency_loc_cld, double * __restrict__ tendency_loc_q, double * __restrict__ tendency_loc_t, double * __restrict__ tendency_tmp_a, double * __restrict__ tendency_tmp_cld, double * __restrict__ tendency_tmp_q, double * __restrict__ tendency_tmp_t, int kfdia, int kidia, int klev, int klon, double ptsphy, double ydcst_rcpd, double ydcst_rd, double ydcst_retv, double ydcst_rg, double ydcst_rlmlt, double ydcst_rlstt, double ydcst_rlvtt, double ydcst_rtt, double ydcst_rv, double ydthf_r2es, double ydthf_r3ies, double ydthf_r3les, double ydthf_r4ies, double ydthf_r4les, double ydthf_r5alscp, double ydthf_r5alvcp, double ydthf_r5ies, double ydthf_r5les, double ydthf_ralfdcp, double ydthf_ralsdcp, double ydthf_ralvdcp, double ydthf_rkoop1, double ydthf_rkoop2, double ydthf_rtice, double ydthf_rticecu, double ydthf_rtwat, double ydthf_rtwat_rtice_r, double ydthf_rtwat_rticecu_r, int yrecldp_laericeauto, int yrecldp_laericesed, int yrecldp_laerliqautolsp, int yrecldp_laerliqcoll, int yrecldp_ncldtop, int yrecldp_nssopt, double yrecldp_ramid, double yrecldp_ramin, double yrecldp_rccn, double yrecldp_rcl_apb1, double yrecldp_rcl_apb2, double yrecldp_rcl_apb3, double yrecldp_rcl_cdenom1, double yrecldp_rcl_cdenom2, double yrecldp_rcl_cdenom3, double yrecldp_rcl_const1i, double yrecldp_rcl_const1r, double yrecldp_rcl_const1s, double yrecldp_rcl_const2i, double yrecldp_rcl_const2r, double yrecldp_rcl_const2s, double yrecldp_rcl_const3i, double yrecldp_rcl_const3r, double yrecldp_rcl_const3s, double yrecldp_rcl_const4i, double yrecldp_rcl_const4r, double yrecldp_rcl_const4s, double yrecldp_rcl_const5i, double yrecldp_rcl_const5r, double yrecldp_rcl_const5s, double yrecldp_rcl_const6i, double yrecldp_rcl_const6r, double yrecldp_rcl_const6s, double yrecldp_rcl_const7s, double yrecldp_rcl_const8s, double yrecldp_rcl_fac1, double yrecldp_rcl_fac2, double yrecldp_rcl_fzrab, double yrecldp_rcl_ka273, double yrecldp_rcl_kk_cloud_num_land, double yrecldp_rcl_kk_cloud_num_sea, double yrecldp_rcl_kkaac, double yrecldp_rcl_kkaau, double yrecldp_rcl_kkbac, double yrecldp_rcl_kkbaun, double yrecldp_rcl_kkbauq, double yrecldp_rclcrit_land, double yrecldp_rclcrit_sea, double yrecldp_rcldiff, double yrecldp_rcldiff_convi, double yrecldp_rcldtopcf, double yrecldp_rcovpmin, double yrecldp_rdensref, double yrecldp_rdepliqrefdepth, double yrecldp_rdepliqrefrate, double yrecldp_riceinit, double yrecldp_rkconv, double yrecldp_rkooptau, double yrecldp_rlcritsnow, double yrecldp_rlmin, double yrecldp_rnice, double yrecldp_rpecons, double yrecldp_rprc1, double yrecldp_rprecrhmax, double yrecldp_rsnowlin1, double yrecldp_rsnowlin2, double yrecldp_rtaumel, double yrecldp_rthomo, double yrecldp_rvice, double yrecldp_rvrain, double yrecldp_rvrfactor, double yrecldp_rvsnow)
{
    double *zlcond1;
    zlcond1 = new double DACE_ALIGN(64)[klon];
    double *zlcond2;
    zlcond2 = new double DACE_ALIGN(64)[klon];
    double *zrainaut;
    zrainaut = new double DACE_ALIGN(64)[klon];
    double *zsnowaut;
    zsnowaut = new double DACE_ALIGN(64)[klon];
    double *zliqcld;
    zliqcld = new double DACE_ALIGN(64)[klon];
    double *zicecld;
    zicecld = new double DACE_ALIGN(64)[klon];
    double *zfokoop;
    zfokoop = new double DACE_ALIGN(64)[klon];
    double *zlicld;
    zlicld = new double DACE_ALIGN(64)[klon];
    double *zlfinalsum;
    zlfinalsum = new double DACE_ALIGN(64)[klon];
    double *zdqs;
    zdqs = new double DACE_ALIGN(64)[klon];
    double *ztold;
    ztold = new double DACE_ALIGN(64)[klon];
    double *zqold;
    zqold = new double DACE_ALIGN(64)[klon];
    double *zdtgdp;
    zdtgdp = new double DACE_ALIGN(64)[klon];
    double *zrdtgdp;
    zrdtgdp = new double DACE_ALIGN(64)[klon];
    double *zcovpclr;
    zcovpclr = new double DACE_ALIGN(64)[klon];
    double *zcovptot;
    zcovptot = new double DACE_ALIGN(64)[klon];
    double *zcovpmax;
    zcovpmax = new double DACE_ALIGN(64)[klon];
    double *zqpretot;
    zqpretot = new double DACE_ALIGN(64)[klon];
    double *zldefr;
    zldefr = new double DACE_ALIGN(64)[klon];
    double *zldifdt;
    zldifdt = new double DACE_ALIGN(64)[klon];
    double *zacust;
    zacust = new double DACE_ALIGN(64)[klon];
    double *zmf;
    zmf = new double DACE_ALIGN(64)[klon];
    double *zrho;
    zrho = new double DACE_ALIGN(64)[klon];
    double *zsolab;
    zsolab = new double DACE_ALIGN(64)[klon];
    double *zsolac;
    zsolac = new double DACE_ALIGN(64)[klon];
    double *zanewm1;
    zanewm1 = new double DACE_ALIGN(64)[klon];
    double *zda;
    zda = new double DACE_ALIGN(64)[klon];
    double *zdp;
    zdp = new double DACE_ALIGN(64)[klon];
    double *zsupsat;
    zsupsat = new double DACE_ALIGN(64)[klon];
    double *zmeltmax;
    zmeltmax = new double DACE_ALIGN(64)[klon];
    double *zfrzmax;
    zfrzmax = new double DACE_ALIGN(64)[klon];
    double *zicetot;
    zicetot = new double DACE_ALIGN(64)[klon];
    double *zdqsmixdt;
    zdqsmixdt = new double DACE_ALIGN(64)[klon];
    double *zcorqsice;
    zcorqsice = new double DACE_ALIGN(64)[klon];
    double *zevaplimmix;
    zevaplimmix = new double DACE_ALIGN(64)[klon];
    double *zcldtopdist;
    zcldtopdist = new double DACE_ALIGN(64)[klon];
    double *zrainacc;
    zrainacc = new double DACE_ALIGN(64)[klon];
    double *zraincld;
    zraincld = new double DACE_ALIGN(64)[klon];
    double *zsnowcld;
    zsnowcld = new double DACE_ALIGN(64)[klon];
    double *psum_solqa;
    psum_solqa = new double DACE_ALIGN(64)[klon];
    int *iphase;
    iphase = new int DACE_ALIGN(64)[5];
    int *imelt;
    imelt = new int DACE_ALIGN(64)[5];
    int *llfall;
    llfall = new int DACE_ALIGN(64)[5];
    double *zvqx;
    zvqx = new double DACE_ALIGN(64)[5];
    double *zfoealfa;
    zfoealfa = new double DACE_ALIGN(64)[(klon * (klev + 1))];
    double *ztp1;
    ztp1 = new double DACE_ALIGN(64)[(klev * klon)];
    double *zlcust;
    zlcust = new double DACE_ALIGN(64)[(5 * klon)];
    double *zli;
    zli = new double DACE_ALIGN(64)[(klev * klon)];
    double *za;
    za = new double DACE_ALIGN(64)[(klev * klon)];
    double *zaorig;
    zaorig = new double DACE_ALIGN(64)[(klev * klon)];
    double *zliqfrac;
    zliqfrac = new double DACE_ALIGN(64)[(klev * klon)];
    double *zicefrac;
    zicefrac = new double DACE_ALIGN(64)[(klev * klon)];
    double *zqx;
    zqx = new double DACE_ALIGN(64)[((5 * klev) * klon)];
    double *zqx0;
    zqx0 = new double DACE_ALIGN(64)[((5 * klev) * klon)];
    double *zqxn;
    zqxn = new double DACE_ALIGN(64)[(5 * klon)];
    double *zqxfg;
    zqxfg = new double DACE_ALIGN(64)[(5 * klon)];
    double *zqxnm1;
    zqxnm1 = new double DACE_ALIGN(64)[(5 * klon)];
    double *zfluxq;
    zfluxq = new double DACE_ALIGN(64)[(5 * klon)];
    double *zpfplsx;
    zpfplsx = new double DACE_ALIGN(64)[((5 * klon) * (klev + 1))];
    double *zlneg;
    zlneg = new double DACE_ALIGN(64)[((5 * klev) * klon)];
    double *zqxn2d;
    zqxn2d = new double DACE_ALIGN(64)[((5 * klev) * klon)];
    double *zqsmix;
    zqsmix = new double DACE_ALIGN(64)[(klev * klon)];
    double *zqsliq;
    zqsliq = new double DACE_ALIGN(64)[(klev * klon)];
    double *zqsice;
    zqsice = new double DACE_ALIGN(64)[(klev * klon)];
    double *zfoeewmt;
    zfoeewmt = new double DACE_ALIGN(64)[(klev * klon)];
    double *zfoeew;
    zfoeew = new double DACE_ALIGN(64)[(klev * klon)];
    double *zfoeeliqt;
    zfoeeliqt = new double DACE_ALIGN(64)[(klev * klon)];
    double *zsolqa;
    zsolqa = new double DACE_ALIGN(64)[(25 * klon)];
    double *zsolqb;
    zsolqb = new double DACE_ALIGN(64)[(25 * klon)];
    double *zqlhs;
    zqlhs = new double DACE_ALIGN(64)[(25 * klon)];
    double *zratio;
    zratio = new double DACE_ALIGN(64)[(5 * klon)];
    double *zsinksum;
    zsinksum = new double DACE_ALIGN(64)[(5 * klon)];
    double *zfallsink;
    zfallsink = new double DACE_ALIGN(64)[(5 * klon)];
    double *zfallsrce;
    zfallsrce = new double DACE_ALIGN(64)[(5 * klon)];
    double *zconvsrce;
    zconvsrce = new double DACE_ALIGN(64)[(5 * klon)];
    double *zconvsink;
    zconvsink = new double DACE_ALIGN(64)[(5 * klon)];
    double *zpsupsatsrce;
    zpsupsatsrce = new double DACE_ALIGN(64)[(5 * klon)];
    double ztw1;
    double ztw2;
    double ztw3;
    double ztw4;
    double ztw5;
    double zepsilon;
    double zqtmst;
    double zrdcp;
    double zrg_r;
    double zrldcp;
    double zqx_slice_plus_zqx_slice;
    double zqadj;
    double tendency_loc_t_slice_minus_ydthf_ralvdcp_zqadj;
    double tendency_loc_t_index_0;
    double zqsliq_index;
    double zli_index;
    double zfacw;
    double zcor;
    double zfaci;
    double zalfaw;
    double zfac;
    double __tmp26;
    double ztmpa;
    double zsolqa_slice_0;
    double zsolqa_slice_2;
    double zconvsrce_slice_0;
    double zdtdp;
    double zdtforc;
    double zevap;
    double zmfdn;
    double zdpmxdt;
    double zqp;
    double zqsat;
    double zcor_2;
    double zcor_3;
    double zdqs_slice_0;
    double zcdmax;
    double zlcond1_index_1;
    double zrhc;
    double zqe;
    double zrhc_zqsice_slice_times_zfac;
    double zqsice_slice_times_zfac;
    double zacond;
    double zfac_times_zdqs_slice;
    double neg_zzdl;
    double zlcond2_index_1;
    double __tmp120;
    double zzco;
    double zrainaut_index_2;
    double zrainacc_index_0;
    double zfallcorr;
    double zicetot_index;
    double zsubsat;
    double zcons1;
    double zmelt;
    double zlambda;
    double zfrz;
    double zzrh;
    double zpreclr;
    double zbeta;
    double zdenom;
    double zdpevap;
    double zzrh_times_zqsliq_slice_0;
    double zcovptot_slice_minus_expr_0;
    double zzrh_times_zqsice_slice;
    double zqxfg_slice_minus_zevap_1;
    double za_slice_plus_zsolac_slice;
    double __tmp264;
    double zzratio;
    double zexplicit;
    double zgdph_r;
    double pfcqsng_slice;
    double zconst;
    int kfdia_plus_1_6;
    int kfdia_plus_1_13;
    int kfdia_plus_1_85;
    int64_t jk;
    int kfdia_plus_1;
    int64_t jl;
    int64_t jm;
    int kfdia_plus_1_0;
    int kfdia_plus_1_1;
    double zvqx_index;
    int kfdia_plus_1_2;
    int kfdia_plus_1_3;
    int kfdia_plus_1_4;
    int kfdia_plus_1_5;
    int kfdia_plus_1_7;
    bool __tmp1;
    double za_index;
    bool __tmp2;
    int kfdia_plus_1_8;
    double zqx_index_9;
    bool __tmp4;
    int iphase_index;
    int iphase_index_0;
    int kfdia_plus_1_9;
    int kfdia_plus_1_10;
    bool __tmp15;
    int kfdia_plus_1_15;
    int kfdia_plus_1_18;
    int kfdia_plus_1_19;
    int kfdia_plus_1_20;
    int kfdia_plus_1_21;
    int kfdia_plus_1_28;
    int kfdia_plus_1_29;
    int kfdia_plus_1_30;
    int kfdia_plus_1_31;
    int kfdia_plus_1_32;
    int kfdia_plus_1_33;
    int kfdia_plus_1_34;
    int kfdia_plus_1_35;
    int kfdia_plus_1_36;
    int kfdia_plus_1_37;
    int kfdia_plus_1_39;
    int kfdia_plus_1_41;
    int kfdia_plus_1_42;
    int kfdia_plus_1_43;
    int kfdia_plus_1_44;
    int kfdia_plus_1_46;
    int kfdia_plus_1_47;
    int kfdia_plus_1_48;
    int kfdia_plus_1_50;
    int kfdia_plus_1_51;
    int kfdia_plus_1_54;
    int kfdia_plus_1_70;
    int kfdia_plus_1_76;
    int kfdia_plus_1_77;
    int kfdia_plus_1_82;
    int kfdia_plus_1_83;
    int kfdia_plus_1_14;
    int64_t jn;
    int kfdia_plus_1_16;
    int kfdia_plus_1_17;
    double zqx_index_22;
    bool __tmp29;
    double zqx_index_24;
    bool __tmp30;
    double ztp1_index_25;
    bool __tmp31;
    bool __tmp32;
    double za_index_4;
    double __tmp35;
    double zsupsat_index;
    double psupsat_index;
    double ztp1_index_26;
    bool __tmp42;
    double ztp1_index_27;
    bool __tmp45;
    int kfdia_plus_1_22;
    int ldcum_index;
    double plude_index_0;
    bool __tmp47;
    double plu_index;
    int ldcum_index_0;
    int kfdia_plus_1_23;
    int kfdia_plus_1_25;
    int kfdia_plus_1_27;
    int llfall_index;
    int iphase_index_1;
    int kfdia_plus_1_24;
    int llfall_index_0;
    int iphase_index_2;
    int kfdia_plus_1_26;
    double zlfinalsum_index_0;
    int ktype_index;
    double plude_index_4;
    double zli_index_1;
    double zdqs_index_0;
    double za_index_12;
    double zdqs_index_2;
    double neg_yrecldp_rlmin;
    double za_index_13;
    bool __tmp97;
    double ztp1_index_63;
    bool __tmp98;
    double zdqs_index_4;
    double neg_yrecldp_rlmin_0;
    double za_index_17;
    double pap_index_10;
    double paph_index_3;
    double ztp1_index_64;
    bool __tmp107;
    bool __tmp108;
    bool __tmp110;
    bool __tmp111;
    bool __tmp118;
    bool __tmp119;
    bool __tmp121;
    double zlcond2_index_2;
    double ztp1_index_65;
    bool __tmp124;
    int llfall_index_1;
    int kfdia_plus_1_40;
    double zqpretot_index_0;
    double ztp1_index_80;
    bool __tmp168;
    double zliqcld_index_1;
    double zicecld_index_3;
    double zlcrit;
    double plsm_index_0;
    double zliqcld_index_4;
    double ztp1_index_83;
    bool __tmp188;
    bool __tmp183;
    bool __tmp187;
    bool __tmp199;
    double ztp1_index_85;
    bool __tmp200;
    int iphase_index_3;
    int kfdia_plus_1_45;
    double zmeltmax_index;
    double zicetot_index_0;
    int imelt_index;
    int imelt_index_0;
    int imelt_index_1;
    int imelt_index_2;
    int imelt_index_3;
    int imelt_index_4;
    double zqx_index_37;
    double ztp1_index_88;
    bool __tmp209;
    double ztp1_index_89;
    bool __tmp210;
    double ztp1_index_90;
    bool __tmp213;
    double prainfrac_toprfz_index_0;
    double zfrzmax_index;
    double zfrzmax_index_1;
    double zqxfg_index_29;
    int imelt_index_5;
    int imelt_index_6;
    int imelt_index_7;
    int imelt_index_8;
    int llfall_index_2;
    int kfdia_plus_1_53;
    double zqxfg_index_47;
    bool __tmp263;
    double zanew;
    int kfdia_plus_1_56;
    int kfdia_plus_1_57;
    int kfdia_plus_1_58;
    int kfdia_plus_1_59;
    int kfdia_plus_1_61;
    int kfdia_plus_1_62;
    int kfdia_plus_1_63;
    int kfdia_plus_1_60;
    double zsolqa_index_58;
    int kfdia_plus_1_64;
    int64_t jo;
    int kfdia_plus_1_65;
    int kfdia_plus_1_66;
    int kfdia_plus_1_67;
    int64_t ik;
    int kfdia_plus_1_68;
    int kfdia_plus_1_69;
    int kfdia_plus_1_72;
    int kfdia_plus_1_71;
    int kfdia_plus_1_73;
    double zqxn_index_5;
    int kfdia_plus_1_74;
    int kfdia_plus_1_75;
    double zqpretot_index_2;
    int kfdia_plus_1_78;
    int iphase_index_4;
    int iphase_index_5;
    int kfdia_plus_1_81;
    int kfdia_plus_1_79;
    int kfdia_plus_1_80;
    double zcovpclr_index_5;
    double zqxfg_index_41;
    bool __tmp250;
    double zcovpclr_index_3;
    double zqxfg_index_36;
    bool __tmp235;
    double ztp1_index_84;
    bool __tmp189;
    double zliqcld_index_7;
    double zsnowcld_index;
    double zcovptot_index_6;
    double za_index_27;
    bool __tmp125;
    double za_index_28;
    bool __tmp126;
    double ztp1_index_66;
    bool __tmp128;
    double zqxfg_index_7;
    bool __tmp129;
    int kfdia_plus_1_84;
    int kfdia_plus_1_86;
    int kfdia_plus_1_87;

    {
        double ydthf_ralsdcp_minus_ydthf_ralvdcp;

        {
            double __out;

            ///////////////////
            // Tasklet code (assign_293_4)
            __out = 1329.31;
            ///////////////////

            ztw1 = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_294_4)
            __out = 0.0074615;
            ///////////////////

            ztw2 = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_295_4)
            __out = 85000.0;
            ///////////////////

            ztw3 = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_296_4)
            __out = 40.637;
            ///////////////////

            ztw4 = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_297_4)
            __out = 275.0;
            ///////////////////

            ztw5 = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_298_4)
            __out = 1e-14;
            ///////////////////

            zepsilon = __out;
        }
        {
            double __in2 = ptsphy;
            double __out;

            ///////////////////
            // Tasklet code (_Div_)
            __out = (1.0 / __in2);
            ///////////////////

            zqtmst = __out;
        }
        {
            double __in2 = ydcst_rg;
            double __out;

            ///////////////////
            // Tasklet code (_Div_)
            __out = (1.0 / __in2);
            ///////////////////

            zrg_r = __out;
        }
        {
            double __in1 = ydcst_rd;
            double __in2 = ydcst_rcpd;
            double __out;

            ///////////////////
            // Tasklet code (_Div_)
            __out = (__in1 / __in2);
            ///////////////////

            zrdcp = __out;
        }
        {
            double __in1 = ydthf_ralsdcp;
            double __in2 = ydthf_ralvdcp;
            double __out;

            ///////////////////
            // Tasklet code (_Sub_)
            __out = (__in1 - __in2);
            ///////////////////

            ydthf_ralsdcp_minus_ydthf_ralvdcp = __out;
        }
        {
            double __in2 = ydthf_ralsdcp_minus_ydthf_ralvdcp;
            double __out;

            ///////////////////
            // Tasklet code (_Div_)
            __out = (1.0 / __in2);
            ///////////////////

            zrldcp = __out;
        }
        {
            int __out;

            ///////////////////
            // Tasklet code (assign_310_4)
            __out = 0;
            ///////////////////

            iphase[4] = __out;
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_311_4)
            __out = 1;
            ///////////////////

            iphase[0] = __out;
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_312_4)
            __out = 1;
            ///////////////////

            iphase[2] = __out;
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_313_4)
            __out = 2;
            ///////////////////

            iphase[1] = __out;
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_314_4)
            __out = 2;
            ///////////////////

            iphase[3] = __out;
        }
        {
            int __out;

            ///////////////////
            // Tasklet code (assign_315_4)
            __out = -99;
            ///////////////////

            imelt[4] = __out;
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_316_4)
            __out = 2;
            ///////////////////

            imelt[0] = __out;
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_317_4)
            __out = 4;
            ///////////////////

            imelt[2] = __out;
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_318_4)
            __out = 3;
            ///////////////////

            imelt[1] = __out;
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_319_4)
            __out = 3;
            ///////////////////

            imelt[3] = __out;
        }

    }
    for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

        kfdia_plus_1 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1); jl = (jl + 1)) {
            {

                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_322_12)
                    __out = 0.0;
                    ///////////////////

                    tendency_loc_t[((jl + (klon * (jk - 1))) - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_323_12)
                    __out = 0.0;
                    ///////////////////

                    tendency_loc_q[((jl + (klon * (jk - 1))) - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_324_12)
                    __out = 0.0;
                    ///////////////////

                    tendency_loc_a[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }

    }
    for (jm = 1; (jm < 5); jm = (jm + 1)) {
        for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

            kfdia_plus_1_0 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_0); jl = (jl + 1)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_328_16)
                        __out = 0.0;
                        ///////////////////

                        tendency_loc_cld[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                    }

                }

            }

        }

    }
    for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

        kfdia_plus_1_1 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_1); jl = (jl + 1)) {
            {

                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_331_12)
                    __out = 0.0;
                    ///////////////////

                    pcovptot[((jl + (klon * (jk - 1))) - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_332_12)
                    __out = 0.0;
                    ///////////////////

                    tendency_loc_cld[(((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }

    }
    {

        {
            double __out;

            ///////////////////
            // Tasklet code (assign_333_4)
            __out = 0.0;
            ///////////////////

            zvqx[4] = __out;
        }

    }
    {

        {
            double __out;

            ///////////////////
            // Tasklet code (assign_334_4)
            __out = 0.0;
            ///////////////////

            zvqx[0] = __out;
        }

    }
    {

        {
            double __inp = yrecldp_rvice;
            double __out;

            ///////////////////
            // Tasklet code (assign_335_4)
            __out = __inp;
            ///////////////////

            zvqx[1] = __out;
        }

    }
    {

        {
            double __inp = yrecldp_rvrain;
            double __out;

            ///////////////////
            // Tasklet code (assign_336_4)
            __out = __inp;
            ///////////////////

            zvqx[2] = __out;
        }

    }
    {

        {
            double __inp = yrecldp_rvsnow;
            double __out;

            ///////////////////
            // Tasklet code (assign_337_4)
            __out = __inp;
            ///////////////////

            zvqx[3] = __out;
        }
        {
            #pragma omp parallel for
            for (auto __i0 = 0; __i0 < 5; __i0 += 1) {
                {
                    int __out;

                    ///////////////////
                    // Tasklet code (assign_338_4)
                    __out = false;
                    ///////////////////

                    llfall[__i0] = __out;
                }
            }
        }

    }
    for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

        zvqx_index = zvqx[(jm - 1)];
        if ((zvqx_index > 0.0)) {
            {

                {
                    int __out;

                    ///////////////////
                    // Tasklet code (assign_341_12)
                    __out = true;
                    ///////////////////

                    llfall[(jm - 1)] = __out;
                }

            }
        }

    }
    {

        {
            int __out;

            ///////////////////
            // Tasklet code (assign_342_4)
            __out = false;
            ///////////////////

            llfall[1] = __out;
        }

    }
    for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

        kfdia_plus_1_2 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_2); jl = (jl + 1)) {
            {
                double pt_index;
                double tendency_tmp_t_index;
                double ptsphy_times_tendency_tmp_t_slice;
                double ztp1_slice;
                double pq_index;
                double tendency_tmp_q_index;
                double ptsphy_times_tendency_tmp_q_slice;
                double zqx_slice;
                double pq_index_0;
                double tendency_tmp_q_index_0;
                double ptsphy_times_tendency_tmp_q_slice_0;
                double zqx0_slice;
                double pa_index;
                double tendency_tmp_a_index;
                double ptsphy_times_tendency_tmp_a_slice;
                double za_slice;
                double pa_index_0;
                double tendency_tmp_a_index_0;
                double ptsphy_times_tendency_tmp_a_slice_0;
                double zaorig_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pt + ((jl + (klon * (jk - 1))) - 1), &pt_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                tendency_tmp_t + ((jl + (klon * (jk - 1))) - 1), &tendency_tmp_t_index, 1);
                {
                    double __in1 = ptsphy;
                    double __in2 = tendency_tmp_t_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ptsphy_times_tendency_tmp_t_slice = __out;
                }
                {
                    double __in1 = pt_index;
                    double __in2 = ptsphy_times_tendency_tmp_t_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    ztp1_slice = __out;
                }
                {
                    double __inp = ztp1_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_345_12)
                    __out = __inp;
                    ///////////////////

                    ztp1[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pq + ((jl + (klon * (jk - 1))) - 1), &pq_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                tendency_tmp_q + ((jl + (klon * (jk - 1))) - 1), &tendency_tmp_q_index, 1);
                {
                    double __in1 = ptsphy;
                    double __in2 = tendency_tmp_q_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ptsphy_times_tendency_tmp_q_slice = __out;
                }
                {
                    double __in1 = pq_index;
                    double __in2 = ptsphy_times_tendency_tmp_q_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zqx_slice = __out;
                }
                {
                    double __inp = zqx_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_346_12)
                    __out = __inp;
                    ///////////////////

                    zqx[(((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pq + ((jl + (klon * (jk - 1))) - 1), &pq_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                tendency_tmp_q + ((jl + (klon * (jk - 1))) - 1), &tendency_tmp_q_index_0, 1);
                {
                    double __in1 = ptsphy;
                    double __in2 = tendency_tmp_q_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ptsphy_times_tendency_tmp_q_slice_0 = __out;
                }
                {
                    double __in1 = pq_index_0;
                    double __in2 = ptsphy_times_tendency_tmp_q_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zqx0_slice = __out;
                }
                {
                    double __inp = zqx0_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_347_12)
                    __out = __inp;
                    ///////////////////

                    zqx0[(((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pa + ((jl + (klon * (jk - 1))) - 1), &pa_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                tendency_tmp_a + ((jl + (klon * (jk - 1))) - 1), &tendency_tmp_a_index, 1);
                {
                    double __in1 = ptsphy;
                    double __in2 = tendency_tmp_a_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ptsphy_times_tendency_tmp_a_slice = __out;
                }
                {
                    double __in1 = pa_index;
                    double __in2 = ptsphy_times_tendency_tmp_a_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    za_slice = __out;
                }
                {
                    double __inp = za_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_348_12)
                    __out = __inp;
                    ///////////////////

                    za[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pa + ((jl + (klon * (jk - 1))) - 1), &pa_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                tendency_tmp_a + ((jl + (klon * (jk - 1))) - 1), &tendency_tmp_a_index_0, 1);
                {
                    double __in1 = ptsphy;
                    double __in2 = tendency_tmp_a_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ptsphy_times_tendency_tmp_a_slice_0 = __out;
                }
                {
                    double __in1 = pa_index_0;
                    double __in2 = ptsphy_times_tendency_tmp_a_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zaorig_slice = __out;
                }
                {
                    double __inp = zaorig_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_349_12)
                    __out = __inp;
                    ///////////////////

                    zaorig[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }

    }
    for (jm = 1; (jm < 5); jm = (jm + 1)) {
        for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

            kfdia_plus_1_3 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_3); jl = (jl + 1)) {
                {
                    double pclv_index;
                    double tendency_tmp_cld_index;
                    double ptsphy_times_tendency_tmp_cld_slice;
                    double zqx_slice_0;
                    double pclv_index_0;
                    double tendency_tmp_cld_index_0;
                    double ptsphy_times_tendency_tmp_cld_slice_0;
                    double zqx0_slice_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pclv + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &pclv_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    tendency_tmp_cld + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &tendency_tmp_cld_index, 1);
                    {
                        double __in1 = ptsphy;
                        double __in2 = tendency_tmp_cld_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ptsphy_times_tendency_tmp_cld_slice = __out;
                    }
                    {
                        double __in1 = pclv_index;
                        double __in2 = ptsphy_times_tendency_tmp_cld_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zqx_slice_0 = __out;
                    }
                    {
                        double __inp = zqx_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_353_16)
                        __out = __inp;
                        ///////////////////

                        zqx[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pclv + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &pclv_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    tendency_tmp_cld + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &tendency_tmp_cld_index_0, 1);
                    {
                        double __in1 = ptsphy;
                        double __in2 = tendency_tmp_cld_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ptsphy_times_tendency_tmp_cld_slice_0 = __out;
                    }
                    {
                        double __in1 = pclv_index_0;
                        double __in2 = ptsphy_times_tendency_tmp_cld_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zqx0_slice_0 = __out;
                    }
                    {
                        double __inp = zqx0_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_354_16)
                        __out = __inp;
                        ///////////////////

                        zqx0[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                    }

                }

            }

        }

    }
    for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {
        for (jk = 1; (jk < (klev + 2)); jk = (jk + 1)) {

            kfdia_plus_1_4 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_4); jl = (jl + 1)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_358_16)
                        __out = 0.0;
                        ///////////////////

                        zpfplsx[(((jl + (klon * (jk - 1))) + ((klon * (jm - 1)) * (klev + 1))) - 1)] = __out;
                    }

                }

            }

        }

    }
    for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {
        for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

            kfdia_plus_1_5 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_5); jl = (jl + 1)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_362_16)
                        __out = 0.0;
                        ///////////////////

                        zqxn2d[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_363_16)
                        __out = 0.0;
                        ///////////////////

                        zlneg[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                    }

                }

            }

        }

    }
    kfdia_plus_1_6 = (kfdia + 1);
    for (jl = kidia; (jl < kfdia_plus_1_6); jl = (jl + 1)) {
        {

            {
                double __out;

                ///////////////////
                // Tasklet code (assign_365_8)
                __out = 0.0;
                ///////////////////

                prainfrac_toprfz[(jl - 1)] = __out;
            }

        }

    }
    for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

        kfdia_plus_1_7 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_7); jl = (jl + 1)) {
            {
                double zqx_index;
                double zqx_index_0;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + ((jl + (klon * (jk - 1))) - 1), &zqx_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqx_index_0, 1);
                {
                    double __in1 = zqx_index;
                    double __in2 = zqx_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zqx_slice_plus_zqx_slice = __out;
                }

            }
            __tmp1 = (zqx_slice_plus_zqx_slice < yrecldp_rlmin);
            za_index = za[((jl + (klon * (jk - 1))) - 1)];

            __tmp2 = (za_index < yrecldp_ramin);
            if ((__tmp1 || __tmp2)) {
                {
                    double zlneg_index;
                    double zqx_index_1;
                    double zlneg_slice_plus_zqx_slice;
                    double zqx_index_2;
                    double tendency_loc_q_index;
                    double tendency_loc_q_slice_plus_zqadj;
                    double tendency_loc_t_index;
                    double ydthf_ralvdcp_times_zqadj;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zlneg + ((jl + (klon * (jk - 1))) - 1), &zlneg_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + ((jl + (klon * (jk - 1))) - 1), &zqx_index_1, 1);
                    {
                        double __in1 = zlneg_index;
                        double __in2 = zqx_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zlneg_slice_plus_zqx_slice = __out;
                    }
                    {
                        double __inp = zlneg_slice_plus_zqx_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_370_16)
                        __out = __inp;
                        ///////////////////

                        zlneg[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + ((jl + (klon * (jk - 1))) - 1), &zqx_index_2, 1);
                    {
                        double __in1 = zqx_index_2;
                        double __in2 = zqtmst;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zqadj = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    tendency_loc_q + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_q_index, 1);
                    {
                        double __in1 = tendency_loc_q_index;
                        double __in2 = zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        tendency_loc_q_slice_plus_zqadj = __out;
                    }
                    {
                        double __inp = tendency_loc_q_slice_plus_zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_372_16)
                        __out = __inp;
                        ///////////////////

                        tendency_loc_q[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    tendency_loc_t + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_t_index, 1);
                    {
                        double __in1 = ydthf_ralvdcp;
                        double __in2 = zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydthf_ralvdcp_times_zqadj = __out;
                    }
                    {
                        double __in1 = tendency_loc_t_index;
                        double __in2 = ydthf_ralvdcp_times_zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        tendency_loc_t_slice_minus_ydthf_ralvdcp_zqadj = __out;
                    }

                }
                {
                    double zqx_index_3;
                    double zqx_index_4;
                    double zqx_slice_plus_zqx_slice_0;

                    {
                        double __inp = tendency_loc_t_slice_minus_ydthf_ralvdcp_zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_373_16)
                        __out = __inp;
                        ///////////////////

                        tendency_loc_t[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_3, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + ((jl + (klon * (jk - 1))) - 1), &zqx_index_4, 1);
                    {
                        double __in1 = zqx_index_3;
                        double __in2 = zqx_index_4;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zqx_slice_plus_zqx_slice_0 = __out;
                    }
                    {
                        double __inp = zqx_slice_plus_zqx_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_374_16)
                        __out = __inp;
                        ///////////////////

                        zqx[(((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1)] = __out;
                    }

                }
                {
                    double zlneg_index_0;
                    double zqx_index_5;
                    double zlneg_slice_plus_zqx_slice_0;
                    double zqx_index_6;
                    double zqadj_0;
                    double tendency_loc_q_index_0;
                    double tendency_loc_q_slice_plus_zqadj_0;

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_375_16)
                        __out = 0.0;
                        ///////////////////

                        zqx[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqx_index_5, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqx_index_6, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zlneg + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zlneg_index_0, 1);
                    {
                        double __in1 = zlneg_index_0;
                        double __in2 = zqx_index_5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zlneg_slice_plus_zqx_slice_0 = __out;
                    }
                    {
                        double __inp = zlneg_slice_plus_zqx_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_376_16)
                        __out = __inp;
                        ///////////////////

                        zlneg[(((jl + (klev * klon)) + (klon * (jk - 1))) - 1)] = __out;
                    }
                    {
                        double __in1 = zqx_index_6;
                        double __in2 = zqtmst;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zqadj_0 = __out;
                    }
                    {
                        double __inp = zqadj_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_377_16)
                        __out = __inp;
                        ///////////////////

                        zqadj = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    tendency_loc_q + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_q_index_0, 1);
                    {
                        double __in1 = tendency_loc_q_index_0;
                        double __in2 = zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        tendency_loc_q_slice_plus_zqadj_0 = __out;
                    }
                    {
                        double __inp = tendency_loc_q_slice_plus_zqadj_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_378_16)
                        __out = __inp;
                        ///////////////////

                        tendency_loc_q[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    tendency_loc_t + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_t_index_0, 1);

                }
                {
                    double ydthf_ralsdcp_times_zqadj;
                    double tendency_loc_t_slice_minus_ydthf_ralsdcp_zqadj;
                    double zqx_index_7;
                    double zqx_index_8;
                    double zqx_slice_plus_zqx_slice_1;

                    {
                        double __in1 = ydthf_ralsdcp;
                        double __in2 = zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydthf_ralsdcp_times_zqadj = __out;
                    }
                    {
                        double __in1 = tendency_loc_t_index_0;
                        double __in2 = ydthf_ralsdcp_times_zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        tendency_loc_t_slice_minus_ydthf_ralsdcp_zqadj = __out;
                    }
                    {
                        double __inp = tendency_loc_t_slice_minus_ydthf_ralsdcp_zqadj;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_379_16)
                        __out = __inp;
                        ///////////////////

                        tendency_loc_t[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_7, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqx_index_8, 1);
                    {
                        double __in1 = zqx_index_7;
                        double __in2 = zqx_index_8;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zqx_slice_plus_zqx_slice_1 = __out;
                    }
                    {
                        double __inp = zqx_slice_plus_zqx_slice_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_380_16)
                        __out = __inp;
                        ///////////////////

                        zqx[(((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1)] = __out;
                    }

                }
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_381_16)
                        __out = 0.0;
                        ///////////////////

                        zqx[(((jl + (klev * klon)) + (klon * (jk - 1))) - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_382_16)
                        __out = 0.0;
                        ///////////////////

                        za[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                }
            }

        }

    }
    for (jm = 1; (jm < 5); jm = (jm + 1)) {
        for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

            kfdia_plus_1_8 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_8); jl = (jl + 1)) {

                zqx_index_9 = zqx[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)];

                __tmp4 = (zqx_index_9 < yrecldp_rlmin);
                if (__tmp4) {
                    {
                        double zlneg_index_1;
                        double zqx_index_10;
                        double zlneg_slice_plus_zqx_slice_1;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlneg + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zlneg_index_1, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx_index_10, 1);
                        {
                            double __in1 = zlneg_index_1;
                            double __in2 = zqx_index_10;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zlneg_slice_plus_zqx_slice_1 = __out;
                        }
                        {
                            double __inp = zlneg_slice_plus_zqx_slice_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_387_20)
                            __out = __inp;
                            ///////////////////

                            zlneg[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                        }

                    }
                    iphase_index = iphase[(jm - 1)];
                    {
                        double zqx_index_11;
                        double zqadj_1;
                        double tendency_loc_q_index_1;
                        double tendency_loc_q_slice_plus_zqadj_1;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx_index_11, 1);
                        {
                            double __in1 = zqx_index_11;
                            double __in2 = zqtmst;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zqadj_1 = __out;
                        }
                        {
                            double __inp = zqadj_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_388_20)
                            __out = __inp;
                            ///////////////////

                            zqadj = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        tendency_loc_q + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_q_index_1, 1);
                        {
                            double __in1 = tendency_loc_q_index_1;
                            double __in2 = zqadj;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            tendency_loc_q_slice_plus_zqadj_1 = __out;
                        }
                        {
                            double __inp = tendency_loc_q_slice_plus_zqadj_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_389_20)
                            __out = __inp;
                            ///////////////////

                            tendency_loc_q[((jl + (klon * (jk - 1))) - 1)] = __out;
                        }

                    }
                    if ((iphase_index == 1)) {
                        {
                            double tendency_loc_t_index_1;
                            double ydthf_ralvdcp_times_zqadj_0;
                            double tendency_loc_t_slice_minus_ydthf_ralvdcp_zqadj_0;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            tendency_loc_t + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_t_index_1, 1);
                            {
                                double __in1 = ydthf_ralvdcp;
                                double __in2 = zqadj;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                ydthf_ralvdcp_times_zqadj_0 = __out;
                            }
                            {
                                double __in1 = tendency_loc_t_index_1;
                                double __in2 = ydthf_ralvdcp_times_zqadj_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                tendency_loc_t_slice_minus_ydthf_ralvdcp_zqadj_0 = __out;
                            }
                            {
                                double __inp = tendency_loc_t_slice_minus_ydthf_ralvdcp_zqadj_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_391_24)
                                __out = __inp;
                                ///////////////////

                                tendency_loc_t[((jl + (klon * (jk - 1))) - 1)] = __out;
                            }

                        }
                    }
                    iphase_index_0 = iphase[(jm - 1)];
                    if ((iphase_index_0 == 2)) {
                        {
                            double tendency_loc_t_index_2;
                            double ydthf_ralsdcp_times_zqadj_0;
                            double tendency_loc_t_slice_minus_ydthf_ralsdcp_zqadj_0;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            tendency_loc_t + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_t_index_2, 1);
                            {
                                double __in1 = ydthf_ralsdcp;
                                double __in2 = zqadj;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                ydthf_ralsdcp_times_zqadj_0 = __out;
                            }
                            {
                                double __in1 = tendency_loc_t_index_2;
                                double __in2 = ydthf_ralsdcp_times_zqadj_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                tendency_loc_t_slice_minus_ydthf_ralsdcp_zqadj_0 = __out;
                            }
                            {
                                double __inp = tendency_loc_t_slice_minus_ydthf_ralsdcp_zqadj_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_393_24)
                                __out = __inp;
                                ///////////////////

                                tendency_loc_t[((jl + (klon * (jk - 1))) - 1)] = __out;
                            }

                        }
                    }
                    {
                        double zqx_index_12;
                        double zqx_index_13;
                        double zqx_slice_plus_zqx_slice_2;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_12, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx_index_13, 1);
                        {
                            double __in1 = zqx_index_12;
                            double __in2 = zqx_index_13;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqx_slice_plus_zqx_slice_2 = __out;
                        }
                        {
                            double __inp = zqx_slice_plus_zqx_slice_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_394_20)
                            __out = __inp;
                            ///////////////////

                            zqx[(((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1)] = __out;
                        }

                    }
                    {

                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_395_20)
                            __out = 0.0;
                            ///////////////////

                            zqx[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                        }

                    }
                }

            }

        }

    }
    for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

        kfdia_plus_1_9 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_9); jl = (jl + 1)) {
            {
                double ztp1_index;
                double min_ydthf_rtwat_ztp1_slice;
                double max_ydthf_rtice_expr;
                double expr_minus_ydthf_rtice;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2;
                double zfoealfa_slice;
                double ztp1_index_0;
                double min_ydthf_rtwat_ztp1_slice_0;
                double max_ydthf_rtice_expr_0;
                double expr_minus_ydthf_rtice_0;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_0;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_0;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2;
                double ztp1_index_1;
                double ztp1_slice_minus_ydcst_rtt;
                double ydthf_r3les_times_ztp1_slice_ydcst_rtt;
                double ztp1_index_2;
                double ztp1_slice_minus_ydthf_r4les;
                double ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les;
                double exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les;
                double expr_times_expr;
                double ztp1_index_3;
                double min_ydthf_rtwat_ztp1_slice_1;
                double max_ydthf_rtice_expr_1;
                double expr_minus_ydthf_rtice_1;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_1;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_1;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_0;
                double __tmp7;
                double ztp1_index_4;
                double ztp1_slice_minus_ydcst_rtt_0;
                double ydthf_r3ies_times_ztp1_slice_ydcst_rtt;
                double ztp1_index_5;
                double ztp1_slice_minus_ydthf_r4ies;
                double ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies;
                double exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies;
                double __tmp8;
                double expr_expr_plus_1_0_expr_expr;
                double ydthf_r2es_times_expr_expr_1_0_expr_expr;
                double pap_index;
                double ydthf_r2es_expr_expr_1_0_expr_expr_div_pap_slice;
                double zfoeewmt_slice;
                double zqsmix_slice;
                double zqsmix_index;
                double zqsmix_index_0;
                double ydcst_retv_times_zqsmix_slice;
                double __tmp9;
                double zqsmix_slice_div_1_0_ydcst_retv_zqsmix_slice;
                double ztp1_index_6;
                double ztp1_slice_minus_ydcst_rtt_1;
                double sign_ztp1_slice_ydcst_rtt;
                double __tmp10;
                double zalfa;
                double ztp1_index_7;
                double ztp1_slice_minus_ydcst_rtt_2;
                double ydthf_r3les_times_ztp1_slice_ydcst_rtt_0;
                double ztp1_index_8;
                double ztp1_slice_minus_ydthf_r4les_0;
                double ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_0;
                double exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_0;
                double ydthf_r2es_times_expr;
                double zalfa_times_ydthf_r2es_expr;
                double __tmp11;
                double ztp1_index_9;
                double ztp1_slice_minus_ydcst_rtt_3;
                double ydthf_r3ies_times_ztp1_slice_ydcst_rtt_0;
                double ztp1_index_10;
                double ztp1_slice_minus_ydthf_r4ies_0;
                double ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_0;
                double exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_0;
                double ydthf_r2es_times_expr_0;
                double __tmp12;
                double zalfa_ydthf_r2es_expr_plus_1_0_zalfa_ydthf_r2es_expr;
                double pap_index_0;
                double zalfa_ydthf_r2es_expr_1_0_zalfa_ydthf_r2es_expr_div_pap_slice;
                double zfoeew_slice;
                double zfoeew_index;
                double min_0_5_zfoeew_slice;
                double zfoeew_index_0;
                double zfoeew_index_1;
                double ydcst_retv_times_zfoeew_slice;
                double __tmp13;
                double zqsice_slice;
                double ztp1_index_11;
                double ztp1_slice_minus_ydcst_rtt_4;
                double ydthf_r3les_times_ztp1_slice_ydcst_rtt_1;
                double ztp1_index_12;
                double ztp1_slice_minus_ydthf_r4les_1;
                double ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_1;
                double exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_1;
                double ydthf_r2es_times_expr_1;
                double pap_index_1;
                double ydthf_r2es_expr_div_pap_slice;
                double zfoeeliqt_slice;
                double zqsliq_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index, 1);
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    zfoealfa_slice = __out;
                }
                {
                    double __inp = zfoealfa_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_398_12)
                    __out = __inp;
                    ///////////////////

                    zfoealfa[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_11, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_12, 1);
                {
                    double __in1 = ztp1_index_11;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_4 = __out;
                }
                {
                    double __in1 = ydthf_r3les;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3les_times_ztp1_slice_ydcst_rtt_1 = __out;
                }
                {
                    double __in1 = ztp1_index_12;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les_1 = __out;
                }
                {
                    double __in1 = ydthf_r3les_times_ztp1_slice_ydcst_rtt_1;
                    double __in2 = ztp1_slice_minus_ydthf_r4les_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_1 = __out;
                }
                {
                    double __in1 = ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_1 = __out;
                }
                {
                    double __in1 = ydthf_r2es;
                    double __in2 = exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r2es_times_expr_1 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_1, 1);
                {
                    double __in1 = ydthf_r2es_times_expr_1;
                    double __in2 = pap_index_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r2es_expr_div_pap_slice = __out;
                }
                {
                    double __in_a = ydthf_r2es_expr_div_pap_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, 0.5);
                    ///////////////////

                    zfoeeliqt_slice = __out;
                }
                {
                    double __inp = zfoeeliqt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_406_12)
                    __out = __inp;
                    ///////////////////

                    zfoeeliqt[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoeeliqt + ((jl + (klon * (jk - 1))) - 1), &zqsliq_slice, 1);
                {
                    double __inp = zqsliq_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_407_12)
                    __out = __inp;
                    ///////////////////

                    zqsliq[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsliq + ((jl + (klon * (jk - 1))) - 1), &zqsliq_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_1, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_2, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_3, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_4, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_5, 1);
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_0 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_1 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_0 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_0;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_0 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_1 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_1;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_1 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_0;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_0 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_0 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_1;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_1 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_1 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_0 = __out;
                }
                {
                    double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp7 = __out;
                }
                {
                    double __in1 = ztp1_index_1;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt = __out;
                }
                {
                    double __in1 = ztp1_index_4;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_0 = __out;
                }
                {
                    double __in1 = ydthf_r3les;
                    double __in2 = ztp1_slice_minus_ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3les_times_ztp1_slice_ydcst_rtt = __out;
                }
                {
                    double __in1 = ztp1_index_2;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les = __out;
                }
                {
                    double __in1 = ydthf_r3les_times_ztp1_slice_ydcst_rtt;
                    double __in2 = ztp1_slice_minus_ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les = __out;
                }
                {
                    double __in1 = ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les = __out;
                }
                {
                    double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2;
                    double __in2 = exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_expr = __out;
                }
                {
                    double __in1 = ydthf_r3ies;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3ies_times_ztp1_slice_ydcst_rtt = __out;
                }
                {
                    double __in1 = ztp1_index_5;
                    double __in2 = ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4ies = __out;
                }
                {
                    double __in1 = ydthf_r3ies_times_ztp1_slice_ydcst_rtt;
                    double __in2 = ztp1_slice_minus_ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies = __out;
                }
                {
                    double __in1 = ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies = __out;
                }
                {
                    double __in1 = __tmp7;
                    double __in2 = exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp8 = __out;
                }
                {
                    double __in1 = expr_times_expr;
                    double __in2 = __tmp8;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    expr_expr_plus_1_0_expr_expr = __out;
                }
                {
                    double __in1 = ydthf_r2es;
                    double __in2 = expr_expr_plus_1_0_expr_expr;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r2es_times_expr_expr_1_0_expr_expr = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pap + ((jl + (klon * (jk - 1))) - 1), &pap_index, 1);
                {
                    double __in1 = ydthf_r2es_times_expr_expr_1_0_expr_expr;
                    double __in2 = pap_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r2es_expr_expr_1_0_expr_expr_div_pap_slice = __out;
                }
                {
                    double __in_a = ydthf_r2es_expr_expr_1_0_expr_expr_div_pap_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, 0.5);
                    ///////////////////

                    zfoeewmt_slice = __out;
                }
                {
                    double __inp = zfoeewmt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_399_12)
                    __out = __inp;
                    ///////////////////

                    zfoeewmt[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoeewmt + ((jl + (klon * (jk - 1))) - 1), &zqsmix_slice, 1);
                {
                    double __inp = zqsmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_400_12)
                    __out = __inp;
                    ///////////////////

                    zqsmix[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_0, 1);
                {
                    double __in1 = ydcst_retv;
                    double __in2 = zqsmix_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_retv_times_zqsmix_slice = __out;
                }
                {
                    double __in2 = ydcst_retv_times_zqsmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp9 = __out;
                }
                {
                    double __in1 = zqsmix_index;
                    double __in2 = __tmp9;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zqsmix_slice_div_1_0_ydcst_retv_zqsmix_slice = __out;
                }
                {
                    double __inp = zqsmix_slice_div_1_0_ydcst_retv_zqsmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_401_12)
                    __out = __inp;
                    ///////////////////

                    zqsmix[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_6, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_7, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_8, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_9, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_10, 1);
                {
                    double __in1 = ztp1_index_6;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_1 = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_ydcst_rtt_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_sign_)
                    __out = sign_numpy_2(__in1);
                    ///////////////////

                    sign_ztp1_slice_ydcst_rtt = __out;
                }
                {
                    double __in2 = sign_ztp1_slice_ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (1.0 * __in2);
                    ///////////////////

                    __tmp10 = __out;
                }
                {
                    double __in_b = __tmp10;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(0.0, __in_b);
                    ///////////////////

                    zalfa = __out;
                }
                {
                    double __in2 = zalfa;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp11 = __out;
                }
                {
                    double __in1 = ztp1_index_7;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_2 = __out;
                }
                {
                    double __in1 = ztp1_index_9;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_3 = __out;
                }
                {
                    double __in1 = ydthf_r3les;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3les_times_ztp1_slice_ydcst_rtt_0 = __out;
                }
                {
                    double __in1 = ztp1_index_8;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les_0 = __out;
                }
                {
                    double __in1 = ydthf_r3les_times_ztp1_slice_ydcst_rtt_0;
                    double __in2 = ztp1_slice_minus_ydthf_r4les_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_0 = __out;
                }
                {
                    double __in1 = ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_0 = __out;
                }
                {
                    double __in1 = ydthf_r3ies;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3ies_times_ztp1_slice_ydcst_rtt_0 = __out;
                }
                {
                    double __in1 = ztp1_index_10;
                    double __in2 = ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4ies_0 = __out;
                }
                {
                    double __in1 = ydthf_r3ies_times_ztp1_slice_ydcst_rtt_0;
                    double __in2 = ztp1_slice_minus_ydthf_r4ies_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_0 = __out;
                }
                {
                    double __in1 = ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_0 = __out;
                }
                {
                    double __in1 = ydthf_r2es;
                    double __in2 = exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r2es_times_expr = __out;
                }
                {
                    double __in1 = zalfa;
                    double __in2 = ydthf_r2es_times_expr;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zalfa_times_ydthf_r2es_expr = __out;
                }
                {
                    double __in1 = ydthf_r2es;
                    double __in2 = exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r2es_times_expr_0 = __out;
                }
                {
                    double __in1 = __tmp11;
                    double __in2 = ydthf_r2es_times_expr_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp12 = __out;
                }
                {
                    double __in1 = zalfa_times_ydthf_r2es_expr;
                    double __in2 = __tmp12;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zalfa_ydthf_r2es_expr_plus_1_0_zalfa_ydthf_r2es_expr = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_0, 1);
                {
                    double __in1 = zalfa_ydthf_r2es_expr_plus_1_0_zalfa_ydthf_r2es_expr;
                    double __in2 = pap_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zalfa_ydthf_r2es_expr_1_0_zalfa_ydthf_r2es_expr_div_pap_slice = __out;
                }
                {
                    double __in_a = zalfa_ydthf_r2es_expr_1_0_zalfa_ydthf_r2es_expr_div_pap_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, 0.5);
                    ///////////////////

                    zfoeew_slice = __out;
                }
                {
                    double __inp = zfoeew_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_403_12)
                    __out = __inp;
                    ///////////////////

                    zfoeew[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoeew + ((jl + (klon * (jk - 1))) - 1), &zfoeew_index, 1);
                {
                    double __in_b = zfoeew_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(0.5, __in_b);
                    ///////////////////

                    min_0_5_zfoeew_slice = __out;
                }
                {
                    double __inp = min_0_5_zfoeew_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_404_12)
                    __out = __inp;
                    ///////////////////

                    zfoeew[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoeew + ((jl + (klon * (jk - 1))) - 1), &zfoeew_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoeew + ((jl + (klon * (jk - 1))) - 1), &zfoeew_index_1, 1);
                {
                    double __in1 = ydcst_retv;
                    double __in2 = zfoeew_index_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_retv_times_zfoeew_slice = __out;
                }
                {
                    double __in2 = ydcst_retv_times_zfoeew_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp13 = __out;
                }
                {
                    double __in1 = zfoeew_index_0;
                    double __in2 = __tmp13;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zqsice_slice = __out;
                }
                {
                    double __inp = zqsice_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_405_12)
                    __out = __inp;
                    ///////////////////

                    zqsice[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }
            {
                double zqsliq_index_0;
                double ydcst_retv_times_zqsliq_slice;
                double __tmp14;
                double zqsliq_slice_div_1_0_ydcst_retv_zqsliq_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsliq + ((jl + (klon * (jk - 1))) - 1), &zqsliq_index_0, 1);
                {
                    double __in1 = ydcst_retv;
                    double __in2 = zqsliq_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_retv_times_zqsliq_slice = __out;
                }
                {
                    double __in2 = ydcst_retv_times_zqsliq_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp14 = __out;
                }
                {
                    double __in1 = zqsliq_index;
                    double __in2 = __tmp14;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zqsliq_slice_div_1_0_ydcst_retv_zqsliq_slice = __out;
                }
                {
                    double __inp = zqsliq_slice_div_1_0_ydcst_retv_zqsliq_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_408_12)
                    __out = __inp;
                    ///////////////////

                    zqsliq[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }

    }
    for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

        kfdia_plus_1_10 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_10); jl = (jl + 1)) {
            {
                double za_index_0;
                double min_1_0_za_slice;
                double max_0_0_expr;
                double zqx_index_14;
                double zqx_index_15;
                double zli_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                za + ((jl + (klon * (jk - 1))) - 1), &za_index_0, 1);
                {
                    double __in_b = za_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_za_slice = __out;
                }
                {
                    double __in_b = min_1_0_za_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(0.0, __in_b);
                    ///////////////////

                    max_0_0_expr = __out;
                }
                {
                    double __inp = max_0_0_expr;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_411_12)
                    __out = __inp;
                    ///////////////////

                    za[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + ((jl + (klon * (jk - 1))) - 1), &zqx_index_14, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqx_index_15, 1);
                {
                    double __in1 = zqx_index_14;
                    double __in2 = zqx_index_15;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zli_slice = __out;
                }
                {
                    double __inp = zli_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_412_12)
                    __out = __inp;
                    ///////////////////

                    zli[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zli + ((jl + (klon * (jk - 1))) - 1), &zli_index, 1);

            }
            __tmp15 = (zli_index > yrecldp_rlmin);
            if (__tmp15) {
                {
                    double zqx_index_16;
                    double zli_index_0;
                    double zliqfrac_slice;
                    double zliqfrac_index;
                    double zicefrac_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + ((jl + (klon * (jk - 1))) - 1), &zqx_index_16, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zli + ((jl + (klon * (jk - 1))) - 1), &zli_index_0, 1);
                    {
                        double __in1 = zqx_index_16;
                        double __in2 = zli_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zliqfrac_slice = __out;
                    }
                    {
                        double __inp = zliqfrac_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_414_16)
                        __out = __inp;
                        ///////////////////

                        zliqfrac[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zliqfrac + ((jl + (klon * (jk - 1))) - 1), &zliqfrac_index, 1);
                    {
                        double __in2 = zliqfrac_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        zicefrac_slice = __out;
                    }
                    {
                        double __inp = zicefrac_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_415_16)
                        __out = __inp;
                        ///////////////////

                        zicefrac[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                }
            } else {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_417_16)
                        __out = 0.0;
                        ///////////////////

                        zliqfrac[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_418_16)
                        __out = 0.0;
                        ///////////////////

                        zicefrac[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                }
            }

        }

    }
    kfdia_plus_1_13 = (kfdia + 1);
    for (jl = kidia; (jl < kfdia_plus_1_13); jl = (jl + 1)) {
        {

            {
                double __out;

                ///////////////////
                // Tasklet code (assign_428_8)
                __out = 0.0;
                ///////////////////

                zanewm1[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_429_8)
                __out = 0.0;
                ///////////////////

                zda[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_430_8)
                __out = 0.0;
                ///////////////////

                zcovpclr[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_431_8)
                __out = 0.0;
                ///////////////////

                zcovpmax[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_432_8)
                __out = 0.0;
                ///////////////////

                zcovptot[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_433_8)
                __out = 0.0;
                ///////////////////

                zcldtopdist[(jl - 1)] = __out;
            }

        }

    }
    for (jk = yrecldp_ncldtop; (jk < (klev + 1)); jk = (jk + 1)) {
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            kfdia_plus_1_14 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_14); jl = (jl + 1)) {
                {
                    double zqxfg_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqxfg_slice, 1);
                    {
                        double __inp = zqxfg_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_437_16)
                        __out = __inp;
                        ///////////////////

                        zqxfg[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }

        }
        kfdia_plus_1_15 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_15); jl = (jl + 1)) {
            {

                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_439_12)
                    __out = 0.0;
                    ///////////////////

                    zlicld[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_440_12)
                    __out = 0.0;
                    ///////////////////

                    zrainaut[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_441_12)
                    __out = 0.0;
                    ///////////////////

                    zrainacc[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_442_12)
                    __out = 0.0;
                    ///////////////////

                    zsnowaut[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_443_12)
                    __out = 0.0;
                    ///////////////////

                    zldefr[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_444_12)
                    __out = 0.0;
                    ///////////////////

                    zacust[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_445_12)
                    __out = 0.0;
                    ///////////////////

                    zqpretot[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_446_12)
                    __out = 0.0;
                    ///////////////////

                    zlfinalsum[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_447_12)
                    __out = 0.0;
                    ///////////////////

                    zlcond1[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_448_12)
                    __out = 0.0;
                    ///////////////////

                    zlcond2[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_449_12)
                    __out = 0.0;
                    ///////////////////

                    zsupsat[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_452_12)
                    __out = 0.0;
                    ///////////////////

                    zsolab[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_453_12)
                    __out = 0.0;
                    ///////////////////

                    zsolac[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_454_12)
                    __out = 0.0;
                    ///////////////////

                    zicetot[(jl - 1)] = __out;
                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {
            for (jn = 1; (jn < (5 + 1)); jn = (jn + 1)) {

                kfdia_plus_1_16 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_16); jl = (jl + 1)) {
                    {

                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_458_20)
                            __out = 0.0;
                            ///////////////////

                            zsolqb[(((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1)] = __out;
                        }
                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_459_20)
                            __out = 0.0;
                            ///////////////////

                            zsolqa[(((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1)] = __out;
                        }

                    }

                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            kfdia_plus_1_17 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_17); jl = (jl + 1)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_462_16)
                        __out = 0.0;
                        ///////////////////

                        zfallsrce[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_463_16)
                        __out = 0.0;
                        ///////////////////

                        zfallsink[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_464_16)
                        __out = 0.0;
                        ///////////////////

                        zconvsrce[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_465_16)
                        __out = 0.0;
                        ///////////////////

                        zconvsink[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_466_16)
                        __out = 0.0;
                        ///////////////////

                        zpsupsatsrce[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_467_16)
                        __out = 0.0;
                        ///////////////////

                        zratio[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }

        }
        kfdia_plus_1_18 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_18); jl = (jl + 1)) {
            {
                double *zgdp;
                zgdp = new double DACE_ALIGN(64)[klon];
                double paph_index_0;
                double paph_index_1;
                double zdp_slice;
                double zdp_index;
                double zgdp_slice;
                double pap_index_3;
                double ztp1_index_15;
                double ydcst_rd_times_ztp1_slice;
                double zrho_slice;
                double zgdp_index;
                double zdtgdp_slice;
                double zdp_index_0;
                double ptsphy_times_ydcst_rg;
                double __tmp21;
                double zrdtgdp_slice;
                double ztp1_index_16;
                double ztp1_slice_minus_ydthf_r4les_2;
                double ztp1_slice_ydthf_r4les_pow_2;
                double zfoeeliqt_index;
                double ydcst_retv_times_zfoeeliqt_slice;
                double __tmp22;
                double ztp1_index_17;
                double ztp1_slice_minus_ydthf_r4ies_1;
                double ztp1_slice_ydthf_r4ies_pow_2;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_3, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_15, 1);
                {
                    double __in1 = ydcst_rd;
                    double __in2 = ztp1_index_15;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_rd_times_ztp1_slice = __out;
                }
                {
                    double __in1 = pap_index_3;
                    double __in2 = ydcst_rd_times_ztp1_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zrho_slice = __out;
                }
                {
                    double __inp = zrho_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_471_12)
                    __out = __inp;
                    ///////////////////

                    zrho[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                paph + (((jk * klon) + jl) - 1), &paph_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                paph + ((jl + (klon * (jk - 1))) - 1), &paph_index_1, 1);
                {
                    double __in1 = paph_index_0;
                    double __in2 = paph_index_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zdp_slice = __out;
                }
                {
                    double __inp = zdp_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_469_12)
                    __out = __inp;
                    ///////////////////

                    zdp[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zdp + (jl - 1), &zdp_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zdp + (jl - 1), &zdp_index_0, 1);
                {
                    double __in1 = ydcst_rg;
                    double __in2 = zdp_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zgdp_slice = __out;
                }
                {
                    double __inp = zgdp_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_470_12)
                    __out = __inp;
                    ///////////////////

                    zgdp[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zgdp + (jl - 1), &zgdp_index, 1);
                {
                    double __in1 = ptsphy;
                    double __in2 = zgdp_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zdtgdp_slice = __out;
                }
                {
                    double __inp = zdtgdp_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_472_12)
                    __out = __inp;
                    ///////////////////

                    zdtgdp[(jl - 1)] = __out;
                }
                {
                    double __in1 = ptsphy;
                    double __in2 = ydcst_rg;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ptsphy_times_ydcst_rg = __out;
                }
                {
                    double __in2 = ptsphy_times_ydcst_rg;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    __tmp21 = __out;
                }
                {
                    double __in1 = zdp_index_0;
                    double __in2 = __tmp21;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zrdtgdp_slice = __out;
                }
                {
                    double __inp = zrdtgdp_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_473_12)
                    __out = __inp;
                    ///////////////////

                    zrdtgdp[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_16, 1);
                {
                    double __in1 = ztp1_index_16;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les_2 = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_ydthf_r4les_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    ztp1_slice_ydthf_r4les_pow_2 = __out;
                }
                {
                    double __in1 = ydthf_r5les;
                    double __in2 = ztp1_slice_ydthf_r4les_pow_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zfacw = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoeeliqt + ((jl + (klon * (jk - 1))) - 1), &zfoeeliqt_index, 1);
                {
                    double __in1 = ydcst_retv;
                    double __in2 = zfoeeliqt_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_retv_times_zfoeeliqt_slice = __out;
                }
                {
                    double __in2 = ydcst_retv_times_zfoeeliqt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp22 = __out;
                }
                {
                    double __in2 = __tmp22;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    zcor = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_17, 1);
                {
                    double __in1 = ztp1_index_17;
                    double __in2 = ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4ies_1 = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_ydthf_r4ies_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    ztp1_slice_ydthf_r4ies_pow_2 = __out;
                }
                {
                    double __in1 = ydthf_r5ies;
                    double __in2 = ztp1_slice_ydthf_r4ies_pow_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zfaci = __out;
                }
                delete[] zgdp;

            }
            {
                double *zdqsicedt;
                zdqsicedt = new double DACE_ALIGN(64)[klon];
                double zfoeew_index_2;
                double ydcst_retv_times_zfoeew_slice_0;
                double __tmp23;
                double zcor_0;
                double zfaci_times_zcor;
                double zqsice_index;
                double zdqsicedt_slice;
                double zdqsicedt_index;
                double ydthf_ralsdcp_times_zdqsicedt_slice;
                double zcorqsice_slice;
                double zalfaw_times_zfacw;
                double __tmp24;
                double __tmp25;
                double zfoeewmt_index;
                double ydcst_retv_times_zfoeewmt_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoeew + ((jl + (klon * (jk - 1))) - 1), &zfoeew_index_2, 1);
                {
                    double __in1 = ydcst_retv;
                    double __in2 = zfoeew_index_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_retv_times_zfoeew_slice_0 = __out;
                }
                {
                    double __in2 = ydcst_retv_times_zfoeew_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp23 = __out;
                }
                {
                    double __in2 = __tmp23;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    zcor_0 = __out;
                }
                {
                    double __inp = zcor_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_481_12)
                    __out = __inp;
                    ///////////////////

                    zcor = __out;
                }
                {
                    double __in1 = zfaci;
                    double __in2 = zcor;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zfaci_times_zcor = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index, 1);
                {
                    double __in1 = zfaci_times_zcor;
                    double __in2 = zqsice_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zdqsicedt_slice = __out;
                }
                {
                    double __inp = zdqsicedt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_482_12)
                    __out = __inp;
                    ///////////////////

                    zdqsicedt[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zdqsicedt + (jl - 1), &zdqsicedt_index, 1);
                {
                    double __in1 = ydthf_ralsdcp;
                    double __in2 = zdqsicedt_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_ralsdcp_times_zdqsicedt_slice = __out;
                }
                {
                    double __in2 = ydthf_ralsdcp_times_zdqsicedt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (1.0 + __in2);
                    ///////////////////

                    zcorqsice_slice = __out;
                }
                {
                    double __inp = zcorqsice_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_483_12)
                    __out = __inp;
                    ///////////////////

                    zcorqsice[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoealfa + ((jl + (klon * (jk - 1))) - 1), &zalfaw, 1);
                {
                    double __in2 = zalfaw;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp24 = __out;
                }
                {
                    double __in1 = __tmp24;
                    double __in2 = zfaci;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp25 = __out;
                }
                {
                    double __in1 = zalfaw;
                    double __in2 = zfacw;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zalfaw_times_zfacw = __out;
                }
                {
                    double __in1 = zalfaw_times_zfacw;
                    double __in2 = __tmp25;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zfac = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoeewmt + ((jl + (klon * (jk - 1))) - 1), &zfoeewmt_index, 1);
                {
                    double __in1 = ydcst_retv;
                    double __in2 = zfoeewmt_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_retv_times_zfoeewmt_slice = __out;
                }
                {
                    double __in2 = ydcst_retv_times_zfoeewmt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp26 = __out;
                }
                delete[] zdqsicedt;

            }
            {
                double *zcorqsmix;
                zcorqsmix = new double DACE_ALIGN(64)[klon];
                double zcor_1;
                double zfac_times_zcor;
                double zqsmix_index_1;
                double zdqsmixdt_slice;
                double ztp1_index_18;
                double min_ydthf_rtwat_ztp1_slice_2;
                double max_ydthf_rtice_expr_2;
                double expr_minus_ydthf_rtice_2;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_2;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_2;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_1;
                double expr_times_ydthf_ralvdcp;
                double ztp1_index_19;
                double min_ydthf_rtwat_ztp1_slice_3;
                double max_ydthf_rtice_expr_3;
                double expr_minus_ydthf_rtice_3;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_3;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_3;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_2;
                double __tmp27;
                double __tmp28;
                double expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp;
                double zdqsmixdt_index;
                double expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zdqsmixdt_slice;
                double zcorqsmix_slice;
                double zqsmix_index_2;
                double zqx_index_17;
                double zqsmix_slice_minus_zqx_slice;
                double zcorqsmix_index;
                double zqsmix_slice_zqx_slice_div_zcorqsmix_slice;
                double zevaplimmix_slice;
                double za_index_1;
                double max_za_slice_zepsec;
                double zqx_index_20;
                double zliqcld_slice;
                double zqx_index_21;
                double zicecld_slice;
                double zliqcld_index;
                double zicecld_index;
                double zlicld_slice;

                {
                    double __in2 = __tmp26;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    zcor_1 = __out;
                }
                {
                    double __inp = zcor_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_487_12)
                    __out = __inp;
                    ///////////////////

                    zcor = __out;
                }
                {
                    double __in1 = zfac;
                    double __in2 = zcor;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zfac_times_zcor = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_1, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_2, 1);
                {
                    double __in1 = zfac_times_zcor;
                    double __in2 = zqsmix_index_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zdqsmixdt_slice = __out;
                }
                {
                    double __inp = zdqsmixdt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_488_12)
                    __out = __inp;
                    ///////////////////

                    zdqsmixdt[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zdqsmixdt + (jl - 1), &zdqsmixdt_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_18, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_19, 1);
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_18;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_2 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_19;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_3 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_2 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_2;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_2 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_3 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_3;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_3 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_2;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_2 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_2 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_1 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_3;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_3 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_3 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_2 = __out;
                }
                {
                    double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp27 = __out;
                }
                {
                    double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_1;
                    double __in2 = ydthf_ralvdcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_ydthf_ralvdcp = __out;
                }
                {
                    double __in1 = __tmp27;
                    double __in2 = ydthf_ralsdcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp28 = __out;
                }
                {
                    double __in1 = expr_times_ydthf_ralvdcp;
                    double __in2 = __tmp28;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp = __out;
                }
                {
                    double __in1 = expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp;
                    double __in2 = zdqsmixdt_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zdqsmixdt_slice = __out;
                }
                {
                    double __in2 = expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zdqsmixdt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (1.0 + __in2);
                    ///////////////////

                    zcorqsmix_slice = __out;
                }
                {
                    double __inp = zcorqsmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_489_12)
                    __out = __inp;
                    ///////////////////

                    zcorqsmix[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zcorqsmix + (jl - 1), &zcorqsmix_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_17, 1);
                {
                    double __in1 = zqsmix_index_2;
                    double __in2 = zqx_index_17;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqsmix_slice_minus_zqx_slice = __out;
                }
                {
                    double __in1 = zqsmix_slice_minus_zqx_slice;
                    double __in2 = zcorqsmix_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zqsmix_slice_zqx_slice_div_zcorqsmix_slice = __out;
                }
                {
                    double __in_a = zqsmix_slice_zqx_slice_div_zcorqsmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, 0.0);
                    ///////////////////

                    zevaplimmix_slice = __out;
                }
                {
                    double __inp = zevaplimmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_490_12)
                    __out = __inp;
                    ///////////////////

                    zevaplimmix[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + ((jl + (klon * (jk - 1))) - 1), &zqx_index_20, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqx_index_21, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                za + ((jl + (klon * (jk - 1))) - 1), &za_index_1, 1);
                {
                    double __in_a = za_index_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, 1e-14);
                    ///////////////////

                    max_za_slice_zepsec = __out;
                }
                {
                    double __in2 = max_za_slice_zepsec;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    ztmpa = __out;
                }
                {
                    double __in1 = zqx_index_20;
                    double __in2 = ztmpa;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zliqcld_slice = __out;
                }
                {
                    double __inp = zliqcld_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_494_12)
                    __out = __inp;
                    ///////////////////

                    zliqcld[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zliqcld + (jl - 1), &zliqcld_index, 1);
                {
                    double __in1 = zqx_index_21;
                    double __in2 = ztmpa;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zicecld_slice = __out;
                }
                {
                    double __inp = zicecld_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_495_12)
                    __out = __inp;
                    ///////////////////

                    zicecld[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zicecld + (jl - 1), &zicecld_index, 1);
                {
                    double __in1 = zliqcld_index;
                    double __in2 = zicecld_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zlicld_slice = __out;
                }
                {
                    double __inp = zlicld_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_496_12)
                    __out = __inp;
                    ///////////////////

                    zlicld[(jl - 1)] = __out;
                }
                delete[] zcorqsmix;

            }

        }
        kfdia_plus_1_19 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_19); jl = (jl + 1)) {

            zqx_index_22 = zqx[((jl + (klon * (jk - 1))) - 1)];

            __tmp29 = (zqx_index_22 < yrecldp_rlmin);
            if (__tmp29) {
                {
                    double zsolqa_slice;
                    double zqx_index_23;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + ((jl + (klon * (jk - 1))) - 1), &zsolqa_slice, 1);
                    {
                        double __inp = zsolqa_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_499_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (4 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + ((jl + (klon * (jk - 1))) - 1), &zqx_index_23, 1);
                    {
                        double __in = zqx_index_23;
                        double __out;

                        ///////////////////
                        // Tasklet code (_USub_)
                        __out = (- __in);
                        ///////////////////

                        zsolqa_slice_0 = __out;
                    }

                }
                {

                    {
                        double __inp = zsolqa_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_500_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (20 * klon)) - 1)] = __out;
                    }

                }
            }
            zqx_index_24 = zqx[(((jl + (klev * klon)) + (klon * (jk - 1))) - 1)];

            __tmp30 = (zqx_index_24 < yrecldp_rlmin);
            if (__tmp30) {
                {
                    double zsolqa_slice_1;
                    double zqx_index_25;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zsolqa_slice_1, 1);
                    {
                        double __inp = zsolqa_slice_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_502_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (9 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqx_index_25, 1);
                    {
                        double __in = zqx_index_25;
                        double __out;

                        ///////////////////
                        // Tasklet code (_USub_)
                        __out = (- __in);
                        ///////////////////

                        zsolqa_slice_2 = __out;
                    }

                }
                {

                    {
                        double __inp = zsolqa_slice_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_503_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (21 * klon)) - 1)] = __out;
                    }

                }
            }

        }
        kfdia_plus_1_20 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_20); jl = (jl + 1)) {
            {
                double ztp1_index_20;
                double ydthf_rkoop2_times_ztp1_slice;
                double ydthf_rkoop1_minus_ydthf_rkoop2_ztp1_slice;
                double ztp1_index_21;
                double ztp1_slice_minus_ydcst_rtt_5;
                double ydthf_r3les_times_ztp1_slice_ydcst_rtt_2;
                double ztp1_index_22;
                double ztp1_slice_minus_ydthf_r4les_3;
                double ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_2;
                double exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_2;
                double ydthf_r2es_times_expr_2;
                double ztp1_index_23;
                double ztp1_slice_minus_ydcst_rtt_6;
                double ydthf_r3ies_times_ztp1_slice_ydcst_rtt_1;
                double ztp1_index_24;
                double ztp1_slice_minus_ydthf_r4ies_2;
                double ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_1;
                double exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_1;
                double ydthf_r2es_times_expr_3;
                double ydthf_r2es_expr_div_ydthf_r2es_expr;
                double zfokoop_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_20, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_21, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_22, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_23, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_24, 1);
                {
                    double __in1 = ydthf_rkoop2;
                    double __in2 = ztp1_index_20;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_rkoop2_times_ztp1_slice = __out;
                }
                {
                    double __in1 = ydthf_rkoop1;
                    double __in2 = ydthf_rkoop2_times_ztp1_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ydthf_rkoop1_minus_ydthf_rkoop2_ztp1_slice = __out;
                }
                {
                    double __in1 = ztp1_index_21;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_5 = __out;
                }
                {
                    double __in1 = ztp1_index_23;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_6 = __out;
                }
                {
                    double __in1 = ydthf_r3les;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3les_times_ztp1_slice_ydcst_rtt_2 = __out;
                }
                {
                    double __in1 = ztp1_index_22;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les_3 = __out;
                }
                {
                    double __in1 = ydthf_r3les_times_ztp1_slice_ydcst_rtt_2;
                    double __in2 = ztp1_slice_minus_ydthf_r4les_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_2 = __out;
                }
                {
                    double __in1 = ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_2 = __out;
                }
                {
                    double __in1 = ydthf_r2es;
                    double __in2 = exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r2es_times_expr_2 = __out;
                }
                {
                    double __in1 = ydthf_r3ies;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3ies_times_ztp1_slice_ydcst_rtt_1 = __out;
                }
                {
                    double __in1 = ztp1_index_24;
                    double __in2 = ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4ies_2 = __out;
                }
                {
                    double __in1 = ydthf_r3ies_times_ztp1_slice_ydcst_rtt_1;
                    double __in2 = ztp1_slice_minus_ydthf_r4ies_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_1 = __out;
                }
                {
                    double __in1 = ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_1 = __out;
                }
                {
                    double __in1 = ydthf_r2es;
                    double __in2 = exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r2es_times_expr_3 = __out;
                }
                {
                    double __in1 = ydthf_r2es_times_expr_2;
                    double __in2 = ydthf_r2es_times_expr_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r2es_expr_div_ydthf_r2es_expr = __out;
                }
                {
                    double __in_a = ydthf_rkoop1_minus_ydthf_rkoop2_ztp1_slice;
                    double __in_b = ydthf_r2es_expr_div_ydthf_r2es_expr;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    zfokoop_slice = __out;
                }
                {
                    double __inp = zfokoop_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_505_12)
                    __out = __inp;
                    ///////////////////

                    zfokoop[(jl - 1)] = __out;
                }

            }

        }
        kfdia_plus_1_21 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_21); jl = (jl + 1)) {

            ztp1_index_25 = ztp1[((jl + (klon * (jk - 1))) - 1)];

            __tmp31 = (ztp1_index_25 >= ydcst_rtt);
            __tmp32 = (yrecldp_nssopt == 0);
            if ((__tmp31 || __tmp32)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_508_16)
                        __out = 1.0;
                        ///////////////////

                        zfac = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_509_16)
                        __out = 1.0;
                        ///////////////////

                        zfaci = __out;
                    }

                }
            } else {
                {
                    double za_index_2;
                    double zfokoop_index;
                    double za_index_3;
                    double __tmp34;
                    double zfokoop_slice_times_1_0_za_slice;
                    double zfac_0;
                    double zfaci_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_2, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_3, 1);
                    {
                        double __in2 = za_index_3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp34 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zfokoop + (jl - 1), &zfokoop_index, 1);
                    {
                        double __in1 = zfokoop_index;
                        double __in2 = __tmp34;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zfokoop_slice_times_1_0_za_slice = __out;
                    }
                    {
                        double __in1 = za_index_2;
                        double __in2 = zfokoop_slice_times_1_0_za_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zfac_0 = __out;
                    }
                    {
                        double __inp = zfac_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_511_16)
                        __out = __inp;
                        ///////////////////

                        zfac = __out;
                    }
                    {
                        double __in1 = ptsphy;
                        double __in2 = yrecldp_rkooptau;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zfaci_0 = __out;
                    }
                    {
                        double __inp = zfaci_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_512_16)
                        __out = __inp;
                        ///////////////////

                        zfaci = __out;
                    }

                }
            }
            za_index_4 = za[((jl + (klon * (jk - 1))) - 1)];
            __tmp35 = (1.0 - yrecldp_ramin);
            if ((za_index_4 > __tmp35)) {
                {
                    double zqx_index_26;
                    double zqsice_index_1;
                    double zfac_times_zqsice_slice;
                    double zqx_slice_minus_zfac_zqsice_slice;
                    double zcorqsice_index_0;
                    double zqx_slice_zfac_zqsice_slice_div_zcorqsice_slice;
                    double zsupsat_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_26, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_1, 1);
                    {
                        double __in1 = zfac;
                        double __in2 = zqsice_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zfac_times_zqsice_slice = __out;
                    }
                    {
                        double __in1 = zqx_index_26;
                        double __in2 = zfac_times_zqsice_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqx_slice_minus_zfac_zqsice_slice = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcorqsice + (jl - 1), &zcorqsice_index_0, 1);
                    {
                        double __in1 = zqx_slice_minus_zfac_zqsice_slice;
                        double __in2 = zcorqsice_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zqx_slice_zfac_zqsice_slice_div_zcorqsice_slice = __out;
                    }
                    {
                        double __in_a = zqx_slice_zfac_zqsice_slice_div_zcorqsice_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        zsupsat_slice = __out;
                    }
                    {
                        double __inp = zsupsat_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_514_16)
                        __out = __inp;
                        ///////////////////

                        zsupsat[(jl - 1)] = __out;
                    }

                }
            } else {
                {
                    double zqx_index_27;
                    double za_index_5;
                    double zqsice_index_2;
                    double za_slice_times_zqsice_slice;
                    double zqx_slice_minus_za_slice_zqsice_slice;
                    double za_index_6;
                    double __tmp37;
                    double max_1_0_za_slice_zepsilon;
                    double zqp1env;
                    double za_index_7;
                    double __tmp38;
                    double zqsice_index_3;
                    double zfac_times_zqsice_slice_0;
                    double zqp1env_minus_zfac_zqsice_slice;
                    double __tmp39;
                    double zcorqsice_index_1;
                    double __tmp40;
                    double zsupsat_slice_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_27, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_5, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_6, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_7, 1);
                    {
                        double __in2 = za_index_6;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp37 = __out;
                    }
                    {
                        double __in2 = za_index_7;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp38 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_2, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_3, 1);
                    {
                        double __in1 = za_index_5;
                        double __in2 = zqsice_index_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        za_slice_times_zqsice_slice = __out;
                    }
                    {
                        double __in1 = zqx_index_27;
                        double __in2 = za_slice_times_zqsice_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqx_slice_minus_za_slice_zqsice_slice = __out;
                    }
                    {
                        double __in_a = __tmp37;
                        double __in_b = zepsilon;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        max_1_0_za_slice_zepsilon = __out;
                    }
                    {
                        double __in1 = zqx_slice_minus_za_slice_zqsice_slice;
                        double __in2 = max_1_0_za_slice_zepsilon;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zqp1env = __out;
                    }
                    {
                        double __in1 = zfac;
                        double __in2 = zqsice_index_3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zfac_times_zqsice_slice_0 = __out;
                    }
                    {
                        double __in1 = zqp1env;
                        double __in2 = zfac_times_zqsice_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqp1env_minus_zfac_zqsice_slice = __out;
                    }
                    {
                        double __in1 = __tmp38;
                        double __in2 = zqp1env_minus_zfac_zqsice_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp39 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcorqsice + (jl - 1), &zcorqsice_index_1, 1);
                    {
                        double __in1 = __tmp39;
                        double __in2 = zcorqsice_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        __tmp40 = __out;
                    }
                    {
                        double __in_a = __tmp40;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        zsupsat_slice_0 = __out;
                    }
                    {
                        double __inp = zsupsat_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_517_16)
                        __out = __inp;
                        ///////////////////

                        zsupsat[(jl - 1)] = __out;
                    }

                }
            }
            zsupsat_index = zsupsat[(jl - 1)];
            if ((zsupsat_index > 1e-14)) {

                ztp1_index_26 = ztp1[((jl + (klon * (jk - 1))) - 1)];

                __tmp42 = (ztp1_index_26 > yrecldp_rthomo);
                if (__tmp42) {
                    {
                        double zsolqa_index;
                        double zsupsat_index_0;
                        double zsolqa_slice_plus_zsupsat_slice;
                        double zsolqa_index_0;
                        double zsupsat_index_1;
                        double zsolqa_slice_minus_zsupsat_slice;
                        double zqxfg_index;
                        double zsupsat_index_2;
                        double zqxfg_slice_plus_zsupsat_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (20 * klon)) - 1), &zsolqa_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsupsat + (jl - 1), &zsupsat_index_0, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsupsat + (jl - 1), &zsupsat_index_1, 1);
                        {
                            double __in1 = zsolqa_index;
                            double __in2 = zsupsat_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zsupsat_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zsupsat_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_520_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (20 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (4 * klon)) - 1), &zsolqa_index_0, 1);
                        {
                            double __in1 = zsolqa_index_0;
                            double __in2 = zsupsat_index_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsolqa_slice_minus_zsupsat_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_zsupsat_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_521_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (4 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsupsat + (jl - 1), &zsupsat_index_2, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxfg + (jl - 1), &zqxfg_index, 1);
                        {
                            double __in1 = zqxfg_index;
                            double __in2 = zsupsat_index_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqxfg_slice_plus_zsupsat_slice = __out;
                        }
                        {
                            double __inp = zqxfg_slice_plus_zsupsat_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_522_20)
                            __out = __inp;
                            ///////////////////

                            zqxfg[(jl - 1)] = __out;
                        }

                    }
                } else {
                    {
                        double zsolqa_index_1;
                        double zsupsat_index_3;
                        double zsolqa_slice_plus_zsupsat_slice_0;
                        double zsolqa_index_2;
                        double zsupsat_index_4;
                        double zsolqa_slice_minus_zsupsat_slice_0;
                        double zqxfg_index_0;
                        double zsupsat_index_5;
                        double zqxfg_slice_plus_zsupsat_slice_0;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (21 * klon)) - 1), &zsolqa_index_1, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsupsat + (jl - 1), &zsupsat_index_3, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsupsat + (jl - 1), &zsupsat_index_4, 1);
                        {
                            double __in1 = zsolqa_index_1;
                            double __in2 = zsupsat_index_3;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zsupsat_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zsupsat_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_524_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (21 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (9 * klon)) - 1), &zsolqa_index_2, 1);
                        {
                            double __in1 = zsolqa_index_2;
                            double __in2 = zsupsat_index_4;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsolqa_slice_minus_zsupsat_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_zsupsat_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_525_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (9 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsupsat + (jl - 1), &zsupsat_index_5, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxfg + ((jl + klon) - 1), &zqxfg_index_0, 1);
                        {
                            double __in1 = zqxfg_index_0;
                            double __in2 = zsupsat_index_5;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqxfg_slice_plus_zsupsat_slice_0 = __out;
                        }
                        {
                            double __inp = zqxfg_slice_plus_zsupsat_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_526_20)
                            __out = __inp;
                            ///////////////////

                            zqxfg[((jl + klon) - 1)] = __out;
                        }

                    }
                }
                {
                    double za_index_8;
                    double __tmp43;
                    double zsolac_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_8, 1);
                    {
                        double __in2 = za_index_8;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp43 = __out;
                    }
                    {
                        double __in1 = __tmp43;
                        double __in2 = zfaci;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zsolac_slice = __out;
                    }
                    {
                        double __inp = zsolac_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_527_16)
                        __out = __inp;
                        ///////////////////

                        zsolac[(jl - 1)] = __out;
                    }

                }
            }
            psupsat_index = psupsat[((jl + (klon * (jk - 1))) - 1)];
            if ((psupsat_index > 1e-14)) {

                ztp1_index_27 = ztp1[((jl + (klon * (jk - 1))) - 1)];

                __tmp45 = (ztp1_index_27 > yrecldp_rthomo);
                if (__tmp45) {
                    {
                        double zsolqa_index_3;
                        double psupsat_index_0;
                        double zsolqa_slice_plus_psupsat_slice;
                        double zpsupsatsrce_slice;
                        double zqxfg_index_1;
                        double psupsat_index_1;
                        double zqxfg_slice_plus_psupsat_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + (jl - 1), &zsolqa_index_3, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        psupsat + ((jl + (klon * (jk - 1))) - 1), &psupsat_index_0, 1);
                        {
                            double __in1 = zsolqa_index_3;
                            double __in2 = psupsat_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_psupsat_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_psupsat_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_530_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        psupsat + ((jl + (klon * (jk - 1))) - 1), &zpsupsatsrce_slice, 1);
                        {
                            double __inp = zpsupsatsrce_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_531_20)
                            __out = __inp;
                            ///////////////////

                            zpsupsatsrce[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        psupsat + ((jl + (klon * (jk - 1))) - 1), &psupsat_index_1, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxfg + (jl - 1), &zqxfg_index_1, 1);
                        {
                            double __in1 = zqxfg_index_1;
                            double __in2 = psupsat_index_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqxfg_slice_plus_psupsat_slice = __out;
                        }
                        {
                            double __inp = zqxfg_slice_plus_psupsat_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_532_20)
                            __out = __inp;
                            ///////////////////

                            zqxfg[(jl - 1)] = __out;
                        }

                    }
                } else {
                    {
                        double zsolqa_index_4;
                        double psupsat_index_2;
                        double zsolqa_slice_plus_psupsat_slice_0;
                        double zpsupsatsrce_slice_0;
                        double zqxfg_index_2;
                        double psupsat_index_3;
                        double zqxfg_slice_plus_psupsat_slice_0;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (6 * klon)) - 1), &zsolqa_index_4, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        psupsat + ((jl + (klon * (jk - 1))) - 1), &psupsat_index_2, 1);
                        {
                            double __in1 = zsolqa_index_4;
                            double __in2 = psupsat_index_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_psupsat_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_psupsat_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_534_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (6 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        psupsat + ((jl + (klon * (jk - 1))) - 1), &zpsupsatsrce_slice_0, 1);
                        {
                            double __inp = zpsupsatsrce_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_535_20)
                            __out = __inp;
                            ///////////////////

                            zpsupsatsrce[((jl + klon) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        psupsat + ((jl + (klon * (jk - 1))) - 1), &psupsat_index_3, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxfg + ((jl + klon) - 1), &zqxfg_index_2, 1);
                        {
                            double __in1 = zqxfg_index_2;
                            double __in2 = psupsat_index_3;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqxfg_slice_plus_psupsat_slice_0 = __out;
                        }
                        {
                            double __inp = zqxfg_slice_plus_psupsat_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_536_20)
                            __out = __inp;
                            ///////////////////

                            zqxfg[((jl + klon) - 1)] = __out;
                        }

                    }
                }
                {
                    double za_index_9;
                    double __tmp46;
                    double zsolac_slice_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_9, 1);
                    {
                        double __in2 = za_index_9;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp46 = __out;
                    }
                    {
                        double __in1 = __tmp46;
                        double __in2 = zfaci;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zsolac_slice_0 = __out;
                    }
                    {
                        double __inp = zsolac_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_537_16)
                        __out = __inp;
                        ///////////////////

                        zsolac[(jl - 1)] = __out;
                    }

                }
            }

        }
        if (((jk < klev) && (jk >= yrecldp_ncldtop))) {

            kfdia_plus_1_22 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_22); jl = (jl + 1)) {
                {
                    double plude_index;
                    double zdtgdp_index;
                    double plude_slice_times_zdtgdp_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    plude + ((jl + (klon * (jk - 1))) - 1), &plude_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdtgdp + (jl - 1), &zdtgdp_index, 1);
                    {
                        double __in1 = plude_index;
                        double __in2 = zdtgdp_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        plude_slice_times_zdtgdp_slice = __out;
                    }
                    {
                        double __inp = plude_slice_times_zdtgdp_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_540_16)
                        __out = __inp;
                        ///////////////////

                        plude[((jl + (klon * (jk - 1))) - 1)] = __out;
                    }

                }
                ldcum_index = ldcum[(jl - 1)];
                plude_index_0 = plude[((jl + (klon * (jk - 1))) - 1)];

                __tmp47 = (plude_index_0 > yrecldp_rlmin);
                plu_index = plu[(((jk * klon) + jl) - 1)];
                if (((ldcum_index && __tmp47) && (plu_index > 1e-14))) {
                    {
                        double zsolac_index;
                        double plude_index_1;
                        double plu_index_0;
                        double plude_slice_div_plu_slice;
                        double zsolac_slice_plus_plude_slice_plu_slice;
                        double zalfaw_0;
                        double plude_index_2;
                        double zconvsrce_slice;
                        double __tmp51;
                        double plude_index_3;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolac + (jl - 1), &zsolac_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        plude + ((jl + (klon * (jk - 1))) - 1), &plude_index_1, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        plu + (((jk * klon) + jl) - 1), &plu_index_0, 1);
                        {
                            double __in1 = plude_index_1;
                            double __in2 = plu_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / __in2);
                            ///////////////////

                            plude_slice_div_plu_slice = __out;
                        }
                        {
                            double __in1 = zsolac_index;
                            double __in2 = plude_slice_div_plu_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolac_slice_plus_plude_slice_plu_slice = __out;
                        }
                        {
                            double __inp = zsolac_slice_plus_plude_slice_plu_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_542_20)
                            __out = __inp;
                            ///////////////////

                            zsolac[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        plude + ((jl + (klon * (jk - 1))) - 1), &plude_index_2, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        plude + ((jl + (klon * (jk - 1))) - 1), &plude_index_3, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zfoealfa + ((jl + (klon * (jk - 1))) - 1), &zalfaw_0, 1);
                        {
                            double __inp = zalfaw_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_543_20)
                            __out = __inp;
                            ///////////////////

                            zalfaw = __out;
                        }
                        {
                            double __in1 = zalfaw;
                            double __in2 = plude_index_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zconvsrce_slice = __out;
                        }
                        {
                            double __inp = zconvsrce_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_544_20)
                            __out = __inp;
                            ///////////////////

                            zconvsrce[(jl - 1)] = __out;
                        }
                        {
                            double __in2 = zalfaw;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp51 = __out;
                        }
                        {
                            double __in1 = __tmp51;
                            double __in2 = plude_index_3;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zconvsrce_slice_0 = __out;
                        }

                    }
                    {
                        double zsolqa_index_5;
                        double zconvsrce_index;
                        double zsolqa_slice_plus_zconvsrce_slice;
                        double zsolqa_index_6;
                        double zconvsrce_index_0;
                        double zsolqa_slice_plus_zconvsrce_slice_0;

                        {
                            double __inp = zconvsrce_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_545_20)
                            __out = __inp;
                            ///////////////////

                            zconvsrce[((jl + klon) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zconvsrce + (jl - 1), &zconvsrce_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zconvsrce + ((jl + klon) - 1), &zconvsrce_index_0, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + (jl - 1), &zsolqa_index_5, 1);
                        {
                            double __in1 = zsolqa_index_5;
                            double __in2 = zconvsrce_index;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zconvsrce_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zconvsrce_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_546_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (6 * klon)) - 1), &zsolqa_index_6, 1);
                        {
                            double __in1 = zsolqa_index_6;
                            double __in2 = zconvsrce_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zconvsrce_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zconvsrce_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_547_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (6 * klon)) - 1)] = __out;
                        }

                    }
                } else {
                    {

                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_549_20)
                            __out = 0.0;
                            ///////////////////

                            plude[((jl + (klon * (jk - 1))) - 1)] = __out;
                        }

                    }
                }
                ldcum_index_0 = ldcum[(jl - 1)];
                if (ldcum_index_0) {
                    {
                        double zsolqa_index_7;
                        double psnde_index;
                        double zdtgdp_index_0;
                        double psnde_slice_times_zdtgdp_slice;
                        double zsolqa_slice_plus_psnde_slice_zdtgdp_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (18 * klon)) - 1), &zsolqa_index_7, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        psnde + ((jl + (klon * (jk - 1))) - 1), &psnde_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zdtgdp + (jl - 1), &zdtgdp_index_0, 1);
                        {
                            double __in1 = psnde_index;
                            double __in2 = zdtgdp_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            psnde_slice_times_zdtgdp_slice = __out;
                        }
                        {
                            double __in1 = zsolqa_index_7;
                            double __in2 = psnde_slice_times_zdtgdp_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_psnde_slice_zdtgdp_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_psnde_slice_zdtgdp_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_551_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (18 * klon)) - 1)] = __out;
                        }

                    }
                }

            }
        }
        if ((jk > yrecldp_ncldtop)) {

            kfdia_plus_1_23 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_23); jl = (jl + 1)) {
                {
                    double pmfu_index;
                    double pmfd_index;
                    double pmfu_slice_plus_pmfd_slice;
                    double zdtgdp_index_1;
                    double pmfu_slice_pmfd_slice_times_zdtgdp_slice;
                    double zmf_slice;
                    double zmf_index;
                    double zanewm1_index;
                    double zacust_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pmfu + ((jl + (klon * (jk - 1))) - 1), &pmfu_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pmfd + ((jl + (klon * (jk - 1))) - 1), &pmfd_index, 1);
                    {
                        double __in1 = pmfu_index;
                        double __in2 = pmfd_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        pmfu_slice_plus_pmfd_slice = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdtgdp + (jl - 1), &zdtgdp_index_1, 1);
                    {
                        double __in1 = pmfu_slice_plus_pmfd_slice;
                        double __in2 = zdtgdp_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        pmfu_slice_pmfd_slice_times_zdtgdp_slice = __out;
                    }
                    {
                        double __in_b = pmfu_slice_pmfd_slice_times_zdtgdp_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(0.0, __in_b);
                        ///////////////////

                        zmf_slice = __out;
                    }
                    {
                        double __inp = zmf_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_554_16)
                        __out = __inp;
                        ///////////////////

                        zmf[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zmf + (jl - 1), &zmf_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zanewm1 + (jl - 1), &zanewm1_index, 1);
                    {
                        double __in1 = zmf_index;
                        double __in2 = zanewm1_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zacust_slice = __out;
                    }
                    {
                        double __inp = zacust_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_555_16)
                        __out = __inp;
                        ///////////////////

                        zacust[(jl - 1)] = __out;
                    }

                }

            }
            for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

                llfall_index = llfall[(jm - 1)];
                iphase_index_1 = iphase[(jm - 1)];
                if (((! llfall_index) && (iphase_index_1 > 0))) {

                    kfdia_plus_1_24 = (kfdia + 1);
                    for (jl = kidia; (jl < kfdia_plus_1_24); jl = (jl + 1)) {
                        {
                            double zmf_index_0;
                            double zqxnm1_index;
                            double zlcust_slice;
                            double zconvsrce_index_1;
                            double zlcust_index;
                            double zconvsrce_slice_plus_zlcust_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zmf + (jl - 1), &zmf_index_0, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxnm1 + ((jl + (klon * (jm - 1))) - 1), &zqxnm1_index, 1);
                            {
                                double __in1 = zmf_index_0;
                                double __in2 = zqxnm1_index;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zlcust_slice = __out;
                            }
                            {
                                double __inp = zlcust_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_559_24)
                                __out = __inp;
                                ///////////////////

                                zlcust[((jl + (klon * (jm - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcust + ((jl + (klon * (jm - 1))) - 1), &zlcust_index, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zconvsrce + ((jl + (klon * (jm - 1))) - 1), &zconvsrce_index_1, 1);
                            {
                                double __in1 = zconvsrce_index_1;
                                double __in2 = zlcust_index;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zconvsrce_slice_plus_zlcust_slice = __out;
                            }
                            {
                                double __inp = zconvsrce_slice_plus_zlcust_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_560_24)
                                __out = __inp;
                                ///////////////////

                                zconvsrce[((jl + (klon * (jm - 1))) - 1)] = __out;
                            }

                        }

                    }
                }

            }
            kfdia_plus_1_25 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_25); jl = (jl + 1)) {
                {
                    double zrdcp_times_0_5;
                    double ztp1_index_28;
                    double ztp1_index_29;
                    double ztp1_slice_plus_ztp1_slice;
                    double zrdcp_0_5_times_ztp1_slice_ztp1_slice;
                    double paph_index_2;
                    double pap_index_6;
                    double pap_index_7;
                    double pap_slice_minus_pap_slice_0;
                    double zanewm1_index_0;
                    double zanewm1_slice_times_zdtforc;
                    double zdqsmixdt_index_0;
                    double zdqs_slice;

                    {
                        double __in1 = zrdcp;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * 0.5);
                        ///////////////////

                        zrdcp_times_0_5 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 2))) - 1), &ztp1_index_28, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_29, 1);
                    {
                        double __in1 = ztp1_index_28;
                        double __in2 = ztp1_index_29;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        ztp1_slice_plus_ztp1_slice = __out;
                    }
                    {
                        double __in1 = zrdcp_times_0_5;
                        double __in2 = ztp1_slice_plus_ztp1_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zrdcp_0_5_times_ztp1_slice_ztp1_slice = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    paph + ((jl + (klon * (jk - 1))) - 1), &paph_index_2, 1);
                    {
                        double __in1 = zrdcp_0_5_times_ztp1_slice_ztp1_slice;
                        double __in2 = paph_index_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zdtdp = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_6, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pap + ((jl + (klon * (jk - 2))) - 1), &pap_index_7, 1);
                    {
                        double __in1 = pap_index_6;
                        double __in2 = pap_index_7;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        pap_slice_minus_pap_slice_0 = __out;
                    }
                    {
                        double __in1 = zdtdp;
                        double __in2 = pap_slice_minus_pap_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zdtforc = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zanewm1 + (jl - 1), &zanewm1_index_0, 1);
                    {
                        double __in1 = zanewm1_index_0;
                        double __in2 = zdtforc;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zanewm1_slice_times_zdtforc = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdqsmixdt + (jl - 1), &zdqsmixdt_index_0, 1);
                    {
                        double __in1 = zanewm1_slice_times_zdtforc;
                        double __in2 = zdqsmixdt_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zdqs_slice = __out;
                    }
                    {
                        double __inp = zdqs_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_564_16)
                        __out = __inp;
                        ///////////////////

                        zdqs[(jl - 1)] = __out;
                    }

                }

            }
            for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

                llfall_index_0 = llfall[(jm - 1)];
                iphase_index_2 = iphase[(jm - 1)];
                if (((! llfall_index_0) && (iphase_index_2 > 0))) {

                    kfdia_plus_1_26 = (kfdia + 1);
                    for (jl = kidia; (jl < kfdia_plus_1_26); jl = (jl + 1)) {
                        {
                            double zlcust_index_0;
                            double zdqs_index;
                            double zlcust_slice_minus_zdqs_slice;
                            double zlfinal;
                            double zlcust_index_1;
                            double zlcust_slice_minus_zlfinal;
                            double zevaplimmix_index;
                            double zlcust_index_2;
                            double zlfinal_0;
                            double zlfinalsum_index;
                            double zlfinalsum_slice_plus_zlfinal;
                            double zsolqa_index_8;
                            double zlcust_index_3;
                            double zsolqa_slice_plus_zlcust_slice;
                            double zsolqa_index_9;
                            double zsolqa_slice_plus_zevap;
                            double zsolqa_index_10;
                            double zsolqa_slice_minus_zevap;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcust + ((jl + (klon * (jm - 1))) - 1), &zlcust_index_0, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcust + ((jl + (klon * (jm - 1))) - 1), &zlcust_index_1, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcust + ((jl + (klon * (jm - 1))) - 1), &zlcust_index_2, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcust + ((jl + (klon * (jm - 1))) - 1), &zlcust_index_3, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zdqs + (jl - 1), &zdqs_index, 1);
                            {
                                double __in1 = zlcust_index_0;
                                double __in2 = zdqs_index;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zlcust_slice_minus_zdqs_slice = __out;
                            }
                            {
                                double __in_b = zlcust_slice_minus_zdqs_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (__max2)
                                __out = max(0.0, __in_b);
                                ///////////////////

                                zlfinal = __out;
                            }
                            {
                                double __in1 = zlcust_index_1;
                                double __in2 = zlfinal;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zlcust_slice_minus_zlfinal = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zevaplimmix + (jl - 1), &zevaplimmix_index, 1);
                            {
                                double __in_a = zlcust_slice_minus_zlfinal;
                                double __in_b = zevaplimmix_index;
                                double __out;

                                ///////////////////
                                // Tasklet code (__min2)
                                __out = min(__in_a, __in_b);
                                ///////////////////

                                zevap = __out;
                            }
                            {
                                double __in1 = zlcust_index_2;
                                double __in2 = zevap;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zlfinal_0 = __out;
                            }
                            {
                                double __inp = zlfinal_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_570_24)
                                __out = __inp;
                                ///////////////////

                                zlfinal = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlfinalsum + (jl - 1), &zlfinalsum_index, 1);
                            {
                                double __in1 = zlfinalsum_index;
                                double __in2 = zlfinal;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zlfinalsum_slice_plus_zlfinal = __out;
                            }
                            {
                                double __inp = zlfinalsum_slice_plus_zlfinal;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_571_24)
                                __out = __inp;
                                ///////////////////

                                zlfinalsum[(jl - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + ((jl + ((6 * klon) * (jm - 1))) - 1), &zsolqa_index_8, 1);
                            {
                                double __in1 = zsolqa_index_8;
                                double __in2 = zlcust_index_3;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zsolqa_slice_plus_zlcust_slice = __out;
                            }
                            {
                                double __inp = zsolqa_slice_plus_zlcust_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_572_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[((jl + ((6 * klon) * (jm - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + (((jl + ((5 * klon) * (jm - 1))) + (4 * klon)) - 1), &zsolqa_index_9, 1);
                            {
                                double __in1 = zsolqa_index_9;
                                double __in2 = zevap;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zsolqa_slice_plus_zevap = __out;
                            }
                            {
                                double __inp = zsolqa_slice_plus_zevap;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_573_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[(((jl + ((5 * klon) * (jm - 1))) + (4 * klon)) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + (((jl + (klon * (jm - 1))) + (20 * klon)) - 1), &zsolqa_index_10, 1);
                            {
                                double __in1 = zsolqa_index_10;
                                double __in2 = zevap;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zsolqa_slice_minus_zevap = __out;
                            }
                            {
                                double __inp = zsolqa_slice_minus_zevap;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_574_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[(((jl + (klon * (jm - 1))) + (20 * klon)) - 1)] = __out;
                            }

                        }

                    }
                }

            }
            kfdia_plus_1_27 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_27); jl = (jl + 1)) {

                zlfinalsum_index_0 = zlfinalsum[(jl - 1)];
                if ((zlfinalsum_index_0 < 1e-14)) {
                    {

                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_577_20)
                            __out = 0.0;
                            ///////////////////

                            zacust[(jl - 1)] = __out;
                        }

                    }
                }
                {
                    double zsolac_index_0;
                    double zacust_index;
                    double zsolac_slice_plus_zacust_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolac + (jl - 1), &zsolac_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zacust + (jl - 1), &zacust_index, 1);
                    {
                        double __in1 = zsolac_index_0;
                        double __in2 = zacust_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolac_slice_plus_zacust_slice = __out;
                    }
                    {
                        double __inp = zsolac_slice_plus_zacust_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_578_16)
                        __out = __inp;
                        ///////////////////

                        zsolac[(jl - 1)] = __out;
                    }

                }

            }
        }
        kfdia_plus_1_28 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_28); jl = (jl + 1)) {
            if ((jk < klev)) {
                {
                    double pmfu_index_0;
                    double pmfd_index_0;
                    double pmfu_slice_plus_pmfd_slice_0;
                    double zdtgdp_index_2;
                    double pmfu_slice_pmfd_slice_times_zdtgdp_slice_0;
                    double zsolab_index;
                    double zsolab_slice_plus_zmfdn;
                    double zsolqb_index;
                    double zsolqb_slice_plus_zmfdn;
                    double zsolqb_index_0;
                    double zsolqb_slice_plus_zmfdn_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pmfu + (((jk * klon) + jl) - 1), &pmfu_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pmfd + (((jk * klon) + jl) - 1), &pmfd_index_0, 1);
                    {
                        double __in1 = pmfu_index_0;
                        double __in2 = pmfd_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        pmfu_slice_plus_pmfd_slice_0 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdtgdp + (jl - 1), &zdtgdp_index_2, 1);
                    {
                        double __in1 = pmfu_slice_plus_pmfd_slice_0;
                        double __in2 = zdtgdp_index_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        pmfu_slice_pmfd_slice_times_zdtgdp_slice_0 = __out;
                    }
                    {
                        double __in_b = pmfu_slice_pmfd_slice_times_zdtgdp_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(0.0, __in_b);
                        ///////////////////

                        zmfdn = __out;
                    }
                    {
                        double __inp = zmfdn;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_585_16)
                        __out = __inp;
                        ///////////////////

                        zconvsink[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolab + (jl - 1), &zsolab_index, 1);
                    {
                        double __in1 = zsolab_index;
                        double __in2 = zmfdn;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolab_slice_plus_zmfdn = __out;
                    }
                    {
                        double __inp = zsolab_slice_plus_zmfdn;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_582_16)
                        __out = __inp;
                        ///////////////////

                        zsolab[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqb + (jl - 1), &zsolqb_index, 1);
                    {
                        double __in1 = zsolqb_index;
                        double __in2 = zmfdn;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqb_slice_plus_zmfdn = __out;
                    }
                    {
                        double __inp = zsolqb_slice_plus_zmfdn;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_583_16)
                        __out = __inp;
                        ///////////////////

                        zsolqb[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqb + ((jl + (6 * klon)) - 1), &zsolqb_index_0, 1);
                    {
                        double __in1 = zsolqb_index_0;
                        double __in2 = zmfdn;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqb_slice_plus_zmfdn_0 = __out;
                    }
                    {
                        double __inp = zsolqb_slice_plus_zmfdn_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_584_16)
                        __out = __inp;
                        ///////////////////

                        zsolqb[((jl + (6 * klon)) - 1)] = __out;
                    }

                }
                {

                    {
                        double __inp = zmfdn;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_586_16)
                        __out = __inp;
                        ///////////////////

                        zconvsink[((jl + klon) - 1)] = __out;
                    }

                }
            }

        }
        kfdia_plus_1_29 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_29); jl = (jl + 1)) {
            {
                double zldifdt_slice;

                {
                    double __in1 = yrecldp_rcldiff;
                    double __in2 = ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zldifdt_slice = __out;
                }
                {
                    double __inp = zldifdt_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_588_12)
                    __out = __inp;
                    ///////////////////

                    zldifdt[(jl - 1)] = __out;
                }

            }
            ktype_index = ktype[(jl - 1)];
            plude_index_4 = plude[((jl + (klon * (jk - 1))) - 1)];
            if (((ktype_index > 0) && (plude_index_4 > 1e-14))) {
                {
                    double zldifdt_index;
                    double yrecldp_rcldiff_convi_times_zldifdt_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zldifdt + (jl - 1), &zldifdt_index, 1);
                    {
                        double __in1 = yrecldp_rcldiff_convi;
                        double __in2 = zldifdt_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        yrecldp_rcldiff_convi_times_zldifdt_slice = __out;
                    }
                    {
                        double __inp = yrecldp_rcldiff_convi_times_zldifdt_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_590_16)
                        __out = __inp;
                        ///////////////////

                        zldifdt[(jl - 1)] = __out;
                    }

                }
            }

        }
        kfdia_plus_1_30 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_30); jl = (jl + 1)) {

            zli_index_1 = zli[((jl + (klon * (jk - 1))) - 1)];
            if ((zli_index_1 > 1e-14)) {
                {
                    double zldifdt_index_0;
                    double zqsmix_index_3;
                    double zqx_index_28;
                    double zqsmix_slice_minus_zqx_slice_0;
                    double max_zqsmix_slice_zqx_slice_0_0;
                    double ze;
                    double za_index_10;
                    double zleros;
                    double zevaplimmix_index_0;
                    double min_zleros_zevaplimmix_slice;
                    double zli_index_2;
                    double min_zleros_zli_slice;
                    double zlicld_index;
                    double zaeros;
                    double zsolac_index_1;
                    double zsolac_slice_minus_zaeros;
                    double zsolqa_index_11;
                    double zliqfrac_index_0;
                    double zliqfrac_slice_times_zleros;
                    double zsolqa_slice_plus_zliqfrac_slice_zleros;
                    double zsolqa_index_12;
                    double zliqfrac_index_1;
                    double zliqfrac_slice_times_zleros_0;
                    double zsolqa_slice_minus_zliqfrac_slice_zleros;
                    double zsolqa_index_13;
                    double zicefrac_index;
                    double zicefrac_slice_times_zleros;
                    double zsolqa_slice_plus_zicefrac_slice_zleros;
                    double zsolqa_index_14;
                    double zicefrac_index_0;
                    double zicefrac_slice_times_zleros_0;
                    double zsolqa_slice_minus_zicefrac_slice_zleros;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zldifdt + (jl - 1), &zldifdt_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_3, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_28, 1);
                    {
                        double __in1 = zqsmix_index_3;
                        double __in2 = zqx_index_28;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqsmix_slice_minus_zqx_slice_0 = __out;
                    }
                    {
                        double __in_a = zqsmix_slice_minus_zqx_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        max_zqsmix_slice_zqx_slice_0_0 = __out;
                    }
                    {
                        double __in1 = zldifdt_index_0;
                        double __in2 = max_zqsmix_slice_zqx_slice_0_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ze = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_10, 1);
                    {
                        double __in1 = za_index_10;
                        double __in2 = ze;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zleros = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zevaplimmix + (jl - 1), &zevaplimmix_index_0, 1);
                    {
                        double __in_a = zleros;
                        double __in_b = zevaplimmix_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        min_zleros_zevaplimmix_slice = __out;
                    }
                    {
                        double __inp = min_zleros_zevaplimmix_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_595_16)
                        __out = __inp;
                        ///////////////////

                        zleros = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zli + ((jl + (klon * (jk - 1))) - 1), &zli_index_2, 1);
                    {
                        double __in_a = zleros;
                        double __in_b = zli_index_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        min_zleros_zli_slice = __out;
                    }
                    {
                        double __inp = min_zleros_zli_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_596_16)
                        __out = __inp;
                        ///////////////////

                        zleros = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zlicld + (jl - 1), &zlicld_index, 1);
                    {
                        double __in1 = zleros;
                        double __in2 = zlicld_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zaeros = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolac + (jl - 1), &zsolac_index_1, 1);
                    {
                        double __in1 = zsolac_index_1;
                        double __in2 = zaeros;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolac_slice_minus_zaeros = __out;
                    }
                    {
                        double __inp = zsolac_slice_minus_zaeros;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_598_16)
                        __out = __inp;
                        ///////////////////

                        zsolac[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (4 * klon)) - 1), &zsolqa_index_11, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zliqfrac + ((jl + (klon * (jk - 1))) - 1), &zliqfrac_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zliqfrac + ((jl + (klon * (jk - 1))) - 1), &zliqfrac_index_1, 1);
                    {
                        double __in1 = zliqfrac_index_0;
                        double __in2 = zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zliqfrac_slice_times_zleros = __out;
                    }
                    {
                        double __in1 = zsolqa_index_11;
                        double __in2 = zliqfrac_slice_times_zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqa_slice_plus_zliqfrac_slice_zleros = __out;
                    }
                    {
                        double __inp = zsolqa_slice_plus_zliqfrac_slice_zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_599_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (4 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (20 * klon)) - 1), &zsolqa_index_12, 1);
                    {
                        double __in1 = zliqfrac_index_1;
                        double __in2 = zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zliqfrac_slice_times_zleros_0 = __out;
                    }
                    {
                        double __in1 = zsolqa_index_12;
                        double __in2 = zliqfrac_slice_times_zleros_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolqa_slice_minus_zliqfrac_slice_zleros = __out;
                    }
                    {
                        double __inp = zsolqa_slice_minus_zliqfrac_slice_zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_600_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (20 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (9 * klon)) - 1), &zsolqa_index_13, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zicefrac + ((jl + (klon * (jk - 1))) - 1), &zicefrac_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zicefrac + ((jl + (klon * (jk - 1))) - 1), &zicefrac_index_0, 1);
                    {
                        double __in1 = zicefrac_index;
                        double __in2 = zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zicefrac_slice_times_zleros = __out;
                    }
                    {
                        double __in1 = zsolqa_index_13;
                        double __in2 = zicefrac_slice_times_zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqa_slice_plus_zicefrac_slice_zleros = __out;
                    }
                    {
                        double __inp = zsolqa_slice_plus_zicefrac_slice_zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_601_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (9 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (21 * klon)) - 1), &zsolqa_index_14, 1);
                    {
                        double __in1 = zicefrac_index_0;
                        double __in2 = zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zicefrac_slice_times_zleros_0 = __out;
                    }
                    {
                        double __in1 = zsolqa_index_14;
                        double __in2 = zicefrac_slice_times_zleros_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolqa_slice_minus_zicefrac_slice_zleros = __out;
                    }
                    {
                        double __inp = zsolqa_slice_minus_zicefrac_slice_zleros;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_602_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (21 * klon)) - 1)] = __out;
                    }

                }
            }

        }
        kfdia_plus_1_31 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_31); jl = (jl + 1)) {
            {
                double ztp1_index_30;
                double zrdcp_times_ztp1_slice;
                double pap_index_8;
                double zdtdp_0;
                double zdp_index_1;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_30, 1);
                {
                    double __in1 = zrdcp;
                    double __in2 = ztp1_index_30;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zrdcp_times_ztp1_slice = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_8, 1);
                {
                    double __in1 = zrdcp_times_ztp1_slice;
                    double __in2 = pap_index_8;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zdtdp_0 = __out;
                }
                {
                    double __inp = zdtdp_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_604_12)
                    __out = __inp;
                    ///////////////////

                    zdtdp = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zdp + (jl - 1), &zdp_index_1, 1);
                {
                    double __in1 = zdp_index_1;
                    double __in2 = zqtmst;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zdpmxdt = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_606_12)
                    __out = 0.0;
                    ///////////////////

                    zmfdn = __out;
                }

            }
            if ((jk < klev)) {
                {
                    double pmfu_index_1;
                    double pmfd_index_1;
                    double zmfdn_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pmfu + (((jk * klon) + jl) - 1), &pmfu_index_1, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pmfd + (((jk * klon) + jl) - 1), &pmfd_index_1, 1);
                    {
                        double __in1 = pmfu_index_1;
                        double __in2 = pmfd_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zmfdn_0 = __out;
                    }
                    {
                        double __inp = zmfdn_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_608_16)
                        __out = __inp;
                        ///////////////////

                        zmfdn = __out;
                    }

                }
            }
            {
                double pvervel_index;
                double __tmp61;
                double pmfu_index_2;
                double pmfd_index_2;
                double pmfu_slice_plus_pmfd_slice_1;
                double pmfu_slice_pmfd_slice_plus_zmfdn;
                double __tmp62;
                double zwtot;
                double neg_zdpmxdt;
                double max_expr_zwtot;
                double min_zdpmxdt_expr;
                double phrsw_index;
                double phrlw_index;
                double zzzdt;
                double zdpmxdt_times_zdtdp;
                double neg_zdpmxdt_0;
                double expr_times_zdtdp;
                double max_expr_zdtdp_zzzdt;
                double min_zdpmxdt_zdtdp_expr;
                double expr_times_ptsphy;
                double zldefr_index;
                double ydthf_ralfdcp_times_zldefr_slice;
                double zdtdiab;
                double zdtdp_times_zwtot;
                double zdtdp_zwtot_times_ptsphy;
                double zdtforc_0;
                double zqold_slice;
                double ztold_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pvervel + ((jl + (klon * (jk - 1))) - 1), &pvervel_index, 1);
                {
                    double __in2 = ydcst_rg;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (0.5 * __in2);
                    ///////////////////

                    __tmp61 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pmfu + ((jl + (klon * (jk - 1))) - 1), &pmfu_index_2, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pmfd + ((jl + (klon * (jk - 1))) - 1), &pmfd_index_2, 1);
                {
                    double __in1 = pmfu_index_2;
                    double __in2 = pmfd_index_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pmfu_slice_plus_pmfd_slice_1 = __out;
                }
                {
                    double __in1 = pmfu_slice_plus_pmfd_slice_1;
                    double __in2 = zmfdn;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pmfu_slice_pmfd_slice_plus_zmfdn = __out;
                }
                {
                    double __in1 = __tmp61;
                    double __in2 = pmfu_slice_pmfd_slice_plus_zmfdn;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp62 = __out;
                }
                {
                    double __in1 = pvervel_index;
                    double __in2 = __tmp62;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zwtot = __out;
                }
                {
                    double __in = zdpmxdt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_USub_)
                    __out = (- __in);
                    ///////////////////

                    neg_zdpmxdt = __out;
                }
                {
                    double __in_a = neg_zdpmxdt;
                    double __in_b = zwtot;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_expr_zwtot = __out;
                }
                {
                    double __in_a = zdpmxdt;
                    double __in_b = max_expr_zwtot;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_zdpmxdt_expr = __out;
                }
                {
                    double __inp = min_zdpmxdt_expr;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_610_12)
                    __out = __inp;
                    ///////////////////

                    zwtot = __out;
                }
                {
                    double __in = zdpmxdt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_USub_)
                    __out = (- __in);
                    ///////////////////

                    neg_zdpmxdt_0 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                phrsw + ((jl + (klon * (jk - 1))) - 1), &phrsw_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                phrlw + ((jl + (klon * (jk - 1))) - 1), &phrlw_index, 1);
                {
                    double __in1 = phrsw_index;
                    double __in2 = phrlw_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zzzdt = __out;
                }
                {
                    double __in1 = zdtdp;
                    double __in2 = zwtot;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zdtdp_times_zwtot = __out;
                }
                {
                    double __in1 = zdpmxdt;
                    double __in2 = zdtdp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zdpmxdt_times_zdtdp = __out;
                }
                {
                    double __in1 = neg_zdpmxdt_0;
                    double __in2 = zdtdp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_zdtdp = __out;
                }
                {
                    double __in_a = expr_times_zdtdp;
                    double __in_b = zzzdt;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_expr_zdtdp_zzzdt = __out;
                }
                {
                    double __in_a = zdpmxdt_times_zdtdp;
                    double __in_b = max_expr_zdtdp_zzzdt;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_zdpmxdt_zdtdp_expr = __out;
                }
                {
                    double __in1 = min_zdpmxdt_zdtdp_expr;
                    double __in2 = ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_ptsphy = __out;
                }
                {
                    double __in1 = zdtdp_times_zwtot;
                    double __in2 = ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zdtdp_zwtot_times_ptsphy = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zldefr + (jl - 1), &zldefr_index, 1);
                {
                    double __in1 = ydthf_ralfdcp;
                    double __in2 = zldefr_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_ralfdcp_times_zldefr_slice = __out;
                }
                {
                    double __in1 = expr_times_ptsphy;
                    double __in2 = ydthf_ralfdcp_times_zldefr_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zdtdiab = __out;
                }
                {
                    double __in1 = zdtdp_zwtot_times_ptsphy;
                    double __in2 = zdtdiab;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zdtforc_0 = __out;
                }
                {
                    double __inp = zdtforc_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_613_12)
                    __out = __inp;
                    ///////////////////

                    zdtforc = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqold_slice, 1);
                {
                    double __inp = zqold_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_614_12)
                    __out = __inp;
                    ///////////////////

                    zqold[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztold_slice, 1);
                {
                    double __inp = ztold_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_615_12)
                    __out = __inp;
                    ///////////////////

                    ztold[(jl - 1)] = __out;
                }

            }
            {
                double ztp1_index_31;
                double ztp1_slice_plus_zdtforc;
                double ztp1_index_32;
                double max_ztp1_slice_160_0;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_31, 1);
                {
                    double __in1 = ztp1_index_31;
                    double __in2 = zdtforc;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    ztp1_slice_plus_zdtforc = __out;
                }
                {
                    double __inp = ztp1_slice_plus_zdtforc;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_616_12)
                    __out = __inp;
                    ///////////////////

                    ztp1[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_32, 1);
                {
                    double __in_a = ztp1_index_32;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, 160);
                    ///////////////////

                    max_ztp1_slice_160_0 = __out;
                }
                {
                    double __inp = max_ztp1_slice_160_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_617_12)
                    __out = __inp;
                    ///////////////////

                    ztp1[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }
        kfdia_plus_1_32 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_32); jl = (jl + 1)) {
            {
                double pap_index_9;
                double ztp1_index_33;
                double min_ydthf_rtwat_ztp1_slice_4;
                double max_ydthf_rtice_expr_4;
                double expr_minus_ydthf_rtice_4;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_4;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_4;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_3;
                double ztp1_index_34;
                double ztp1_slice_minus_ydcst_rtt_7;
                double ydthf_r3les_times_ztp1_slice_ydcst_rtt_3;
                double ztp1_index_35;
                double ztp1_slice_minus_ydthf_r4les_4;
                double ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_3;
                double exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_3;
                double expr_times_expr_0;
                double ztp1_index_36;
                double min_ydthf_rtwat_ztp1_slice_5;
                double max_ydthf_rtice_expr_5;
                double expr_minus_ydthf_rtice_5;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_5;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_5;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_4;
                double __tmp63;
                double ztp1_index_37;
                double ztp1_slice_minus_ydcst_rtt_8;
                double ydthf_r3ies_times_ztp1_slice_ydcst_rtt_2;
                double ztp1_index_38;
                double ztp1_slice_minus_ydthf_r4ies_3;
                double ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_2;
                double exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_2;
                double __tmp64;
                double expr_expr_plus_1_0_expr_expr_0;
                double ydthf_r2es_times_expr_expr_1_0_expr_expr_0;
                double min_0_5_zqsat;
                double ydcst_retv_times_zqsat;
                double __tmp65;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_9, 1);
                {
                    double __in2 = pap_index_9;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    zqp = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_33, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_34, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_35, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_36, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_37, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_38, 1);
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_33;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_4 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_36;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_5 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_4 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_4;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_4 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_5 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_5;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_5 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_4;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_4 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_4 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_3 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_5;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_5 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_5 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_4 = __out;
                }
                {
                    double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp63 = __out;
                }
                {
                    double __in1 = ztp1_index_34;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_7 = __out;
                }
                {
                    double __in1 = ztp1_index_37;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_8 = __out;
                }
                {
                    double __in1 = ydthf_r3les;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_7;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3les_times_ztp1_slice_ydcst_rtt_3 = __out;
                }
                {
                    double __in1 = ztp1_index_35;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les_4 = __out;
                }
                {
                    double __in1 = ydthf_r3les_times_ztp1_slice_ydcst_rtt_3;
                    double __in2 = ztp1_slice_minus_ydthf_r4les_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_3 = __out;
                }
                {
                    double __in1 = ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_3 = __out;
                }
                {
                    double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_3;
                    double __in2 = exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_expr_0 = __out;
                }
                {
                    double __in1 = ydthf_r3ies;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_8;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3ies_times_ztp1_slice_ydcst_rtt_2 = __out;
                }
                {
                    double __in1 = ztp1_index_38;
                    double __in2 = ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4ies_3 = __out;
                }
                {
                    double __in1 = ydthf_r3ies_times_ztp1_slice_ydcst_rtt_2;
                    double __in2 = ztp1_slice_minus_ydthf_r4ies_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_2 = __out;
                }
                {
                    double __in1 = ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_2 = __out;
                }
                {
                    double __in1 = __tmp63;
                    double __in2 = exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp64 = __out;
                }
                {
                    double __in1 = expr_times_expr_0;
                    double __in2 = __tmp64;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    expr_expr_plus_1_0_expr_expr_0 = __out;
                }
                {
                    double __in1 = ydthf_r2es;
                    double __in2 = expr_expr_plus_1_0_expr_expr_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r2es_times_expr_expr_1_0_expr_expr_0 = __out;
                }
                {
                    double __in1 = ydthf_r2es_times_expr_expr_1_0_expr_expr_0;
                    double __in2 = zqp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqsat = __out;
                }
                {
                    double __in_b = zqsat;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(0.5, __in_b);
                    ///////////////////

                    min_0_5_zqsat = __out;
                }
                {
                    double __inp = min_0_5_zqsat;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_622_12)
                    __out = __inp;
                    ///////////////////

                    zqsat = __out;
                }
                {
                    double __in1 = ydcst_retv;
                    double __in2 = zqsat;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_retv_times_zqsat = __out;
                }
                {
                    double __in2 = ydcst_retv_times_zqsat;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp65 = __out;
                }
                {
                    double __in2 = __tmp65;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    zcor_2 = __out;
                }

            }
            {
                double zqsat_times_zcor;
                double zqsmix_index_4;
                double zqsmix_slice_minus_zqsat;
                double zqsat_times_zcor_0;
                double ztp1_index_39;
                double min_ydthf_rtwat_ztp1_slice_6;
                double max_ydthf_rtice_expr_6;
                double expr_minus_ydthf_rtice_6;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_6;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_6;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_5;
                double expr_times_ydthf_r5alvcp;
                double ztp1_index_40;
                double ztp1_slice_minus_ydthf_r4les_5;
                double ztp1_slice_ydthf_r4les_pow_2_0;
                double __tmp66;
                double expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2;
                double ztp1_index_41;
                double min_ydthf_rtwat_ztp1_slice_7;
                double max_ydthf_rtice_expr_7;
                double expr_minus_ydthf_rtice_7;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_7;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_7;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_6;
                double __tmp67;
                double __tmp68;
                double ztp1_index_42;
                double ztp1_slice_minus_ydthf_r4ies_4;
                double ztp1_slice_ydthf_r4ies_pow_2_0;
                double __tmp69;
                double __tmp70;
                double expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2;
                double zqsat_zcor_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2;
                double __tmp71;
                double zcond;
                double ztp1_index_43;
                double ztp1_index_44;
                double min_ydthf_rtwat_ztp1_slice_8;
                double max_ydthf_rtice_expr_8;
                double expr_minus_ydthf_rtice_8;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_8;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_8;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_7;
                double expr_times_ydthf_ralvdcp_0;
                double ztp1_index_45;
                double min_ydthf_rtwat_ztp1_slice_9;
                double max_ydthf_rtice_expr_9;
                double expr_minus_ydthf_rtice_9;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_9;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_9;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_8;
                double __tmp72;
                double __tmp73;
                double expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp_0;
                double expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zcond;
                double ztp1_slice_plus_expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_zcond;
                double zqsmix_index_5;
                double zqsmix_slice_minus_zcond;
                double ztp1_index_46;
                double min_ydthf_rtwat_ztp1_slice_10;
                double max_ydthf_rtice_expr_10;
                double expr_minus_ydthf_rtice_10;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_10;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_10;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_9;
                double ztp1_index_47;
                double ztp1_slice_minus_ydcst_rtt_9;
                double ydthf_r3les_times_ztp1_slice_ydcst_rtt_4;
                double ztp1_index_48;
                double ztp1_slice_minus_ydthf_r4les_6;
                double ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_4;
                double exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_4;
                double expr_times_expr_1;
                double ztp1_index_49;
                double min_ydthf_rtwat_ztp1_slice_11;
                double max_ydthf_rtice_expr_11;
                double expr_minus_ydthf_rtice_11;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_11;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_11;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_10;
                double __tmp74;
                double ztp1_index_50;
                double ztp1_slice_minus_ydcst_rtt_10;
                double ydthf_r3ies_times_ztp1_slice_ydcst_rtt_3;
                double ztp1_index_51;
                double ztp1_slice_minus_ydthf_r4ies_5;
                double ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_3;
                double exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_3;
                double __tmp75;
                double expr_expr_plus_1_0_expr_expr_1;
                double ydthf_r2es_times_expr_expr_1_0_expr_expr_1;
                double zqsat_0;
                double min_0_5_zqsat_0;
                double ydcst_retv_times_zqsat_0;
                double __tmp76;

                {
                    double __inp = zcor_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_623_12)
                    __out = __inp;
                    ///////////////////

                    zcor = __out;
                }
                {
                    double __in1 = zqsat;
                    double __in2 = zcor;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqsat_times_zcor = __out;
                }
                {
                    double __inp = zqsat_times_zcor;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_624_12)
                    __out = __inp;
                    ///////////////////

                    zqsat = __out;
                }
                {
                    double __in1 = zqsat;
                    double __in2 = zcor;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqsat_times_zcor_0 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_4, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_5, 1);
                {
                    double __in1 = zqsmix_index_4;
                    double __in2 = zqsat;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqsmix_slice_minus_zqsat = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_39, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_40, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_41, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_42, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_43, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_44, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_45, 1);
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_39;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_6 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_41;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_7 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_44;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_8 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_45;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_9 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_6 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_6;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_6 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_7;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_7 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_7;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_7 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_8;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_8 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_8;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_8 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_9;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_9 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_9;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_9 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_6;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_6 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_6 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_5 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_7;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_7 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_7;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_7 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_7;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_6 = __out;
                }
                {
                    double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp67 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_8;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_8 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_8;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_8 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_8;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_7 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_9;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_9 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_9;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_9 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_9;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_8 = __out;
                }
                {
                    double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_8;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp72 = __out;
                }
                {
                    double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_5;
                    double __in2 = ydthf_r5alvcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_ydthf_r5alvcp = __out;
                }
                {
                    double __in1 = ztp1_index_40;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les_5 = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_ydthf_r4les_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    ztp1_slice_ydthf_r4les_pow_2_0 = __out;
                }
                {
                    double __in2 = ztp1_slice_ydthf_r4les_pow_2_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    __tmp66 = __out;
                }
                {
                    double __in1 = expr_times_ydthf_r5alvcp;
                    double __in2 = __tmp66;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2 = __out;
                }
                {
                    double __in1 = __tmp67;
                    double __in2 = ydthf_r5alscp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp68 = __out;
                }
                {
                    double __in1 = ztp1_index_42;
                    double __in2 = ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4ies_4 = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_ydthf_r4ies_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    ztp1_slice_ydthf_r4ies_pow_2_0 = __out;
                }
                {
                    double __in2 = ztp1_slice_ydthf_r4ies_pow_2_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    __tmp69 = __out;
                }
                {
                    double __in1 = __tmp68;
                    double __in2 = __tmp69;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp70 = __out;
                }
                {
                    double __in1 = expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2;
                    double __in2 = __tmp70;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2 = __out;
                }
                {
                    double __in1 = zqsat_times_zcor_0;
                    double __in2 = expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqsat_zcor_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2 = __out;
                }
                {
                    double __in2 = zqsat_zcor_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (1.0 + __in2);
                    ///////////////////

                    __tmp71 = __out;
                }
                {
                    double __in1 = zqsmix_slice_minus_zqsat;
                    double __in2 = __tmp71;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zcond = __out;
                }
                {
                    double __in1 = zqsmix_index_5;
                    double __in2 = zcond;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqsmix_slice_minus_zcond = __out;
                }
                {
                    double __inp = zqsmix_slice_minus_zcond;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_627_12)
                    __out = __inp;
                    ///////////////////

                    zqsmix[((jl + (klon * (jk - 1))) - 1)] = __out;
                }
                {
                    double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_7;
                    double __in2 = ydthf_ralvdcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_ydthf_ralvdcp_0 = __out;
                }
                {
                    double __in1 = __tmp72;
                    double __in2 = ydthf_ralsdcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp73 = __out;
                }
                {
                    double __in1 = expr_times_ydthf_ralvdcp_0;
                    double __in2 = __tmp73;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp_0 = __out;
                }
                {
                    double __in1 = expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp_0;
                    double __in2 = zcond;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zcond = __out;
                }
                {
                    double __in1 = ztp1_index_43;
                    double __in2 = expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zcond;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    ztp1_slice_plus_expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_zcond = __out;
                }
                {
                    double __inp = ztp1_slice_plus_expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_zcond;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_626_12)
                    __out = __inp;
                    ///////////////////

                    ztp1[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_46, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_47, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_48, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_49, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_50, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_51, 1);
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_46;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_10 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_10;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_10 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_10;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_10 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_10;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_10 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_10;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_10 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_10;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_9 = __out;
                }
                {
                    double __in1 = ztp1_index_48;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les_6 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_49;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_11 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_11;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_11 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_11;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_11 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_11;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_11 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_11;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_11 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_11;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_10 = __out;
                }
                {
                    double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_10;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp74 = __out;
                }
                {
                    double __in1 = ztp1_index_51;
                    double __in2 = ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4ies_5 = __out;
                }
                {
                    double __in1 = ztp1_index_47;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_9 = __out;
                }
                {
                    double __in1 = ztp1_index_50;
                    double __in2 = ydcst_rtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydcst_rtt_10 = __out;
                }
                {
                    double __in1 = ydthf_r3les;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_9;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3les_times_ztp1_slice_ydcst_rtt_4 = __out;
                }
                {
                    double __in1 = ydthf_r3les_times_ztp1_slice_ydcst_rtt_4;
                    double __in2 = ztp1_slice_minus_ydthf_r4les_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_4 = __out;
                }
                {
                    double __in1 = ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_4 = __out;
                }
                {
                    double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_9;
                    double __in2 = exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_expr_1 = __out;
                }
                {
                    double __in1 = ydthf_r3ies;
                    double __in2 = ztp1_slice_minus_ydcst_rtt_10;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r3ies_times_ztp1_slice_ydcst_rtt_3 = __out;
                }
                {
                    double __in1 = ydthf_r3ies_times_ztp1_slice_ydcst_rtt_3;
                    double __in2 = ztp1_slice_minus_ydthf_r4ies_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_3 = __out;
                }
                {
                    double __in1 = ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_exp_)
                    __out = exp(__in1);
                    ///////////////////

                    exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_3 = __out;
                }
                {
                    double __in1 = __tmp74;
                    double __in2 = exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp75 = __out;
                }
                {
                    double __in1 = expr_times_expr_1;
                    double __in2 = __tmp75;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    expr_expr_plus_1_0_expr_expr_1 = __out;
                }
                {
                    double __in1 = ydthf_r2es;
                    double __in2 = expr_expr_plus_1_0_expr_expr_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydthf_r2es_times_expr_expr_1_0_expr_expr_1 = __out;
                }
                {
                    double __in1 = ydthf_r2es_times_expr_expr_1_0_expr_expr_1;
                    double __in2 = zqp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqsat_0 = __out;
                }
                {
                    double __inp = zqsat_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_628_12)
                    __out = __inp;
                    ///////////////////

                    zqsat = __out;
                }
                {
                    double __in_b = zqsat;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(0.5, __in_b);
                    ///////////////////

                    min_0_5_zqsat_0 = __out;
                }
                {
                    double __inp = min_0_5_zqsat_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_629_12)
                    __out = __inp;
                    ///////////////////

                    zqsat = __out;
                }
                {
                    double __in1 = ydcst_retv;
                    double __in2 = zqsat;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ydcst_retv_times_zqsat_0 = __out;
                }
                {
                    double __in2 = ydcst_retv_times_zqsat_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp76 = __out;
                }
                {
                    double __in2 = __tmp76;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    zcor_3 = __out;
                }

            }
            {
                double zqsat_times_zcor_1;
                double zqsmix_index_6;
                double zqsmix_slice_minus_zqsat_0;
                double zqsat_times_zcor_2;
                double ztp1_index_52;
                double min_ydthf_rtwat_ztp1_slice_12;
                double max_ydthf_rtice_expr_12;
                double expr_minus_ydthf_rtice_12;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_12;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_12;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_11;
                double expr_times_ydthf_r5alvcp_0;
                double ztp1_index_53;
                double ztp1_slice_minus_ydthf_r4les_7;
                double ztp1_slice_ydthf_r4les_pow_2_1;
                double __tmp77;
                double expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2_0;
                double ztp1_index_54;
                double min_ydthf_rtwat_ztp1_slice_13;
                double max_ydthf_rtice_expr_13;
                double expr_minus_ydthf_rtice_13;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_13;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_13;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_12;
                double __tmp78;
                double __tmp79;
                double ztp1_index_55;
                double ztp1_slice_minus_ydthf_r4ies_6;
                double ztp1_slice_ydthf_r4ies_pow_2_1;
                double __tmp80;
                double __tmp81;
                double expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_0;
                double zqsat_zcor_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_0;
                double __tmp82;
                double zcond1;
                double ztp1_index_56;
                double ztp1_index_57;
                double min_ydthf_rtwat_ztp1_slice_14;
                double max_ydthf_rtice_expr_14;
                double expr_minus_ydthf_rtice_14;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_14;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_14;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_13;
                double expr_times_ydthf_ralvdcp_1;
                double ztp1_index_58;
                double min_ydthf_rtwat_ztp1_slice_15;
                double max_ydthf_rtice_expr_15;
                double expr_minus_ydthf_rtice_15;
                double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_15;
                double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_15;
                double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_14;
                double __tmp83;
                double __tmp84;
                double expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp_1;
                double expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zcond1;
                double ztp1_slice_plus_expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_zcond1;
                double zqsmix_index_7;
                double zqsmix_slice_minus_zcond1;

                {
                    double __inp = zcor_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_630_12)
                    __out = __inp;
                    ///////////////////

                    zcor = __out;
                }
                {
                    double __in1 = zqsat;
                    double __in2 = zcor;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqsat_times_zcor_1 = __out;
                }
                {
                    double __inp = zqsat_times_zcor_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_631_12)
                    __out = __inp;
                    ///////////////////

                    zqsat = __out;
                }
                {
                    double __in1 = zqsat;
                    double __in2 = zcor;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqsat_times_zcor_2 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_6, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_7, 1);
                {
                    double __in1 = zqsmix_index_6;
                    double __in2 = zqsat;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqsmix_slice_minus_zqsat_0 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_52, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_53, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_54, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_55, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_56, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_57, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_58, 1);
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_52;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_12 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_54;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_13 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_57;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_14 = __out;
                }
                {
                    double __in_a = ydthf_rtwat;
                    double __in_b = ztp1_index_58;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_ydthf_rtwat_ztp1_slice_15 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_12;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_12 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_12;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_12 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_13;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_13 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_13;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_13 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_14;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_14 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_14;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_14 = __out;
                }
                {
                    double __in_a = ydthf_rtice;
                    double __in_b = min_ydthf_rtwat_ztp1_slice_15;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_ydthf_rtice_expr_15 = __out;
                }
                {
                    double __in1 = max_ydthf_rtice_expr_15;
                    double __in2 = ydthf_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_ydthf_rtice_15 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_12;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_12 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_12;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_12 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_12;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_11 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_13;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_13 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_13;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_13 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_13;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_12 = __out;
                }
                {
                    double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_12;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp78 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_14;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_14 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_14;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_14 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_14;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_13 = __out;
                }
                {
                    double __in1 = expr_minus_ydthf_rtice_15;
                    double __in2 = ydthf_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_15 = __out;
                }
                {
                    double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_15;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_15 = __out;
                }
                {
                    double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_15;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_14 = __out;
                }
                {
                    double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_14;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp83 = __out;
                }
                {
                    double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_11;
                    double __in2 = ydthf_r5alvcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_ydthf_r5alvcp_0 = __out;
                }
                {
                    double __in1 = ztp1_index_53;
                    double __in2 = ydthf_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4les_7 = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_ydthf_r4les_7;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    ztp1_slice_ydthf_r4les_pow_2_1 = __out;
                }
                {
                    double __in2 = ztp1_slice_ydthf_r4les_pow_2_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    __tmp77 = __out;
                }
                {
                    double __in1 = expr_times_ydthf_r5alvcp_0;
                    double __in2 = __tmp77;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2_0 = __out;
                }
                {
                    double __in1 = __tmp78;
                    double __in2 = ydthf_r5alscp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp79 = __out;
                }
                {
                    double __in1 = ztp1_index_55;
                    double __in2 = ydthf_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_ydthf_r4ies_6 = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_ydthf_r4ies_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    ztp1_slice_ydthf_r4ies_pow_2_1 = __out;
                }
                {
                    double __in2 = ztp1_slice_ydthf_r4ies_pow_2_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    __tmp80 = __out;
                }
                {
                    double __in1 = __tmp79;
                    double __in2 = __tmp80;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp81 = __out;
                }
                {
                    double __in1 = expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2_0;
                    double __in2 = __tmp81;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_0 = __out;
                }
                {
                    double __in1 = zqsat_times_zcor_2;
                    double __in2 = expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqsat_zcor_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_0 = __out;
                }
                {
                    double __in2 = zqsat_zcor_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (1.0 + __in2);
                    ///////////////////

                    __tmp82 = __out;
                }
                {
                    double __in1 = zqsmix_slice_minus_zqsat_0;
                    double __in2 = __tmp82;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zcond1 = __out;
                }
                {
                    double __in1 = zqsmix_index_7;
                    double __in2 = zcond1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqsmix_slice_minus_zcond1 = __out;
                }
                {
                    double __inp = zqsmix_slice_minus_zcond1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_634_12)
                    __out = __inp;
                    ///////////////////

                    zqsmix[((jl + (klon * (jk - 1))) - 1)] = __out;
                }
                {
                    double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_13;
                    double __in2 = ydthf_ralvdcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_ydthf_ralvdcp_1 = __out;
                }
                {
                    double __in1 = __tmp83;
                    double __in2 = ydthf_ralsdcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp84 = __out;
                }
                {
                    double __in1 = expr_times_ydthf_ralvdcp_1;
                    double __in2 = __tmp84;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp_1 = __out;
                }
                {
                    double __in1 = expr_ydthf_ralvdcp_plus_1_0_expr_ydthf_ralsdcp_1;
                    double __in2 = zcond1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zcond1 = __out;
                }
                {
                    double __in1 = ztp1_index_56;
                    double __in2 = expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_times_zcond1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    ztp1_slice_plus_expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_zcond1 = __out;
                }
                {
                    double __inp = ztp1_slice_plus_expr_ydthf_ralvdcp_1_0_expr_ydthf_ralsdcp_zcond1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_633_12)
                    __out = __inp;
                    ///////////////////

                    ztp1[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }
        kfdia_plus_1_33 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_33); jl = (jl + 1)) {
            {
                double zqsmix_index_8;
                double zqold_index;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_8, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqold + (jl - 1), &zqold_index, 1);
                {
                    double __in1 = zqsmix_index_8;
                    double __in2 = zqold_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zdqs_slice_0 = __out;
                }

            }
            {
                double zqsmix_slice_0;
                double ztp1_slice_0;

                {
                    double __inp = zdqs_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_636_12)
                    __out = __inp;
                    ///////////////////

                    zdqs[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqold + (jl - 1), &zqsmix_slice_0, 1);
                {
                    double __inp = zqsmix_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_637_12)
                    __out = __inp;
                    ///////////////////

                    zqsmix[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztold + (jl - 1), &ztp1_slice_0, 1);
                {
                    double __inp = ztp1_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_638_12)
                    __out = __inp;
                    ///////////////////

                    ztp1[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }
        kfdia_plus_1_34 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_34); jl = (jl + 1)) {

            zdqs_index_0 = zdqs[(jl - 1)];
            if ((zdqs_index_0 > 0.0)) {
                {
                    double za_index_11;
                    double zdqs_index_1;
                    double zlicld_index_0;
                    double min_zdqs_slice_zlicld_slice;
                    double zlevap;
                    double zevaplimmix_index_1;
                    double min_zlevap_zevaplimmix_slice;
                    double zqsmix_index_9;
                    double zqx_index_29;
                    double zqsmix_slice_minus_zqx_slice_1;
                    double max_zqsmix_slice_zqx_slice_0_0_0;
                    double min_zlevap_expr;
                    double zsolqa_index_15;
                    double zliqfrac_index_3;
                    double zliqfrac_slice_times_zlevap;
                    double zsolqa_slice_plus_zliqfrac_slice_zlevap;
                    double zsolqa_index_16;
                    double zliqfrac_index_4;
                    double zliqfrac_slice_times_zlevap_0;
                    double zsolqa_slice_minus_zliqfrac_slice_zlevap;
                    double zsolqa_index_17;
                    double zicefrac_index_2;
                    double zicefrac_slice_times_zlevap;
                    double zsolqa_slice_plus_zicefrac_slice_zlevap;
                    double zsolqa_index_18;
                    double zicefrac_index_3;
                    double zicefrac_slice_times_zlevap_0;
                    double zsolqa_slice_minus_zicefrac_slice_zlevap;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_11, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdqs + (jl - 1), &zdqs_index_1, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zlicld + (jl - 1), &zlicld_index_0, 1);
                    {
                        double __in_a = zdqs_index_1;
                        double __in_b = zlicld_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        min_zdqs_slice_zlicld_slice = __out;
                    }
                    {
                        double __in1 = za_index_11;
                        double __in2 = min_zdqs_slice_zlicld_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zlevap = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zevaplimmix + (jl - 1), &zevaplimmix_index_1, 1);
                    {
                        double __in_a = zlevap;
                        double __in_b = zevaplimmix_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        min_zlevap_zevaplimmix_slice = __out;
                    }
                    {
                        double __inp = min_zlevap_zevaplimmix_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_642_16)
                        __out = __inp;
                        ///////////////////

                        zlevap = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_9, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_29, 1);
                    {
                        double __in1 = zqsmix_index_9;
                        double __in2 = zqx_index_29;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqsmix_slice_minus_zqx_slice_1 = __out;
                    }
                    {
                        double __in_a = zqsmix_slice_minus_zqx_slice_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        max_zqsmix_slice_zqx_slice_0_0_0 = __out;
                    }
                    {
                        double __in_a = zlevap;
                        double __in_b = max_zqsmix_slice_zqx_slice_0_0_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        min_zlevap_expr = __out;
                    }
                    {
                        double __inp = min_zlevap_expr;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_643_16)
                        __out = __inp;
                        ///////////////////

                        zlevap = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zliqfrac + ((jl + (klon * (jk - 1))) - 1), &zliqfrac_index_3, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zliqfrac + ((jl + (klon * (jk - 1))) - 1), &zliqfrac_index_4, 1);
                    {
                        double __in1 = zliqfrac_index_3;
                        double __in2 = zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zliqfrac_slice_times_zlevap = __out;
                    }
                    {
                        double __in1 = zliqfrac_index_4;
                        double __in2 = zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zliqfrac_slice_times_zlevap_0 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zicefrac + ((jl + (klon * (jk - 1))) - 1), &zicefrac_index_2, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zicefrac + ((jl + (klon * (jk - 1))) - 1), &zicefrac_index_3, 1);
                    {
                        double __in1 = zicefrac_index_2;
                        double __in2 = zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zicefrac_slice_times_zlevap = __out;
                    }
                    {
                        double __in1 = zicefrac_index_3;
                        double __in2 = zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zicefrac_slice_times_zlevap_0 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (4 * klon)) - 1), &zsolqa_index_15, 1);
                    {
                        double __in1 = zsolqa_index_15;
                        double __in2 = zliqfrac_slice_times_zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqa_slice_plus_zliqfrac_slice_zlevap = __out;
                    }
                    {
                        double __inp = zsolqa_slice_plus_zliqfrac_slice_zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_646_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (4 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (20 * klon)) - 1), &zsolqa_index_16, 1);
                    {
                        double __in1 = zsolqa_index_16;
                        double __in2 = zliqfrac_slice_times_zlevap_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolqa_slice_minus_zliqfrac_slice_zlevap = __out;
                    }
                    {
                        double __inp = zsolqa_slice_minus_zliqfrac_slice_zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_647_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (20 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (9 * klon)) - 1), &zsolqa_index_17, 1);
                    {
                        double __in1 = zsolqa_index_17;
                        double __in2 = zicefrac_slice_times_zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqa_slice_plus_zicefrac_slice_zlevap = __out;
                    }
                    {
                        double __inp = zsolqa_slice_plus_zicefrac_slice_zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_648_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (9 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (21 * klon)) - 1), &zsolqa_index_18, 1);
                    {
                        double __in1 = zsolqa_index_18;
                        double __in2 = zicefrac_slice_times_zlevap_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolqa_slice_minus_zicefrac_slice_zlevap = __out;
                    }
                    {
                        double __inp = zsolqa_slice_minus_zicefrac_slice_zlevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_649_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (21 * klon)) - 1)] = __out;
                    }

                }
            }

        }
        kfdia_plus_1_35 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_35); jl = (jl + 1)) {

            za_index_12 = za[((jl + (klon * (jk - 1))) - 1)];
            zdqs_index_2 = zdqs[(jl - 1)];
            neg_yrecldp_rlmin = (- yrecldp_rlmin);
            if (((za_index_12 > 1e-14) && (zdqs_index_2 <= neg_yrecldp_rlmin))) {

                za_index_13 = za[((jl + (klon * (jk - 1))) - 1)];
                {
                    double zdqs_index_3;
                    double neg_zdqs_slice;
                    double zlcond1_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdqs + (jl - 1), &zdqs_index_3, 1);
                    {
                        double __in = zdqs_index_3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_USub_)
                        __out = (- __in);
                        ///////////////////

                        neg_zdqs_slice = __out;
                    }
                    {
                        double __in_a = neg_zdqs_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        zlcond1_slice = __out;
                    }
                    {
                        double __inp = zlcond1_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_652_16)
                        __out = __inp;
                        ///////////////////

                        zlcond1[(jl - 1)] = __out;
                    }

                }
                if ((za_index_13 > 0.99)) {
                    {
                        double zqsmix_index_10;
                        double ydcst_retv_times_zqsmix_slice_0;
                        double __tmp90;
                        double zcor_4;
                        double zqx_index_30;
                        double zqsmix_index_11;
                        double zqx_slice_minus_zqsmix_slice;
                        double zqsmix_index_12;
                        double zcor_times_zqsmix_slice;
                        double ztp1_index_59;
                        double min_ydthf_rtwat_ztp1_slice_16;
                        double max_ydthf_rtice_expr_16;
                        double expr_minus_ydthf_rtice_16;
                        double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_16;
                        double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_16;
                        double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_15;
                        double expr_times_ydthf_r5alvcp_1;
                        double ztp1_index_60;
                        double ztp1_slice_minus_ydthf_r4les_8;
                        double ztp1_slice_ydthf_r4les_pow_2_2;
                        double __tmp91;
                        double expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2_1;
                        double ztp1_index_61;
                        double min_ydthf_rtwat_ztp1_slice_17;
                        double max_ydthf_rtice_expr_17;
                        double expr_minus_ydthf_rtice_17;
                        double expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_17;
                        double expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_17;
                        double min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_16;
                        double __tmp92;
                        double __tmp93;
                        double ztp1_index_62;
                        double ztp1_slice_minus_ydthf_r4ies_7;
                        double ztp1_slice_ydthf_r4ies_pow_2_2;
                        double __tmp94;
                        double __tmp95;
                        double expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_1;
                        double zcor_zqsmix_slice_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2;
                        double __tmp96;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_10, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_11, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_12, 1);
                        {
                            double __in1 = ydcst_retv;
                            double __in2 = zqsmix_index_10;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            ydcst_retv_times_zqsmix_slice_0 = __out;
                        }
                        {
                            double __in2 = ydcst_retv_times_zqsmix_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp90 = __out;
                        }
                        {
                            double __in2 = __tmp90;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (1.0 / __in2);
                            ///////////////////

                            zcor_4 = __out;
                        }
                        {
                            double __inp = zcor_4;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_654_20)
                            __out = __inp;
                            ///////////////////

                            zcor = __out;
                        }
                        {
                            double __in1 = zcor;
                            double __in2 = zqsmix_index_12;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zcor_times_zqsmix_slice = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_30, 1);
                        {
                            double __in1 = zqx_index_30;
                            double __in2 = zqsmix_index_11;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqx_slice_minus_zqsmix_slice = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_59, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_60, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_61, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_62, 1);
                        {
                            double __in_a = ydthf_rtwat;
                            double __in_b = ztp1_index_59;
                            double __out;

                            ///////////////////
                            // Tasklet code (__min2)
                            __out = min(__in_a, __in_b);
                            ///////////////////

                            min_ydthf_rtwat_ztp1_slice_16 = __out;
                        }
                        {
                            double __in_a = ydthf_rtwat;
                            double __in_b = ztp1_index_61;
                            double __out;

                            ///////////////////
                            // Tasklet code (__min2)
                            __out = min(__in_a, __in_b);
                            ///////////////////

                            min_ydthf_rtwat_ztp1_slice_17 = __out;
                        }
                        {
                            double __in_a = ydthf_rtice;
                            double __in_b = min_ydthf_rtwat_ztp1_slice_16;
                            double __out;

                            ///////////////////
                            // Tasklet code (__max2)
                            __out = max(__in_a, __in_b);
                            ///////////////////

                            max_ydthf_rtice_expr_16 = __out;
                        }
                        {
                            double __in1 = max_ydthf_rtice_expr_16;
                            double __in2 = ydthf_rtice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            expr_minus_ydthf_rtice_16 = __out;
                        }
                        {
                            double __in_a = ydthf_rtice;
                            double __in_b = min_ydthf_rtwat_ztp1_slice_17;
                            double __out;

                            ///////////////////
                            // Tasklet code (__max2)
                            __out = max(__in_a, __in_b);
                            ///////////////////

                            max_ydthf_rtice_expr_17 = __out;
                        }
                        {
                            double __in1 = max_ydthf_rtice_expr_17;
                            double __in2 = ydthf_rtice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            expr_minus_ydthf_rtice_17 = __out;
                        }
                        {
                            double __in1 = expr_minus_ydthf_rtice_16;
                            double __in2 = ydthf_rtwat_rtice_r;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_16 = __out;
                        }
                        {
                            double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_16;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = (dace::math::ipow(__in1, 2));
                            ///////////////////

                            expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_16 = __out;
                        }
                        {
                            double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_16;
                            double __out;

                            ///////////////////
                            // Tasklet code (__min2)
                            __out = min(1, __in_b);
                            ///////////////////

                            min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_15 = __out;
                        }
                        {
                            double __in1 = expr_minus_ydthf_rtice_17;
                            double __in2 = ydthf_rtwat_rtice_r;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_17 = __out;
                        }
                        {
                            double __in1 = expr_ydthf_rtice_times_ydthf_rtwat_rtice_r_17;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = (dace::math::ipow(__in1, 2));
                            ///////////////////

                            expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_17 = __out;
                        }
                        {
                            double __in_b = expr_ydthf_rtice_ydthf_rtwat_rtice_r_pow_2_17;
                            double __out;

                            ///////////////////
                            // Tasklet code (__min2)
                            __out = min(1, __in_b);
                            ///////////////////

                            min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_16 = __out;
                        }
                        {
                            double __in2 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_16;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp92 = __out;
                        }
                        {
                            double __in1 = min_1_0_expr_ydthf_rtice_ydthf_rtwat_rtice_r_2_15;
                            double __in2 = ydthf_r5alvcp;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            expr_times_ydthf_r5alvcp_1 = __out;
                        }
                        {
                            double __in1 = ztp1_index_60;
                            double __in2 = ydthf_r4les;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            ztp1_slice_minus_ydthf_r4les_8 = __out;
                        }
                        {
                            double __in1 = ztp1_slice_minus_ydthf_r4les_8;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = (dace::math::ipow(__in1, 2));
                            ///////////////////

                            ztp1_slice_ydthf_r4les_pow_2_2 = __out;
                        }
                        {
                            double __in2 = ztp1_slice_ydthf_r4les_pow_2_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (1.0 / __in2);
                            ///////////////////

                            __tmp91 = __out;
                        }
                        {
                            double __in1 = expr_times_ydthf_r5alvcp_1;
                            double __in2 = __tmp91;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2_1 = __out;
                        }
                        {
                            double __in1 = __tmp92;
                            double __in2 = ydthf_r5alscp;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp93 = __out;
                        }
                        {
                            double __in1 = ztp1_index_62;
                            double __in2 = ydthf_r4ies;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            ztp1_slice_minus_ydthf_r4ies_7 = __out;
                        }
                        {
                            double __in1 = ztp1_slice_minus_ydthf_r4ies_7;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = (dace::math::ipow(__in1, 2));
                            ///////////////////

                            ztp1_slice_ydthf_r4ies_pow_2_2 = __out;
                        }
                        {
                            double __in2 = ztp1_slice_ydthf_r4ies_pow_2_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (1.0 / __in2);
                            ///////////////////

                            __tmp94 = __out;
                        }
                        {
                            double __in1 = __tmp93;
                            double __in2 = __tmp94;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp95 = __out;
                        }
                        {
                            double __in1 = expr_ydthf_r5alvcp_times_1_0_ztp1_slice_ydthf_r4les_2_1;
                            double __in2 = __tmp95;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_1 = __out;
                        }
                        {
                            double __in1 = zcor_times_zqsmix_slice;
                            double __in2 = expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_plus_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zcor_zqsmix_slice_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2 = __out;
                        }
                        {
                            double __in2 = zcor_zqsmix_slice_times_expr_ydthf_r5alvcp_1_0_ztp1_slice_ydthf_r4les_2_1_0_expr_ydthf_r5alscp_1_0_ztp1_slice_ydthf_r4ies_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (1.0 + __in2);
                            ///////////////////

                            __tmp96 = __out;
                        }
                        {
                            double __in1 = zqx_slice_minus_zqsmix_slice;
                            double __in2 = __tmp96;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / __in2);
                            ///////////////////

                            zcdmax = __out;
                        }

                    }
                } else {
                    {
                        double zqx_index_31;
                        double za_index_14;
                        double zqsmix_index_13;
                        double za_slice_times_zqsmix_slice;
                        double zqx_slice_minus_za_slice_zqsmix_slice;
                        double za_index_15;
                        double zcdmax_0;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_31, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_14, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_15, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqsmix + ((jl + (klon * (jk - 1))) - 1), &zqsmix_index_13, 1);
                        {
                            double __in1 = za_index_14;
                            double __in2 = zqsmix_index_13;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            za_slice_times_zqsmix_slice = __out;
                        }
                        {
                            double __in1 = zqx_index_31;
                            double __in2 = za_slice_times_zqsmix_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqx_slice_minus_za_slice_zqsmix_slice = __out;
                        }
                        {
                            double __in1 = zqx_slice_minus_za_slice_zqsmix_slice;
                            double __in2 = za_index_15;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / __in2);
                            ///////////////////

                            zcdmax_0 = __out;
                        }
                        {
                            double __inp = zcdmax_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_657_20)
                            __out = __inp;
                            ///////////////////

                            zcdmax = __out;
                        }

                    }
                }
                {
                    double zlcond1_index;
                    double min_zlcond1_slice_zcdmax;
                    double max_expr_0_0;
                    double za_index_16;
                    double zlcond1_index_0;
                    double za_slice_times_zlcond1_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zlcond1 + (jl - 1), &zlcond1_index, 1);
                    {
                        double __in_a = zlcond1_index;
                        double __in_b = zcdmax;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        min_zlcond1_slice_zcdmax = __out;
                    }
                    {
                        double __in_a = min_zlcond1_slice_zcdmax;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        max_expr_0_0 = __out;
                    }
                    {
                        double __inp = max_expr_0_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_658_16)
                        __out = __inp;
                        ///////////////////

                        zlcond1[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zlcond1 + (jl - 1), &zlcond1_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_16, 1);
                    {
                        double __in1 = za_index_16;
                        double __in2 = zlcond1_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        za_slice_times_zlcond1_slice = __out;
                    }
                    {
                        double __inp = za_slice_times_zlcond1_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_659_16)
                        __out = __inp;
                        ///////////////////

                        zlcond1[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zlcond1 + (jl - 1), &zlcond1_index_1, 1);

                }
                __tmp97 = (zlcond1_index_1 < yrecldp_rlmin);
                if (__tmp97) {
                    {

                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_661_20)
                            __out = 0.0;
                            ///////////////////

                            zlcond1[(jl - 1)] = __out;
                        }

                    }
                }
                ztp1_index_63 = ztp1[((jl + (klon * (jk - 1))) - 1)];

                __tmp98 = (ztp1_index_63 > yrecldp_rthomo);
                if (__tmp98) {
                    {
                        double zsolqa_index_19;
                        double zlcond1_index_2;
                        double zsolqa_slice_plus_zlcond1_slice;
                        double zsolqa_index_20;
                        double zlcond1_index_3;
                        double zsolqa_slice_minus_zlcond1_slice;
                        double zqxfg_index_3;
                        double zlcond1_index_4;
                        double zqxfg_slice_plus_zlcond1_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (20 * klon)) - 1), &zsolqa_index_19, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlcond1 + (jl - 1), &zlcond1_index_2, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlcond1 + (jl - 1), &zlcond1_index_3, 1);
                        {
                            double __in1 = zsolqa_index_19;
                            double __in2 = zlcond1_index_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zlcond1_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zlcond1_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_663_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (20 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (4 * klon)) - 1), &zsolqa_index_20, 1);
                        {
                            double __in1 = zsolqa_index_20;
                            double __in2 = zlcond1_index_3;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsolqa_slice_minus_zlcond1_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_zlcond1_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_664_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (4 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlcond1 + (jl - 1), &zlcond1_index_4, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxfg + (jl - 1), &zqxfg_index_3, 1);
                        {
                            double __in1 = zqxfg_index_3;
                            double __in2 = zlcond1_index_4;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqxfg_slice_plus_zlcond1_slice = __out;
                        }
                        {
                            double __inp = zqxfg_slice_plus_zlcond1_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_665_20)
                            __out = __inp;
                            ///////////////////

                            zqxfg[(jl - 1)] = __out;
                        }

                    }
                } else {
                    {
                        double zsolqa_index_21;
                        double zlcond1_index_5;
                        double zsolqa_slice_plus_zlcond1_slice_0;
                        double zsolqa_index_22;
                        double zlcond1_index_6;
                        double zsolqa_slice_minus_zlcond1_slice_0;
                        double zqxfg_index_4;
                        double zlcond1_index_7;
                        double zqxfg_slice_plus_zlcond1_slice_0;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (21 * klon)) - 1), &zsolqa_index_21, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlcond1 + (jl - 1), &zlcond1_index_5, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlcond1 + (jl - 1), &zlcond1_index_6, 1);
                        {
                            double __in1 = zsolqa_index_21;
                            double __in2 = zlcond1_index_5;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zlcond1_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zlcond1_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_667_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (21 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (9 * klon)) - 1), &zsolqa_index_22, 1);
                        {
                            double __in1 = zsolqa_index_22;
                            double __in2 = zlcond1_index_6;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsolqa_slice_minus_zlcond1_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_zlcond1_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_668_20)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (9 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlcond1 + (jl - 1), &zlcond1_index_7, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxfg + ((jl + klon) - 1), &zqxfg_index_4, 1);
                        {
                            double __in1 = zqxfg_index_4;
                            double __in2 = zlcond1_index_7;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqxfg_slice_plus_zlcond1_slice_0 = __out;
                        }
                        {
                            double __inp = zqxfg_slice_plus_zlcond1_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_669_20)
                            __out = __inp;
                            ///////////////////

                            zqxfg[((jl + klon) - 1)] = __out;
                        }

                    }
                }
            }

        }
        kfdia_plus_1_36 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_36); jl = (jl + 1)) {

            zdqs_index_4 = zdqs[(jl - 1)];
            neg_yrecldp_rlmin_0 = (- yrecldp_rlmin);
            za_index_17 = za[((jl + (klon * (jk - 1))) - 1)];
            if (((zdqs_index_4 <= neg_yrecldp_rlmin_0) && (za_index_17 < (1.0 - 1e-14)))) {

                pap_index_10 = pap[((jl + (klon * (jk - 1))) - 1)];
                paph_index_3 = paph[((jl + (klev * klon)) - 1)];
                if (((pap_index_10 / paph_index_3) > 0.8)) {
                    {
                        double __tmp103;
                        double zsigk_minus_0_8;
                        double zsigk_0_8_div_0_2;
                        double zsigk_0_8_0_2_pow_2;
                        double __tmp104;

                        {
                            double __in2 = yrecldp_ramid;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp103 = __out;
                        }
                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = ((pap_index_10 / paph_index_3) - 0.8);
                            ///////////////////

                            zsigk_minus_0_8 = __out;
                        }
                        {
                            double __in1 = zsigk_minus_0_8;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / 0.2);
                            ///////////////////

                            zsigk_0_8_div_0_2 = __out;
                        }
                        {
                            double __in1 = zsigk_0_8_div_0_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = (dace::math::ipow(__in1, 2));
                            ///////////////////

                            zsigk_0_8_0_2_pow_2 = __out;
                        }
                        {
                            double __in1 = __tmp103;
                            double __in2 = zsigk_0_8_0_2_pow_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp104 = __out;
                        }
                        {
                            double __in1 = yrecldp_ramid;
                            double __in2 = __tmp104;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zrhc = __out;
                        }

                    }
                } else {
                    {

                        {
                            double __inp = yrecldp_ramid;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_676_20)
                            __out = __inp;
                            ///////////////////

                            zrhc = __out;
                        }

                    }
                }
                if ((yrecldp_nssopt == 0)) {
                    {
                        double zqx_index_32;
                        double za_index_18;
                        double zqsice_index_4;
                        double za_slice_times_zqsice_slice_0;
                        double zqx_slice_minus_za_slice_zqsice_slice_0;
                        double za_index_19;
                        double __tmp105;
                        double max_zepsec_1_0_za_slice;
                        double max_0_0_zqe;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_32, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_18, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_19, 1);
                        {
                            double __in2 = za_index_19;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp105 = __out;
                        }
                        {
                            double __in_b = __tmp105;
                            double __out;

                            ///////////////////
                            // Tasklet code (__max2)
                            __out = max(1e-14, __in_b);
                            ///////////////////

                            max_zepsec_1_0_za_slice = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_4, 1);
                        {
                            double __in1 = za_index_18;
                            double __in2 = zqsice_index_4;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            za_slice_times_zqsice_slice_0 = __out;
                        }
                        {
                            double __in1 = zqx_index_32;
                            double __in2 = za_slice_times_zqsice_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqx_slice_minus_za_slice_zqsice_slice_0 = __out;
                        }
                        {
                            double __in1 = zqx_slice_minus_za_slice_zqsice_slice_0;
                            double __in2 = max_zepsec_1_0_za_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / __in2);
                            ///////////////////

                            zqe = __out;
                        }
                        {
                            double __in_b = zqe;
                            double __out;

                            ///////////////////
                            // Tasklet code (__max2)
                            __out = max(0.0, __in_b);
                            ///////////////////

                            max_0_0_zqe = __out;
                        }
                        {
                            double __inp = max_0_0_zqe;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_679_20)
                            __out = __inp;
                            ///////////////////

                            zqe = __out;
                        }

                    }
                } else {
                    if ((yrecldp_nssopt == 1)) {
                        {
                            double zqx_index_33;
                            double za_index_20;
                            double zqsice_index_5;
                            double za_slice_times_zqsice_slice_1;
                            double zqx_slice_minus_za_slice_zqsice_slice_1;
                            double za_index_21;
                            double __tmp106;
                            double max_zepsec_1_0_za_slice_0;
                            double zqe_0;
                            double max_0_0_zqe_0;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_33, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            za + ((jl + (klon * (jk - 1))) - 1), &za_index_20, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            za + ((jl + (klon * (jk - 1))) - 1), &za_index_21, 1);
                            {
                                double __in2 = za_index_21;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (1.0 - __in2);
                                ///////////////////

                                __tmp106 = __out;
                            }
                            {
                                double __in_b = __tmp106;
                                double __out;

                                ///////////////////
                                // Tasklet code (__max2)
                                __out = max(1e-14, __in_b);
                                ///////////////////

                                max_zepsec_1_0_za_slice_0 = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_5, 1);
                            {
                                double __in1 = za_index_20;
                                double __in2 = zqsice_index_5;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                za_slice_times_zqsice_slice_1 = __out;
                            }
                            {
                                double __in1 = zqx_index_33;
                                double __in2 = za_slice_times_zqsice_slice_1;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zqx_slice_minus_za_slice_zqsice_slice_1 = __out;
                            }
                            {
                                double __in1 = zqx_slice_minus_za_slice_zqsice_slice_1;
                                double __in2 = max_zepsec_1_0_za_slice_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Div_)
                                __out = (__in1 / __in2);
                                ///////////////////

                                zqe_0 = __out;
                            }
                            {
                                double __inp = zqe_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_681_20)
                                __out = __inp;
                                ///////////////////

                                zqe = __out;
                            }
                            {
                                double __in_b = zqe;
                                double __out;

                                ///////////////////
                                // Tasklet code (__max2)
                                __out = max(0.0, __in_b);
                                ///////////////////

                                max_0_0_zqe_0 = __out;
                            }
                            {
                                double __inp = max_0_0_zqe_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_682_20)
                                __out = __inp;
                                ///////////////////

                                zqe = __out;
                            }

                        }
                    } else {
                        if ((yrecldp_nssopt == 2)) {
                            {
                                double zqe_1;


                                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                                zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqe_1, 1);
                                {
                                    double __inp = zqe_1;
                                    double __out;

                                    ///////////////////
                                    // Tasklet code (assign_684_20)
                                    __out = __inp;
                                    ///////////////////

                                    zqe = __out;
                                }

                            }
                        } else {
                            if ((yrecldp_nssopt == 3)) {
                                {
                                    double zqx_index_34;
                                    double zli_index_3;
                                    double zqe_2;


                                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                                    zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_34, 1);

                                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                                    zli + ((jl + (klon * (jk - 1))) - 1), &zli_index_3, 1);
                                    {
                                        double __in1 = zqx_index_34;
                                        double __in2 = zli_index_3;
                                        double __out;

                                        ///////////////////
                                        // Tasklet code (_Add_)
                                        __out = (__in1 + __in2);
                                        ///////////////////

                                        zqe_2 = __out;
                                    }
                                    {
                                        double __inp = zqe_2;
                                        double __out;

                                        ///////////////////
                                        // Tasklet code (assign_686_20)
                                        __out = __inp;
                                        ///////////////////

                                        zqe = __out;
                                    }

                                }
                            }
                        }
                    }
                }
                ztp1_index_64 = ztp1[((jl + (klon * (jk - 1))) - 1)];

                __tmp107 = (ztp1_index_64 >= ydcst_rtt);
                __tmp108 = (yrecldp_nssopt == 0);
                if ((__tmp107 || __tmp108)) {
                    {

                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_688_20)
                            __out = 1.0;
                            ///////////////////

                            zfac = __out;
                        }

                    }
                } else {
                    {
                        double zfac_1;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zfokoop + (jl - 1), &zfac_1, 1);
                        {
                            double __inp = zfac_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_690_20)
                            __out = __inp;
                            ///////////////////

                            zfac = __out;
                        }

                    }
                }
                {
                    double zqsice_index_6;
                    double zrhc_times_zqsice_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_6, 1);
                    {
                        double __in1 = zrhc;
                        double __in2 = zqsice_index_6;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zrhc_times_zqsice_slice = __out;
                    }
                    {
                        double __in1 = zrhc_times_zqsice_slice;
                        double __in2 = zfac;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zrhc_zqsice_slice_times_zfac = __out;
                    }

                }
                __tmp110 = (zqe >= zrhc_zqsice_slice_times_zfac);
                {
                    double zqsice_index_7;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_7, 1);
                    {
                        double __in1 = zqsice_index_7;
                        double __in2 = zfac;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zqsice_slice_times_zfac = __out;
                    }

                }
                __tmp111 = (zqe < zqsice_slice_times_zfac);
                if ((__tmp110 && __tmp111)) {
                    {
                        double za_index_22;
                        double __tmp113;
                        double neg_1_0_za_slice;
                        double expr_times_zfac;
                        double zdqs_index_5;
                        double expr_zfac_times_zdqs_slice;
                        double zqsice_index_8;
                        double zfac_times_zqsice_slice_1;
                        double zfac_zqsice_slice_minus_zqe;
                        double __tmp114;
                        double max_2_0_zfac_zqsice_slice_zqe_zepsec;
                        double za_index_23;
                        double __tmp115;
                        double min_zacond_1_0_za_slice;
                        double neg_zfac;
                        double zdqs_index_6;
                        double expr_times_zdqs_slice;
                        double expr_zdqs_slice_times_0_5;
                        double zlcond2_slice;
                        double zqsice_index_9;
                        double zfac_times_zqsice_slice_2;
                        double zfac_zqsice_slice_minus_zqe_0;
                        double __tmp116;
                        double za_index_24;
                        double __tmp117;
                        double max_zepsec_1_0_za_slice_1;
                        double zzdl;
                        double zdqs_index_7;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zdqs + (jl - 1), &zdqs_index_7, 1);
                        {
                            double __in1 = zfac;
                            double __in2 = zdqs_index_7;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zfac_times_zdqs_slice = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_22, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_23, 1);
                        {
                            double __in2 = za_index_22;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp113 = __out;
                        }
                        {
                            double __in = __tmp113;
                            double __out;

                            ///////////////////
                            // Tasklet code (_USub_)
                            __out = (- __in);
                            ///////////////////

                            neg_1_0_za_slice = __out;
                        }
                        {
                            double __in2 = za_index_23;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp115 = __out;
                        }
                        {
                            double __in = zfac;
                            double __out;

                            ///////////////////
                            // Tasklet code (_USub_)
                            __out = (- __in);
                            ///////////////////

                            neg_zfac = __out;
                        }
                        {
                            double __in1 = neg_1_0_za_slice;
                            double __in2 = zfac;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            expr_times_zfac = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zdqs + (jl - 1), &zdqs_index_5, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zdqs + (jl - 1), &zdqs_index_6, 1);
                        {
                            double __in1 = expr_times_zfac;
                            double __in2 = zdqs_index_5;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            expr_zfac_times_zdqs_slice = __out;
                        }
                        {
                            double __in1 = neg_zfac;
                            double __in2 = zdqs_index_6;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            expr_times_zdqs_slice = __out;
                        }
                        {
                            double __in1 = expr_times_zdqs_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * 0.5);
                            ///////////////////

                            expr_zdqs_slice_times_0_5 = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_8, 1);
                        {
                            double __in1 = zfac;
                            double __in2 = zqsice_index_8;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zfac_times_zqsice_slice_1 = __out;
                        }
                        {
                            double __in1 = zfac_times_zqsice_slice_1;
                            double __in2 = zqe;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zfac_zqsice_slice_minus_zqe = __out;
                        }
                        {
                            double __in2 = zfac_zqsice_slice_minus_zqe;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (2.0 * __in2);
                            ///////////////////

                            __tmp114 = __out;
                        }
                        {
                            double __in_a = __tmp114;
                            double __out;

                            ///////////////////
                            // Tasklet code (__max2)
                            __out = max(__in_a, 1e-14);
                            ///////////////////

                            max_2_0_zfac_zqsice_slice_zqe_zepsec = __out;
                        }
                        {
                            double __in1 = expr_zfac_times_zdqs_slice;
                            double __in2 = max_2_0_zfac_zqsice_slice_zqe_zepsec;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / __in2);
                            ///////////////////

                            zacond = __out;
                        }
                        {
                            double __in_a = zacond;
                            double __in_b = __tmp115;
                            double __out;

                            ///////////////////
                            // Tasklet code (__min2)
                            __out = min(__in_a, __in_b);
                            ///////////////////

                            min_zacond_1_0_za_slice = __out;
                        }
                        {
                            double __inp = min_zacond_1_0_za_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_693_20)
                            __out = __inp;
                            ///////////////////

                            zacond = __out;
                        }
                        {
                            double __in1 = expr_zdqs_slice_times_0_5;
                            double __in2 = zacond;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zlcond2_slice = __out;
                        }
                        {
                            double __inp = zlcond2_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_694_20)
                            __out = __inp;
                            ///////////////////

                            zlcond2[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_24, 1);
                        {
                            double __in2 = za_index_24;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp117 = __out;
                        }
                        {
                            double __in_b = __tmp117;
                            double __out;

                            ///////////////////
                            // Tasklet code (__max2)
                            __out = max(1e-14, __in_b);
                            ///////////////////

                            max_zepsec_1_0_za_slice_1 = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_9, 1);
                        {
                            double __in1 = zfac;
                            double __in2 = zqsice_index_9;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zfac_times_zqsice_slice_2 = __out;
                        }
                        {
                            double __in1 = zfac_times_zqsice_slice_2;
                            double __in2 = zqe;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zfac_zqsice_slice_minus_zqe_0 = __out;
                        }
                        {
                            double __in2 = zfac_zqsice_slice_minus_zqe_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (2.0 * __in2);
                            ///////////////////

                            __tmp116 = __out;
                        }
                        {
                            double __in1 = __tmp116;
                            double __in2 = max_zepsec_1_0_za_slice_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / __in2);
                            ///////////////////

                            zzdl = __out;
                        }
                        {
                            double __in = zzdl;
                            double __out;

                            ///////////////////
                            // Tasklet code (_USub_)
                            __out = (- __in);
                            ///////////////////

                            neg_zzdl = __out;
                        }

                    }
                    __tmp118 = (zfac_times_zdqs_slice < neg_zzdl);
                    if (__tmp118) {
                        {
                            double za_index_25;
                            double za_slice_minus_1_0;
                            double za_slice_1_0_times_zfac;
                            double zdqs_index_8;
                            double za_slice_1_0_zfac_times_zdqs_slice;
                            double zqsice_index_10;
                            double zfac_times_zqsice_slice_3;
                            double za_slice_1_0_zfac_zdqs_slice_minus_zfac_zqsice_slice;
                            double zqx_index_35;
                            double zlcondlim;
                            double zlcond2_index;
                            double min_zlcond2_slice_zlcondlim;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            za + ((jl + (klon * (jk - 1))) - 1), &za_index_25, 1);
                            {
                                double __in1 = za_index_25;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - 1.0);
                                ///////////////////

                                za_slice_minus_1_0 = __out;
                            }
                            {
                                double __in1 = za_slice_minus_1_0;
                                double __in2 = zfac;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                za_slice_1_0_times_zfac = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zdqs + (jl - 1), &zdqs_index_8, 1);
                            {
                                double __in1 = za_slice_1_0_times_zfac;
                                double __in2 = zdqs_index_8;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                za_slice_1_0_zfac_times_zdqs_slice = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_10, 1);
                            {
                                double __in1 = zfac;
                                double __in2 = zqsice_index_10;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zfac_times_zqsice_slice_3 = __out;
                            }
                            {
                                double __in1 = za_slice_1_0_zfac_times_zdqs_slice;
                                double __in2 = zfac_times_zqsice_slice_3;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                za_slice_1_0_zfac_zdqs_slice_minus_zfac_zqsice_slice = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_35, 1);
                            {
                                double __in1 = za_slice_1_0_zfac_zdqs_slice_minus_zfac_zqsice_slice;
                                double __in2 = zqx_index_35;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zlcondlim = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcond2 + (jl - 1), &zlcond2_index, 1);
                            {
                                double __in_a = zlcond2_index;
                                double __in_b = zlcondlim;
                                double __out;

                                ///////////////////
                                // Tasklet code (__min2)
                                __out = min(__in_a, __in_b);
                                ///////////////////

                                min_zlcond2_slice_zlcondlim = __out;
                            }
                            {
                                double __inp = min_zlcond2_slice_zlcondlim;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_698_24)
                                __out = __inp;
                                ///////////////////

                                zlcond2[(jl - 1)] = __out;
                            }

                        }
                    }
                    {
                        double zlcond2_index_0;
                        double max_zlcond2_slice_0_0;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlcond2 + (jl - 1), &zlcond2_index_0, 1);
                        {
                            double __in_a = zlcond2_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (__max2)
                            __out = max(__in_a, 0.0);
                            ///////////////////

                            max_zlcond2_slice_0_0 = __out;
                        }
                        {
                            double __inp = max_zlcond2_slice_0_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_699_20)
                            __out = __inp;
                            ///////////////////

                            zlcond2[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zlcond2 + (jl - 1), &zlcond2_index_1, 1);

                    }
                    __tmp119 = (zlcond2_index_1 < yrecldp_rlmin);
                    {
                        double za_index_26;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_26, 1);
                        {
                            double __in2 = za_index_26;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp120 = __out;
                        }

                    }
                    __tmp121 = (__tmp120 < 1e-14);
                    if ((__tmp119 || __tmp121)) {
                        {

                            {
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_701_24)
                                __out = 0.0;
                                ///////////////////

                                zlcond2[(jl - 1)] = __out;
                            }
                            {
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_702_24)
                                __out = 0.0;
                                ///////////////////

                                zacond = __out;
                            }

                        }
                    }
                    zlcond2_index_2 = zlcond2[(jl - 1)];
                    if ((zlcond2_index_2 == 0.0)) {
                        {

                            {
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_704_24)
                                __out = 0.0;
                                ///////////////////

                                zacond = __out;
                            }

                        }
                    }
                    ztp1_index_65 = ztp1[((jl + (klon * (jk - 1))) - 1)];
                    {
                        double zsolac_index_2;
                        double zsolac_slice_plus_zacond;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolac + (jl - 1), &zsolac_index_2, 1);
                        {
                            double __in1 = zsolac_index_2;
                            double __in2 = zacond;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolac_slice_plus_zacond = __out;
                        }
                        {
                            double __inp = zsolac_slice_plus_zacond;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_705_20)
                            __out = __inp;
                            ///////////////////

                            zsolac[(jl - 1)] = __out;
                        }

                    }
                    __tmp124 = (ztp1_index_65 > yrecldp_rthomo);
                    if (__tmp124) {
                        {
                            double zsolqa_index_23;
                            double zlcond2_index_3;
                            double zsolqa_slice_plus_zlcond2_slice;
                            double zsolqa_index_24;
                            double zlcond2_index_4;
                            double zsolqa_slice_minus_zlcond2_slice;
                            double zqxfg_index_5;
                            double zlcond2_index_5;
                            double zqxfg_slice_plus_zlcond2_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + ((jl + (20 * klon)) - 1), &zsolqa_index_23, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcond2 + (jl - 1), &zlcond2_index_3, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcond2 + (jl - 1), &zlcond2_index_4, 1);
                            {
                                double __in1 = zsolqa_index_23;
                                double __in2 = zlcond2_index_3;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zsolqa_slice_plus_zlcond2_slice = __out;
                            }
                            {
                                double __inp = zsolqa_slice_plus_zlcond2_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_707_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[((jl + (20 * klon)) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + ((jl + (4 * klon)) - 1), &zsolqa_index_24, 1);
                            {
                                double __in1 = zsolqa_index_24;
                                double __in2 = zlcond2_index_4;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zsolqa_slice_minus_zlcond2_slice = __out;
                            }
                            {
                                double __inp = zsolqa_slice_minus_zlcond2_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_708_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[((jl + (4 * klon)) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcond2 + (jl - 1), &zlcond2_index_5, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + (jl - 1), &zqxfg_index_5, 1);
                            {
                                double __in1 = zqxfg_index_5;
                                double __in2 = zlcond2_index_5;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zqxfg_slice_plus_zlcond2_slice = __out;
                            }
                            {
                                double __inp = zqxfg_slice_plus_zlcond2_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_709_24)
                                __out = __inp;
                                ///////////////////

                                zqxfg[(jl - 1)] = __out;
                            }

                        }
                    } else {
                        {
                            double zsolqa_index_25;
                            double zlcond2_index_6;
                            double zsolqa_slice_plus_zlcond2_slice_0;
                            double zsolqa_index_26;
                            double zlcond2_index_7;
                            double zsolqa_slice_minus_zlcond2_slice_0;
                            double zqxfg_index_6;
                            double zlcond2_index_8;
                            double zqxfg_slice_plus_zlcond2_slice_0;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + ((jl + (21 * klon)) - 1), &zsolqa_index_25, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcond2 + (jl - 1), &zlcond2_index_6, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcond2 + (jl - 1), &zlcond2_index_7, 1);
                            {
                                double __in1 = zsolqa_index_25;
                                double __in2 = zlcond2_index_6;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zsolqa_slice_plus_zlcond2_slice_0 = __out;
                            }
                            {
                                double __inp = zsolqa_slice_plus_zlcond2_slice_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_711_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[((jl + (21 * klon)) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + ((jl + (9 * klon)) - 1), &zsolqa_index_26, 1);
                            {
                                double __in1 = zsolqa_index_26;
                                double __in2 = zlcond2_index_7;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zsolqa_slice_minus_zlcond2_slice_0 = __out;
                            }
                            {
                                double __inp = zsolqa_slice_minus_zlcond2_slice_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_712_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[((jl + (9 * klon)) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zlcond2 + (jl - 1), &zlcond2_index_8, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + klon) - 1), &zqxfg_index_6, 1);
                            {
                                double __in1 = zqxfg_index_6;
                                double __in2 = zlcond2_index_8;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zqxfg_slice_plus_zlcond2_slice_0 = __out;
                            }
                            {
                                double __inp = zqxfg_slice_plus_zlcond2_slice_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_713_24)
                                __out = __inp;
                                ///////////////////

                                zqxfg[((jl + klon) - 1)] = __out;
                            }

                        }
                    }
                }
            }

        }
        kfdia_plus_1_37 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_37); jl = (jl + 1)) {

            za_index_27 = za[((jl + (klon * (jk - 2))) - 1)];

            __tmp125 = (za_index_27 < yrecldp_rcldtopcf);
            za_index_28 = za[((jl + (klon * (jk - 1))) - 1)];

            __tmp126 = (za_index_28 >= yrecldp_rcldtopcf);
            if ((__tmp125 && __tmp126)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_717_20)
                        __out = 0.0;
                        ///////////////////

                        zcldtopdist[(jl - 1)] = __out;
                    }

                }
            } else {
                {
                    double zcldtopdist_index;
                    double zdp_index_2;
                    double zrho_index;
                    double zrho_slice_times_ydcst_rg;
                    double zdp_slice_div_zrho_slice_ydcst_rg;
                    double zcldtopdist_slice_plus_zdp_slice_zrho_slice_ydcst_rg;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcldtopdist + (jl - 1), &zcldtopdist_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdp + (jl - 1), &zdp_index_2, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zrho + (jl - 1), &zrho_index, 1);
                    {
                        double __in1 = zrho_index;
                        double __in2 = ydcst_rg;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zrho_slice_times_ydcst_rg = __out;
                    }
                    {
                        double __in1 = zdp_index_2;
                        double __in2 = zrho_slice_times_ydcst_rg;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zdp_slice_div_zrho_slice_ydcst_rg = __out;
                    }
                    {
                        double __in1 = zcldtopdist_index;
                        double __in2 = zdp_slice_div_zrho_slice_ydcst_rg;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zcldtopdist_slice_plus_zdp_slice_zrho_slice_ydcst_rg = __out;
                    }
                    {
                        double __inp = zcldtopdist_slice_plus_zdp_slice_zrho_slice_ydcst_rg;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_719_20)
                        __out = __inp;
                        ///////////////////

                        zcldtopdist[(jl - 1)] = __out;
                    }

                }
            }
            ztp1_index_66 = ztp1[((jl + (klon * (jk - 1))) - 1)];

            __tmp128 = (ztp1_index_66 < ydcst_rtt);
            zqxfg_index_7 = zqxfg[(jl - 1)];

            __tmp129 = (zqxfg_index_7 > yrecldp_rlmin);
            if ((__tmp128 && __tmp129)) {
                {
                    double *zicenuclei;
                    zicenuclei = new double DACE_ALIGN(64)[klon];
                    double ztp1_index_67;
                    double ztp1_slice_minus_ydcst_rtt_11;
                    double ydthf_r3ies_times_ztp1_slice_ydcst_rtt_4;
                    double ztp1_index_68;
                    double ztp1_slice_minus_ydthf_r4ies_8;
                    double ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_4;
                    double exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_4;
                    double ydthf_r2es_times_expr_4;
                    double ydthf_r2es_expr_times_ydcst_rv;
                    double zvpice;
                    double zfokoop_index_0;
                    double zvpliq;
                    double zvpliq_minus_zvpice;
                    double __tmp131;
                    double __tmp132;
                    double __tmp133;
                    double exp_12_96_zvpliq_zvpice_zvpliq_0_639;
                    double zicenuclei_slice;
                    double ztp1_index_69;
                    double ydcst_rv_times_ztp1_slice;
                    double ydcst_rlstt_div_ydcst_rv_ztp1_slice;
                    double ydcst_rlstt_ydcst_rv_ztp1_slice_minus_1_0;
                    double ydcst_rlstt_times_ydcst_rlstt_ydcst_rv_ztp1_slice_1_0;
                    double ztp1_index_70;
                    double __tmp134;
                    double zadd;
                    double ztp1_index_71;
                    double ydcst_rv_times_ztp1_slice_0;
                    double pap_index_11;
                    double ydcst_rv_ztp1_slice_times_pap_slice;
                    double __tmp135;
                    double zbdd;
                    double zicenuclei_index;
                    double zrho_index_0;
                    double zicenuclei_slice_div_zrho_slice;
                    double zicenuclei_slice_zrho_slice_pow_0_666;
                    double __tmp136;
                    double zvpliq_minus_zvpice_0;
                    double __tmp137;
                    double zadd_plus_zbdd;
                    double __tmp138;
                    double __tmp139;
                    double zcvds;
                    double zicecld_index_0;
                    double zicenuclei_index_0;
                    double zicenuclei_slice_times_yrecldp_riceinit;
                    double zrho_index_1;
                    double zicenuclei_slice_yrecldp_riceinit_div_zrho_slice;
                    double zice0;
                    double __tmp140;
                    double __tmp141;
                    double zice0_pow_0_666;
                    double __tmp142;
                    double zinew;
                    double za_index_29;
                    double zinew_minus_zice0;
                    double za_slice_times_zinew_zice0;
                    double zdepos;
                    double zqxfg_index_8;
                    double min_zdepos_zqxfg_slice;
                    double zicenuclei_index_1;
                    double zicenuclei_slice_div_15000_0;
                    double zinfactor;
                    double __tmp143;
                    double zcldtopdist_index_0;
                    double zcldtopdist_slice_div_yrecldp_rdepliqrefdepth;
                    double yrecldp_rdepliqrefrate_plus_zcldtopdist_slice_yrecldp_rdepliqrefdepth;
                    double __tmp144;
                    double zinfactor_plus_1_0_zinfactor_yrecldp_rdepliqrefrate_zcldtopdist_slice_yrecldp_rdepliqrefdepth;
                    double min_zinfactor_1_0_zinfactor_yrecldp_rdepliqrefrate_zcldtopdist_slice_yrecldp_rdepliqrefdepth_1_0;
                    double zdepos_times_expr;
                    double zsolqa_index_27;
                    double zsolqa_slice_plus_zdepos;
                    double zsolqa_index_28;
                    double zsolqa_slice_minus_zdepos;
                    double zqxfg_index_9;
                    double zqxfg_slice_plus_zdepos;
                    double zqxfg_index_10;
                    double zqxfg_slice_minus_zdepos;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_67, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_68, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_69, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_70, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_71, 1);
                    {
                        double __in2 = ztp1_index_70;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (0.024 * __in2);
                        ///////////////////

                        __tmp134 = __out;
                    }
                    {
                        double __in1 = ztp1_index_67;
                        double __in2 = ydcst_rtt;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        ztp1_slice_minus_ydcst_rtt_11 = __out;
                    }
                    {
                        double __in1 = ydthf_r3ies;
                        double __in2 = ztp1_slice_minus_ydcst_rtt_11;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydthf_r3ies_times_ztp1_slice_ydcst_rtt_4 = __out;
                    }
                    {
                        double __in1 = ztp1_index_68;
                        double __in2 = ydthf_r4ies;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        ztp1_slice_minus_ydthf_r4ies_8 = __out;
                    }
                    {
                        double __in1 = ydthf_r3ies_times_ztp1_slice_ydcst_rtt_4;
                        double __in2 = ztp1_slice_minus_ydthf_r4ies_8;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_4 = __out;
                    }
                    {
                        double __in1 = ydthf_r3ies_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4ies_4;
                        double __out;

                        ///////////////////
                        // Tasklet code (_numpy_exp_)
                        __out = exp(__in1);
                        ///////////////////

                        exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_4 = __out;
                    }
                    {
                        double __in1 = ydthf_r2es;
                        double __in2 = exp_ydthf_r3ies_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4ies_4;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydthf_r2es_times_expr_4 = __out;
                    }
                    {
                        double __in1 = ydcst_rv;
                        double __in2 = ztp1_index_69;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydcst_rv_times_ztp1_slice = __out;
                    }
                    {
                        double __in1 = ydcst_rv;
                        double __in2 = ztp1_index_71;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydcst_rv_times_ztp1_slice_0 = __out;
                    }
                    {
                        double __in1 = ydthf_r2es_times_expr_4;
                        double __in2 = ydcst_rv;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydthf_r2es_expr_times_ydcst_rv = __out;
                    }
                    {
                        double __in1 = ydthf_r2es_expr_times_ydcst_rv;
                        double __in2 = ydcst_rd;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zvpice = __out;
                    }
                    {
                        double __in2 = zvpice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (2.21 * __in2);
                        ///////////////////

                        __tmp135 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zfokoop + (jl - 1), &zfokoop_index_0, 1);
                    {
                        double __in1 = zvpice;
                        double __in2 = zfokoop_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zvpliq = __out;
                    }
                    {
                        double __in1 = zvpliq;
                        double __in2 = zvpice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zvpliq_minus_zvpice = __out;
                    }
                    {
                        double __in2 = zvpliq_minus_zvpice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (12.96 * __in2);
                        ///////////////////

                        __tmp131 = __out;
                    }
                    {
                        double __in1 = __tmp131;
                        double __in2 = zvpliq;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        __tmp132 = __out;
                    }
                    {
                        double __in1 = __tmp132;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - 0.639);
                        ///////////////////

                        __tmp133 = __out;
                    }
                    {
                        double __in1 = __tmp133;
                        double __out;

                        ///////////////////
                        // Tasklet code (_numpy_exp_)
                        __out = exp(__in1);
                        ///////////////////

                        exp_12_96_zvpliq_zvpice_zvpliq_0_639 = __out;
                    }
                    {
                        double __in2 = exp_12_96_zvpliq_zvpice_zvpliq_0_639;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (1000.0 * __in2);
                        ///////////////////

                        zicenuclei_slice = __out;
                    }
                    {
                        double __inp = zicenuclei_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_723_20)
                        __out = __inp;
                        ///////////////////

                        zicenuclei[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zicenuclei + (jl - 1), &zicenuclei_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zicenuclei + (jl - 1), &zicenuclei_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zicenuclei + (jl - 1), &zicenuclei_index_1, 1);
                    {
                        double __in1 = zicenuclei_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / 15000.0);
                        ///////////////////

                        zicenuclei_slice_div_15000_0 = __out;
                    }
                    {
                        double __in_a = zicenuclei_slice_div_15000_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, 1);
                        ///////////////////

                        zinfactor = __out;
                    }
                    {
                        double __in2 = zinfactor;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp143 = __out;
                    }
                    {
                        double __in1 = zvpliq;
                        double __in2 = zvpice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zvpliq_minus_zvpice_0 = __out;
                    }
                    {
                        double __in1 = ydcst_rlstt;
                        double __in2 = ydcst_rv_times_ztp1_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        ydcst_rlstt_div_ydcst_rv_ztp1_slice = __out;
                    }
                    {
                        double __in1 = ydcst_rlstt_div_ydcst_rv_ztp1_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - 1.0);
                        ///////////////////

                        ydcst_rlstt_ydcst_rv_ztp1_slice_minus_1_0 = __out;
                    }
                    {
                        double __in1 = ydcst_rlstt;
                        double __in2 = ydcst_rlstt_ydcst_rv_ztp1_slice_minus_1_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydcst_rlstt_times_ydcst_rlstt_ydcst_rv_ztp1_slice_1_0 = __out;
                    }
                    {
                        double __in1 = ydcst_rlstt_times_ydcst_rlstt_ydcst_rv_ztp1_slice_1_0;
                        double __in2 = __tmp134;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zadd = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_11, 1);
                    {
                        double __in1 = ydcst_rv_times_ztp1_slice_0;
                        double __in2 = pap_index_11;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydcst_rv_ztp1_slice_times_pap_slice = __out;
                    }
                    {
                        double __in1 = ydcst_rv_ztp1_slice_times_pap_slice;
                        double __in2 = __tmp135;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zbdd = __out;
                    }
                    {
                        double __in1 = zadd;
                        double __in2 = zbdd;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zadd_plus_zbdd = __out;
                    }
                    {
                        double __in2 = zadd_plus_zbdd;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (8.87 * __in2);
                        ///////////////////

                        __tmp138 = __out;
                    }
                    {
                        double __in1 = __tmp138;
                        double __in2 = zvpice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp139 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zrho + (jl - 1), &zrho_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zrho + (jl - 1), &zrho_index_1, 1);
                    {
                        double __in1 = zicenuclei_index;
                        double __in2 = zrho_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zicenuclei_slice_div_zrho_slice = __out;
                    }
                    {
                        double __in1 = zicenuclei_slice_div_zrho_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, 0.666);
                        ///////////////////

                        zicenuclei_slice_zrho_slice_pow_0_666 = __out;
                    }
                    {
                        double __in2 = zicenuclei_slice_zrho_slice_pow_0_666;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (7.8 * __in2);
                        ///////////////////

                        __tmp136 = __out;
                    }
                    {
                        double __in1 = __tmp136;
                        double __in2 = zvpliq_minus_zvpice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp137 = __out;
                    }
                    {
                        double __in1 = __tmp137;
                        double __in2 = __tmp139;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zcvds = __out;
                    }
                    {
                        double __in2 = zcvds;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (0.666 * __in2);
                        ///////////////////

                        __tmp140 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zicecld + (jl - 1), &zicecld_index_0, 1);
                    {
                        double __in1 = zicenuclei_index_0;
                        double __in2 = yrecldp_riceinit;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zicenuclei_slice_times_yrecldp_riceinit = __out;
                    }
                    {
                        double __in1 = zicenuclei_slice_times_yrecldp_riceinit;
                        double __in2 = zrho_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zicenuclei_slice_yrecldp_riceinit_div_zrho_slice = __out;
                    }
                    {
                        double __in_a = zicecld_index_0;
                        double __in_b = zicenuclei_slice_yrecldp_riceinit_div_zrho_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        zice0 = __out;
                    }
                    {
                        double __in1 = zice0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, 0.666);
                        ///////////////////

                        zice0_pow_0_666 = __out;
                    }
                    {
                        double __in1 = __tmp140;
                        double __in2 = ptsphy;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp141 = __out;
                    }
                    {
                        double __in1 = __tmp141;
                        double __in2 = zice0_pow_0_666;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        __tmp142 = __out;
                    }
                    {
                        double __in1 = __tmp142;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, 1.5);
                        ///////////////////

                        zinew = __out;
                    }
                    {
                        double __in1 = zinew;
                        double __in2 = zice0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zinew_minus_zice0 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_29, 1);
                    {
                        double __in1 = za_index_29;
                        double __in2 = zinew_minus_zice0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        za_slice_times_zinew_zice0 = __out;
                    }
                    {
                        double __in_a = za_slice_times_zinew_zice0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        zdepos = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + (jl - 1), &zqxfg_index_8, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + klon) - 1), &zqxfg_index_9, 1);
                    {
                        double __in_a = zdepos;
                        double __in_b = zqxfg_index_8;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        min_zdepos_zqxfg_slice = __out;
                    }
                    {
                        double __inp = min_zdepos_zqxfg_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_730_20)
                        __out = __inp;
                        ///////////////////

                        zdepos = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcldtopdist + (jl - 1), &zcldtopdist_index_0, 1);
                    {
                        double __in1 = zcldtopdist_index_0;
                        double __in2 = yrecldp_rdepliqrefdepth;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zcldtopdist_slice_div_yrecldp_rdepliqrefdepth = __out;
                    }
                    {
                        double __in1 = yrecldp_rdepliqrefrate;
                        double __in2 = zcldtopdist_slice_div_yrecldp_rdepliqrefdepth;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        yrecldp_rdepliqrefrate_plus_zcldtopdist_slice_yrecldp_rdepliqrefdepth = __out;
                    }
                    {
                        double __in1 = __tmp143;
                        double __in2 = yrecldp_rdepliqrefrate_plus_zcldtopdist_slice_yrecldp_rdepliqrefdepth;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp144 = __out;
                    }
                    {
                        double __in1 = zinfactor;
                        double __in2 = __tmp144;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zinfactor_plus_1_0_zinfactor_yrecldp_rdepliqrefrate_zcldtopdist_slice_yrecldp_rdepliqrefdepth = __out;
                    }
                    {
                        double __in_a = zinfactor_plus_1_0_zinfactor_yrecldp_rdepliqrefrate_zcldtopdist_slice_yrecldp_rdepliqrefdepth;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, 1);
                        ///////////////////

                        min_zinfactor_1_0_zinfactor_yrecldp_rdepliqrefrate_zcldtopdist_slice_yrecldp_rdepliqrefdepth_1_0 = __out;
                    }
                    {
                        double __in1 = zdepos;
                        double __in2 = min_zinfactor_1_0_zinfactor_yrecldp_rdepliqrefrate_zcldtopdist_slice_yrecldp_rdepliqrefdepth_1_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zdepos_times_expr = __out;
                    }
                    {
                        double __inp = zdepos_times_expr;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_732_20)
                        __out = __inp;
                        ///////////////////

                        zdepos = __out;
                    }
                    {
                        double __in1 = zqxfg_index_9;
                        double __in2 = zdepos;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zqxfg_slice_plus_zdepos = __out;
                    }
                    {
                        double __inp = zqxfg_slice_plus_zdepos;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_735_20)
                        __out = __inp;
                        ///////////////////

                        zqxfg[((jl + klon) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + (jl - 1), &zqxfg_index_10, 1);
                    {
                        double __in1 = zqxfg_index_10;
                        double __in2 = zdepos;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqxfg_slice_minus_zdepos = __out;
                    }
                    {
                        double __inp = zqxfg_slice_minus_zdepos;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_736_20)
                        __out = __inp;
                        ///////////////////

                        zqxfg[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + klon) - 1), &zsolqa_index_27, 1);
                    {
                        double __in1 = zsolqa_index_27;
                        double __in2 = zdepos;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqa_slice_plus_zdepos = __out;
                    }
                    {
                        double __inp = zsolqa_slice_plus_zdepos;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_733_20)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + klon) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (5 * klon)) - 1), &zsolqa_index_28, 1);
                    {
                        double __in1 = zsolqa_index_28;
                        double __in2 = zdepos;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolqa_slice_minus_zdepos = __out;
                    }
                    {
                        double __inp = zsolqa_slice_minus_zdepos;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_734_20)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (5 * klon)) - 1)] = __out;
                    }
                    delete[] zicenuclei;

                }
            }

        }
        kfdia_plus_1_39 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_39); jl = (jl + 1)) {
            {
                double za_index_33;
                double max_za_slice_zepsec_0;
                double ztmpa_0;
                double zqxfg_index_15;
                double zliqcld_slice_0;
                double zqxfg_index_16;
                double zicecld_slice_0;
                double zliqcld_index_0;
                double zicecld_index_2;
                double zlicld_slice_0;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                za + ((jl + (klon * (jk - 1))) - 1), &za_index_33, 1);
                {
                    double __in_a = za_index_33;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, 1e-14);
                    ///////////////////

                    max_za_slice_zepsec_0 = __out;
                }
                {
                    double __in2 = max_za_slice_zepsec_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    ztmpa_0 = __out;
                }
                {
                    double __inp = ztmpa_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_765_12)
                    __out = __inp;
                    ///////////////////

                    ztmpa = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxfg + (jl - 1), &zqxfg_index_15, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxfg + ((jl + klon) - 1), &zqxfg_index_16, 1);
                {
                    double __in1 = zqxfg_index_15;
                    double __in2 = ztmpa;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zliqcld_slice_0 = __out;
                }
                {
                    double __inp = zliqcld_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_766_12)
                    __out = __inp;
                    ///////////////////

                    zliqcld[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zliqcld + (jl - 1), &zliqcld_index_0, 1);
                {
                    double __in1 = zqxfg_index_16;
                    double __in2 = ztmpa;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zicecld_slice_0 = __out;
                }
                {
                    double __inp = zicecld_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_767_12)
                    __out = __inp;
                    ///////////////////

                    zicecld[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zicecld + (jl - 1), &zicecld_index_2, 1);
                {
                    double __in1 = zliqcld_index_0;
                    double __in2 = zicecld_index_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zlicld_slice_0 = __out;
                }
                {
                    double __inp = zlicld_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_768_12)
                    __out = __inp;
                    ///////////////////

                    zlicld[(jl - 1)] = __out;
                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            llfall_index_1 = llfall[(jm - 1)];
            if ((llfall_index_1 || (jm == 2))) {

                kfdia_plus_1_40 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_40); jl = (jl + 1)) {
                    if ((jk > yrecldp_ncldtop)) {
                        {
                            double zpfplsx_index;
                            double zdtgdp_index_3;
                            double zfallsrce_slice;
                            double zsolqa_index_31;
                            double zfallsrce_index;
                            double zsolqa_slice_plus_zfallsrce_slice;
                            double zqxfg_index_17;
                            double zfallsrce_index_0;
                            double zqxfg_slice_plus_zfallsrce_slice;
                            double zqpretot_index;
                            double zqxfg_index_18;
                            double zqpretot_slice_plus_zqxfg_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zpfplsx + (((jl + (klon * (jk - 1))) + ((klon * (jm - 1)) * (klev + 1))) - 1), &zpfplsx_index, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zdtgdp + (jl - 1), &zdtgdp_index_3, 1);
                            {
                                double __in1 = zpfplsx_index;
                                double __in2 = zdtgdp_index_3;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zfallsrce_slice = __out;
                            }
                            {
                                double __inp = zfallsrce_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_773_24)
                                __out = __inp;
                                ///////////////////

                                zfallsrce[((jl + (klon * (jm - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zfallsrce + ((jl + (klon * (jm - 1))) - 1), &zfallsrce_index, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zfallsrce + ((jl + (klon * (jm - 1))) - 1), &zfallsrce_index_0, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + ((jl + ((6 * klon) * (jm - 1))) - 1), &zsolqa_index_31, 1);
                            {
                                double __in1 = zsolqa_index_31;
                                double __in2 = zfallsrce_index;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zsolqa_slice_plus_zfallsrce_slice = __out;
                            }
                            {
                                double __inp = zsolqa_slice_plus_zfallsrce_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_774_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[((jl + ((6 * klon) * (jm - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + (klon * (jm - 1))) - 1), &zqxfg_index_17, 1);
                            {
                                double __in1 = zqxfg_index_17;
                                double __in2 = zfallsrce_index_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zqxfg_slice_plus_zfallsrce_slice = __out;
                            }
                            {
                                double __inp = zqxfg_slice_plus_zfallsrce_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_775_24)
                                __out = __inp;
                                ///////////////////

                                zqxfg[((jl + (klon * (jm - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + (klon * (jm - 1))) - 1), &zqxfg_index_18, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqpretot + (jl - 1), &zqpretot_index, 1);
                            {
                                double __in1 = zqpretot_index;
                                double __in2 = zqxfg_index_18;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zqpretot_slice_plus_zqxfg_slice = __out;
                            }
                            {
                                double __inp = zqpretot_slice_plus_zqxfg_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_776_24)
                                __out = __inp;
                                ///////////////////

                                zqpretot[(jl - 1)] = __out;
                            }

                        }
                    }
                    if ((yrecldp_laericesed && (jm == 2))) {
                        {
                            double zre_ice;
                            double zre_ice_pow_1_0;
                            double zvqx_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            pre_ice + ((jl + (klon * (jk - 1))) - 1), &zre_ice, 1);
                            {
                                double __in1 = zre_ice;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Pow_)
                                __out = (dace::math::ipow(__in1, 1));
                                ///////////////////

                                zre_ice_pow_1_0 = __out;
                            }
                            {
                                double __in2 = zre_ice_pow_1_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (0.002 * __in2);
                                ///////////////////

                                zvqx_slice = __out;
                            }
                            {
                                double __inp = zvqx_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_779_24)
                                __out = __inp;
                                ///////////////////

                                zvqx[1] = __out;
                            }

                        }
                    }
                    {
                        double zvqx_index_0;
                        double zrho_index_8;
                        double zfall;
                        double zdtgdp_index_4;
                        double zfallsink_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zvqx + (jm - 1), &zvqx_index_0, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrho + (jl - 1), &zrho_index_8, 1);
                        {
                            double __in1 = zvqx_index_0;
                            double __in2 = zrho_index_8;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zfall = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zdtgdp + (jl - 1), &zdtgdp_index_4, 1);
                        {
                            double __in1 = zdtgdp_index_4;
                            double __in2 = zfall;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zfallsink_slice = __out;
                        }
                        {
                            double __inp = zfallsink_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_781_20)
                            __out = __inp;
                            ///////////////////

                            zfallsink[((jl + (klon * (jm - 1))) - 1)] = __out;
                        }

                    }

                }
            }

        }
        kfdia_plus_1_41 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_41); jl = (jl + 1)) {

            zqpretot_index_0 = zqpretot[(jl - 1)];
            if ((zqpretot_index_0 > 1e-14)) {
                {
                    double zcovptot_index;
                    double __tmp162;
                    double za_index_34;
                    double za_index_35;
                    double max_za_slice_za_slice;
                    double __tmp163;
                    double __tmp164;
                    double za_index_36;
                    double min_za_slice_1_0_1eneg06;
                    double __tmp165;
                    double __tmp166;
                    double __tmp167;
                    double zcovptot_index_0;
                    double max_zcovptot_slice_yrecldp_rcovpmin;
                    double zcovptot_index_1;
                    double za_index_37;
                    double zcovptot_slice_minus_za_slice;
                    double zcovpclr_slice;
                    double zqxfg_index_19;
                    double zcovptot_index_2;
                    double zraincld_slice;
                    double zqxfg_index_20;
                    double zcovptot_index_3;
                    double zsnowcld_slice;
                    double zcovptot_index_4;
                    double zcovpmax_index;
                    double max_zcovptot_slice_zcovpmax_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index, 1);
                    {
                        double __in2 = zcovptot_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp162 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_34, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 2))) - 1), &za_index_35, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 2))) - 1), &za_index_36, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_37, 1);
                    {
                        double __in_a = za_index_34;
                        double __in_b = za_index_35;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        max_za_slice_za_slice = __out;
                    }
                    {
                        double __in2 = max_za_slice_za_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp163 = __out;
                    }
                    {
                        double __in1 = __tmp162;
                        double __in2 = __tmp163;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp164 = __out;
                    }
                    {
                        double __in_a = za_index_36;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, 0.999999);
                        ///////////////////

                        min_za_slice_1_0_1eneg06 = __out;
                    }
                    {
                        double __in2 = min_za_slice_1_0_1eneg06;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp165 = __out;
                    }
                    {
                        double __in1 = __tmp164;
                        double __in2 = __tmp165;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        __tmp166 = __out;
                    }
                    {
                        double __in2 = __tmp166;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (1.0 - __in2);
                        ///////////////////

                        __tmp167 = __out;
                    }
                    {
                        double __inp = __tmp167;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_784_16)
                        __out = __inp;
                        ///////////////////

                        zcovptot[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_0, 1);
                    {
                        double __in_a = zcovptot_index_0;
                        double __in_b = yrecldp_rcovpmin;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        max_zcovptot_slice_yrecldp_rcovpmin = __out;
                    }
                    {
                        double __inp = max_zcovptot_slice_yrecldp_rcovpmin;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_785_16)
                        __out = __inp;
                        ///////////////////

                        zcovptot[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_1, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_2, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_3, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_4, 1);
                    {
                        double __in1 = zcovptot_index_1;
                        double __in2 = za_index_37;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zcovptot_slice_minus_za_slice = __out;
                    }
                    {
                        double __in_b = zcovptot_slice_minus_za_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(0.0, __in_b);
                        ///////////////////

                        zcovpclr_slice = __out;
                    }
                    {
                        double __inp = zcovpclr_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_786_16)
                        __out = __inp;
                        ///////////////////

                        zcovpclr[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (2 * klon)) - 1), &zqxfg_index_19, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (3 * klon)) - 1), &zqxfg_index_20, 1);
                    {
                        double __in1 = zqxfg_index_19;
                        double __in2 = zcovptot_index_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zraincld_slice = __out;
                    }
                    {
                        double __inp = zraincld_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_787_16)
                        __out = __inp;
                        ///////////////////

                        zraincld[(jl - 1)] = __out;
                    }
                    {
                        double __in1 = zqxfg_index_20;
                        double __in2 = zcovptot_index_3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zsnowcld_slice = __out;
                    }
                    {
                        double __inp = zsnowcld_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_788_16)
                        __out = __inp;
                        ///////////////////

                        zsnowcld[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovpmax + (jl - 1), &zcovpmax_index, 1);
                    {
                        double __in_a = zcovptot_index_4;
                        double __in_b = zcovpmax_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        max_zcovptot_slice_zcovpmax_slice = __out;
                    }
                    {
                        double __inp = max_zcovptot_slice_zcovpmax_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_789_16)
                        __out = __inp;
                        ///////////////////

                        zcovpmax[(jl - 1)] = __out;
                    }

                }
            } else {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_791_16)
                        __out = 0.0;
                        ///////////////////

                        zraincld[(jl - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_792_16)
                        __out = 0.0;
                        ///////////////////

                        zsnowcld[(jl - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_793_16)
                        __out = 0.0;
                        ///////////////////

                        zcovptot[(jl - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_794_16)
                        __out = 0.0;
                        ///////////////////

                        zcovpclr[(jl - 1)] = __out;
                    }
                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_795_16)
                        __out = 0.0;
                        ///////////////////

                        zcovpmax[(jl - 1)] = __out;
                    }

                }
            }

        }
        kfdia_plus_1_42 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_42); jl = (jl + 1)) {

            ztp1_index_80 = ztp1[((jl + (klon * (jk - 1))) - 1)];

            __tmp168 = (ztp1_index_80 <= ydcst_rtt);
            if (__tmp168) {

                zicecld_index_3 = zicecld[(jl - 1)];
                if ((zicecld_index_3 > 1e-14)) {
                    {
                        double ptsphy_times_yrecldp_rsnowlin1;
                        double ztp1_index_81;
                        double ztp1_slice_minus_ydcst_rtt_13;
                        double yrecldp_rsnowlin2_times_ztp1_slice_ydcst_rtt;
                        double exp_yrecldp_rsnowlin2_ztp1_slice_ydcst_rtt;

                        {
                            double __in1 = ptsphy;
                            double __in2 = yrecldp_rsnowlin1;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            ptsphy_times_yrecldp_rsnowlin1 = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_81, 1);
                        {
                            double __in1 = ztp1_index_81;
                            double __in2 = ydcst_rtt;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            ztp1_slice_minus_ydcst_rtt_13 = __out;
                        }
                        {
                            double __in1 = yrecldp_rsnowlin2;
                            double __in2 = ztp1_slice_minus_ydcst_rtt_13;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            yrecldp_rsnowlin2_times_ztp1_slice_ydcst_rtt = __out;
                        }
                        {
                            double __in1 = yrecldp_rsnowlin2_times_ztp1_slice_ydcst_rtt;
                            double __out;

                            ///////////////////
                            // Tasklet code (_numpy_exp_)
                            __out = exp(__in1);
                            ///////////////////

                            exp_yrecldp_rsnowlin2_ztp1_slice_ydcst_rtt = __out;
                        }
                        {
                            double __in1 = ptsphy_times_yrecldp_rsnowlin1;
                            double __in2 = exp_yrecldp_rsnowlin2_ztp1_slice_ydcst_rtt;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zzco = __out;
                        }

                    }
                    if (yrecldp_laericeauto) {

                        zlcrit = picrit_aer[((jl + (klon * (jk - 1))) - 1)];
                        {
                            double pnice_index;
                            double yrecldp_rnice_div_pnice_slice;
                            double yrecldp_rnice_pnice_slice_pow_0_333;
                            double zzco_times_yrecldp_rnice_pnice_slice_0_333;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            pnice + ((jl + (klon * (jk - 1))) - 1), &pnice_index, 1);
                            {
                                double __in1 = yrecldp_rnice;
                                double __in2 = pnice_index;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Div_)
                                __out = (__in1 / __in2);
                                ///////////////////

                                yrecldp_rnice_div_pnice_slice = __out;
                            }
                            {
                                double __in1 = yrecldp_rnice_div_pnice_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Pow_)
                                __out = dace::math::pow(__in1, 0.333);
                                ///////////////////

                                yrecldp_rnice_pnice_slice_pow_0_333 = __out;
                            }
                            {
                                double __in1 = zzco;
                                double __in2 = yrecldp_rnice_pnice_slice_pow_0_333;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zzco_times_yrecldp_rnice_pnice_slice_0_333 = __out;
                            }
                            {
                                double __inp = zzco_times_yrecldp_rnice_pnice_slice_0_333;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_802_24)
                                __out = __inp;
                                ///////////////////

                                zzco = __out;
                            }

                        }
                    } else {

                        zlcrit = yrecldp_rlcritsnow;

                    }
                    {
                        double zicecld_index_4;
                        double zicecld_slice_div_zlcrit;
                        double zicecld_slice_zlcrit_pow_2;
                        double neg_zicecld_slice_zlcrit_2;
                        double exp_expr;
                        double __tmp170;
                        double zsnowaut_slice;
                        double zsolqb_index_1;
                        double zsnowaut_index;
                        double zsolqb_slice_plus_zsnowaut_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zicecld + (jl - 1), &zicecld_index_4, 1);
                        {
                            double __in1 = zicecld_index_4;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / zlcrit);
                            ///////////////////

                            zicecld_slice_div_zlcrit = __out;
                        }
                        {
                            double __in1 = zicecld_slice_div_zlcrit;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = (dace::math::ipow(__in1, 2));
                            ///////////////////

                            zicecld_slice_zlcrit_pow_2 = __out;
                        }
                        {
                            double __in = zicecld_slice_zlcrit_pow_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_USub_)
                            __out = (- __in);
                            ///////////////////

                            neg_zicecld_slice_zlcrit_2 = __out;
                        }
                        {
                            double __in1 = neg_zicecld_slice_zlcrit_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_numpy_exp_)
                            __out = exp(__in1);
                            ///////////////////

                            exp_expr = __out;
                        }
                        {
                            double __in2 = exp_expr;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (1.0 - __in2);
                            ///////////////////

                            __tmp170 = __out;
                        }
                        {
                            double __in1 = zzco;
                            double __in2 = __tmp170;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zsnowaut_slice = __out;
                        }
                        {
                            double __inp = zsnowaut_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_805_20)
                            __out = __inp;
                            ///////////////////

                            zsnowaut[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsnowaut + (jl - 1), &zsnowaut_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqb + ((jl + (8 * klon)) - 1), &zsolqb_index_1, 1);
                        {
                            double __in1 = zsolqb_index_1;
                            double __in2 = zsnowaut_index;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqb_slice_plus_zsnowaut_slice = __out;
                        }
                        {
                            double __inp = zsolqb_slice_plus_zsnowaut_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_806_20)
                            __out = __inp;
                            ///////////////////

                            zsolqb[((jl + (8 * klon)) - 1)] = __out;
                        }

                    }
                }
            }
            zliqcld_index_1 = zliqcld[(jl - 1)];
            if ((zliqcld_index_1 > 1e-14)) {

                plsm_index_0 = plsm[(jl - 1)];
                if ((plsm_index_0 > 0.5)) {

                    zlcrit = yrecldp_rclcrit_land;
                    {

                        {
                            double __inp = yrecldp_rcl_kk_cloud_num_land;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_833_24)
                            __out = __inp;
                            ///////////////////

                            zconst = __out;
                        }

                    }
                } else {

                    zlcrit = yrecldp_rclcrit_sea;
                    {

                        {
                            double __inp = yrecldp_rcl_kk_cloud_num_sea;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_836_24)
                            __out = __inp;
                            ///////////////////

                            zconst = __out;
                        }

                    }
                }
                zliqcld_index_4 = zliqcld[(jl - 1)];
                if ((zliqcld_index_4 > zlcrit)) {
                    {
                        double za_index_38;
                        double __tmp179;
                        double __tmp180;
                        double __tmp181;
                        double zliqcld_index_5;
                        double zliqcld_slice_pow_yrecldp_rcl_kkbauq;
                        double __tmp182;
                        double zconst_pow_yrecldp_rcl_kkbaun;
                        double zrainaut_slice_0;
                        double zrainaut_index_1;
                        double zqxfg_index_21;
                        double min_zrainaut_slice_zqxfg_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_38, 1);
                        {
                            double __in2 = za_index_38;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (1.5 * __in2);
                            ///////////////////

                            __tmp179 = __out;
                        }
                        {
                            double __in1 = __tmp179;
                            double __in2 = ptsphy;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp180 = __out;
                        }
                        {
                            double __in1 = __tmp180;
                            double __in2 = yrecldp_rcl_kkaau;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp181 = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zliqcld + (jl - 1), &zliqcld_index_5, 1);
                        {
                            double __in1 = zliqcld_index_5;
                            double __in2 = yrecldp_rcl_kkbauq;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = dace::math::pow(__in1, __in2);
                            ///////////////////

                            zliqcld_slice_pow_yrecldp_rcl_kkbauq = __out;
                        }
                        {
                            double __in1 = __tmp181;
                            double __in2 = zliqcld_slice_pow_yrecldp_rcl_kkbauq;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp182 = __out;
                        }
                        {
                            double __in1 = zconst;
                            double __in2 = yrecldp_rcl_kkbaun;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = dace::math::pow(__in1, __in2);
                            ///////////////////

                            zconst_pow_yrecldp_rcl_kkbaun = __out;
                        }
                        {
                            double __in1 = __tmp182;
                            double __in2 = zconst_pow_yrecldp_rcl_kkbaun;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zrainaut_slice_0 = __out;
                        }
                        {
                            double __inp = zrainaut_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_839_24)
                            __out = __inp;
                            ///////////////////

                            zrainaut[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainaut + (jl - 1), &zrainaut_index_1, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxfg + (jl - 1), &zqxfg_index_21, 1);
                        {
                            double __in_a = zrainaut_index_1;
                            double __in_b = zqxfg_index_21;
                            double __out;

                            ///////////////////
                            // Tasklet code (__min2)
                            __out = min(__in_a, __in_b);
                            ///////////////////

                            min_zrainaut_slice_zqxfg_slice = __out;
                        }
                        {
                            double __inp = min_zrainaut_slice_zqxfg_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_840_24)
                            __out = __inp;
                            ///////////////////

                            zrainaut[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainaut + (jl - 1), &zrainaut_index_2, 1);

                    }
                    __tmp183 = (zrainaut_index_2 < 1e-14);
                    if (__tmp183) {
                        {

                            {
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_842_28)
                                __out = 0.0;
                                ///////////////////

                                zrainaut[(jl - 1)] = __out;
                            }

                        }
                    }
                    {
                        double za_index_39;
                        double __tmp184;
                        double __tmp185;
                        double __tmp186;
                        double zliqcld_index_6;
                        double zraincld_index;
                        double zliqcld_slice_times_zraincld_slice;
                        double zliqcld_slice_zraincld_slice_pow_yrecldp_rcl_kkbac;
                        double zrainacc_slice;
                        double zrainacc_index;
                        double zqxfg_index_22;
                        double min_zrainacc_slice_zqxfg_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        za + ((jl + (klon * (jk - 1))) - 1), &za_index_39, 1);
                        {
                            double __in2 = za_index_39;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (2.0 * __in2);
                            ///////////////////

                            __tmp184 = __out;
                        }
                        {
                            double __in1 = __tmp184;
                            double __in2 = ptsphy;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp185 = __out;
                        }
                        {
                            double __in1 = __tmp185;
                            double __in2 = yrecldp_rcl_kkaac;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp186 = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zliqcld + (jl - 1), &zliqcld_index_6, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zraincld + (jl - 1), &zraincld_index, 1);
                        {
                            double __in1 = zliqcld_index_6;
                            double __in2 = zraincld_index;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zliqcld_slice_times_zraincld_slice = __out;
                        }
                        {
                            double __in1 = zliqcld_slice_times_zraincld_slice;
                            double __in2 = yrecldp_rcl_kkbac;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = dace::math::pow(__in1, __in2);
                            ///////////////////

                            zliqcld_slice_zraincld_slice_pow_yrecldp_rcl_kkbac = __out;
                        }
                        {
                            double __in1 = __tmp186;
                            double __in2 = zliqcld_slice_zraincld_slice_pow_yrecldp_rcl_kkbac;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zrainacc_slice = __out;
                        }
                        {
                            double __inp = zrainacc_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_843_24)
                            __out = __inp;
                            ///////////////////

                            zrainacc[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainacc + (jl - 1), &zrainacc_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxfg + (jl - 1), &zqxfg_index_22, 1);
                        {
                            double __in_a = zrainacc_index;
                            double __in_b = zqxfg_index_22;
                            double __out;

                            ///////////////////
                            // Tasklet code (__min2)
                            __out = min(__in_a, __in_b);
                            ///////////////////

                            min_zrainacc_slice_zqxfg_slice = __out;
                        }
                        {
                            double __inp = min_zrainacc_slice_zqxfg_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_844_24)
                            __out = __inp;
                            ///////////////////

                            zrainacc[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainacc + (jl - 1), &zrainacc_index_0, 1);

                    }
                    __tmp187 = (zrainacc_index_0 < 1e-14);
                    if (__tmp187) {
                        {

                            {
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_846_28)
                                __out = 0.0;
                                ///////////////////

                                zrainacc[(jl - 1)] = __out;
                            }

                        }
                    }
                } else {
                    {

                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_848_24)
                            __out = 0.0;
                            ///////////////////

                            zrainaut[(jl - 1)] = __out;
                        }
                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_849_24)
                            __out = 0.0;
                            ///////////////////

                            zrainacc[(jl - 1)] = __out;
                        }

                    }
                }
                ztp1_index_83 = ztp1[((jl + (klon * (jk - 1))) - 1)];

                __tmp188 = (ztp1_index_83 <= ydcst_rtt);
                if (__tmp188) {
                    {
                        double zsolqa_index_32;
                        double zrainaut_index_3;
                        double zsolqa_slice_plus_zrainaut_slice;
                        double zsolqa_index_33;
                        double zrainacc_index_1;
                        double zsolqa_slice_plus_zrainacc_slice;
                        double zsolqa_index_34;
                        double zrainaut_index_4;
                        double zsolqa_slice_minus_zrainaut_slice;
                        double zsolqa_index_35;
                        double zrainacc_index_2;
                        double zsolqa_slice_minus_zrainacc_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (3 * klon)) - 1), &zsolqa_index_32, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainaut + (jl - 1), &zrainaut_index_3, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainaut + (jl - 1), &zrainaut_index_4, 1);
                        {
                            double __in1 = zsolqa_index_32;
                            double __in2 = zrainaut_index_3;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zrainaut_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zrainaut_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_851_24)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (3 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (3 * klon)) - 1), &zsolqa_index_33, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainacc + (jl - 1), &zrainacc_index_1, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainacc + (jl - 1), &zrainacc_index_2, 1);
                        {
                            double __in1 = zsolqa_index_33;
                            double __in2 = zrainacc_index_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zrainacc_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zrainacc_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_852_24)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (3 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (15 * klon)) - 1), &zsolqa_index_34, 1);
                        {
                            double __in1 = zsolqa_index_34;
                            double __in2 = zrainaut_index_4;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsolqa_slice_minus_zrainaut_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_zrainaut_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_853_24)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (15 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (15 * klon)) - 1), &zsolqa_index_35, 1);
                        {
                            double __in1 = zsolqa_index_35;
                            double __in2 = zrainacc_index_2;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsolqa_slice_minus_zrainacc_slice = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_zrainacc_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_854_24)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (15 * klon)) - 1)] = __out;
                        }

                    }
                } else {
                    {
                        double zsolqa_index_36;
                        double zrainaut_index_5;
                        double zsolqa_slice_plus_zrainaut_slice_0;
                        double zsolqa_index_37;
                        double zrainacc_index_3;
                        double zsolqa_slice_plus_zrainacc_slice_0;
                        double zsolqa_index_38;
                        double zrainaut_index_6;
                        double zsolqa_slice_minus_zrainaut_slice_0;
                        double zsolqa_index_39;
                        double zrainacc_index_4;
                        double zsolqa_slice_minus_zrainacc_slice_0;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (2 * klon)) - 1), &zsolqa_index_36, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainaut + (jl - 1), &zrainaut_index_5, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainaut + (jl - 1), &zrainaut_index_6, 1);
                        {
                            double __in1 = zsolqa_index_36;
                            double __in2 = zrainaut_index_5;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zrainaut_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zrainaut_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_856_24)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (2 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (2 * klon)) - 1), &zsolqa_index_37, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainacc + (jl - 1), &zrainacc_index_3, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrainacc + (jl - 1), &zrainacc_index_4, 1);
                        {
                            double __in1 = zsolqa_index_37;
                            double __in2 = zrainacc_index_3;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqa_slice_plus_zrainacc_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_zrainacc_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_857_24)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (2 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (10 * klon)) - 1), &zsolqa_index_38, 1);
                        {
                            double __in1 = zsolqa_index_38;
                            double __in2 = zrainaut_index_6;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsolqa_slice_minus_zrainaut_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_zrainaut_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_858_24)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (10 * klon)) - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + ((jl + (10 * klon)) - 1), &zsolqa_index_39, 1);
                        {
                            double __in1 = zsolqa_index_39;
                            double __in2 = zrainacc_index_4;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsolqa_slice_minus_zrainacc_slice_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_zrainacc_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_859_24)
                            __out = __inp;
                            ///////////////////

                            zsolqa[((jl + (10 * klon)) - 1)] = __out;
                        }

                    }
                }
            }

        }
        kfdia_plus_1_43 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_43); jl = (jl + 1)) {

            ztp1_index_84 = ztp1[((jl + (klon * (jk - 1))) - 1)];

            __tmp189 = (ztp1_index_84 <= ydcst_rtt);
            zliqcld_index_7 = zliqcld[(jl - 1)];
            if ((__tmp189 && (zliqcld_index_7 > 1e-14))) {

                zsnowcld_index = zsnowcld[(jl - 1)];
                zcovptot_index_6 = zcovptot[(jl - 1)];
                {
                    double zrho_index_9;
                    double yrecldp_rdensref_div_zrho_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zrho + (jl - 1), &zrho_index_9, 1);
                    {
                        double __in1 = yrecldp_rdensref;
                        double __in2 = zrho_index_9;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        yrecldp_rdensref_div_zrho_slice = __out;
                    }
                    {
                        double __in1 = yrecldp_rdensref_div_zrho_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, 0.4);
                        ///////////////////

                        zfallcorr = __out;
                    }

                }
                if (((zsnowcld_index > 1e-14) && (zcovptot_index_6 > 0.01))) {
                    {
                        double *zsnowrime;
                        zsnowrime = new double DACE_ALIGN(64)[klon];
                        double zcovptot_index_7;
                        double __tmp195;
                        double __tmp196;
                        double __tmp197;
                        double __tmp198;
                        double zrho_index_10;
                        double zsnowcld_index_0;
                        double zrho_slice_times_zsnowcld_slice;
                        double zrho_slice_zsnowcld_slice_times_yrecldp_rcl_const1s;
                        double zrho_slice_zsnowcld_slice_yrecldp_rcl_const1s_pow_yrecldp_rcl_const8s;
                        double zsnowrime_slice;
                        double zsnowrime_index;
                        double min_zsnowrime_slice_1_0;
                        double zsolqb_index_4;
                        double zsnowrime_index_0;
                        double zsolqb_slice_plus_zsnowrime_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zcovptot + (jl - 1), &zcovptot_index_7, 1);
                        {
                            double __in2 = zcovptot_index_7;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (0.3 * __in2);
                            ///////////////////

                            __tmp195 = __out;
                        }
                        {
                            double __in1 = __tmp195;
                            double __in2 = ptsphy;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp196 = __out;
                        }
                        {
                            double __in1 = __tmp196;
                            double __in2 = yrecldp_rcl_const7s;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp197 = __out;
                        }
                        {
                            double __in1 = __tmp197;
                            double __in2 = zfallcorr;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            __tmp198 = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zrho + (jl - 1), &zrho_index_10, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsnowcld + (jl - 1), &zsnowcld_index_0, 1);
                        {
                            double __in1 = zrho_index_10;
                            double __in2 = zsnowcld_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zrho_slice_times_zsnowcld_slice = __out;
                        }
                        {
                            double __in1 = zrho_slice_times_zsnowcld_slice;
                            double __in2 = yrecldp_rcl_const1s;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zrho_slice_zsnowcld_slice_times_yrecldp_rcl_const1s = __out;
                        }
                        {
                            double __in1 = zrho_slice_zsnowcld_slice_times_yrecldp_rcl_const1s;
                            double __in2 = yrecldp_rcl_const8s;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Pow_)
                            __out = dace::math::pow(__in1, __in2);
                            ///////////////////

                            zrho_slice_zsnowcld_slice_yrecldp_rcl_const1s_pow_yrecldp_rcl_const8s = __out;
                        }
                        {
                            double __in1 = __tmp198;
                            double __in2 = zrho_slice_zsnowcld_slice_yrecldp_rcl_const1s_pow_yrecldp_rcl_const8s;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zsnowrime_slice = __out;
                        }
                        {
                            double __inp = zsnowrime_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_865_24)
                            __out = __inp;
                            ///////////////////

                            zsnowrime[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsnowrime + (jl - 1), &zsnowrime_index, 1);
                        {
                            double __in_a = zsnowrime_index;
                            double __out;

                            ///////////////////
                            // Tasklet code (__min2)
                            __out = min(__in_a, 1);
                            ///////////////////

                            min_zsnowrime_slice_1_0 = __out;
                        }
                        {
                            double __inp = min_zsnowrime_slice_1_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_866_24)
                            __out = __inp;
                            ///////////////////

                            zsnowrime[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsnowrime + (jl - 1), &zsnowrime_index_0, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqb + ((jl + (3 * klon)) - 1), &zsolqb_index_4, 1);
                        {
                            double __in1 = zsolqb_index_4;
                            double __in2 = zsnowrime_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zsolqb_slice_plus_zsnowrime_slice = __out;
                        }
                        {
                            double __inp = zsolqb_slice_plus_zsnowrime_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_867_24)
                            __out = __inp;
                            ///////////////////

                            zsolqb[((jl + (3 * klon)) - 1)] = __out;
                        }
                        delete[] zsnowrime;

                    }
                }
            }

        }
        kfdia_plus_1_44 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_44); jl = (jl + 1)) {
            {
                double zqxfg_index_23;
                double zqxfg_index_24;
                double zicetot_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxfg + ((jl + klon) - 1), &zqxfg_index_23, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxfg + ((jl + (3 * klon)) - 1), &zqxfg_index_24, 1);
                {
                    double __in1 = zqxfg_index_23;
                    double __in2 = zqxfg_index_24;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zicetot_slice = __out;
                }
                {
                    double __inp = zicetot_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_869_12)
                    __out = __inp;
                    ///////////////////

                    zicetot[(jl - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zicetot + (jl - 1), &zicetot_index, 1);
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_870_12)
                    __out = 0.0;
                    ///////////////////

                    zmeltmax[(jl - 1)] = __out;
                }

            }
            __tmp199 = (zicetot_index > 1e-14);
            ztp1_index_85 = ztp1[((jl + (klon * (jk - 1))) - 1)];

            __tmp200 = (ztp1_index_85 > ydcst_rtt);
            if ((__tmp199 && __tmp200)) {
                {
                    double zqsice_index_11;
                    double zqx_index_36;
                    double zqsice_slice_minus_zqx_slice_0;
                    double ztp1_index_86;
                    double ztp1_slice_minus_ydcst_rtt_14;
                    double pap_index_13;
                    double pap_slice_minus_ztw3;
                    double ztw2_times_pap_slice_ztw3;
                    double ztw1_plus_ztw2_pap_slice_ztw3;
                    double ztp1_index_87;
                    double ztp1_slice_minus_ztw5;
                    double ztw4_times_ztp1_slice_ztw5;
                    double ztw1_ztw2_pap_slice_ztw3_minus_ztw4_ztp1_slice_ztw5;
                    double zsubsat_times_ztw1_ztw2_pap_slice_ztw3_ztw4_ztp1_slice_ztw5;
                    double ztdmtw0;
                    double __tmp202;
                    double __tmp203;
                    double ptsphy_times_1_0_0_5_ztdmtw0;
                    double ptsphy_1_0_0_5_ztdmtw0_div_yrecldp_rtaumel;
                    double ztdmtw0_times_zcons1;
                    double ztdmtw0_zcons1_times_zrldcp;
                    double zmeltmax_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_11, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_36, 1);
                    {
                        double __in1 = zqsice_index_11;
                        double __in2 = zqx_index_36;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqsice_slice_minus_zqx_slice_0 = __out;
                    }
                    {
                        double __in_a = zqsice_slice_minus_zqx_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        zsubsat = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_86, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_87, 1);
                    {
                        double __in1 = ztp1_index_86;
                        double __in2 = ydcst_rtt;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        ztp1_slice_minus_ydcst_rtt_14 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_13, 1);
                    {
                        double __in1 = pap_index_13;
                        double __in2 = ztw3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        pap_slice_minus_ztw3 = __out;
                    }
                    {
                        double __in1 = ztw2;
                        double __in2 = pap_slice_minus_ztw3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ztw2_times_pap_slice_ztw3 = __out;
                    }
                    {
                        double __in1 = ztw1;
                        double __in2 = ztw2_times_pap_slice_ztw3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        ztw1_plus_ztw2_pap_slice_ztw3 = __out;
                    }
                    {
                        double __in1 = ztp1_index_87;
                        double __in2 = ztw5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        ztp1_slice_minus_ztw5 = __out;
                    }
                    {
                        double __in1 = ztw4;
                        double __in2 = ztp1_slice_minus_ztw5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ztw4_times_ztp1_slice_ztw5 = __out;
                    }
                    {
                        double __in1 = ztw1_plus_ztw2_pap_slice_ztw3;
                        double __in2 = ztw4_times_ztp1_slice_ztw5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        ztw1_ztw2_pap_slice_ztw3_minus_ztw4_ztp1_slice_ztw5 = __out;
                    }
                    {
                        double __in1 = zsubsat;
                        double __in2 = ztw1_ztw2_pap_slice_ztw3_minus_ztw4_ztp1_slice_ztw5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zsubsat_times_ztw1_ztw2_pap_slice_ztw3_ztw4_ztp1_slice_ztw5 = __out;
                    }
                    {
                        double __in1 = ztp1_slice_minus_ydcst_rtt_14;
                        double __in2 = zsubsat_times_ztw1_ztw2_pap_slice_ztw3_ztw4_ztp1_slice_ztw5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        ztdmtw0 = __out;
                    }
                    {
                        double __in2 = ztdmtw0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (0.5 * __in2);
                        ///////////////////

                        __tmp202 = __out;
                    }
                    {
                        double __in2 = __tmp202;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (1.0 + __in2);
                        ///////////////////

                        __tmp203 = __out;
                    }
                    {
                        double __in1 = ptsphy;
                        double __in2 = __tmp203;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ptsphy_times_1_0_0_5_ztdmtw0 = __out;
                    }
                    {
                        double __in1 = ptsphy_times_1_0_0_5_ztdmtw0;
                        double __in2 = yrecldp_rtaumel;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        ptsphy_1_0_0_5_ztdmtw0_div_yrecldp_rtaumel = __out;
                    }
                    {
                        double __inp = ptsphy_1_0_0_5_ztdmtw0_div_yrecldp_rtaumel;
                        double __out;

                        ///////////////////
                        // Tasklet code (abs)
                        __out = abs(__inp);
                        ///////////////////

                        zcons1 = __out;
                    }
                    {
                        double __in1 = ztdmtw0;
                        double __in2 = zcons1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ztdmtw0_times_zcons1 = __out;
                    }
                    {
                        double __in1 = ztdmtw0_times_zcons1;
                        double __in2 = zrldcp;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ztdmtw0_zcons1_times_zrldcp = __out;
                    }
                    {
                        double __in_a = ztdmtw0_zcons1_times_zrldcp;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        zmeltmax_slice = __out;
                    }
                    {
                        double __inp = zmeltmax_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_875_16)
                        __out = __inp;
                        ///////////////////

                        zmeltmax[(jl - 1)] = __out;
                    }

                }
            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            iphase_index_3 = iphase[(jm - 1)];
            if ((iphase_index_3 == 2)) {

                kfdia_plus_1_45 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_45); jl = (jl + 1)) {

                    zmeltmax_index = zmeltmax[(jl - 1)];
                    zicetot_index_0 = zicetot[(jl - 1)];
                    if (((zmeltmax_index > 1e-14) && (zicetot_index_0 > 1e-14))) {
                        {
                            double zqxfg_index_25;
                            double zicetot_index_1;
                            double zalfa2;
                            double zqxfg_index_26;
                            double zmeltmax_index_0;
                            double zalfa2_times_zmeltmax_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + (klon * (jm - 1))) - 1), &zqxfg_index_25, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + (klon * (jm - 1))) - 1), &zqxfg_index_26, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zicetot + (jl - 1), &zicetot_index_1, 1);
                            {
                                double __in1 = zqxfg_index_25;
                                double __in2 = zicetot_index_1;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Div_)
                                __out = (__in1 / __in2);
                                ///////////////////

                                zalfa2 = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zmeltmax + (jl - 1), &zmeltmax_index_0, 1);
                            {
                                double __in1 = zalfa2;
                                double __in2 = zmeltmax_index_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zalfa2_times_zmeltmax_slice = __out;
                            }
                            {
                                double __in_a = zqxfg_index_26;
                                double __in_b = zalfa2_times_zmeltmax_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (__min2)
                                __out = min(__in_a, __in_b);
                                ///////////////////

                                zmelt = __out;
                            }

                        }
                        imelt_index = imelt[(jm - 1)];
                        imelt_index_0 = imelt[(jm - 1)];
                        imelt_index_1 = imelt[(jm - 1)];
                        imelt_index_2 = imelt[(jm - 1)];
                        imelt_index_3 = imelt[(jm - 1)];
                        imelt_index_4 = imelt[(jm - 1)];
                        {
                            double zqxfg_index_27;
                            double zqxfg_slice_minus_zmelt;
                            double zqxfg_index_28;
                            double zqxfg_slice_plus_zmelt;
                            double zsolqa_index_40;
                            double zsolqa_slice_plus_zmelt;
                            double zsolqa_index_41;
                            double zsolqa_slice_minus_zmelt;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + (klon * (jm - 1))) - 1), &zqxfg_index_27, 1);
                            {
                                double __in1 = zqxfg_index_27;
                                double __in2 = zmelt;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zqxfg_slice_minus_zmelt = __out;
                            }
                            {
                                double __inp = zqxfg_slice_minus_zmelt;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_882_24)
                                __out = __inp;
                                ///////////////////

                                zqxfg[((jl + (klon * (jm - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + (klon * (imelt_index - 1))) - 1), &zqxfg_index_28, 1);
                            {
                                double __in1 = zqxfg_index_28;
                                double __in2 = zmelt;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zqxfg_slice_plus_zmelt = __out;
                            }
                            {
                                double __inp = zqxfg_slice_plus_zmelt;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_883_24)
                                __out = __inp;
                                ///////////////////

                                zqxfg[((jl + (klon * (imelt_index_0 - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + (((jl + (klon * (imelt_index_1 - 1))) + ((5 * klon) * (jm - 1))) - 1), &zsolqa_index_40, 1);
                            {
                                double __in1 = zsolqa_index_40;
                                double __in2 = zmelt;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zsolqa_slice_plus_zmelt = __out;
                            }
                            {
                                double __inp = zsolqa_slice_plus_zmelt;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_884_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[(((jl + (klon * (imelt_index_2 - 1))) + ((5 * klon) * (jm - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + (((jl + ((5 * klon) * (imelt_index_3 - 1))) + (klon * (jm - 1))) - 1), &zsolqa_index_41, 1);
                            {
                                double __in1 = zsolqa_index_41;
                                double __in2 = zmelt;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zsolqa_slice_minus_zmelt = __out;
                            }
                            {
                                double __inp = zsolqa_slice_minus_zmelt;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_885_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[(((jl + ((5 * klon) * (imelt_index_4 - 1))) + (klon * (jm - 1))) - 1)] = __out;
                            }

                        }
                    }

                }
            }

        }
        kfdia_plus_1_46 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_46); jl = (jl + 1)) {

            zqx_index_37 = zqx[(((jl + ((2 * klev) * klon)) + (klon * (jk - 1))) - 1)];
            if ((zqx_index_37 > 1e-14)) {

                ztp1_index_88 = ztp1[((jl + (klon * (jk - 1))) - 1)];

                __tmp209 = (ztp1_index_88 <= ydcst_rtt);
                ztp1_index_89 = ztp1[((jl + (klon * (jk - 2))) - 1)];

                __tmp210 = (ztp1_index_89 > ydcst_rtt);
                if ((__tmp209 && __tmp210)) {
                    {
                        double zqx_index_38;
                        double zqx_index_39;
                        double zqx_slice_plus_zqx_slice_3;
                        double zqpretot_slice;
                        double zqx_index_40;
                        double zqpretot_index_1;
                        double prainfrac_toprfz_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((3 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_38, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((2 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_39, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((2 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_40, 1);
                        {
                            double __in1 = zqx_index_38;
                            double __in2 = zqx_index_39;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqx_slice_plus_zqx_slice_3 = __out;
                        }
                        {
                            double __in_a = zqx_slice_plus_zqx_slice_3;
                            double __out;

                            ///////////////////
                            // Tasklet code (__max2)
                            __out = max(__in_a, 1e-14);
                            ///////////////////

                            zqpretot_slice = __out;
                        }
                        {
                            double __inp = zqpretot_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_889_20)
                            __out = __inp;
                            ///////////////////

                            zqpretot[(jl - 1)] = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqpretot + (jl - 1), &zqpretot_index_1, 1);
                        {
                            double __in1 = zqx_index_40;
                            double __in2 = zqpretot_index_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / __in2);
                            ///////////////////

                            prainfrac_toprfz_slice = __out;
                        }
                        {
                            double __inp = prainfrac_toprfz_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_890_20)
                            __out = __inp;
                            ///////////////////

                            prainfrac_toprfz[(jl - 1)] = __out;
                        }

                    }
                }
                ztp1_index_90 = ztp1[((jl + (klon * (jk - 1))) - 1)];

                __tmp213 = (ztp1_index_90 < ydcst_rtt);
                if (__tmp213) {

                    prainfrac_toprfz_index_0 = prainfrac_toprfz[(jl - 1)];
                    if ((prainfrac_toprfz_index_0 > 0.8)) {
                        {
                            double zrho_index_11;
                            double zqx_index_41;
                            double zrho_slice_times_zqx_slice;
                            double yrecldp_rcl_fac1_div_zrho_slice_zqx_slice;
                            double ztp1_index_91;
                            double ztp1_slice_minus_ydcst_rtt_15;
                            double ztemp;
                            double zrho_index_12;
                            double yrecldp_rcl_const5r_div_zrho_slice;
                            double ptsphy_times_yrecldp_rcl_const5r_zrho_slice;
                            double exp_ztemp;
                            double expr_minus_1_0;
                            double ptsphy_yrecldp_rcl_const5r_zrho_slice_times_expr_1_0;
                            double zlambda_pow_yrecldp_rcl_const6r;
                            double zfrzmax_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zrho + (jl - 1), &zrho_index_11, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zrho + (jl - 1), &zrho_index_12, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqx + (((jl + ((2 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_41, 1);
                            {
                                double __in1 = zrho_index_11;
                                double __in2 = zqx_index_41;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zrho_slice_times_zqx_slice = __out;
                            }
                            {
                                double __in1 = yrecldp_rcl_fac1;
                                double __in2 = zrho_slice_times_zqx_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Div_)
                                __out = (__in1 / __in2);
                                ///////////////////

                                yrecldp_rcl_fac1_div_zrho_slice_zqx_slice = __out;
                            }
                            {
                                double __in1 = yrecldp_rcl_fac1_div_zrho_slice_zqx_slice;
                                double __in2 = yrecldp_rcl_fac2;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Pow_)
                                __out = dace::math::pow(__in1, __in2);
                                ///////////////////

                                zlambda = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_91, 1);
                            {
                                double __in1 = ztp1_index_91;
                                double __in2 = ydcst_rtt;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                ztp1_slice_minus_ydcst_rtt_15 = __out;
                            }
                            {
                                double __in1 = yrecldp_rcl_fzrab;
                                double __in2 = ztp1_slice_minus_ydcst_rtt_15;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                ztemp = __out;
                            }
                            {
                                double __in1 = ztemp;
                                double __out;

                                ///////////////////
                                // Tasklet code (_numpy_exp_)
                                __out = exp(__in1);
                                ///////////////////

                                exp_ztemp = __out;
                            }
                            {
                                double __in1 = exp_ztemp;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - 1.0);
                                ///////////////////

                                expr_minus_1_0 = __out;
                            }
                            {
                                double __in1 = yrecldp_rcl_const5r;
                                double __in2 = zrho_index_12;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Div_)
                                __out = (__in1 / __in2);
                                ///////////////////

                                yrecldp_rcl_const5r_div_zrho_slice = __out;
                            }
                            {
                                double __in1 = ptsphy;
                                double __in2 = yrecldp_rcl_const5r_div_zrho_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                ptsphy_times_yrecldp_rcl_const5r_zrho_slice = __out;
                            }
                            {
                                double __in1 = ptsphy_times_yrecldp_rcl_const5r_zrho_slice;
                                double __in2 = expr_minus_1_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                ptsphy_yrecldp_rcl_const5r_zrho_slice_times_expr_1_0 = __out;
                            }
                            {
                                double __in1 = zlambda;
                                double __in2 = yrecldp_rcl_const6r;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Pow_)
                                __out = dace::math::pow(__in1, __in2);
                                ///////////////////

                                zlambda_pow_yrecldp_rcl_const6r = __out;
                            }
                            {
                                double __in1 = ptsphy_yrecldp_rcl_const5r_zrho_slice_times_expr_1_0;
                                double __in2 = zlambda_pow_yrecldp_rcl_const6r;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zfrz = __out;
                            }
                            {
                                double __in_a = zfrz;
                                double __out;

                                ///////////////////
                                // Tasklet code (__max2)
                                __out = max(__in_a, 0.0);
                                ///////////////////

                                zfrzmax_slice = __out;
                            }
                            {
                                double __inp = zfrzmax_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_900_24)
                                __out = __inp;
                                ///////////////////

                                zfrzmax[(jl - 1)] = __out;
                            }

                        }
                    } else {
                        {
                            double ztp1_index_92;
                            double ydcst_rtt_minus_ztp1_slice;
                            double __tmp215;
                            double __tmp216;
                            double ptsphy_times_1_0_0_5_ydcst_rtt_ztp1_slice;
                            double ptsphy_1_0_0_5_ydcst_rtt_ztp1_slice_div_yrecldp_rtaumel;
                            double zcons1_0;
                            double ztp1_index_93;
                            double ydcst_rtt_minus_ztp1_slice_0;
                            double ydcst_rtt_ztp1_slice_times_zcons1;
                            double ydcst_rtt_ztp1_slice_zcons1_times_zrldcp;
                            double zfrzmax_slice_0;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_92, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_93, 1);
                            {
                                double __in1 = ydcst_rtt;
                                double __in2 = ztp1_index_92;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                ydcst_rtt_minus_ztp1_slice = __out;
                            }
                            {
                                double __in2 = ydcst_rtt_minus_ztp1_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (0.5 * __in2);
                                ///////////////////

                                __tmp215 = __out;
                            }
                            {
                                double __in2 = __tmp215;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (1.0 + __in2);
                                ///////////////////

                                __tmp216 = __out;
                            }
                            {
                                double __in1 = ydcst_rtt;
                                double __in2 = ztp1_index_93;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                ydcst_rtt_minus_ztp1_slice_0 = __out;
                            }
                            {
                                double __in1 = ptsphy;
                                double __in2 = __tmp216;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                ptsphy_times_1_0_0_5_ydcst_rtt_ztp1_slice = __out;
                            }
                            {
                                double __in1 = ptsphy_times_1_0_0_5_ydcst_rtt_ztp1_slice;
                                double __in2 = yrecldp_rtaumel;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Div_)
                                __out = (__in1 / __in2);
                                ///////////////////

                                ptsphy_1_0_0_5_ydcst_rtt_ztp1_slice_div_yrecldp_rtaumel = __out;
                            }
                            {
                                double __inp = ptsphy_1_0_0_5_ydcst_rtt_ztp1_slice_div_yrecldp_rtaumel;
                                double __out;

                                ///////////////////
                                // Tasklet code (abs)
                                __out = abs(__inp);
                                ///////////////////

                                zcons1_0 = __out;
                            }
                            {
                                double __inp = zcons1_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_902_24)
                                __out = __inp;
                                ///////////////////

                                zcons1 = __out;
                            }
                            {
                                double __in1 = ydcst_rtt_minus_ztp1_slice_0;
                                double __in2 = zcons1;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                ydcst_rtt_ztp1_slice_times_zcons1 = __out;
                            }
                            {
                                double __in1 = ydcst_rtt_ztp1_slice_times_zcons1;
                                double __in2 = zrldcp;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                ydcst_rtt_ztp1_slice_zcons1_times_zrldcp = __out;
                            }
                            {
                                double __in_a = ydcst_rtt_ztp1_slice_zcons1_times_zrldcp;
                                double __out;

                                ///////////////////
                                // Tasklet code (__max2)
                                __out = max(__in_a, 0.0);
                                ///////////////////

                                zfrzmax_slice_0 = __out;
                            }
                            {
                                double __inp = zfrzmax_slice_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_903_24)
                                __out = __inp;
                                ///////////////////

                                zfrzmax[(jl - 1)] = __out;
                            }

                        }
                    }
                    zfrzmax_index = zfrzmax[(jl - 1)];
                    if ((zfrzmax_index > 1e-14)) {
                        {
                            double zqx_index_42;
                            double zfrzmax_index_0;
                            double zfrz_0;
                            double zsolqa_index_42;
                            double zsolqa_slice_plus_zfrz;
                            double zsolqa_index_43;
                            double zsolqa_slice_minus_zfrz;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqx + (((jl + ((2 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_42, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zfrzmax + (jl - 1), &zfrzmax_index_0, 1);
                            {
                                double __in_a = zqx_index_42;
                                double __in_b = zfrzmax_index_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (__min2)
                                __out = min(__in_a, __in_b);
                                ///////////////////

                                zfrz_0 = __out;
                            }
                            {
                                double __inp = zfrz_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_905_24)
                                __out = __inp;
                                ///////////////////

                                zfrz = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + ((jl + (13 * klon)) - 1), &zsolqa_index_42, 1);
                            {
                                double __in1 = zsolqa_index_42;
                                double __in2 = zfrz;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zsolqa_slice_plus_zfrz = __out;
                            }
                            {
                                double __inp = zsolqa_slice_plus_zfrz;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_906_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[((jl + (13 * klon)) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + ((jl + (17 * klon)) - 1), &zsolqa_index_43, 1);
                            {
                                double __in1 = zsolqa_index_43;
                                double __in2 = zfrz;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zsolqa_slice_minus_zfrz = __out;
                            }
                            {
                                double __inp = zsolqa_slice_minus_zfrz;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_907_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[((jl + (17 * klon)) - 1)] = __out;
                            }

                        }
                    }
                }
            }

        }
        kfdia_plus_1_47 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_47); jl = (jl + 1)) {
            {
                double ztp1_index_94;
                double yrecldp_rthomo_minus_ztp1_slice;
                double yrecldp_rthomo_ztp1_slice_times_zrldcp;
                double zfrzmax_slice_1;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_94, 1);
                {
                    double __in1 = yrecldp_rthomo;
                    double __in2 = ztp1_index_94;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    yrecldp_rthomo_minus_ztp1_slice = __out;
                }
                {
                    double __in1 = yrecldp_rthomo_minus_ztp1_slice;
                    double __in2 = zrldcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    yrecldp_rthomo_ztp1_slice_times_zrldcp = __out;
                }
                {
                    double __in_a = yrecldp_rthomo_ztp1_slice_times_zrldcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, 0.0);
                    ///////////////////

                    zfrzmax_slice_1 = __out;
                }
                {
                    double __inp = zfrzmax_slice_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_909_12)
                    __out = __inp;
                    ///////////////////

                    zfrzmax[(jl - 1)] = __out;
                }

            }

        }
        kfdia_plus_1_48 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_48); jl = (jl + 1)) {

            zfrzmax_index_1 = zfrzmax[(jl - 1)];
            zqxfg_index_29 = zqxfg[(jl - 1)];
            if (((zfrzmax_index_1 > 1e-14) && (zqxfg_index_29 > 1e-14))) {
                {
                    double zqxfg_index_30;
                    double zfrzmax_index_2;
                    double zfrz_1;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + (jl - 1), &zqxfg_index_30, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zfrzmax + (jl - 1), &zfrzmax_index_2, 1);
                    {
                        double __in_a = zqxfg_index_30;
                        double __in_b = zfrzmax_index_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        zfrz_1 = __out;
                    }
                    {
                        double __inp = zfrz_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_913_16)
                        __out = __inp;
                        ///////////////////

                        zfrz = __out;
                    }

                }
                imelt_index_5 = imelt[0];
                imelt_index_6 = imelt[0];
                imelt_index_7 = imelt[0];
                imelt_index_8 = imelt[0];
                {
                    double zsolqa_index_44;
                    double zsolqa_slice_plus_zfrz_0;
                    double zsolqa_index_45;
                    double zsolqa_slice_minus_zfrz_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (klon * (imelt_index_5 - 1))) - 1), &zsolqa_index_44, 1);
                    {
                        double __in1 = zsolqa_index_44;
                        double __in2 = zfrz;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqa_slice_plus_zfrz_0 = __out;
                    }
                    {
                        double __inp = zsolqa_slice_plus_zfrz_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_914_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (klon * (imelt_index_6 - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + ((5 * klon) * (imelt_index_7 - 1))) - 1), &zsolqa_index_45, 1);
                    {
                        double __in1 = zsolqa_index_45;
                        double __in2 = zfrz;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolqa_slice_minus_zfrz_0 = __out;
                    }
                    {
                        double __inp = zsolqa_slice_minus_zfrz_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_915_16)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + ((5 * klon) * (imelt_index_8 - 1))) - 1)] = __out;
                    }

                }
            }

        }
        kfdia_plus_1_50 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_50); jl = (jl + 1)) {

            zcovpclr_index_3 = zcovpclr[(jl - 1)];
            zqxfg_index_36 = zqxfg[((jl + (2 * klon)) - 1)];
            {
                double __tmp229;
                double zcovpmax_index_1;
                double __tmp230;
                double za_index_44;
                double __tmp231;
                double max_zepsec_1_0_za_slice_4;
                double __tmp232;
                double zzrh_0;
                double max_zzrh_yrecldp_rprecrhmax_0;
                double min_expr_1_0_0;
                double min_0_8_zzrh;
                double zqx_index_44;
                double zqsliq_index_7;
                double min_zqx_slice_zqsliq_slice;
                double zqe_4;
                double zqsliq_index_8;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_44, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsliq + ((jl + (klon * (jk - 1))) - 1), &zqsliq_index_7, 1);
                {
                    double __in_a = zqx_index_44;
                    double __in_b = zqsliq_index_7;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_zqx_slice_zqsliq_slice = __out;
                }
                {
                    double __in_b = min_zqx_slice_zqsliq_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(0.0, __in_b);
                    ///////////////////

                    zqe_4 = __out;
                }
                {
                    double __inp = zqe_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_940_16)
                    __out = __inp;
                    ///////////////////

                    zqe = __out;
                }
                {
                    double __in2 = yrecldp_rprecrhmax;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp229 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zcovpmax + (jl - 1), &zcovpmax_index_1, 1);
                {
                    double __in1 = __tmp229;
                    double __in2 = zcovpmax_index_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp230 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                za + ((jl + (klon * (jk - 1))) - 1), &za_index_44, 1);
                {
                    double __in2 = za_index_44;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp231 = __out;
                }
                {
                    double __in_b = __tmp231;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(1e-14, __in_b);
                    ///////////////////

                    max_zepsec_1_0_za_slice_4 = __out;
                }
                {
                    double __in1 = __tmp230;
                    double __in2 = max_zepsec_1_0_za_slice_4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    __tmp232 = __out;
                }
                {
                    double __in1 = yrecldp_rprecrhmax;
                    double __in2 = __tmp232;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zzrh_0 = __out;
                }
                {
                    double __inp = zzrh_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_937_16)
                    __out = __inp;
                    ///////////////////

                    zzrh = __out;
                }
                {
                    double __in_a = zzrh;
                    double __in_b = yrecldp_rprecrhmax;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_zzrh_yrecldp_rprecrhmax_0 = __out;
                }
                {
                    double __in_a = max_zzrh_yrecldp_rprecrhmax_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, 1);
                    ///////////////////

                    min_expr_1_0_0 = __out;
                }
                {
                    double __inp = min_expr_1_0_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_938_16)
                    __out = __inp;
                    ///////////////////

                    zzrh = __out;
                }
                {
                    double __in_b = zzrh;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(0.8, __in_b);
                    ///////////////////

                    min_0_8_zzrh = __out;
                }
                {
                    double __inp = min_0_8_zzrh;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_939_16)
                    __out = __inp;
                    ///////////////////

                    zzrh = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsliq + ((jl + (klon * (jk - 1))) - 1), &zqsliq_index_8, 1);
                {
                    double __in1 = zzrh;
                    double __in2 = zqsliq_index_8;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zzrh_times_zqsliq_slice_0 = __out;
                }

            }
            __tmp235 = (zqe < zzrh_times_zqsliq_slice_0);
            if ((((zcovpclr_index_3 > 1e-14) && (zqxfg_index_36 > 1e-14)) && __tmp235)) {
                {
                    double zqxfg_index_37;
                    double zcovptot_index_12;
                    double zpreclr_0;
                    double zrho_index_13;
                    double yrecldp_rdensref_div_zrho_slice_0;
                    double zfallcorr_0;
                    double ydcst_rv_div_ydcst_rd;
                    double ztp1_index_95;
                    double ztp1_slice_minus_ydcst_rtt_16;
                    double ydthf_r3les_times_ztp1_slice_ydcst_rtt_5;
                    double ztp1_index_96;
                    double ztp1_slice_minus_ydthf_r4les_9;
                    double ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_5;
                    double exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_5;
                    double ydthf_r2es_times_expr_6;
                    double zesatliq;
                    double zrho_index_14;
                    double zrho_slice_times_zpreclr;
                    double yrecldp_rcl_fac1_div_zrho_slice_zpreclr;
                    double zlambda_0;
                    double yrecldp_rcl_cdenom1_times_zesatliq;
                    double ztp1_index_97;
                    double yrecldp_rcl_cdenom2_times_ztp1_slice;
                    double yrecldp_rcl_cdenom2_ztp1_slice_times_zesatliq;
                    double yrecldp_rcl_cdenom1_zesatliq_minus_yrecldp_rcl_cdenom2_ztp1_slice_zesatliq;
                    double ztp1_index_98;
                    double ztp1_slice_pow_3_0_0;
                    double yrecldp_rcl_cdenom3_times_ztp1_slice_3_0;
                    double pap_index_15;
                    double yrecldp_rcl_cdenom3_ztp1_slice_3_0_times_pap_slice;
                    double zevap_denom;
                    double ztp1_index_99;
                    double ztp1_slice_div_273_0_0;
                    double ztp1_slice_273_0_pow_1_5_0;
                    double ztp1_slice_273_0_1_5_times_393_0;
                    double ztp1_index_100;
                    double ztp1_slice_plus_120_0_0;
                    double zcorr2;
                    double zqsliq_index_9;
                    double zzrh_times_zqsliq_slice_1;
                    double zzrh_zqsliq_slice_minus_zqe;
                    double zsubsat_0;
                    double zqsliq_index_10;
                    double __tmp236;
                    double ztp1_index_101;
                    double ztp1_slice_pow_2_0_0;
                    double __tmp237;
                    double __tmp238;
                    double __tmp239;
                    double zcorr2_div_zevap_denom;
                    double __tmp240;
                    double zlambda_pow_yrecldp_rcl_const4r;
                    double __tmp241;
                    double zrho_index_15;
                    double zrho_slice_times_zfallcorr;
                    double zrho_slice_zfallcorr_pow_0_5;
                    double yrecldp_rcl_const2r_times_zrho_slice_zfallcorr_0_5;
                    double zcorr2_pow_0_5;
                    double zlambda_pow_yrecldp_rcl_const3r;
                    double zcorr2_0_5_times_zlambda_yrecldp_rcl_const3r;
                    double yrecldp_rcl_const2r_zrho_slice_zfallcorr_0_5_div_zcorr2_0_5_zlambda_yrecldp_rcl_const3r;
                    double __tmp242;
                    double zbeta_0;
                    double zbeta_times_ptsphy_0;
                    double zdenom_0;
                    double zcovpclr_index_4;
                    double zcovpclr_slice_times_zbeta_0;
                    double zcovpclr_slice_zbeta_times_ptsphy;
                    double zcovpclr_slice_zbeta_ptsphy_times_zsubsat;
                    double zdpevap_0;
                    double zqxfg_index_38;
                    double zevap_1;
                    double zsolqa_index_48;
                    double zsolqa_slice_plus_zevap_1;
                    double zsolqa_index_49;
                    double zsolqa_slice_minus_zevap_1;
                    double zcovptot_index_13;
                    double zcovptot_index_14;
                    double za_index_45;
                    double zcovptot_slice_minus_za_slice_1;
                    double zcovptot_slice_za_slice_times_zevap_0;
                    double zqxfg_index_39;
                    double zcovptot_slice_za_slice_zevap_div_zqxfg_slice_0;
                    double max_0_0_zcovptot_slice_za_slice_zevap_zqxfg_slice_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (2 * klon)) - 1), &zqxfg_index_37, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (2 * klon)) - 1), &zqxfg_index_38, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (2 * klon)) - 1), &zqxfg_index_39, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_12, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_13, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_14, 1);
                    {
                        double __in1 = zqxfg_index_37;
                        double __in2 = zcovptot_index_12;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zpreclr_0 = __out;
                    }
                    {
                        double __inp = zpreclr_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_943_20)
                        __out = __inp;
                        ///////////////////

                        zpreclr = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zrho + (jl - 1), &zrho_index_13, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zrho + (jl - 1), &zrho_index_14, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zrho + (jl - 1), &zrho_index_15, 1);
                    {
                        double __in1 = zrho_index_14;
                        double __in2 = zpreclr;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zrho_slice_times_zpreclr = __out;
                    }
                    {
                        double __in1 = yrecldp_rdensref;
                        double __in2 = zrho_index_13;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        yrecldp_rdensref_div_zrho_slice_0 = __out;
                    }
                    {
                        double __in1 = yrecldp_rdensref_div_zrho_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, 0.4);
                        ///////////////////

                        zfallcorr_0 = __out;
                    }
                    {
                        double __inp = zfallcorr_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_944_20)
                        __out = __inp;
                        ///////////////////

                        zfallcorr = __out;
                    }
                    {
                        double __in1 = zrho_index_15;
                        double __in2 = zfallcorr;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zrho_slice_times_zfallcorr = __out;
                    }
                    {
                        double __in1 = zrho_slice_times_zfallcorr;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::sqrt(__in1);
                        ///////////////////

                        zrho_slice_zfallcorr_pow_0_5 = __out;
                    }
                    {
                        double __in1 = ydcst_rv;
                        double __in2 = ydcst_rd;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        ydcst_rv_div_ydcst_rd = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_95, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_96, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_97, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_98, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_99, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_100, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    ztp1 + ((jl + (klon * (jk - 1))) - 1), &ztp1_index_101, 1);
                    {
                        double __in1 = ztp1_index_98;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = (dace::math::ipow(__in1, 3));
                        ///////////////////

                        ztp1_slice_pow_3_0_0 = __out;
                    }
                    {
                        double __in1 = ztp1_index_99;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / 273.0);
                        ///////////////////

                        ztp1_slice_div_273_0_0 = __out;
                    }
                    {
                        double __in1 = ztp1_slice_div_273_0_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, 1.5);
                        ///////////////////

                        ztp1_slice_273_0_pow_1_5_0 = __out;
                    }
                    {
                        double __in1 = ztp1_slice_273_0_pow_1_5_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * 393.0);
                        ///////////////////

                        ztp1_slice_273_0_1_5_times_393_0 = __out;
                    }
                    {
                        double __in1 = ztp1_index_100;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + 120.0);
                        ///////////////////

                        ztp1_slice_plus_120_0_0 = __out;
                    }
                    {
                        double __in1 = ztp1_slice_273_0_1_5_times_393_0;
                        double __in2 = ztp1_slice_plus_120_0_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zcorr2 = __out;
                    }
                    {
                        double __in1 = zcorr2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::sqrt(__in1);
                        ///////////////////

                        zcorr2_pow_0_5 = __out;
                    }
                    {
                        double __in1 = ztp1_index_101;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = (dace::math::ipow(__in1, 2));
                        ///////////////////

                        ztp1_slice_pow_2_0_0 = __out;
                    }
                    {
                        double __in1 = ztp1_index_95;
                        double __in2 = ydcst_rtt;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        ztp1_slice_minus_ydcst_rtt_16 = __out;
                    }
                    {
                        double __in1 = ydthf_r3les;
                        double __in2 = ztp1_slice_minus_ydcst_rtt_16;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydthf_r3les_times_ztp1_slice_ydcst_rtt_5 = __out;
                    }
                    {
                        double __in1 = ztp1_index_96;
                        double __in2 = ydthf_r4les;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        ztp1_slice_minus_ydthf_r4les_9 = __out;
                    }
                    {
                        double __in1 = ydthf_r3les_times_ztp1_slice_ydcst_rtt_5;
                        double __in2 = ztp1_slice_minus_ydthf_r4les_9;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_5 = __out;
                    }
                    {
                        double __in1 = ydthf_r3les_ztp1_slice_ydcst_rtt_div_ztp1_slice_ydthf_r4les_5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_numpy_exp_)
                        __out = exp(__in1);
                        ///////////////////

                        exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_5 = __out;
                    }
                    {
                        double __in1 = ydthf_r2es;
                        double __in2 = exp_ydthf_r3les_ztp1_slice_ydcst_rtt_ztp1_slice_ydthf_r4les_5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydthf_r2es_times_expr_6 = __out;
                    }
                    {
                        double __in1 = ydcst_rv_div_ydcst_rd;
                        double __in2 = ydthf_r2es_times_expr_6;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zesatliq = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_fac1;
                        double __in2 = zrho_slice_times_zpreclr;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        yrecldp_rcl_fac1_div_zrho_slice_zpreclr = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_fac1_div_zrho_slice_zpreclr;
                        double __in2 = yrecldp_rcl_fac2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, __in2);
                        ///////////////////

                        zlambda_0 = __out;
                    }
                    {
                        double __inp = zlambda_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_946_20)
                        __out = __inp;
                        ///////////////////

                        zlambda = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_cdenom1;
                        double __in2 = zesatliq;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        yrecldp_rcl_cdenom1_times_zesatliq = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_cdenom2;
                        double __in2 = ztp1_index_97;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        yrecldp_rcl_cdenom2_times_ztp1_slice = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_cdenom2_times_ztp1_slice;
                        double __in2 = zesatliq;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        yrecldp_rcl_cdenom2_ztp1_slice_times_zesatliq = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_cdenom1_times_zesatliq;
                        double __in2 = yrecldp_rcl_cdenom2_ztp1_slice_times_zesatliq;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        yrecldp_rcl_cdenom1_zesatliq_minus_yrecldp_rcl_cdenom2_ztp1_slice_zesatliq = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_cdenom3;
                        double __in2 = ztp1_slice_pow_3_0_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        yrecldp_rcl_cdenom3_times_ztp1_slice_3_0 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_15, 1);
                    {
                        double __in1 = yrecldp_rcl_cdenom3_times_ztp1_slice_3_0;
                        double __in2 = pap_index_15;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        yrecldp_rcl_cdenom3_ztp1_slice_3_0_times_pap_slice = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_cdenom1_zesatliq_minus_yrecldp_rcl_cdenom2_ztp1_slice_zesatliq;
                        double __in2 = yrecldp_rcl_cdenom3_ztp1_slice_3_0_times_pap_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zevap_denom = __out;
                    }
                    {
                        double __in1 = zcorr2;
                        double __in2 = zevap_denom;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zcorr2_div_zevap_denom = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsliq + ((jl + (klon * (jk - 1))) - 1), &zqsliq_index_9, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsliq + ((jl + (klon * (jk - 1))) - 1), &zqsliq_index_10, 1);
                    {
                        double __in2 = zqsliq_index_10;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (0.5 / __in2);
                        ///////////////////

                        __tmp236 = __out;
                    }
                    {
                        double __in1 = __tmp236;
                        double __in2 = ztp1_slice_pow_2_0_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp237 = __out;
                    }
                    {
                        double __in1 = __tmp237;
                        double __in2 = zesatliq;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp238 = __out;
                    }
                    {
                        double __in1 = zzrh;
                        double __in2 = zqsliq_index_9;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zzrh_times_zqsliq_slice_1 = __out;
                    }
                    {
                        double __in1 = zzrh_times_zqsliq_slice_1;
                        double __in2 = zqe;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zzrh_zqsliq_slice_minus_zqe = __out;
                    }
                    {
                        double __in_a = zzrh_zqsliq_slice_minus_zqe;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 0.0);
                        ///////////////////

                        zsubsat_0 = __out;
                    }
                    {
                        double __inp = zsubsat_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_950_20)
                        __out = __inp;
                        ///////////////////

                        zsubsat = __out;
                    }
                    {
                        double __in1 = __tmp238;
                        double __in2 = yrecldp_rcl_const1r;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp239 = __out;
                    }
                    {
                        double __in1 = __tmp239;
                        double __in2 = zcorr2_div_zevap_denom;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        __tmp240 = __out;
                    }
                    {
                        double __in1 = zlambda;
                        double __in2 = yrecldp_rcl_const4r;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, __in2);
                        ///////////////////

                        zlambda_pow_yrecldp_rcl_const4r = __out;
                    }
                    {
                        double __in2 = zlambda_pow_yrecldp_rcl_const4r;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (0.78 / __in2);
                        ///////////////////

                        __tmp241 = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_const2r;
                        double __in2 = zrho_slice_zfallcorr_pow_0_5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        yrecldp_rcl_const2r_times_zrho_slice_zfallcorr_0_5 = __out;
                    }
                    {
                        double __in1 = zlambda;
                        double __in2 = yrecldp_rcl_const3r;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, __in2);
                        ///////////////////

                        zlambda_pow_yrecldp_rcl_const3r = __out;
                    }
                    {
                        double __in1 = zcorr2_pow_0_5;
                        double __in2 = zlambda_pow_yrecldp_rcl_const3r;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcorr2_0_5_times_zlambda_yrecldp_rcl_const3r = __out;
                    }
                    {
                        double __in1 = yrecldp_rcl_const2r_times_zrho_slice_zfallcorr_0_5;
                        double __in2 = zcorr2_0_5_times_zlambda_yrecldp_rcl_const3r;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        yrecldp_rcl_const2r_zrho_slice_zfallcorr_0_5_div_zcorr2_0_5_zlambda_yrecldp_rcl_const3r = __out;
                    }
                    {
                        double __in1 = __tmp241;
                        double __in2 = yrecldp_rcl_const2r_zrho_slice_zfallcorr_0_5_div_zcorr2_0_5_zlambda_yrecldp_rcl_const3r;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        __tmp242 = __out;
                    }
                    {
                        double __in1 = __tmp240;
                        double __in2 = __tmp242;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zbeta_0 = __out;
                    }
                    {
                        double __inp = zbeta_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_951_20)
                        __out = __inp;
                        ///////////////////

                        zbeta = __out;
                    }
                    {
                        double __in1 = zbeta;
                        double __in2 = ptsphy;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zbeta_times_ptsphy_0 = __out;
                    }
                    {
                        double __in2 = zbeta_times_ptsphy_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (1.0 + __in2);
                        ///////////////////

                        zdenom_0 = __out;
                    }
                    {
                        double __inp = zdenom_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_952_20)
                        __out = __inp;
                        ///////////////////

                        zdenom = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovpclr + (jl - 1), &zcovpclr_index_4, 1);
                    {
                        double __in1 = zcovpclr_index_4;
                        double __in2 = zbeta;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovpclr_slice_times_zbeta_0 = __out;
                    }
                    {
                        double __in1 = zcovpclr_slice_times_zbeta_0;
                        double __in2 = ptsphy;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovpclr_slice_zbeta_times_ptsphy = __out;
                    }
                    {
                        double __in1 = zcovpclr_slice_zbeta_times_ptsphy;
                        double __in2 = zsubsat;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovpclr_slice_zbeta_ptsphy_times_zsubsat = __out;
                    }
                    {
                        double __in1 = zcovpclr_slice_zbeta_ptsphy_times_zsubsat;
                        double __in2 = zdenom;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zdpevap_0 = __out;
                    }
                    {
                        double __inp = zdpevap_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_953_20)
                        __out = __inp;
                        ///////////////////

                        zdpevap = __out;
                    }
                    {
                        double __in_a = zdpevap;
                        double __in_b = zqxfg_index_38;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        zevap_1 = __out;
                    }
                    {
                        double __inp = zevap_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_954_20)
                        __out = __inp;
                        ///////////////////

                        zevap = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (14 * klon)) - 1), &zsolqa_index_48, 1);
                    {
                        double __in1 = zsolqa_index_48;
                        double __in2 = zevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqa_slice_plus_zevap_1 = __out;
                    }
                    {
                        double __inp = zsolqa_slice_plus_zevap_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_955_20)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (14 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (22 * klon)) - 1), &zsolqa_index_49, 1);
                    {
                        double __in1 = zsolqa_index_49;
                        double __in2 = zevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolqa_slice_minus_zevap_1 = __out;
                    }
                    {
                        double __inp = zsolqa_slice_minus_zevap_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_956_20)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (22 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_45, 1);
                    {
                        double __in1 = zcovptot_index_14;
                        double __in2 = za_index_45;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zcovptot_slice_minus_za_slice_1 = __out;
                    }
                    {
                        double __in1 = zcovptot_slice_minus_za_slice_1;
                        double __in2 = zevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovptot_slice_za_slice_times_zevap_0 = __out;
                    }
                    {
                        double __in1 = zcovptot_slice_za_slice_times_zevap_0;
                        double __in2 = zqxfg_index_39;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zcovptot_slice_za_slice_zevap_div_zqxfg_slice_0 = __out;
                    }
                    {
                        double __in_b = zcovptot_slice_za_slice_zevap_div_zqxfg_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(0.0, __in_b);
                        ///////////////////

                        max_0_0_zcovptot_slice_za_slice_zevap_zqxfg_slice_0 = __out;
                    }
                    {
                        double __in1 = zcovptot_index_13;
                        double __in2 = max_0_0_zcovptot_slice_za_slice_zevap_zqxfg_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zcovptot_slice_minus_expr_0 = __out;
                    }

                }
                {
                    double max_yrecldp_rcovpmin_zcovptot_slice_expr_0;
                    double zqxfg_index_40;
                    double zqxfg_slice_minus_zevap_0;

                    {
                        double __in_a = yrecldp_rcovpmin;
                        double __in_b = zcovptot_slice_minus_expr_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        max_yrecldp_rcovpmin_zcovptot_slice_expr_0 = __out;
                    }
                    {
                        double __inp = max_yrecldp_rcovpmin_zcovptot_slice_expr_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_957_20)
                        __out = __inp;
                        ///////////////////

                        zcovptot[(jl - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (2 * klon)) - 1), &zqxfg_index_40, 1);
                    {
                        double __in1 = zqxfg_index_40;
                        double __in2 = zevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqxfg_slice_minus_zevap_0 = __out;
                    }
                    {
                        double __inp = zqxfg_slice_minus_zevap_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_958_20)
                        __out = __inp;
                        ///////////////////

                        zqxfg[((jl + (2 * klon)) - 1)] = __out;
                    }

                }
            }

        }
        kfdia_plus_1_51 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_51); jl = (jl + 1)) {

            zcovpclr_index_5 = zcovpclr[(jl - 1)];
            zqxfg_index_41 = zqxfg[((jl + (3 * klon)) - 1)];
            {
                double __tmp243;
                double zcovpmax_index_2;
                double __tmp244;
                double za_index_46;
                double __tmp245;
                double max_zepsec_1_0_za_slice_5;
                double __tmp246;
                double zzrh_1;
                double max_zzrh_yrecldp_rprecrhmax_1;
                double min_expr_1_0_1;
                double zqx_index_45;
                double za_index_47;
                double zqsice_index_12;
                double za_slice_times_zqsice_slice_2;
                double zqx_slice_minus_za_slice_zqsice_slice_2;
                double za_index_48;
                double __tmp247;
                double max_zepsec_1_0_za_slice_6;
                double zqe_5;
                double zqsice_index_13;
                double min_zqe_zqsice_slice;
                double max_0_0_expr_1;
                double zqsice_index_14;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                za + ((jl + (klon * (jk - 1))) - 1), &za_index_47, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                za + ((jl + (klon * (jk - 1))) - 1), &za_index_48, 1);
                {
                    double __in2 = za_index_48;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp247 = __out;
                }
                {
                    double __in_b = __tmp247;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(1e-14, __in_b);
                    ///////////////////

                    max_zepsec_1_0_za_slice_6 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_45, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_12, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_13, 1);
                {
                    double __in1 = za_index_47;
                    double __in2 = zqsice_index_12;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    za_slice_times_zqsice_slice_2 = __out;
                }
                {
                    double __in1 = zqx_index_45;
                    double __in2 = za_slice_times_zqsice_slice_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqx_slice_minus_za_slice_zqsice_slice_2 = __out;
                }
                {
                    double __in1 = zqx_slice_minus_za_slice_zqsice_slice_2;
                    double __in2 = max_zepsec_1_0_za_slice_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zqe_5 = __out;
                }
                {
                    double __inp = zqe_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_963_16)
                    __out = __inp;
                    ///////////////////

                    zqe = __out;
                }
                {
                    double __in_a = zqe;
                    double __in_b = zqsice_index_13;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_zqe_zqsice_slice = __out;
                }
                {
                    double __in_b = min_zqe_zqsice_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(0.0, __in_b);
                    ///////////////////

                    max_0_0_expr_1 = __out;
                }
                {
                    double __inp = max_0_0_expr_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_964_16)
                    __out = __inp;
                    ///////////////////

                    zqe = __out;
                }
                {
                    double __in2 = yrecldp_rprecrhmax;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp243 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zcovpmax + (jl - 1), &zcovpmax_index_2, 1);
                {
                    double __in1 = __tmp243;
                    double __in2 = zcovpmax_index_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp244 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                za + ((jl + (klon * (jk - 1))) - 1), &za_index_46, 1);
                {
                    double __in2 = za_index_46;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp245 = __out;
                }
                {
                    double __in_b = __tmp245;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(1e-14, __in_b);
                    ///////////////////

                    max_zepsec_1_0_za_slice_5 = __out;
                }
                {
                    double __in1 = __tmp244;
                    double __in2 = max_zepsec_1_0_za_slice_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    __tmp246 = __out;
                }
                {
                    double __in1 = yrecldp_rprecrhmax;
                    double __in2 = __tmp246;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zzrh_1 = __out;
                }
                {
                    double __inp = zzrh_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_961_16)
                    __out = __inp;
                    ///////////////////

                    zzrh = __out;
                }
                {
                    double __in_a = zzrh;
                    double __in_b = yrecldp_rprecrhmax;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_zzrh_yrecldp_rprecrhmax_1 = __out;
                }
                {
                    double __in_a = max_zzrh_yrecldp_rprecrhmax_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, 1);
                    ///////////////////

                    min_expr_1_0_1 = __out;
                }
                {
                    double __inp = min_expr_1_0_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_962_16)
                    __out = __inp;
                    ///////////////////

                    zzrh = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_14, 1);
                {
                    double __in1 = zzrh;
                    double __in2 = zqsice_index_14;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zzrh_times_zqsice_slice = __out;
                }

            }
            __tmp250 = (zqe < zzrh_times_zqsice_slice);
            if ((((zcovpclr_index_5 > 1e-14) && (zqxfg_index_41 > 1e-14)) && __tmp250)) {
                {
                    double zbeta1;
                    double zdpr;
                    double zqxfg_index_42;
                    double zcovpclr_index_6;
                    double zqxfg_slice_times_zcovpclr_slice_0;
                    double zcovptot_index_15;
                    double zdtgdp_index_8;
                    double zcovptot_slice_times_zdtgdp_slice_1;
                    double abs_zcovptot_slice_zdtgdp_slice_0;
                    double max_expr_zepsilon_0;
                    double zcovptot_index_16;
                    double zdtgdp_index_9;
                    double zcovptot_slice_times_zdtgdp_slice_2;
                    double sign_zcovptot_slice_zdtgdp_slice_0;
                    double expr_times_expr_3;
                    double zpreclr_1;
                    double pap_index_16;
                    double paph_index_5;
                    double pap_slice_div_paph_slice_0;
                    double sqrt_pap_slice_paph_slice_0;
                    double expr_div_yrecldp_rvrfactor_0;
                    double expr_yrecldp_rvrfactor_times_zpreclr_0;
                    double zcovpclr_index_7;
                    double max_zcovpclr_slice_zepsec_0;
                    double zbeta1_0;
                    double ydcst_rg_times_yrecldp_rpecons_0;
                    double zbeta1_pow_0_5777_0;
                    double zbeta_1;
                    double zbeta_times_ptsphy_1;
                    double zcorqsice_index_2;
                    double zbeta_ptsphy_times_zcorqsice_slice;
                    double zdenom_1;
                    double zcovpclr_index_8;
                    double zcovpclr_slice_times_zbeta_1;
                    double zqsice_index_15;
                    double zqsice_slice_minus_zqe;
                    double zcovpclr_slice_zbeta_times_zqsice_slice_zqe;
                    double zcovpclr_slice_zbeta_zqsice_slice_zqe_div_zdenom;
                    double zdp_index_5;
                    double zcovpclr_slice_zbeta_zqsice_slice_zqe_zdenom_times_zdp_slice;
                    double zdpr_0;
                    double zdtgdp_index_10;
                    double zdpevap_1;
                    double zqxfg_index_43;
                    double zevap_2;
                    double zsolqa_index_50;
                    double zsolqa_slice_plus_zevap_2;
                    double zsolqa_index_51;
                    double zsolqa_slice_minus_zevap_2;
                    double zcovptot_index_17;
                    double zcovptot_index_18;
                    double za_index_49;
                    double zcovptot_slice_minus_za_slice_2;
                    double zcovptot_slice_za_slice_times_zevap_1;
                    double zqxfg_index_44;
                    double zcovptot_slice_za_slice_zevap_div_zqxfg_slice_1;
                    double max_0_0_zcovptot_slice_za_slice_zevap_zqxfg_slice_1;
                    double zcovptot_slice_minus_expr_1;
                    double max_yrecldp_rcovpmin_zcovptot_slice_expr_1;
                    double zqxfg_index_45;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (3 * klon)) - 1), &zqxfg_index_42, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (3 * klon)) - 1), &zqxfg_index_43, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (3 * klon)) - 1), &zqxfg_index_44, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxfg + ((jl + (3 * klon)) - 1), &zqxfg_index_45, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovpclr + (jl - 1), &zcovpclr_index_6, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovpclr + (jl - 1), &zcovpclr_index_7, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovpclr + (jl - 1), &zcovpclr_index_8, 1);
                    {
                        double __in1 = zqxfg_index_42;
                        double __in2 = zcovpclr_index_6;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zqxfg_slice_times_zcovpclr_slice_0 = __out;
                    }
                    {
                        double __in_a = zcovpclr_index_7;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 1e-14);
                        ///////////////////

                        max_zcovpclr_slice_zepsec_0 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_15, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_16, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_17, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcovptot + (jl - 1), &zcovptot_index_18, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdtgdp + (jl - 1), &zdtgdp_index_8, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdtgdp + (jl - 1), &zdtgdp_index_9, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdtgdp + (jl - 1), &zdtgdp_index_10, 1);
                    {
                        double __in1 = zcovptot_index_15;
                        double __in2 = zdtgdp_index_8;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovptot_slice_times_zdtgdp_slice_1 = __out;
                    }
                    {
                        double __inp = zcovptot_slice_times_zdtgdp_slice_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (abs)
                        __out = abs(__inp);
                        ///////////////////

                        abs_zcovptot_slice_zdtgdp_slice_0 = __out;
                    }
                    {
                        double __in1 = zcovptot_index_16;
                        double __in2 = zdtgdp_index_9;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovptot_slice_times_zdtgdp_slice_2 = __out;
                    }
                    {
                        double __in1 = zcovptot_slice_times_zdtgdp_slice_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_numpy_sign_)
                        __out = sign_numpy_2(__in1);
                        ///////////////////

                        sign_zcovptot_slice_zdtgdp_slice_0 = __out;
                    }
                    {
                        double __in_a = abs_zcovptot_slice_zdtgdp_slice_0;
                        double __in_b = zepsilon;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        max_expr_zepsilon_0 = __out;
                    }
                    {
                        double __in1 = max_expr_zepsilon_0;
                        double __in2 = sign_zcovptot_slice_zdtgdp_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        expr_times_expr_3 = __out;
                    }
                    {
                        double __in1 = zqxfg_slice_times_zcovpclr_slice_0;
                        double __in2 = expr_times_expr_3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zpreclr_1 = __out;
                    }
                    {
                        double __inp = zpreclr_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_967_20)
                        __out = __inp;
                        ///////////////////

                        zpreclr = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    pap + ((jl + (klon * (jk - 1))) - 1), &pap_index_16, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    paph + ((jl + (klev * klon)) - 1), &paph_index_5, 1);
                    {
                        double __in1 = pap_index_16;
                        double __in2 = paph_index_5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        pap_slice_div_paph_slice_0 = __out;
                    }
                    {
                        double __in1 = pap_slice_div_paph_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_numpy_sqrt_)
                        __out = sqrt(__in1);
                        ///////////////////

                        sqrt_pap_slice_paph_slice_0 = __out;
                    }
                    {
                        double __in1 = sqrt_pap_slice_paph_slice_0;
                        double __in2 = yrecldp_rvrfactor;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        expr_div_yrecldp_rvrfactor_0 = __out;
                    }
                    {
                        double __in1 = expr_div_yrecldp_rvrfactor_0;
                        double __in2 = zpreclr;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        expr_yrecldp_rvrfactor_times_zpreclr_0 = __out;
                    }
                    {
                        double __in1 = expr_yrecldp_rvrfactor_times_zpreclr_0;
                        double __in2 = max_zcovpclr_slice_zepsec_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zbeta1_0 = __out;
                    }
                    {
                        double __inp = zbeta1_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_968_20)
                        __out = __inp;
                        ///////////////////

                        zbeta1 = __out;
                    }
                    {
                        double __in1 = zbeta1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Pow_)
                        __out = dace::math::pow(__in1, 0.5777);
                        ///////////////////

                        zbeta1_pow_0_5777_0 = __out;
                    }
                    {
                        double __in1 = ydcst_rg;
                        double __in2 = yrecldp_rpecons;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        ydcst_rg_times_yrecldp_rpecons_0 = __out;
                    }
                    {
                        double __in1 = ydcst_rg_times_yrecldp_rpecons_0;
                        double __in2 = zbeta1_pow_0_5777_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zbeta_1 = __out;
                    }
                    {
                        double __inp = zbeta_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_969_20)
                        __out = __inp;
                        ///////////////////

                        zbeta = __out;
                    }
                    {
                        double __in1 = zcovpclr_index_8;
                        double __in2 = zbeta;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovpclr_slice_times_zbeta_1 = __out;
                    }
                    {
                        double __in1 = zbeta;
                        double __in2 = ptsphy;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zbeta_times_ptsphy_1 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zcorqsice + (jl - 1), &zcorqsice_index_2, 1);
                    {
                        double __in1 = zbeta_times_ptsphy_1;
                        double __in2 = zcorqsice_index_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zbeta_ptsphy_times_zcorqsice_slice = __out;
                    }
                    {
                        double __in2 = zbeta_ptsphy_times_zcorqsice_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (1.0 + __in2);
                        ///////////////////

                        zdenom_1 = __out;
                    }
                    {
                        double __inp = zdenom_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_970_20)
                        __out = __inp;
                        ///////////////////

                        zdenom = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqsice + ((jl + (klon * (jk - 1))) - 1), &zqsice_index_15, 1);
                    {
                        double __in1 = zqsice_index_15;
                        double __in2 = zqe;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqsice_slice_minus_zqe = __out;
                    }
                    {
                        double __in1 = zcovpclr_slice_times_zbeta_1;
                        double __in2 = zqsice_slice_minus_zqe;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovpclr_slice_zbeta_times_zqsice_slice_zqe = __out;
                    }
                    {
                        double __in1 = zcovpclr_slice_zbeta_times_zqsice_slice_zqe;
                        double __in2 = zdenom;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zcovpclr_slice_zbeta_zqsice_slice_zqe_div_zdenom = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zdp + (jl - 1), &zdp_index_5, 1);
                    {
                        double __in1 = zcovpclr_slice_zbeta_zqsice_slice_zqe_div_zdenom;
                        double __in2 = zdp_index_5;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovpclr_slice_zbeta_zqsice_slice_zqe_zdenom_times_zdp_slice = __out;
                    }
                    {
                        double __in1 = zcovpclr_slice_zbeta_zqsice_slice_zqe_zdenom_times_zdp_slice;
                        double __in2 = zrg_r;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zdpr_0 = __out;
                    }
                    {
                        double __inp = zdpr_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_971_20)
                        __out = __inp;
                        ///////////////////

                        zdpr = __out;
                    }
                    {
                        double __in1 = zdpr;
                        double __in2 = zdtgdp_index_10;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zdpevap_1 = __out;
                    }
                    {
                        double __inp = zdpevap_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_972_20)
                        __out = __inp;
                        ///////////////////

                        zdpevap = __out;
                    }
                    {
                        double __in_a = zdpevap;
                        double __in_b = zqxfg_index_43;
                        double __out;

                        ///////////////////
                        // Tasklet code (__min2)
                        __out = min(__in_a, __in_b);
                        ///////////////////

                        zevap_2 = __out;
                    }
                    {
                        double __inp = zevap_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_973_20)
                        __out = __inp;
                        ///////////////////

                        zevap = __out;
                    }
                    {
                        double __in1 = zqxfg_index_45;
                        double __in2 = zevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqxfg_slice_minus_zevap_1 = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (19 * klon)) - 1), &zsolqa_index_50, 1);
                    {
                        double __in1 = zsolqa_index_50;
                        double __in2 = zevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zsolqa_slice_plus_zevap_2 = __out;
                    }
                    {
                        double __inp = zsolqa_slice_plus_zevap_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_974_20)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (19 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsolqa + ((jl + (23 * klon)) - 1), &zsolqa_index_51, 1);
                    {
                        double __in1 = zsolqa_index_51;
                        double __in2 = zevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsolqa_slice_minus_zevap_2 = __out;
                    }
                    {
                        double __inp = zsolqa_slice_minus_zevap_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_975_20)
                        __out = __inp;
                        ///////////////////

                        zsolqa[((jl + (23 * klon)) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    za + ((jl + (klon * (jk - 1))) - 1), &za_index_49, 1);
                    {
                        double __in1 = zcovptot_index_18;
                        double __in2 = za_index_49;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zcovptot_slice_minus_za_slice_2 = __out;
                    }
                    {
                        double __in1 = zcovptot_slice_minus_za_slice_2;
                        double __in2 = zevap;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zcovptot_slice_za_slice_times_zevap_1 = __out;
                    }
                    {
                        double __in1 = zcovptot_slice_za_slice_times_zevap_1;
                        double __in2 = zqxfg_index_44;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zcovptot_slice_za_slice_zevap_div_zqxfg_slice_1 = __out;
                    }
                    {
                        double __in_b = zcovptot_slice_za_slice_zevap_div_zqxfg_slice_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(0.0, __in_b);
                        ///////////////////

                        max_0_0_zcovptot_slice_za_slice_zevap_zqxfg_slice_1 = __out;
                    }
                    {
                        double __in1 = zcovptot_index_17;
                        double __in2 = max_0_0_zcovptot_slice_za_slice_zevap_zqxfg_slice_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zcovptot_slice_minus_expr_1 = __out;
                    }
                    {
                        double __in_a = yrecldp_rcovpmin;
                        double __in_b = zcovptot_slice_minus_expr_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        max_yrecldp_rcovpmin_zcovptot_slice_expr_1 = __out;
                    }
                    {
                        double __inp = max_yrecldp_rcovpmin_zcovptot_slice_expr_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_976_20)
                        __out = __inp;
                        ///////////////////

                        zcovptot[(jl - 1)] = __out;
                    }

                }
                {

                    {
                        double __inp = zqxfg_slice_minus_zevap_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_977_20)
                        __out = __inp;
                        ///////////////////

                        zqxfg[((jl + (3 * klon)) - 1)] = __out;
                    }

                }
            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            llfall_index_2 = llfall[(jm - 1)];
            if (llfall_index_2) {

                kfdia_plus_1_53 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_53); jl = (jl + 1)) {

                    zqxfg_index_47 = zqxfg[((jl + (klon * (jm - 1))) - 1)];

                    __tmp263 = (zqxfg_index_47 < yrecldp_rlmin);
                    if (__tmp263) {
                        {
                            double zsolqa_index_54;
                            double zqxfg_index_48;
                            double zsolqa_slice_plus_zqxfg_slice;
                            double zsolqa_index_55;
                            double zqxfg_index_49;
                            double zsolqa_slice_minus_zqxfg_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + (((jl + ((5 * klon) * (jm - 1))) + (4 * klon)) - 1), &zsolqa_index_54, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + (klon * (jm - 1))) - 1), &zqxfg_index_48, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqxfg + ((jl + (klon * (jm - 1))) - 1), &zqxfg_index_49, 1);
                            {
                                double __in1 = zsolqa_index_54;
                                double __in2 = zqxfg_index_48;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (__in1 + __in2);
                                ///////////////////

                                zsolqa_slice_plus_zqxfg_slice = __out;
                            }
                            {
                                double __inp = zsolqa_slice_plus_zqxfg_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_1007_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[(((jl + ((5 * klon) * (jm - 1))) + (4 * klon)) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + (((jl + (klon * (jm - 1))) + (20 * klon)) - 1), &zsolqa_index_55, 1);
                            {
                                double __in1 = zsolqa_index_55;
                                double __in2 = zqxfg_index_49;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zsolqa_slice_minus_zqxfg_slice = __out;
                            }
                            {
                                double __inp = zsolqa_slice_minus_zqxfg_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_1008_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[(((jl + (klon * (jm - 1))) + (20 * klon)) - 1)] = __out;
                            }

                        }
                    }

                }
            }

        }
        kfdia_plus_1_54 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_54); jl = (jl + 1)) {
            {
                double za_index_54;
                double zsolac_index_3;
                double zsolab_index_0;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                za + ((jl + (klon * (jk - 1))) - 1), &za_index_54, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zsolac + (jl - 1), &zsolac_index_3, 1);
                {
                    double __in1 = za_index_54;
                    double __in2 = zsolac_index_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    za_slice_plus_zsolac_slice = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zsolab + (jl - 1), &zsolab_index_0, 1);
                {
                    double __in2 = zsolab_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (1.0 + __in2);
                    ///////////////////

                    __tmp264 = __out;
                }

            }
            zanew = (za_slice_plus_zsolac_slice / __tmp264);

            zanew = min(zanew, 1);
            if ((zanew < yrecldp_ramin)) {

                zanew = 0.0;

            }
            {
                double zaorig_index;
                double zda_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zaorig + ((jl + (klon * (jk - 1))) - 1), &zaorig_index, 1);
                {
                    double __in2 = zaorig_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (zanew - __in2);
                    ///////////////////

                    zda_slice = __out;
                }
                {
                    double __inp = zda_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1014_12)
                    __out = __inp;
                    ///////////////////

                    zda[(jl - 1)] = __out;
                }
                {
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1015_12)
                    __out = zanew;
                    ///////////////////

                    zanewm1[(jl - 1)] = __out;
                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            kfdia_plus_1_56 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_56); jl = (jl + 1)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1021_16)
                        __out = 0.0;
                        ///////////////////

                        zsinksum[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {
            for (jn = 1; (jn < (5 + 1)); jn = (jn + 1)) {

                kfdia_plus_1_57 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_57); jl = (jl + 1)) {
                    {
                        double zsinksum_index;
                        double zsolqa_index_56;
                        double zsinksum_slice_minus_zsolqa_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsinksum + ((jl + (klon * (jm - 1))) - 1), &zsinksum_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + (((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1), &zsolqa_index_56, 1);
                        {
                            double __in1 = zsinksum_index;
                            double __in2 = zsolqa_index_56;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zsinksum_slice_minus_zsolqa_slice = __out;
                        }
                        {
                            double __inp = zsinksum_slice_minus_zsolqa_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1025_20)
                            __out = __inp;
                            ///////////////////

                            zsinksum[((jl + (klon * (jm - 1))) - 1)] = __out;
                        }

                    }

                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            kfdia_plus_1_58 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_58); jl = (jl + 1)) {
                {
                    double zqx_index_51;
                    double zmax;
                    double zsinksum_index_0;
                    double zrat;
                    double zratio_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx_index_51, 1);
                    {
                        double __in_a = zqx_index_51;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 1e-14);
                        ///////////////////

                        zmax = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsinksum + ((jl + (klon * (jm - 1))) - 1), &zsinksum_index_0, 1);
                    {
                        double __in_a = zsinksum_index_0;
                        double __in_b = zmax;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        zrat = __out;
                    }
                    {
                        double __in1 = zmax;
                        double __in2 = zrat;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zratio_slice = __out;
                    }
                    {
                        double __inp = zratio_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1030_16)
                        __out = __inp;
                        ///////////////////

                        zratio[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            kfdia_plus_1_59 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_59); jl = (jl + 1)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1033_16)
                        __out = 0.0;
                        ///////////////////

                        zsinksum[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {
            {

                {
                    #pragma omp parallel for
                    for (auto __i0 = 0; __i0 < klon; __i0 += 1) {
                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1035_12)
                            __out = 0.0;
                            ///////////////////

                            psum_solqa[__i0] = __out;
                        }
                    }
                }

            }
            for (jn = 1; (jn < (5 + 1)); jn = (jn + 1)) {

                kfdia_plus_1_60 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_60); jl = (jl + 1)) {
                    {
                        double psum_solqa_index;
                        double zsolqa_index_57;
                        double psum_solqa_slice_plus_zsolqa_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        psum_solqa + (jl - 1), &psum_solqa_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + (((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1), &zsolqa_index_57, 1);
                        {
                            double __in1 = psum_solqa_index;
                            double __in2 = zsolqa_index_57;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            psum_solqa_slice_plus_zsolqa_slice = __out;
                        }
                        {
                            double __inp = psum_solqa_slice_plus_zsolqa_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1038_20)
                            __out = __inp;
                            ///////////////////

                            psum_solqa[(jl - 1)] = __out;
                        }

                    }

                }

            }
            kfdia_plus_1_61 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_61); jl = (jl + 1)) {
                {
                    double zsinksum_index_1;
                    double psum_solqa_index_0;
                    double zsinksum_slice_minus_psum_solqa_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsinksum + ((jl + (klon * (jm - 1))) - 1), &zsinksum_index_1, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    psum_solqa + (jl - 1), &psum_solqa_index_0, 1);
                    {
                        double __in1 = zsinksum_index_1;
                        double __in2 = psum_solqa_index_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zsinksum_slice_minus_psum_solqa_slice = __out;
                    }
                    {
                        double __inp = zsinksum_slice_minus_psum_solqa_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1040_16)
                        __out = __inp;
                        ///////////////////

                        zsinksum[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }
            kfdia_plus_1_62 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_62); jl = (jl + 1)) {
                {
                    double zqx_index_52;
                    double zmm;
                    double zsinksum_index_2;
                    double zrr;
                    double zratio_slice_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx_index_52, 1);
                    {
                        double __in_a = zqx_index_52;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, 1e-14);
                        ///////////////////

                        zmm = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zsinksum + ((jl + (klon * (jm - 1))) - 1), &zsinksum_index_2, 1);
                    {
                        double __in_a = zsinksum_index_2;
                        double __in_b = zmm;
                        double __out;

                        ///////////////////
                        // Tasklet code (__max2)
                        __out = max(__in_a, __in_b);
                        ///////////////////

                        zrr = __out;
                    }
                    {
                        double __in1 = zmm;
                        double __in2 = zrr;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zratio_slice_0 = __out;
                    }
                    {
                        double __inp = zratio_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1044_16)
                        __out = __inp;
                        ///////////////////

                        zratio[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }
            kfdia_plus_1_63 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_63); jl = (jl + 1)) {
                {


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zratio + ((jl + (klon * (jm - 1))) - 1), &zzratio, 1);

                }
                for (jn = 1; (jn < (5 + 1)); jn = (jn + 1)) {

                    zsolqa_index_58 = zsolqa[(((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1)];
                    if ((zsolqa_index_58 < 0.0)) {
                        {
                            double zsolqa_index_59;
                            double zsolqa_slice_times_zzratio;
                            double zsolqa_index_60;
                            double zsolqa_slice_times_zzratio_0;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + (((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1), &zsolqa_index_59, 1);
                            {
                                double __in1 = zsolqa_index_59;
                                double __in2 = zzratio;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zsolqa_slice_times_zzratio = __out;
                            }
                            {
                                double __inp = zsolqa_slice_times_zzratio;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_1049_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[(((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1)] = __out;
                            }

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqa + (((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1), &zsolqa_index_60, 1);
                            {
                                double __in1 = zsolqa_index_60;
                                double __in2 = zzratio;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zsolqa_slice_times_zzratio_0 = __out;
                            }
                            {
                                double __inp = zsolqa_slice_times_zzratio_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_1050_24)
                                __out = __inp;
                                ///////////////////

                                zsolqa[(((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1)] = __out;
                            }

                        }
                    }

                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {
            for (jn = 1; (jn < (5 + 1)); jn = (jn + 1)) {
                if ((jn == jm)) {

                    kfdia_plus_1_64 = (kfdia + 1);
                    for (jl = kidia; (jl < kfdia_plus_1_64); jl = (jl + 1)) {
                        {
                            double zfallsink_index;
                            double zqlhs_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zfallsink + ((jl + (klon * (jm - 1))) - 1), &zfallsink_index, 1);
                            {
                                double __in2 = zfallsink_index;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Add_)
                                __out = (1.0 + __in2);
                                ///////////////////

                                zqlhs_slice = __out;
                            }
                            {
                                double __inp = zqlhs_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_1055_24)
                                __out = __inp;
                                ///////////////////

                                zqlhs[(((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1)] = __out;
                            }

                        }
                        for (jo = 1; (jo < (5 + 1)); jo = (jo + 1)) {
                            {
                                double zqlhs_index;
                                double zsolqb_index_5;
                                double zqlhs_slice_plus_zsolqb_slice;


                                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                                zqlhs + (((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1), &zqlhs_index, 1);

                                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                                zsolqb + (((jl + ((5 * klon) * (jn - 1))) + (klon * (jo - 1))) - 1), &zsolqb_index_5, 1);
                                {
                                    double __in1 = zqlhs_index;
                                    double __in2 = zsolqb_index_5;
                                    double __out;

                                    ///////////////////
                                    // Tasklet code (_Add_)
                                    __out = (__in1 + __in2);
                                    ///////////////////

                                    zqlhs_slice_plus_zsolqb_slice = __out;
                                }
                                {
                                    double __inp = zqlhs_slice_plus_zsolqb_slice;
                                    double __out;

                                    ///////////////////
                                    // Tasklet code (assign_1057_28)
                                    __out = __inp;
                                    ///////////////////

                                    zqlhs[(((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1)] = __out;
                                }

                            }

                        }

                    }
                } else {

                    kfdia_plus_1_65 = (kfdia + 1);
                    for (jl = kidia; (jl < kfdia_plus_1_65); jl = (jl + 1)) {
                        {
                            double zsolqb_index_6;
                            double zqlhs_slice_0;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zsolqb + (((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1), &zsolqb_index_6, 1);
                            {
                                double __in = zsolqb_index_6;
                                double __out;

                                ///////////////////
                                // Tasklet code (_USub_)
                                __out = (- __in);
                                ///////////////////

                                zqlhs_slice_0 = __out;
                            }
                            {
                                double __inp = zqlhs_slice_0;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_1060_24)
                                __out = __inp;
                                ///////////////////

                                zqlhs[(((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1)] = __out;
                            }

                        }

                    }
                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            kfdia_plus_1_66 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_66); jl = (jl + 1)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1063_16)
                        __out = 0.0;
                        ///////////////////

                        zexplicit = __out;
                    }

                }
                for (jn = 1; (jn < (5 + 1)); jn = (jn + 1)) {
                    {
                        double zsolqa_index_61;
                        double zexplicit_plus_zsolqa_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zsolqa + (((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1), &zsolqa_index_61, 1);
                        {
                            double __in1 = zexplicit;
                            double __in2 = zsolqa_index_61;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zexplicit_plus_zsolqa_slice = __out;
                        }
                        {
                            double __inp = zexplicit_plus_zsolqa_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1065_20)
                            __out = __inp;
                            ///////////////////

                            zexplicit = __out;
                        }

                    }

                }
                {
                    double zqx_index_53;
                    double zqxn_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx_index_53, 1);
                    {
                        double __in1 = zqx_index_53;
                        double __in2 = zexplicit;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zqxn_slice = __out;
                    }
                    {
                        double __inp = zqxn_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1066_16)
                        __out = __inp;
                        ///////////////////

                        zqxn[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }

        }
        for (jn = 1; (jn < 5); jn = (jn + 1)) {
            for (jm = (jn + 1); (jm < (5 + 1)); jm = (jm + 1)) {

                kfdia_plus_1_67 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_67); jl = (jl + 1)) {
                    {
                        double zqlhs_index_0;
                        double zqlhs_index_1;
                        double zqlhs_slice_div_zqlhs_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqlhs + (((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1), &zqlhs_index_0, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqlhs + ((jl + ((6 * klon) * (jn - 1))) - 1), &zqlhs_index_1, 1);
                        {
                            double __in1 = zqlhs_index_0;
                            double __in2 = zqlhs_index_1;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Div_)
                            __out = (__in1 / __in2);
                            ///////////////////

                            zqlhs_slice_div_zqlhs_slice = __out;
                        }
                        {
                            double __inp = zqlhs_slice_div_zqlhs_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1070_20)
                            __out = __inp;
                            ///////////////////

                            zqlhs[(((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1)] = __out;
                        }

                    }

                }

            }

        }
        for (jn = 1; (jn < 5); jn = (jn + 1)) {
            for (jm = (jn + 1); (jm < (5 + 1)); jm = (jm + 1)) {
                for (ik = (jn + 1); (ik < (5 + 1)); ik = (ik + 1)) {

                    kfdia_plus_1_68 = (kfdia + 1);
                    for (jl = kidia; (jl < kfdia_plus_1_68); jl = (jl + 1)) {
                        {
                            double zqlhs_index_2;
                            double zqlhs_index_3;
                            double zqlhs_index_4;
                            double zqlhs_slice_times_zqlhs_slice;
                            double zqlhs_slice_minus_zqlhs_slice_zqlhs_slice;


                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqlhs + (((jl + ((5 * klon) * (ik - 1))) + (klon * (jm - 1))) - 1), &zqlhs_index_2, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqlhs + (((jl + (klon * (jm - 1))) + ((5 * klon) * (jn - 1))) - 1), &zqlhs_index_3, 1);

                            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                            zqlhs + (((jl + ((5 * klon) * (ik - 1))) + (klon * (jn - 1))) - 1), &zqlhs_index_4, 1);
                            {
                                double __in1 = zqlhs_index_3;
                                double __in2 = zqlhs_index_4;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Mult_)
                                __out = (__in1 * __in2);
                                ///////////////////

                                zqlhs_slice_times_zqlhs_slice = __out;
                            }
                            {
                                double __in1 = zqlhs_index_2;
                                double __in2 = zqlhs_slice_times_zqlhs_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (_Sub_)
                                __out = (__in1 - __in2);
                                ///////////////////

                                zqlhs_slice_minus_zqlhs_slice_zqlhs_slice = __out;
                            }
                            {
                                double __inp = zqlhs_slice_minus_zqlhs_slice_zqlhs_slice;
                                double __out;

                                ///////////////////
                                // Tasklet code (assign_1075_24)
                                __out = __inp;
                                ///////////////////

                                zqlhs[(((jl + ((5 * klon) * (ik - 1))) + (klon * (jm - 1))) - 1)] = __out;
                            }

                        }

                    }

                }

            }

        }
        for (jn = 2; (jn < (5 + 1)); jn = (jn + 1)) {
            for (jm = 1; (jm < jn); jm = (jm + 1)) {

                kfdia_plus_1_69 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_69); jl = (jl + 1)) {
                    {
                        double zqxn_index;
                        double zqlhs_index_5;
                        double zqxn_index_0;
                        double zqlhs_slice_times_zqxn_slice;
                        double zqxn_slice_minus_zqlhs_slice_zqxn_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxn + ((jl + (klon * (jn - 1))) - 1), &zqxn_index, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxn_index_0, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqlhs + (((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1), &zqlhs_index_5, 1);
                        {
                            double __in1 = zqlhs_index_5;
                            double __in2 = zqxn_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zqlhs_slice_times_zqxn_slice = __out;
                        }
                        {
                            double __in1 = zqxn_index;
                            double __in2 = zqlhs_slice_times_zqxn_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqxn_slice_minus_zqlhs_slice_zqxn_slice = __out;
                        }
                        {
                            double __inp = zqxn_slice_minus_zqlhs_slice_zqxn_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1079_20)
                            __out = __inp;
                            ///////////////////

                            zqxn[((jl + (klon * (jn - 1))) - 1)] = __out;
                        }

                    }

                }

            }

        }
        kfdia_plus_1_70 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_70); jl = (jl + 1)) {
            {
                double zqxn_index_1;
                double zqlhs_index_6;
                double zqxn_slice_div_zqlhs_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxn + ((jl + (4 * klon)) - 1), &zqxn_index_1, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqlhs + ((jl + (24 * klon)) - 1), &zqlhs_index_6, 1);
                {
                    double __in1 = zqxn_index_1;
                    double __in2 = zqlhs_index_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    zqxn_slice_div_zqlhs_slice = __out;
                }
                {
                    double __inp = zqxn_slice_div_zqlhs_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1081_12)
                    __out = __inp;
                    ///////////////////

                    zqxn[((jl + (4 * klon)) - 1)] = __out;
                }

            }

        }
        for (jn = (5 - 1); (jn > 0); jn = (jn + -1)) {
            for (jm = (jn + 1); (jm < (5 + 1)); jm = (jm + 1)) {

                kfdia_plus_1_71 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_71); jl = (jl + 1)) {
                    {
                        double zqxn_index_2;
                        double zqlhs_index_7;
                        double zqxn_index_3;
                        double zqlhs_slice_times_zqxn_slice_0;
                        double zqxn_slice_minus_zqlhs_slice_zqxn_slice_0;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxn + ((jl + (klon * (jn - 1))) - 1), &zqxn_index_2, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxn_index_3, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqlhs + (((jl + ((5 * klon) * (jm - 1))) + (klon * (jn - 1))) - 1), &zqlhs_index_7, 1);
                        {
                            double __in1 = zqlhs_index_7;
                            double __in2 = zqxn_index_3;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            zqlhs_slice_times_zqxn_slice_0 = __out;
                        }
                        {
                            double __in1 = zqxn_index_2;
                            double __in2 = zqlhs_slice_times_zqxn_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqxn_slice_minus_zqlhs_slice_zqxn_slice_0 = __out;
                        }
                        {
                            double __inp = zqxn_slice_minus_zqlhs_slice_zqxn_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1085_20)
                            __out = __inp;
                            ///////////////////

                            zqxn[((jl + (klon * (jn - 1))) - 1)] = __out;
                        }

                    }

                }

            }
            kfdia_plus_1_72 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_72); jl = (jl + 1)) {
                {
                    double zqxn_index_4;
                    double zqlhs_index_8;
                    double zqxn_slice_div_zqlhs_slice_0;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxn + ((jl + (klon * (jn - 1))) - 1), &zqxn_index_4, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqlhs + ((jl + ((6 * klon) * (jn - 1))) - 1), &zqlhs_index_8, 1);
                    {
                        double __in1 = zqxn_index_4;
                        double __in2 = zqlhs_index_8;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (__in1 / __in2);
                        ///////////////////

                        zqxn_slice_div_zqlhs_slice_0 = __out;
                    }
                    {
                        double __inp = zqxn_slice_div_zqlhs_slice_0;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1087_16)
                        __out = __inp;
                        ///////////////////

                        zqxn[((jl + (klon * (jn - 1))) - 1)] = __out;
                    }

                }

            }

        }
        for (jn = 1; (jn < 5); jn = (jn + 1)) {

            kfdia_plus_1_73 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_73); jl = (jl + 1)) {

                zqxn_index_5 = zqxn[((jl + (klon * (jn - 1))) - 1)];
                if ((zqxn_index_5 < 1e-14)) {
                    {
                        double zqxn_index_6;
                        double zqxn_index_7;
                        double zqxn_slice_plus_zqxn_slice;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxn + ((jl + (4 * klon)) - 1), &zqxn_index_6, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxn + ((jl + (klon * (jn - 1))) - 1), &zqxn_index_7, 1);
                        {
                            double __in1 = zqxn_index_6;
                            double __in2 = zqxn_index_7;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            zqxn_slice_plus_zqxn_slice = __out;
                        }
                        {
                            double __inp = zqxn_slice_plus_zqxn_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1091_20)
                            __out = __inp;
                            ///////////////////

                            zqxn[((jl + (4 * klon)) - 1)] = __out;
                        }

                    }
                    {

                        {
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1092_20)
                            __out = 0.0;
                            ///////////////////

                            zqxn[((jl + (klon * (jn - 1))) - 1)] = __out;
                        }

                    }
                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            kfdia_plus_1_74 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_74); jl = (jl + 1)) {
                {
                    double zqxnm1_slice;
                    double zqxn2d_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxnm1_slice, 1);
                    {
                        double __inp = zqxnm1_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1095_16)
                        __out = __inp;
                        ///////////////////

                        zqxnm1[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxn2d_slice, 1);
                    {
                        double __inp = zqxn2d_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1096_16)
                        __out = __inp;
                        ///////////////////

                        zqxn2d[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                    }

                }

            }

        }
        for (jm = 1; (jm < (5 + 1)); jm = (jm + 1)) {

            kfdia_plus_1_75 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_75); jl = (jl + 1)) {
                {
                    double zfallsink_index_0;
                    double zqxn_index_8;
                    double zfallsink_slice_times_zqxn_slice;
                    double zrdtgdp_index;
                    double zpfplsx_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zfallsink + ((jl + (klon * (jm - 1))) - 1), &zfallsink_index_0, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxn_index_8, 1);
                    {
                        double __in1 = zfallsink_index_0;
                        double __in2 = zqxn_index_8;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zfallsink_slice_times_zqxn_slice = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zrdtgdp + (jl - 1), &zrdtgdp_index, 1);
                    {
                        double __in1 = zfallsink_slice_times_zqxn_slice;
                        double __in2 = zrdtgdp_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zpfplsx_slice = __out;
                    }
                    {
                        double __inp = zpfplsx_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1099_16)
                        __out = __inp;
                        ///////////////////

                        zpfplsx[((((jk * klon) + jl) + ((klon * (jm - 1)) * (klev + 1))) - 1)] = __out;
                    }

                }

            }

        }
        kfdia_plus_1_76 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_76); jl = (jl + 1)) {
            {
                double zpfplsx_index_2;
                double zpfplsx_index_3;
                double zqpretot_slice_0;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zpfplsx + ((((jk * klon) + jl) + ((3 * klon) * (klev + 1))) - 1), &zpfplsx_index_2, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zpfplsx + ((((jk * klon) + jl) + ((2 * klon) * (klev + 1))) - 1), &zpfplsx_index_3, 1);
                {
                    double __in1 = zpfplsx_index_2;
                    double __in2 = zpfplsx_index_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zqpretot_slice_0 = __out;
                }
                {
                    double __inp = zqpretot_slice_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1101_12)
                    __out = __inp;
                    ///////////////////

                    zqpretot[(jl - 1)] = __out;
                }

            }

        }
        kfdia_plus_1_77 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_77); jl = (jl + 1)) {

            zqpretot_index_2 = zqpretot[(jl - 1)];
            if ((zqpretot_index_2 < 1e-14)) {
                {

                    {
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1104_16)
                        __out = 0.0;
                        ///////////////////

                        zcovptot[(jl - 1)] = __out;
                    }

                }
            }

        }
        for (jm = 1; (jm < 5); jm = (jm + 1)) {

            kfdia_plus_1_78 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_78); jl = (jl + 1)) {
                {
                    double zpsupsatsrce_index;
                    double zconvsrce_index_2;
                    double zpsupsatsrce_slice_plus_zconvsrce_slice;
                    double zfallsrce_index_1;
                    double zpsupsatsrce_slice_zconvsrce_slice_plus_zfallsrce_slice;
                    double zfallsink_index_1;
                    double zconvsink_index;
                    double zfallsink_slice_plus_zconvsink_slice;
                    double zqxn_index_9;
                    double zfallsink_slice_zconvsink_slice_times_zqxn_slice;
                    double zfluxq_slice;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zpsupsatsrce + ((jl + (klon * (jm - 1))) - 1), &zpsupsatsrce_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zconvsrce + ((jl + (klon * (jm - 1))) - 1), &zconvsrce_index_2, 1);
                    {
                        double __in1 = zpsupsatsrce_index;
                        double __in2 = zconvsrce_index_2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zpsupsatsrce_slice_plus_zconvsrce_slice = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zfallsrce + ((jl + (klon * (jm - 1))) - 1), &zfallsrce_index_1, 1);
                    {
                        double __in1 = zpsupsatsrce_slice_plus_zconvsrce_slice;
                        double __in2 = zfallsrce_index_1;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zpsupsatsrce_slice_zconvsrce_slice_plus_zfallsrce_slice = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zfallsink + ((jl + (klon * (jm - 1))) - 1), &zfallsink_index_1, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zconvsink + ((jl + (klon * (jm - 1))) - 1), &zconvsink_index, 1);
                    {
                        double __in1 = zfallsink_index_1;
                        double __in2 = zconvsink_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        zfallsink_slice_plus_zconvsink_slice = __out;
                    }

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxn_index_9, 1);
                    {
                        double __in1 = zfallsink_slice_plus_zconvsink_slice;
                        double __in2 = zqxn_index_9;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zfallsink_slice_zconvsink_slice_times_zqxn_slice = __out;
                    }
                    {
                        double __in1 = zpsupsatsrce_slice_zconvsrce_slice_plus_zfallsrce_slice;
                        double __in2 = zfallsink_slice_zconvsink_slice_times_zqxn_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zfluxq_slice = __out;
                    }
                    {
                        double __inp = zfluxq_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1107_16)
                        __out = __inp;
                        ///////////////////

                        zfluxq[((jl + (klon * (jm - 1))) - 1)] = __out;
                    }

                }

            }
            iphase_index_4 = iphase[(jm - 1)];
            if ((iphase_index_4 == 1)) {

                kfdia_plus_1_79 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_79); jl = (jl + 1)) {
                    {
                        double tendency_loc_t_index_3;
                        double zqxn_index_10;
                        double zqx_index_54;
                        double zqxn_slice_minus_zqx_slice;
                        double zfluxq_index;
                        double zqxn_slice_zqx_slice_minus_zfluxq_slice;
                        double ydthf_ralvdcp_times_zqxn_slice_zqx_slice_zfluxq_slice;
                        double ydthf_ralvdcp_zqxn_slice_zqx_slice_zfluxq_slice_times_zqtmst;
                        double tendency_loc_t_slice_plus_ydthf_ralvdcp_zqxn_slice_zqx_slice_zfluxq_slice_zqtmst;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        tendency_loc_t + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_t_index_3, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxn_index_10, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx_index_54, 1);
                        {
                            double __in1 = zqxn_index_10;
                            double __in2 = zqx_index_54;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqxn_slice_minus_zqx_slice = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zfluxq + ((jl + (klon * (jm - 1))) - 1), &zfluxq_index, 1);
                        {
                            double __in1 = zqxn_slice_minus_zqx_slice;
                            double __in2 = zfluxq_index;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqxn_slice_zqx_slice_minus_zfluxq_slice = __out;
                        }
                        {
                            double __in1 = ydthf_ralvdcp;
                            double __in2 = zqxn_slice_zqx_slice_minus_zfluxq_slice;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            ydthf_ralvdcp_times_zqxn_slice_zqx_slice_zfluxq_slice = __out;
                        }
                        {
                            double __in1 = ydthf_ralvdcp_times_zqxn_slice_zqx_slice_zfluxq_slice;
                            double __in2 = zqtmst;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            ydthf_ralvdcp_zqxn_slice_zqx_slice_zfluxq_slice_times_zqtmst = __out;
                        }
                        {
                            double __in1 = tendency_loc_t_index_3;
                            double __in2 = ydthf_ralvdcp_zqxn_slice_zqx_slice_zfluxq_slice_times_zqtmst;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            tendency_loc_t_slice_plus_ydthf_ralvdcp_zqxn_slice_zqx_slice_zfluxq_slice_zqtmst = __out;
                        }
                        {
                            double __inp = tendency_loc_t_slice_plus_ydthf_ralvdcp_zqxn_slice_zqx_slice_zfluxq_slice_zqtmst;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1110_20)
                            __out = __inp;
                            ///////////////////

                            tendency_loc_t[((jl + (klon * (jk - 1))) - 1)] = __out;
                        }

                    }

                }
            }
            iphase_index_5 = iphase[(jm - 1)];
            if ((iphase_index_5 == 2)) {

                kfdia_plus_1_80 = (kfdia + 1);
                for (jl = kidia; (jl < kfdia_plus_1_80); jl = (jl + 1)) {
                    {
                        double tendency_loc_t_index_4;
                        double zqxn_index_11;
                        double zqx_index_55;
                        double zqxn_slice_minus_zqx_slice_0;
                        double zfluxq_index_0;
                        double zqxn_slice_zqx_slice_minus_zfluxq_slice_0;
                        double ydthf_ralsdcp_times_zqxn_slice_zqx_slice_zfluxq_slice;
                        double ydthf_ralsdcp_zqxn_slice_zqx_slice_zfluxq_slice_times_zqtmst;
                        double tendency_loc_t_slice_plus_ydthf_ralsdcp_zqxn_slice_zqx_slice_zfluxq_slice_zqtmst;


                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        tendency_loc_t + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_t_index_4, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxn_index_11, 1);

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zqx + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx_index_55, 1);
                        {
                            double __in1 = zqxn_index_11;
                            double __in2 = zqx_index_55;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqxn_slice_minus_zqx_slice_0 = __out;
                        }

                        dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                        zfluxq + ((jl + (klon * (jm - 1))) - 1), &zfluxq_index_0, 1);
                        {
                            double __in1 = zqxn_slice_minus_zqx_slice_0;
                            double __in2 = zfluxq_index_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - __in2);
                            ///////////////////

                            zqxn_slice_zqx_slice_minus_zfluxq_slice_0 = __out;
                        }
                        {
                            double __in1 = ydthf_ralsdcp;
                            double __in2 = zqxn_slice_zqx_slice_minus_zfluxq_slice_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            ydthf_ralsdcp_times_zqxn_slice_zqx_slice_zfluxq_slice = __out;
                        }
                        {
                            double __in1 = ydthf_ralsdcp_times_zqxn_slice_zqx_slice_zfluxq_slice;
                            double __in2 = zqtmst;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Mult_)
                            __out = (__in1 * __in2);
                            ///////////////////

                            ydthf_ralsdcp_zqxn_slice_zqx_slice_zfluxq_slice_times_zqtmst = __out;
                        }
                        {
                            double __in1 = tendency_loc_t_index_4;
                            double __in2 = ydthf_ralsdcp_zqxn_slice_zqx_slice_zfluxq_slice_times_zqtmst;
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + __in2);
                            ///////////////////

                            tendency_loc_t_slice_plus_ydthf_ralsdcp_zqxn_slice_zqx_slice_zfluxq_slice_zqtmst = __out;
                        }
                        {
                            double __inp = tendency_loc_t_slice_plus_ydthf_ralsdcp_zqxn_slice_zqx_slice_zfluxq_slice_zqtmst;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_1113_20)
                            __out = __inp;
                            ///////////////////

                            tendency_loc_t[((jl + (klon * (jk - 1))) - 1)] = __out;
                        }

                    }

                }
            }
            kfdia_plus_1_81 = (kfdia + 1);
            for (jl = kidia; (jl < kfdia_plus_1_81); jl = (jl + 1)) {
                {
                    double tendency_loc_cld_index;
                    double zqxn_index_12;
                    double zqx0_index;
                    double zqxn_slice_minus_zqx0_slice;
                    double zqxn_slice_zqx0_slice_times_zqtmst;
                    double tendency_loc_cld_slice_plus_zqxn_slice_zqx0_slice_zqtmst;


                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    tendency_loc_cld + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &tendency_loc_cld_index, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqxn + ((jl + (klon * (jm - 1))) - 1), &zqxn_index_12, 1);

                    dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                    zqx0 + (((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1), &zqx0_index, 1);
                    {
                        double __in1 = zqxn_index_12;
                        double __in2 = zqx0_index;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        zqxn_slice_minus_zqx0_slice = __out;
                    }
                    {
                        double __in1 = zqxn_slice_minus_zqx0_slice;
                        double __in2 = zqtmst;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        zqxn_slice_zqx0_slice_times_zqtmst = __out;
                    }
                    {
                        double __in1 = tendency_loc_cld_index;
                        double __in2 = zqxn_slice_zqx0_slice_times_zqtmst;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        tendency_loc_cld_slice_plus_zqxn_slice_zqx0_slice_zqtmst = __out;
                    }
                    {
                        double __inp = tendency_loc_cld_slice_plus_zqxn_slice_zqx0_slice_zqtmst;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_1115_16)
                        __out = __inp;
                        ///////////////////

                        tendency_loc_cld[(((jl + ((klev * klon) * (jm - 1))) + (klon * (jk - 1))) - 1)] = __out;
                    }

                }

            }

        }
        kfdia_plus_1_82 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_82); jl = (jl + 1)) {
            {
                double tendency_loc_q_index_2;
                double zqxn_index_13;
                double zqx_index_56;
                double zqxn_slice_minus_zqx_slice_1;
                double zqxn_slice_zqx_slice_times_zqtmst;
                double tendency_loc_q_slice_plus_zqxn_slice_zqx_slice_zqtmst;
                double tendency_loc_a_index;
                double zda_index;
                double zda_slice_times_zqtmst;
                double tendency_loc_a_slice_plus_zda_slice_zqtmst;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                tendency_loc_q + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_q_index_2, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxn + ((jl + (4 * klon)) - 1), &zqxn_index_13, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx + (((jl + ((4 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx_index_56, 1);
                {
                    double __in1 = zqxn_index_13;
                    double __in2 = zqx_index_56;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqxn_slice_minus_zqx_slice_1 = __out;
                }
                {
                    double __in1 = zqxn_slice_minus_zqx_slice_1;
                    double __in2 = zqtmst;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqxn_slice_zqx_slice_times_zqtmst = __out;
                }
                {
                    double __in1 = tendency_loc_q_index_2;
                    double __in2 = zqxn_slice_zqx_slice_times_zqtmst;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    tendency_loc_q_slice_plus_zqxn_slice_zqx_slice_zqtmst = __out;
                }
                {
                    double __inp = tendency_loc_q_slice_plus_zqxn_slice_zqx_slice_zqtmst;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1117_12)
                    __out = __inp;
                    ///////////////////

                    tendency_loc_q[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                tendency_loc_a + ((jl + (klon * (jk - 1))) - 1), &tendency_loc_a_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zda + (jl - 1), &zda_index, 1);
                {
                    double __in1 = zda_index;
                    double __in2 = zqtmst;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zda_slice_times_zqtmst = __out;
                }
                {
                    double __in1 = tendency_loc_a_index;
                    double __in2 = zda_slice_times_zqtmst;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    tendency_loc_a_slice_plus_zda_slice_zqtmst = __out;
                }
                {
                    double __inp = tendency_loc_a_slice_plus_zda_slice_zqtmst;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1118_12)
                    __out = __inp;
                    ///////////////////

                    tendency_loc_a[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }
        kfdia_plus_1_83 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_83); jl = (jl + 1)) {
            {
                double pcovptot_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zcovptot + (jl - 1), &pcovptot_slice, 1);
                {
                    double __inp = pcovptot_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1120_12)
                    __out = __inp;
                    ///////////////////

                    pcovptot[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }

    }
    for (jk = 1; (jk < (klev + 2)); jk = (jk + 1)) {

        kfdia_plus_1_84 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_84); jl = (jl + 1)) {
            {
                double zpfplsx_index_4;
                double zpfplsx_index_5;
                double pfplsl_slice;
                double zpfplsx_index_6;
                double zpfplsx_index_7;
                double pfplsn_slice;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zpfplsx + (((jl + (klon * (jk - 1))) + ((2 * klon) * (klev + 1))) - 1), &zpfplsx_index_4, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zpfplsx + ((jl + (klon * (jk - 1))) - 1), &zpfplsx_index_5, 1);
                {
                    double __in1 = zpfplsx_index_4;
                    double __in2 = zpfplsx_index_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfplsl_slice = __out;
                }
                {
                    double __inp = pfplsl_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1123_12)
                    __out = __inp;
                    ///////////////////

                    pfplsl[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zpfplsx + (((jl + (klon * (jk - 1))) + ((3 * klon) * (klev + 1))) - 1), &zpfplsx_index_6, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zpfplsx + (((jl + (klon * (jk - 1))) + (klon * (klev + 1))) - 1), &zpfplsx_index_7, 1);
                {
                    double __in1 = zpfplsx_index_6;
                    double __in2 = zpfplsx_index_7;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfplsn_slice = __out;
                }
                {
                    double __inp = pfplsn_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1124_12)
                    __out = __inp;
                    ///////////////////

                    pfplsn[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }

    }
    kfdia_plus_1_85 = (kfdia + 1);
    for (jl = kidia; (jl < kfdia_plus_1_85); jl = (jl + 1)) {
        {

            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1126_8)
                __out = 0.0;
                ///////////////////

                pfsqlf[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1127_8)
                __out = 0.0;
                ///////////////////

                pfsqif[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1128_8)
                __out = 0.0;
                ///////////////////

                pfsqrf[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1129_8)
                __out = 0.0;
                ///////////////////

                pfsqsf[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1130_8)
                __out = 0.0;
                ///////////////////

                pfcqlng[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1131_8)
                __out = 0.0;
                ///////////////////

                pfcqnng[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1132_8)
                __out = 0.0;
                ///////////////////

                pfcqrng[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1133_8)
                __out = 0.0;
                ///////////////////

                pfcqsng[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1134_8)
                __out = 0.0;
                ///////////////////

                pfsqltur[(jl - 1)] = __out;
            }
            {
                double __out;

                ///////////////////
                // Tasklet code (assign_1135_8)
                __out = 0.0;
                ///////////////////

                pfsqitur[(jl - 1)] = __out;
            }

        }

    }
    for (jk = 1; (jk < (klev + 1)); jk = (jk + 1)) {

        kfdia_plus_1_86 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_86); jl = (jl + 1)) {
            {
                double neg_zrg_r;
                double paph_index_6;
                double paph_index_7;
                double paph_slice_minus_paph_slice;
                double expr_times_paph_slice_paph_slice;
                double pfsqlf_index;
                double pfsqif_index;
                double pfsqrf_slice;
                double pfsqsf_slice;
                double pfcqlng_index;
                double pfcqnng_index;
                double pfcqrng_slice;

                {
                    double __in = zrg_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_USub_)
                    __out = (- __in);
                    ///////////////////

                    neg_zrg_r = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                paph + (((jk * klon) + jl) - 1), &paph_index_6, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                paph + ((jl + (klon * (jk - 1))) - 1), &paph_index_7, 1);
                {
                    double __in1 = paph_index_6;
                    double __in2 = paph_index_7;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    paph_slice_minus_paph_slice = __out;
                }
                {
                    double __in1 = neg_zrg_r;
                    double __in2 = paph_slice_minus_paph_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_times_paph_slice_paph_slice = __out;
                }
                {
                    double __in1 = expr_times_paph_slice_paph_slice;
                    double __in2 = zqtmst;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zgdph_r = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqlf + ((jl + (klon * (jk - 1))) - 1), &pfsqlf_index, 1);
                {
                    double __inp = pfsqlf_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1139_12)
                    __out = __inp;
                    ///////////////////

                    pfsqlf[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqlf + ((jl + (klon * (jk - 1))) - 1), &pfsqrf_slice, 1);
                {
                    double __inp = pfsqrf_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1141_12)
                    __out = __inp;
                    ///////////////////

                    pfsqrf[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqif + ((jl + (klon * (jk - 1))) - 1), &pfsqif_index, 1);
                {
                    double __inp = pfsqif_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1140_12)
                    __out = __inp;
                    ///////////////////

                    pfsqif[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqif + ((jl + (klon * (jk - 1))) - 1), &pfsqsf_slice, 1);
                {
                    double __inp = pfsqsf_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1142_12)
                    __out = __inp;
                    ///////////////////

                    pfsqsf[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfcqlng + ((jl + (klon * (jk - 1))) - 1), &pfcqlng_index, 1);
                {
                    double __inp = pfcqlng_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1143_12)
                    __out = __inp;
                    ///////////////////

                    pfcqlng[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfcqlng + ((jl + (klon * (jk - 1))) - 1), &pfcqrng_slice, 1);
                {
                    double __inp = pfcqrng_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1145_12)
                    __out = __inp;
                    ///////////////////

                    pfcqrng[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfcqnng + ((jl + (klon * (jk - 1))) - 1), &pfcqnng_index, 1);
                {
                    double __inp = pfcqnng_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1144_12)
                    __out = __inp;
                    ///////////////////

                    pfcqnng[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfcqnng + ((jl + (klon * (jk - 1))) - 1), &pfcqsng_slice, 1);

            }
            {
                double pfsqltur_index;
                double pfsqitur_index;
                double zalfaw_1;
                double pfsqlf_index_0;
                double zqxn2d_index;
                double zqx0_index_0;
                double zqxn2d_slice_minus_zqx0_slice;
                double pvfl_index;
                double pvfl_slice_times_ptsphy;
                double zqxn2d_slice_zqx0_slice_plus_pvfl_slice_ptsphy;
                double plude_index_5;
                double zalfaw_times_plude_slice;
                double zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_minus_zalfaw_plude_slice;
                double zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_zalfaw_plude_slice_times_zgdph_r;
                double pfsqlf_slice_plus_zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_zalfaw_plude_slice_zgdph_r;
                double pfcqlng_index_0;
                double zlneg_index_2;
                double zlneg_slice_times_zgdph_r;
                double pfcqlng_slice_plus_zlneg_slice_zgdph_r;
                double pfsqltur_index_0;
                double pvfl_index_0;
                double pvfl_slice_times_ptsphy_0;
                double pvfl_slice_ptsphy_times_zgdph_r;
                double pfsqltur_slice_plus_pvfl_slice_ptsphy_zgdph_r;
                double pfsqrf_index;
                double zqxn2d_index_0;
                double zqx0_index_1;
                double zqxn2d_slice_minus_zqx0_slice_0;
                double zqxn2d_slice_zqx0_slice_times_zgdph_r;
                double pfsqrf_slice_plus_zqxn2d_slice_zqx0_slice_zgdph_r;
                double pfcqrng_index;
                double zlneg_index_3;
                double zlneg_slice_times_zgdph_r_0;
                double pfcqrng_slice_plus_zlneg_slice_zgdph_r;
                double pfsqif_index_0;
                double zqxn2d_index_1;
                double zqx0_index_2;
                double zqxn2d_slice_minus_zqx0_slice_1;
                double pvfi_index;
                double pvfi_slice_times_ptsphy;
                double zqxn2d_slice_zqx0_slice_plus_pvfi_slice_ptsphy;
                double __tmp270;
                double plude_index_6;
                double __tmp271;
                double zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_minus_1_0_zalfaw_plude_slice;
                double zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_1_0_zalfaw_plude_slice_times_zgdph_r;
                double pfsqif_slice_plus_zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_1_0_zalfaw_plude_slice_zgdph_r;
                double pfcqnng_index_0;
                double zlneg_index_4;
                double zlneg_slice_times_zgdph_r_1;
                double pfcqnng_slice_plus_zlneg_slice_zgdph_r;
                double pfsqitur_index_0;
                double pvfi_index_0;
                double pvfi_slice_times_ptsphy_0;
                double pvfi_slice_ptsphy_times_zgdph_r;
                double pfsqitur_slice_plus_pvfi_slice_ptsphy_zgdph_r;
                double pfsqsf_index;
                double zqxn2d_index_2;
                double zqx0_index_3;
                double zqxn2d_slice_minus_zqx0_slice_2;
                double zqxn2d_slice_zqx0_slice_times_zgdph_r_0;
                double pfsqsf_slice_plus_zqxn2d_slice_zqx0_slice_zgdph_r;
                double pfcqsng_index;
                double zlneg_index_5;
                double zlneg_slice_times_zgdph_r_2;
                double pfcqsng_slice_plus_zlneg_slice_zgdph_r;


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqltur + ((jl + (klon * (jk - 1))) - 1), &pfsqltur_index, 1);
                {
                    double __inp = pfsqltur_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1147_12)
                    __out = __inp;
                    ///////////////////

                    pfsqltur[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqltur + (((jk * klon) + jl) - 1), &pfsqltur_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pvfl + ((jl + (klon * (jk - 1))) - 1), &pvfl_index_0, 1);
                {
                    double __in1 = pvfl_index_0;
                    double __in2 = ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    pvfl_slice_times_ptsphy_0 = __out;
                }
                {
                    double __in1 = pvfl_slice_times_ptsphy_0;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    pvfl_slice_ptsphy_times_zgdph_r = __out;
                }
                {
                    double __in1 = pfsqltur_index_0;
                    double __in2 = pvfl_slice_ptsphy_times_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfsqltur_slice_plus_pvfl_slice_ptsphy_zgdph_r = __out;
                }
                {
                    double __inp = pfsqltur_slice_plus_pvfl_slice_ptsphy_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1152_12)
                    __out = __inp;
                    ///////////////////

                    pfsqltur[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfcqlng + (((jk * klon) + jl) - 1), &pfcqlng_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zlneg + ((jl + (klon * (jk - 1))) - 1), &zlneg_index_2, 1);
                {
                    double __in1 = zlneg_index_2;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zlneg_slice_times_zgdph_r = __out;
                }
                {
                    double __in1 = pfcqlng_index_0;
                    double __in2 = zlneg_slice_times_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfcqlng_slice_plus_zlneg_slice_zgdph_r = __out;
                }
                {
                    double __inp = pfcqlng_slice_plus_zlneg_slice_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1151_12)
                    __out = __inp;
                    ///////////////////

                    pfcqlng[(((jk * klon) + jl) - 1)] = __out;
                }
                {
                    double __inp = pfcqsng_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1146_12)
                    __out = __inp;
                    ///////////////////

                    pfcqsng[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfcqsng + (((jk * klon) + jl) - 1), &pfcqsng_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zlneg + (((jl + ((3 * klev) * klon)) + (klon * (jk - 1))) - 1), &zlneg_index_5, 1);
                {
                    double __in1 = zlneg_index_5;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zlneg_slice_times_zgdph_r_2 = __out;
                }
                {
                    double __in1 = pfcqsng_index;
                    double __in2 = zlneg_slice_times_zgdph_r_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfcqsng_slice_plus_zlneg_slice_zgdph_r = __out;
                }
                {
                    double __inp = pfcqsng_slice_plus_zlneg_slice_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1159_12)
                    __out = __inp;
                    ///////////////////

                    pfcqsng[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxn2d + (((jl + ((2 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqxn2d_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx0 + (((jl + ((2 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx0_index_1, 1);
                {
                    double __in1 = zqxn2d_index_0;
                    double __in2 = zqx0_index_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqxn2d_slice_minus_zqx0_slice_0 = __out;
                }
                {
                    double __in1 = zqxn2d_slice_minus_zqx0_slice_0;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqxn2d_slice_zqx0_slice_times_zgdph_r = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqrf + (((jk * klon) + jl) - 1), &pfsqrf_index, 1);
                {
                    double __in1 = pfsqrf_index;
                    double __in2 = zqxn2d_slice_zqx0_slice_times_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfsqrf_slice_plus_zqxn2d_slice_zqx0_slice_zgdph_r = __out;
                }
                {
                    double __inp = pfsqrf_slice_plus_zqxn2d_slice_zqx0_slice_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1153_12)
                    __out = __inp;
                    ///////////////////

                    pfsqrf[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zlneg + (((jl + ((2 * klev) * klon)) + (klon * (jk - 1))) - 1), &zlneg_index_3, 1);
                {
                    double __in1 = zlneg_index_3;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zlneg_slice_times_zgdph_r_0 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfcqrng + (((jk * klon) + jl) - 1), &pfcqrng_index, 1);
                {
                    double __in1 = pfcqrng_index;
                    double __in2 = zlneg_slice_times_zgdph_r_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfcqrng_slice_plus_zlneg_slice_zgdph_r = __out;
                }
                {
                    double __inp = pfcqrng_slice_plus_zlneg_slice_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1154_12)
                    __out = __inp;
                    ///////////////////

                    pfcqrng[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zfoealfa + ((jl + (klon * (jk - 1))) - 1), &zalfaw_1, 1);
                {
                    double __inp = zalfaw_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1149_12)
                    __out = __inp;
                    ///////////////////

                    zalfaw = __out;
                }
                {
                    double __in2 = zalfaw;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp270 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqlf + (((jk * klon) + jl) - 1), &pfsqlf_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxn2d + ((jl + (klon * (jk - 1))) - 1), &zqxn2d_index, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxn2d + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqxn2d_index_1, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx0 + ((jl + (klon * (jk - 1))) - 1), &zqx0_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx0 + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zqx0_index_2, 1);
                {
                    double __in1 = zqxn2d_index;
                    double __in2 = zqx0_index_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqxn2d_slice_minus_zqx0_slice = __out;
                }
                {
                    double __in1 = zqxn2d_index_1;
                    double __in2 = zqx0_index_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqxn2d_slice_minus_zqx0_slice_1 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pvfl + ((jl + (klon * (jk - 1))) - 1), &pvfl_index, 1);
                {
                    double __in1 = pvfl_index;
                    double __in2 = ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    pvfl_slice_times_ptsphy = __out;
                }
                {
                    double __in1 = zqxn2d_slice_minus_zqx0_slice;
                    double __in2 = pvfl_slice_times_ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zqxn2d_slice_zqx0_slice_plus_pvfl_slice_ptsphy = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                plude + ((jl + (klon * (jk - 1))) - 1), &plude_index_5, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                plude + ((jl + (klon * (jk - 1))) - 1), &plude_index_6, 1);
                {
                    double __in1 = zalfaw;
                    double __in2 = plude_index_5;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zalfaw_times_plude_slice = __out;
                }
                {
                    double __in1 = zqxn2d_slice_zqx0_slice_plus_pvfl_slice_ptsphy;
                    double __in2 = zalfaw_times_plude_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_minus_zalfaw_plude_slice = __out;
                }
                {
                    double __in1 = __tmp270;
                    double __in2 = plude_index_6;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp271 = __out;
                }
                {
                    double __in1 = zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_minus_zalfaw_plude_slice;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_zalfaw_plude_slice_times_zgdph_r = __out;
                }
                {
                    double __in1 = pfsqlf_index_0;
                    double __in2 = zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_zalfaw_plude_slice_times_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfsqlf_slice_plus_zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_zalfaw_plude_slice_zgdph_r = __out;
                }
                {
                    double __inp = pfsqlf_slice_plus_zqxn2d_slice_zqx0_slice_pvfl_slice_ptsphy_zalfaw_plude_slice_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1150_12)
                    __out = __inp;
                    ///////////////////

                    pfsqlf[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqif + (((jk * klon) + jl) - 1), &pfsqif_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pvfi + ((jl + (klon * (jk - 1))) - 1), &pvfi_index, 1);
                {
                    double __in1 = pvfi_index;
                    double __in2 = ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    pvfi_slice_times_ptsphy = __out;
                }
                {
                    double __in1 = zqxn2d_slice_minus_zqx0_slice_1;
                    double __in2 = pvfi_slice_times_ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    zqxn2d_slice_zqx0_slice_plus_pvfi_slice_ptsphy = __out;
                }
                {
                    double __in1 = zqxn2d_slice_zqx0_slice_plus_pvfi_slice_ptsphy;
                    double __in2 = __tmp271;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_minus_1_0_zalfaw_plude_slice = __out;
                }
                {
                    double __in1 = zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_minus_1_0_zalfaw_plude_slice;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_1_0_zalfaw_plude_slice_times_zgdph_r = __out;
                }
                {
                    double __in1 = pfsqif_index_0;
                    double __in2 = zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_1_0_zalfaw_plude_slice_times_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfsqif_slice_plus_zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_1_0_zalfaw_plude_slice_zgdph_r = __out;
                }
                {
                    double __inp = pfsqif_slice_plus_zqxn2d_slice_zqx0_slice_pvfi_slice_ptsphy_1_0_zalfaw_plude_slice_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1155_12)
                    __out = __inp;
                    ///////////////////

                    pfsqif[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqitur + ((jl + (klon * (jk - 1))) - 1), &pfsqitur_index, 1);
                {
                    double __inp = pfsqitur_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1148_12)
                    __out = __inp;
                    ///////////////////

                    pfsqitur[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqitur + (((jk * klon) + jl) - 1), &pfsqitur_index_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pvfi + ((jl + (klon * (jk - 1))) - 1), &pvfi_index_0, 1);
                {
                    double __in1 = pvfi_index_0;
                    double __in2 = ptsphy;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    pvfi_slice_times_ptsphy_0 = __out;
                }
                {
                    double __in1 = pvfi_slice_times_ptsphy_0;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    pvfi_slice_ptsphy_times_zgdph_r = __out;
                }
                {
                    double __in1 = pfsqitur_index_0;
                    double __in2 = pvfi_slice_ptsphy_times_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfsqitur_slice_plus_pvfi_slice_ptsphy_zgdph_r = __out;
                }
                {
                    double __inp = pfsqitur_slice_plus_pvfi_slice_ptsphy_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1157_12)
                    __out = __inp;
                    ///////////////////

                    pfsqitur[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zlneg + (((jl + (klev * klon)) + (klon * (jk - 1))) - 1), &zlneg_index_4, 1);
                {
                    double __in1 = zlneg_index_4;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zlneg_slice_times_zgdph_r_1 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfcqnng + (((jk * klon) + jl) - 1), &pfcqnng_index_0, 1);
                {
                    double __in1 = pfcqnng_index_0;
                    double __in2 = zlneg_slice_times_zgdph_r_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfcqnng_slice_plus_zlneg_slice_zgdph_r = __out;
                }
                {
                    double __inp = pfcqnng_slice_plus_zlneg_slice_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1156_12)
                    __out = __inp;
                    ///////////////////

                    pfcqnng[(((jk * klon) + jl) - 1)] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqxn2d + (((jl + ((3 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqxn2d_index_2, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                zqx0 + (((jl + ((3 * klev) * klon)) + (klon * (jk - 1))) - 1), &zqx0_index_3, 1);
                {
                    double __in1 = zqxn2d_index_2;
                    double __in2 = zqx0_index_3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqxn2d_slice_minus_zqx0_slice_2 = __out;
                }
                {
                    double __in1 = zqxn2d_slice_minus_zqx0_slice_2;
                    double __in2 = zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zqxn2d_slice_zqx0_slice_times_zgdph_r_0 = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfsqsf + (((jk * klon) + jl) - 1), &pfsqsf_index, 1);
                {
                    double __in1 = pfsqsf_index;
                    double __in2 = zqxn2d_slice_zqx0_slice_times_zgdph_r_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    pfsqsf_slice_plus_zqxn2d_slice_zqx0_slice_zgdph_r = __out;
                }
                {
                    double __inp = pfsqsf_slice_plus_zqxn2d_slice_zqx0_slice_zgdph_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1158_12)
                    __out = __inp;
                    ///////////////////

                    pfsqsf[(((jk * klon) + jl) - 1)] = __out;
                }

            }

        }

    }
    for (jk = 1; (jk < (klev + 2)); jk = (jk + 1)) {

        kfdia_plus_1_87 = (kfdia + 1);
        for (jl = kidia; (jl < kfdia_plus_1_87); jl = (jl + 1)) {
            {
                double neg_ydcst_rlvtt;
                double pfplsl_index;
                double pfhpsl_slice;
                double neg_ydcst_rlstt;
                double pfplsn_index;
                double pfhpsn_slice;

                {
                    double __in = ydcst_rlvtt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_USub_)
                    __out = (- __in);
                    ///////////////////

                    neg_ydcst_rlvtt = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfplsl + ((jl + (klon * (jk - 1))) - 1), &pfplsl_index, 1);
                {
                    double __in1 = neg_ydcst_rlvtt;
                    double __in2 = pfplsl_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    pfhpsl_slice = __out;
                }
                {
                    double __inp = pfhpsl_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1162_12)
                    __out = __inp;
                    ///////////////////

                    pfhpsl[((jl + (klon * (jk - 1))) - 1)] = __out;
                }
                {
                    double __in = ydcst_rlstt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_USub_)
                    __out = (- __in);
                    ///////////////////

                    neg_ydcst_rlstt = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                pfplsn + ((jl + (klon * (jk - 1))) - 1), &pfplsn_index, 1);
                {
                    double __in1 = neg_ydcst_rlstt;
                    double __in2 = pfplsn_index;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    pfhpsn_slice = __out;
                }
                {
                    double __inp = pfhpsn_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_1163_12)
                    __out = __inp;
                    ///////////////////

                    pfhpsn[((jl + (klon * (jk - 1))) - 1)] = __out;
                }

            }

        }

    }
    delete[] zlcond1;
    delete[] zlcond2;
    delete[] zrainaut;
    delete[] zsnowaut;
    delete[] zliqcld;
    delete[] zicecld;
    delete[] zfokoop;
    delete[] zlicld;
    delete[] zlfinalsum;
    delete[] zdqs;
    delete[] ztold;
    delete[] zqold;
    delete[] zdtgdp;
    delete[] zrdtgdp;
    delete[] zcovpclr;
    delete[] zcovptot;
    delete[] zcovpmax;
    delete[] zqpretot;
    delete[] zldefr;
    delete[] zldifdt;
    delete[] zacust;
    delete[] zmf;
    delete[] zrho;
    delete[] zsolab;
    delete[] zsolac;
    delete[] zanewm1;
    delete[] zda;
    delete[] zdp;
    delete[] zsupsat;
    delete[] zmeltmax;
    delete[] zfrzmax;
    delete[] zicetot;
    delete[] zdqsmixdt;
    delete[] zcorqsice;
    delete[] zevaplimmix;
    delete[] zcldtopdist;
    delete[] zrainacc;
    delete[] zraincld;
    delete[] zsnowcld;
    delete[] psum_solqa;
    delete[] iphase;
    delete[] imelt;
    delete[] llfall;
    delete[] zvqx;
    delete[] zfoealfa;
    delete[] ztp1;
    delete[] zlcust;
    delete[] zli;
    delete[] za;
    delete[] zaorig;
    delete[] zliqfrac;
    delete[] zicefrac;
    delete[] zqx;
    delete[] zqx0;
    delete[] zqxn;
    delete[] zqxfg;
    delete[] zqxnm1;
    delete[] zfluxq;
    delete[] zpfplsx;
    delete[] zlneg;
    delete[] zqxn2d;
    delete[] zqsmix;
    delete[] zqsliq;
    delete[] zqsice;
    delete[] zfoeewmt;
    delete[] zfoeew;
    delete[] zfoeeliqt;
    delete[] zsolqa;
    delete[] zsolqb;
    delete[] zqlhs;
    delete[] zratio;
    delete[] zsinksum;
    delete[] zfallsink;
    delete[] zfallsrce;
    delete[] zconvsrce;
    delete[] zconvsink;
    delete[] zpsupsatsrce;
}

DACE_EXPORTED void __program_cloudsc_py(cloudsc_py_state_t *__state, int * __restrict__ ktype, int * __restrict__ ldcum, double * __restrict__ pa, double * __restrict__ pap, double * __restrict__ paph, double * __restrict__ pccn, double * __restrict__ pclv, double * __restrict__ pcovptot, double * __restrict__ pdyna, double * __restrict__ pdyni, double * __restrict__ pdynl, double * __restrict__ pfcqlng, double * __restrict__ pfcqnng, double * __restrict__ pfcqrng, double * __restrict__ pfcqsng, double * __restrict__ pfhpsl, double * __restrict__ pfhpsn, double * __restrict__ pfplsl, double * __restrict__ pfplsn, double * __restrict__ pfsqif, double * __restrict__ pfsqitur, double * __restrict__ pfsqlf, double * __restrict__ pfsqltur, double * __restrict__ pfsqrf, double * __restrict__ pfsqsf, double * __restrict__ phrlw, double * __restrict__ phrsw, double * __restrict__ picrit_aer, double * __restrict__ plcrit_aer, double * __restrict__ plsm, double * __restrict__ plu, double * __restrict__ plude, double * __restrict__ pmfd, double * __restrict__ pmfu, double * __restrict__ pnice, double * __restrict__ pq, double * __restrict__ prainfrac_toprfz, double * __restrict__ pre_ice, double * __restrict__ psnde, double * __restrict__ psupsat, double * __restrict__ pt, double * __restrict__ pvervel, double * __restrict__ pvfa, double * __restrict__ pvfi, double * __restrict__ pvfl, double * __restrict__ tendency_loc_a, double * __restrict__ tendency_loc_cld, double * __restrict__ tendency_loc_q, double * __restrict__ tendency_loc_t, double * __restrict__ tendency_tmp_a, double * __restrict__ tendency_tmp_cld, double * __restrict__ tendency_tmp_q, double * __restrict__ tendency_tmp_t, int kfdia, int kidia, int klev, int klon, double ptsphy, double ydcst_rcpd, double ydcst_rd, double ydcst_retv, double ydcst_rg, double ydcst_rlmlt, double ydcst_rlstt, double ydcst_rlvtt, double ydcst_rtt, double ydcst_rv, double ydthf_r2es, double ydthf_r3ies, double ydthf_r3les, double ydthf_r4ies, double ydthf_r4les, double ydthf_r5alscp, double ydthf_r5alvcp, double ydthf_r5ies, double ydthf_r5les, double ydthf_ralfdcp, double ydthf_ralsdcp, double ydthf_ralvdcp, double ydthf_rkoop1, double ydthf_rkoop2, double ydthf_rtice, double ydthf_rticecu, double ydthf_rtwat, double ydthf_rtwat_rtice_r, double ydthf_rtwat_rticecu_r, int yrecldp_laericeauto, int yrecldp_laericesed, int yrecldp_laerliqautolsp, int yrecldp_laerliqcoll, int yrecldp_ncldtop, int yrecldp_nssopt, double yrecldp_ramid, double yrecldp_ramin, double yrecldp_rccn, double yrecldp_rcl_apb1, double yrecldp_rcl_apb2, double yrecldp_rcl_apb3, double yrecldp_rcl_cdenom1, double yrecldp_rcl_cdenom2, double yrecldp_rcl_cdenom3, double yrecldp_rcl_const1i, double yrecldp_rcl_const1r, double yrecldp_rcl_const1s, double yrecldp_rcl_const2i, double yrecldp_rcl_const2r, double yrecldp_rcl_const2s, double yrecldp_rcl_const3i, double yrecldp_rcl_const3r, double yrecldp_rcl_const3s, double yrecldp_rcl_const4i, double yrecldp_rcl_const4r, double yrecldp_rcl_const4s, double yrecldp_rcl_const5i, double yrecldp_rcl_const5r, double yrecldp_rcl_const5s, double yrecldp_rcl_const6i, double yrecldp_rcl_const6r, double yrecldp_rcl_const6s, double yrecldp_rcl_const7s, double yrecldp_rcl_const8s, double yrecldp_rcl_fac1, double yrecldp_rcl_fac2, double yrecldp_rcl_fzrab, double yrecldp_rcl_ka273, double yrecldp_rcl_kk_cloud_num_land, double yrecldp_rcl_kk_cloud_num_sea, double yrecldp_rcl_kkaac, double yrecldp_rcl_kkaau, double yrecldp_rcl_kkbac, double yrecldp_rcl_kkbaun, double yrecldp_rcl_kkbauq, double yrecldp_rclcrit_land, double yrecldp_rclcrit_sea, double yrecldp_rcldiff, double yrecldp_rcldiff_convi, double yrecldp_rcldtopcf, double yrecldp_rcovpmin, double yrecldp_rdensref, double yrecldp_rdepliqrefdepth, double yrecldp_rdepliqrefrate, double yrecldp_riceinit, double yrecldp_rkconv, double yrecldp_rkooptau, double yrecldp_rlcritsnow, double yrecldp_rlmin, double yrecldp_rnice, double yrecldp_rpecons, double yrecldp_rprc1, double yrecldp_rprecrhmax, double yrecldp_rsnowlin1, double yrecldp_rsnowlin2, double yrecldp_rtaumel, double yrecldp_rthomo, double yrecldp_rvice, double yrecldp_rvrain, double yrecldp_rvrfactor, double yrecldp_rvsnow)
{
    __program_cloudsc_py_internal(__state, ktype, ldcum, pa, pap, paph, pccn, pclv, pcovptot, pdyna, pdyni, pdynl, pfcqlng, pfcqnng, pfcqrng, pfcqsng, pfhpsl, pfhpsn, pfplsl, pfplsn, pfsqif, pfsqitur, pfsqlf, pfsqltur, pfsqrf, pfsqsf, phrlw, phrsw, picrit_aer, plcrit_aer, plsm, plu, plude, pmfd, pmfu, pnice, pq, prainfrac_toprfz, pre_ice, psnde, psupsat, pt, pvervel, pvfa, pvfi, pvfl, tendency_loc_a, tendency_loc_cld, tendency_loc_q, tendency_loc_t, tendency_tmp_a, tendency_tmp_cld, tendency_tmp_q, tendency_tmp_t, kfdia, kidia, klev, klon, ptsphy, ydcst_rcpd, ydcst_rd, ydcst_retv, ydcst_rg, ydcst_rlmlt, ydcst_rlstt, ydcst_rlvtt, ydcst_rtt, ydcst_rv, ydthf_r2es, ydthf_r3ies, ydthf_r3les, ydthf_r4ies, ydthf_r4les, ydthf_r5alscp, ydthf_r5alvcp, ydthf_r5ies, ydthf_r5les, ydthf_ralfdcp, ydthf_ralsdcp, ydthf_ralvdcp, ydthf_rkoop1, ydthf_rkoop2, ydthf_rtice, ydthf_rticecu, ydthf_rtwat, ydthf_rtwat_rtice_r, ydthf_rtwat_rticecu_r, yrecldp_laericeauto, yrecldp_laericesed, yrecldp_laerliqautolsp, yrecldp_laerliqcoll, yrecldp_ncldtop, yrecldp_nssopt, yrecldp_ramid, yrecldp_ramin, yrecldp_rccn, yrecldp_rcl_apb1, yrecldp_rcl_apb2, yrecldp_rcl_apb3, yrecldp_rcl_cdenom1, yrecldp_rcl_cdenom2, yrecldp_rcl_cdenom3, yrecldp_rcl_const1i, yrecldp_rcl_const1r, yrecldp_rcl_const1s, yrecldp_rcl_const2i, yrecldp_rcl_const2r, yrecldp_rcl_const2s, yrecldp_rcl_const3i, yrecldp_rcl_const3r, yrecldp_rcl_const3s, yrecldp_rcl_const4i, yrecldp_rcl_const4r, yrecldp_rcl_const4s, yrecldp_rcl_const5i, yrecldp_rcl_const5r, yrecldp_rcl_const5s, yrecldp_rcl_const6i, yrecldp_rcl_const6r, yrecldp_rcl_const6s, yrecldp_rcl_const7s, yrecldp_rcl_const8s, yrecldp_rcl_fac1, yrecldp_rcl_fac2, yrecldp_rcl_fzrab, yrecldp_rcl_ka273, yrecldp_rcl_kk_cloud_num_land, yrecldp_rcl_kk_cloud_num_sea, yrecldp_rcl_kkaac, yrecldp_rcl_kkaau, yrecldp_rcl_kkbac, yrecldp_rcl_kkbaun, yrecldp_rcl_kkbauq, yrecldp_rclcrit_land, yrecldp_rclcrit_sea, yrecldp_rcldiff, yrecldp_rcldiff_convi, yrecldp_rcldtopcf, yrecldp_rcovpmin, yrecldp_rdensref, yrecldp_rdepliqrefdepth, yrecldp_rdepliqrefrate, yrecldp_riceinit, yrecldp_rkconv, yrecldp_rkooptau, yrecldp_rlcritsnow, yrecldp_rlmin, yrecldp_rnice, yrecldp_rpecons, yrecldp_rprc1, yrecldp_rprecrhmax, yrecldp_rsnowlin1, yrecldp_rsnowlin2, yrecldp_rtaumel, yrecldp_rthomo, yrecldp_rvice, yrecldp_rvrain, yrecldp_rvrfactor, yrecldp_rvsnow);
}

DACE_EXPORTED cloudsc_py_state_t *__dace_init_cloudsc_py(int klev, int klon)
{

    int __result = 0;
    cloudsc_py_state_t *__state = new cloudsc_py_state_t;

    if (__result) {
        delete __state;
        return nullptr;
    }

    return __state;
}

DACE_EXPORTED int __dace_exit_cloudsc_py(cloudsc_py_state_t *__state)
{

    int __err = 0;
    delete __state;
    return __err;
}
