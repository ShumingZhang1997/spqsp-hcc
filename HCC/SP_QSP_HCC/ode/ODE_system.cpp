#include "ODE_system.h"
    
#define SPVAR(x) NV_DATA_S(y)[x]
#define NSPVAR(x) ptrOde->_nonspecies_var[x]
#define PARAM(x) _class_parameter[x]
#define PFILE(x) param.getVal(x)

namespace CancerVCT{
#define QSP_W ODE_system::_QSP_weight

bool ODE_system::use_steady_state = false;
bool ODE_system::use_resection = false;
double ODE_system::_QSP_weight = 1.0;

ODE_system::ODE_system()
:CVODEBase()
{
    setupVariables();
    setupEvents();
    setupCVODE();
    update_y_other();
}

ODE_system::ODE_system(const ODE_system& c)
{
    setupCVODE();
}

ODE_system::~ODE_system()
{
}

void ODE_system::initSolver(realtype t){

    restore_y();
    int flag;

    flag = CVodeInit(_cvode_mem, f, t, _y);
    check_flag(&flag, "CVodeInit", 1);

    /* Call CVodeRootInit to specify the root function g */
    flag = CVodeRootInit(_cvode_mem, _nroot, g);
    check_flag(&flag, "CVodeRootInit", 1);
    
    	/*Do not do this. Event only trigger when turn from false to true.
	  If this is reset before trigger evaluation at the beginning of simulation,
	  t=0 events might be missed.*/
    //updateTriggerComponentConditionsOnValue(t);
    //resetEventTriggers();

    return;
} 

state_type ODE_system::_class_parameter = state_type(193, 0);

void ODE_system::setup_class_parameters(Param& param){
    //V_C, mwbc732826_b725_481d_9afa_ad349757f63b, index: 0
    //Unit: metre^(3)
    _class_parameter[0] = PFILE(5) * 0.0010000000000000007;
    //V_P, mw93376915_721f_4994_a60c_bfa753cf2b47, index: 1
    //Unit: metre^(3)
    _class_parameter[1] = PFILE(6) * 0.0010000000000000007;
    //V_LN, mw6aaa8bb3_c3b2_4e25_a81e_5fb91c8c01eb, index: 2
    //Unit: metre^(3)
    _class_parameter[2] = PFILE(8) * 1.0000000000000006e-06;
    //V_e, mw780ee5b3_cb6f_423c_b6e5_b240087e9705, index: 3
    //Unit: metre^(3)
    _class_parameter[3] = PFILE(9) * 0.0010000000000000007;
    //A_e, mwa1f28a11_5416_417b_a9cd_3668f1346fad, index: 4
    //Unit: metre^(2)
    _class_parameter[4] = PFILE(10) * 1e-12;
    //A_s, mwd7401aea_a0d9_4b12_8f00_61bf492fd2e1, index: 5
    //Unit: metre^(2)
    _class_parameter[5] = PFILE(11) * 1e-12;
    //syn_CT, mw23b3bab2_cec3_42bf_b415_24b3184045c4, index: 6
    //Unit: metre^(2)
    _class_parameter[6] = PFILE(12) * 1e-12;
    //syn_DN, mw876fbeb3_d6ce_4c37_8aac_d21e10f281e0, index: 7
    //Unit: metre^(2)
    _class_parameter[7] = PFILE(13) * 1e-12;
    //Treg_P, mw39a0ed72_200b_477d_a26b_f191dee7473b, index: 8
    //Unit: metre^(2)
    _class_parameter[8] = PFILE(14) * 1e-12;
    //Treg_T, mw9898b54a_0dcc_4701_81d0_6746f2f1668d, index: 9
    //Unit: metre^(2)
    _class_parameter[9] = PFILE(15) * 1e-12;
    //syn_CT.PDL1_aPDL1, mw4769224a_85bf_4dc3_a5dd_229e3f56c9b0, index: 10
    //Unit: mole^(1)metre^(-2)
    _class_parameter[10] = PFILE(68) * 1.66053872801495e-12;
    //syn_CT.PDL1_aPDL1_PDL1, mw4147800f_beb3_4f47_b558_83e065b456dc, index: 11
    //Unit: mole^(1)metre^(-2)
    _class_parameter[11] = PFILE(69) * 1.66053872801495e-12;
    //syn_DN.CD28_CD80, mwbe693294_0273_4fa9_952a_abd7abefa74f, index: 12
    //Unit: mole^(1)metre^(-2)
    _class_parameter[12] = PFILE(70) * 1.66053872801495e-12;
    //syn_DN.CD28_CD80_CD28, mwa570af6e_3dd5_4fa6_900b_4ac846f13d21, index: 13
    //Unit: mole^(1)metre^(-2)
    _class_parameter[13] = PFILE(71) * 1.66053872801495e-12;
    //syn_DN.CD28_CD86, mwc481dd0a_ef3b_4101_a68c_de28092c1697, index: 14
    //Unit: mole^(1)metre^(-2)
    _class_parameter[14] = PFILE(72) * 1.66053872801495e-12;
    //syn_DN.PDL1_CD80, mw0b7a350f_5d7c_4252_bb37_7fd4eddce557, index: 15
    //Unit: mole^(1)metre^(-2)
    _class_parameter[15] = PFILE(73) * 1.66053872801495e-12;
    //syn_DN.CTLA4_CD80, mw50223b59_1047_42a6_bf47_848b0a86bf1e, index: 16
    //Unit: mole^(1)metre^(-2)
    _class_parameter[16] = PFILE(74) * 1.66053872801495e-12;
    //syn_DN.CTLA4_CD80_CTLA4, mwc9351bd5_2f95_4e42_bd22_ecdc53984d6d, index: 17
    //Unit: mole^(1)metre^(-2)
    _class_parameter[17] = PFILE(75) * 1.66053872801495e-12;
    //syn_DN.CD80_CTLA4_CD80, mw1f313d18_a65d_4201_8c3d_09e0d87676b5, index: 18
    //Unit: mole^(1)metre^(-2)
    _class_parameter[18] = PFILE(76) * 1.66053872801495e-12;
    //syn_DN.CTLA4_CD80_CTLA4_CD80, mw62cde88c_00aa_4d27_93dd_6f2e42bcc710, index: 19
    //Unit: mole^(1)metre^(-2)
    _class_parameter[19] = PFILE(77) * 1.66053872801495e-12;
    //syn_DN.CTLA4_CD86, mw8203fab7_e56c_49ec_baae_94c3362ee52b, index: 20
    //Unit: mole^(1)metre^(-2)
    _class_parameter[20] = PFILE(78) * 1.66053872801495e-12;
    //syn_DN.CD86_CTLA4_CD86, mw7fcd03da_022c_41b0_9141_d54cc37b8720, index: 21
    //Unit: mole^(1)metre^(-2)
    _class_parameter[21] = PFILE(79) * 1.66053872801495e-12;
    //syn_DN.CTLA4_aCTLA4, mw894d7002_ad5e_4ecf_b510_9744e74c661c, index: 22
    //Unit: mole^(1)metre^(-2)
    _class_parameter[22] = PFILE(80) * 1.66053872801495e-12;
    //syn_DN.CTLA4_aCTLA4_CTLA4, mw3997908f_73de_47eb_9591_d8addb9c4f94, index: 23
    //Unit: mole^(1)metre^(-2)
    _class_parameter[23] = PFILE(81) * 1.66053872801495e-12;
    //syn_DN.PDL1_aPDL1, mwba5c5a4f_db68_4d9c_b67b_0335bd02b464, index: 24
    //Unit: mole^(1)metre^(-2)
    _class_parameter[24] = PFILE(82) * 1.66053872801495e-12;
    //syn_DN.PDL1_aPDL1_PDL1, mw48f83962_a996_483a_9f72_ca7454c22dd9, index: 25
    //Unit: mole^(1)metre^(-2)
    _class_parameter[25] = PFILE(83) * 1.66053872801495e-12;
    //Treg_P.CTLA4_aCTLA4, mw5968c89a_8754_44d2_ac32_5fcce49beea0, index: 26
    //Unit: mole^(1)metre^(-2)
    _class_parameter[26] = PFILE(84) * 1.66053872801495e-12;
    //Treg_P.CTLA4_aCTLA4_CTLA4, mwe2ca6e7c_7c55_4f60_9328_7c5775891b97, index: 27
    //Unit: mole^(1)metre^(-2)
    _class_parameter[27] = PFILE(85) * 1.66053872801495e-12;
    //Treg_T.CTLA4_aCTLA4, mwf19e5b25_78ce_456f_9aeb_75949364479a, index: 28
    //Unit: mole^(1)metre^(-2)
    _class_parameter[28] = PFILE(86) * 1.66053872801495e-12;
    //Treg_T.CTLA4_aCTLA4_CTLA4, mwa663e7d9_4e05_421a_8eaa_33b0e09f815d, index: 29
    //Unit: mole^(1)metre^(-2)
    _class_parameter[29] = PFILE(87) * 1.66053872801495e-12;
    //k_cell_clear, mw3736d194_dc0e_4644_bbe1_e5e58b544242, index: 30
    //Unit: second^(-1)
    _class_parameter[30] = PFILE(88) * 1.15740740740741e-05;
    //cell, mw89b2e6bc_1636_41b3_b28c_387dd021a3d6, index: 31
    //Unit: mole^(1)
    _class_parameter[31] = PFILE(89) * 1.66053872801495e-24;
    //day, mwfcf7adad_26aa_4d37_9d95_a395acfba16a, index: 32
    //Unit: second^(1)
    _class_parameter[32] = PFILE(90) * 86400.0;
    //cellmcm2, mw72cae028_7753_43d9_96f8_2a9f5bde3a22, index: 33
    //Unit: mole^(1)
    _class_parameter[33] = PFILE(91) * 1.66053872801495e-24;
    //molarity, mw597a9262_cb9a_46bc_a95b_71ea4cb95d08, index: 34
    //Unit: metre^(-3)mole^(1)
    _class_parameter[34] = PFILE(92) * 999.9999999999994;
    //vol_cell, mw42f51148_3733_4f84_959c_03584162e095, index: 35
    //Unit: metre^(3)mole^(-1)
    _class_parameter[35] = PFILE(93) * 602214.1989999996;
    //vol_Tcell, mw84d4a8ab_09ee_4330_9b2d_3f0d043e687d, index: 36
    //Unit: metre^(3)mole^(-1)
    _class_parameter[36] = PFILE(94) * 602214.1989999996;
    //V_Tmin, mw070b0542_5a65_474a_a091_b1d376fc6de9, index: 37
    //Unit: metre^(3)
    _class_parameter[37] = PFILE(95) * 1.0000000000000013e-09;
    //blood_volume, mw2dc08fce_859d_498c_9c5d_e2f25208bf75, index: 38
    //Unit: metre^(3)
    _class_parameter[38] = PFILE(108) * 0.0010000000000000007;
    //peripheral_volume, mwab68b6bd_73fb_41ed_bcab_c9ad82e5700c, index: 39
    //Unit: metre^(3)
    _class_parameter[39] = PFILE(109) * 0.0010000000000000007;
    //LN_volume, mwf0480d87_c182_4e9e_be7e_9e74c8df9764, index: 40
    //Unit: metre^(3)
    _class_parameter[40] = PFILE(110) * 1.0000000000000006e-06;
    //H_CD28, mwa65730ee_acf3_4370_9e4d_3a37d33d3bb0, index: 41
    //Unit: dimensionless^(1)
    _class_parameter[41] = PFILE(112) * 1.0;
    //k_C1_growth, mweb294556_ca0e_4267_aba9_446f0611b07e, index: 42
    //Unit: second^(-1)
    _class_parameter[42] = PFILE(113) * 1.15740740740741e-05;
    //k_C1_death, mw9563f919_e491_4b16_a20a_252bca40ab81, index: 43
    //Unit: second^(-1)
    _class_parameter[43] = PFILE(115) * 1.15740740740741e-05;
    //initial_tumour_diameter, mwc55a0fd7_9db3_46c3_9b5f_178b5d171bdc, index: 44
    //Unit: metre^(1)
    _class_parameter[44] = PFILE(116) * 0.01;
    //div_T1, mwcd41b01a_bb9f_426a_b44e_3c34f28fe575, index: 45
    //Unit: dimensionless^(1)
    _class_parameter[45] = PFILE(117) * 1.0;
    //n_T1_clones, mwafa56f2d_151f_4c3c_bb9e_b9a6635d9841, index: 46
    //Unit: dimensionless^(1)
    _class_parameter[46] = PFILE(118) * 1.0;
    //q_T1_LN_in, mw4e5b695b_bb0c_4f25_9aa6_ff5aed37ec23, index: 47
    //Unit: metre^(-3)second^(-1)
    _class_parameter[47] = PFILE(119) * 11.574074074074097;
    //q_T1_LN_out, mw45f3ed36_9ec6_4b43_bf4c_557ee7911cdf, index: 48
    //Unit: second^(-1)
    _class_parameter[48] = PFILE(120) * 1.15740740740741e-05;
    //k_T1_act, mwa6f23259_f9a0_4079_b24e_a7114e9cd682, index: 49
    //Unit: second^(-1)
    _class_parameter[49] = PFILE(121) * 1.15740740740741e-05;
    //k_T1_pro, mw75a5db23_de85_4f8c_bfd2_06d30e38bec4, index: 50
    //Unit: second^(-1)
    _class_parameter[50] = PFILE(122) * 1.15740740740741e-05;
    //k_T1_death, mw85cad532_0cc2_4190_8830_ce832a8d9466, index: 51
    //Unit: second^(-1)
    _class_parameter[51] = PFILE(123) * 1.15740740740741e-05;
    //q_T1_P_in, mw27c14556_d5fe_4765_af30_d75d13cc729c, index: 52
    //Unit: metre^(-3)second^(-1)
    _class_parameter[52] = PFILE(124) * 11.574074074074097;
    //q_T1_P_out, mw14e1796f_b7e0_4166_8508_487601fa39f9, index: 53
    //Unit: second^(-1)
    _class_parameter[53] = PFILE(125) * 1.15740740740741e-05;
    //q_T1_T_in, mw21066b8a_a7bd_4e5a_a511_5033654eff1e, index: 54
    //Unit: metre^(-3)second^(-1)
    _class_parameter[54] = PFILE(126) * 11.574074074074097;
    //Q_nT1_thym, mwd61ad4b1_cd7a_4ba0_af68_aa2b1455a97a, index: 55
    //Unit: mole^(1)second^(-1)
    _class_parameter[55] = PFILE(127) * 1.92191982409137e-29;
    //k_nT1_pro, mw186f5aee_e0dd_421a_bf24_461f84edbdc1, index: 56
    //Unit: mole^(1)second^(-1)
    _class_parameter[56] = PFILE(128) * 1.92191982409137e-29;
    //K_nT1_pro, mw01bf65b2_1a94_46df_a7bd_1ea0b6eb7727, index: 57
    //Unit: mole^(1)
    _class_parameter[57] = PFILE(129) * 1.66053872801495e-24;
    //k_nT1_death, mw253f9ad7_a317_43f2_bb51_5ae85c61c641, index: 58
    //Unit: second^(-1)
    _class_parameter[58] = PFILE(130) * 1.15740740740741e-05;
    //k_IL2_deg, mw47f53f24_6205_4871_833d_bcd0bfa94f8e, index: 59
    //Unit: second^(-1)
    _class_parameter[59] = PFILE(131) * 1.0;
    //k_IL2_cons, mw69bf5f42_12df_456c_b95f_1b01b746374a, index: 60
    //Unit: second^(-1)
    _class_parameter[60] = PFILE(132) * 1.15740740740741e-05;
    //k_IL2_sec, mwf076e74c_809f_48a5_bfee_91ccd0a3216e, index: 61
    //Unit: second^(-1)
    _class_parameter[61] = PFILE(133) * 1.0;
    //IL2_50, mwacf95e3f_1404_4842_90ca_e4cff2725b6d, index: 62
    //Unit: metre^(-3)mole^(1)
    _class_parameter[62] = PFILE(134) * 999.9999999999994;
    //N0, mw400a7b2d_8783_4b28_93a7_7566add719a9, index: 63
    //Unit: dimensionless^(1)
    _class_parameter[63] = PFILE(135) * 1.0;
    //N_IL2, mw04cd6136_b750_4f06_826b_c419fa1185db, index: 64
    //Unit: dimensionless^(1)
    _class_parameter[64] = PFILE(136) * 1.0;
    //k_T1, mw209ee6ab_be81_4d49_9124_e8d6ebad85e4, index: 65
    //Unit: second^(-1)
    _class_parameter[65] = PFILE(137) * 1.15740740740741e-05;
    //k_C_T1, mwcf0a8e62_19bc_4b53_8d55_a9be2d57fe7b, index: 66
    //Unit: second^(-1)
    _class_parameter[66] = PFILE(138) * 1.15740740740741e-05;
    //k_Treg, mwa340c013_d826_459a_acd9_9ef5d8ac7967, index: 67
    //Unit: second^(-1)
    _class_parameter[67] = PFILE(139) * 1.15740740740741e-05;
    //div_T0, mw4532b02e_8399_424f_ae29_ecac5033fc39, index: 68
    //Unit: dimensionless^(1)
    _class_parameter[68] = PFILE(143) * 1.0;
    //n_T0_clones, mwdfedf252_ad12_4d3f_a8ad_839e08981f95, index: 69
    //Unit: dimensionless^(1)
    _class_parameter[69] = PFILE(144) * 1.0;
    //q_T0_LN_in, mwce2513a5_019b_4ec9_b7aa_f58ca6141fdb, index: 70
    //Unit: metre^(-3)second^(-1)
    _class_parameter[70] = PFILE(145) * 11.574074074074097;
    //q_T0_LN_out, mw58fec153_e0bd_4da3_bb66_c22fc923f91e, index: 71
    //Unit: second^(-1)
    _class_parameter[71] = PFILE(146) * 1.15740740740741e-05;
    //k_T0_act, mw8c0227b0_0818_46a3_a157_6a4f0145601f, index: 72
    //Unit: second^(-1)
    _class_parameter[72] = PFILE(147) * 1.15740740740741e-05;
    //k_T0_pro, mw2fe03360_be47_4fb6_8e71_6135b4af0715, index: 73
    //Unit: second^(-1)
    _class_parameter[73] = PFILE(148) * 1.15740740740741e-05;
    //k_T0_death, mwdba2f84c_1542_4b7f_960a_7af03036a097, index: 74
    //Unit: second^(-1)
    _class_parameter[74] = PFILE(149) * 1.15740740740741e-05;
    //q_T0_P_in, mw724caa3b_953b_43b4_8061_8a2d8e3e82c6, index: 75
    //Unit: metre^(-3)second^(-1)
    _class_parameter[75] = PFILE(150) * 11.574074074074097;
    //q_T0_P_out, mw2ce112a1_b63b_4481_8441_e112009fa007, index: 76
    //Unit: second^(-1)
    _class_parameter[76] = PFILE(151) * 1.15740740740741e-05;
    //q_T0_T_in, mw251082a5_1bd9_4d64_84a0_ab1f26308116, index: 77
    //Unit: metre^(-3)second^(-1)
    _class_parameter[77] = PFILE(152) * 11.574074074074097;
    //Q_nT0_thym, mw2cdf9c70_21b9_465e_9cf5_64630fe0a6d9, index: 78
    //Unit: mole^(1)second^(-1)
    _class_parameter[78] = PFILE(153) * 1.92191982409137e-29;
    //k_nT0_pro, mwdadf1ca9_cec0_4275_8a5e_e54503c4e59d, index: 79
    //Unit: mole^(1)second^(-1)
    _class_parameter[79] = PFILE(154) * 1.92191982409137e-29;
    //K_nT0_pro, mw922335fa_257c_42eb_850b_e60807de5247, index: 80
    //Unit: mole^(1)
    _class_parameter[80] = PFILE(155) * 1.66053872801495e-24;
    //k_nT0_death, mwb751b880_cb0b_47ea_b498_66dcea9b718e, index: 81
    //Unit: second^(-1)
    _class_parameter[81] = PFILE(156) * 1.15740740740741e-05;
    //k_APC_mat, mw3e1ad180_8c17_4df8_9d07_8c2c6ec6bd2e, index: 82
    //Unit: second^(-1)
    _class_parameter[82] = PFILE(159) * 1.15740740740741e-05;
    //k_APC_mig, mw8160e6d7_ed43_41e0_8f4f_da89d99cd8ad, index: 83
    //Unit: second^(-1)
    _class_parameter[83] = PFILE(160) * 1.15740740740741e-05;
    //k_APC_death, mw8cc3dfe6_0da4_49a4_9386_9247b0d7f2fa, index: 84
    //Unit: second^(-1)
    _class_parameter[84] = PFILE(161) * 1.15740740740741e-05;
    //k_mAPC_death, mw980532ae_7dd3_4941_a8dc_951143298352, index: 85
    //Unit: second^(-1)
    _class_parameter[85] = PFILE(162) * 1.15740740740741e-05;
    //APC0_T, mw3b3f2927_2894_4459_9921_06e8b8b5bc19, index: 86
    //Unit: metre^(-3)mole^(1)
    _class_parameter[86] = PFILE(163) * 1.6605387280149534e-18;
    //APC0_LN, mwb6024d40_bb08_491a_83a2_de5ede07e099, index: 87
    //Unit: metre^(-3)mole^(1)
    _class_parameter[87] = PFILE(164) * 1.6605387280149534e-18;
    //k_c, mw41ca0542_ca03_4ef0_b0e0_7f01aa6fff9f, index: 88
    //Unit: second^(-1)
    _class_parameter[88] = PFILE(165) * 1.15740740740741e-05;
    //c0, mw072bc42a_40bf_45d8_824b_c830cd5ae75b, index: 89
    //Unit: metre^(-3)mole^(1)
    _class_parameter[89] = PFILE(166) * 999.9999999999994;
    //c50, mw7ed927f9_4afa_4423_9ef4_8afdee4bf735, index: 90
    //Unit: metre^(-3)mole^(1)
    _class_parameter[90] = PFILE(167) * 999.9999999999994;
    //DAMPs, mw4592e420_c091_49bb_b208_bd2e38aaf187, index: 91
    //Unit: dimensionless^(1)
    _class_parameter[91] = PFILE(168) * 6.02214199e+23;
    //n_sites_APC, mw8c8266db_44ce_4a86_bb81_93a1e050f9cd, index: 92
    //Unit: dimensionless^(1)
    _class_parameter[92] = PFILE(169) * 1.0;
    //kin, mwb0b80d60_2c25_4829_8c94_878da1be7cf1, index: 93
    //Unit: second^(-1)
    _class_parameter[93] = PFILE(170) * 1.15740740740741e-05;
    //kout, mw642453c1_cf5a_4882_9a07_f1092e4b4cff, index: 94
    //Unit: second^(-1)
    _class_parameter[94] = PFILE(171) * 1.15740740740741e-05;
    //k_P0_up, mw0ce4d373_6c72_4222_8c7f_201108455ebf, index: 95
    //Unit: mole^(-1)second^(-1)
    _class_parameter[95] = PFILE(172) * 6.97007174768519e+18;
    //k_xP0_deg, mw6c549174_67f0_495a_bdb4_ff02901c5470, index: 96
    //Unit: second^(-1)
    _class_parameter[96] = PFILE(173) * 1.15740740740741e-05;
    //k_P0_deg, mw85973cc5_cb66_46c3_896c_9c1afcf4947f, index: 97
    //Unit: second^(-1)
    _class_parameter[97] = PFILE(174) * 1.15740740740741e-05;
    //k_p0_deg, mw184dc020_4a7d_4b9f_b82b_3b97f7f2162f, index: 98
    //Unit: second^(-1)
    _class_parameter[98] = PFILE(175) * 1.15740740740741e-05;
    //k_P0_on, mwb8990525_887d_42ac_b46b_f0c39833210d, index: 99
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[99] = PFILE(176) * 1.1574074074074112e-08;
    //k_P0_d1, mw8d8dd086_b66d_40a1_8a9c_66b122846a4a, index: 100
    //Unit: metre^(-3)mole^(1)
    _class_parameter[100] = PFILE(177) * 999.9999999999994;
    //p0_50, mw40e8d966_ce69_4695_b18f_32e99ed14051, index: 101
    //Unit: metre^(-2)mole^(1)
    _class_parameter[101] = PFILE(178) * 1.66053872801495e-12;
    //P0_C1, mw42241e84_2461_451c_ba85_7a82e38005cd, index: 102
    //Unit: metre^(-3)
    _class_parameter[102] = PFILE(179) * 6.022141989999979e+26;
    //A_syn, mw1dec1cd7_7be1_47ed_8120_f40cd58469f8, index: 103
    //Unit: metre^(2)
    _class_parameter[103] = PFILE(180) * 1e-12;
    //A_Tcell, mw16847761_9bad_4f43_9ff1_c59a128f0906, index: 104
    //Unit: metre^(2)
    _class_parameter[104] = PFILE(181) * 1e-12;
    //A_cell, mwc02f4480_110d_4b28_b9d9_e637906e49cb, index: 105
    //Unit: metre^(2)
    _class_parameter[105] = PFILE(182) * 1e-12;
    //k_M1p0_TCR_on, mw35d7f1c9_231a_4674_92cc_5248430db78b, index: 106
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[106] = PFILE(183) * 602214199000.0;
    //k_M1p0_TCR_off, mwebd45a34_96eb_4fce_94a1_88c4a786499e, index: 107
    //Unit: second^(-1)
    _class_parameter[107] = PFILE(184) * 1.0;
    //TCR_p0_tot, mwce0ec5dc_cf04_49e7_a145_99f037bee232, index: 108
    //Unit: metre^(-2)mole^(1)
    _class_parameter[108] = PFILE(185) * 1.66053872801495e-12;
    //k_P1_up, mwce8e5a59_9118_4088_b2da_55a783e979a4, index: 109
    //Unit: mole^(-1)second^(-1)
    _class_parameter[109] = PFILE(187) * 6.97007174768519e+18;
    //k_xP1_deg, mw88856bf1_0fc4_4408_b136_c8dcf986937b, index: 110
    //Unit: second^(-1)
    _class_parameter[110] = PFILE(188) * 1.15740740740741e-05;
    //k_P1_deg, mw930d56a4_fd83_4e3b_a3ce_b2dfb8eb2803, index: 111
    //Unit: second^(-1)
    _class_parameter[111] = PFILE(189) * 1.15740740740741e-05;
    //k_p1_deg, mw1030dd2e_e3f3_4068_a3a4_5cbccf77efc1, index: 112
    //Unit: second^(-1)
    _class_parameter[112] = PFILE(190) * 1.15740740740741e-05;
    //k_P1_on, mwd1ee8636_a779_400b_bdc9_83f6fd6f9657, index: 113
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[113] = PFILE(191) * 1.1574074074074112e-08;
    //k_P1_d1, mw75a8d04d_653f_4fb9_943c_007065c667ac, index: 114
    //Unit: metre^(-3)mole^(1)
    _class_parameter[114] = PFILE(192) * 999.9999999999994;
    //p1_50, mw692f76fd_d22f_4f38_9c2b_cbbeec3a84b7, index: 115
    //Unit: metre^(-2)mole^(1)
    _class_parameter[115] = PFILE(193) * 1.66053872801495e-12;
    //P1_C1, mw712a01e3_4713_463b_b711_24cc418a17dd, index: 116
    //Unit: metre^(-3)
    _class_parameter[116] = PFILE(194) * 6.022141989999979e+26;
    //k_M1p1_TCR_on, mw505e0384_75a3_407f_8968_3efb54b799f3, index: 117
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[117] = PFILE(195) * 602214199000.0;
    //k_M1p1_TCR_off, mw5ce1633b_566b_4973_a984_a3e2056317aa, index: 118
    //Unit: second^(-1)
    _class_parameter[118] = PFILE(196) * 1.0;
    //k_M1p1_TCR_p, mw9399840e_9583_4188_9e03_b5322ab1bd5a, index: 119
    //Unit: second^(-1)
    _class_parameter[119] = PFILE(197) * 1.0;
    //phi_M1p1_TCR, mwefd10fe0_8b86_4cc1_ac65_df6d1504aa48, index: 120
    //Unit: second^(-1)
    _class_parameter[120] = PFILE(198) * 1.0;
    //N_M1p1_TCR, mw85bfe2bb_117c_4e79_998f_a55e66b01239, index: 121
    //Unit: dimensionless^(1)
    _class_parameter[121] = PFILE(199) * 1.0;
    //TCR_p1_tot, mw88c46664_2c0f_46da_9656_4eaed3f5c2cb, index: 122
    //Unit: metre^(-2)mole^(1)
    _class_parameter[122] = PFILE(200) * 1.66053872801495e-12;
    //q_P_nivolumab, mwabfc765b_667d_42b9_b1da_c000d2ed0988, index: 123
    //Unit: second^(-1)
    _class_parameter[123] = PFILE(202) * 1.0;
    //q_T_nivolumab, mw3f04dc7f_889a_422d_87ee_f6398df72fa9, index: 124
    //Unit: second^(-1)
    _class_parameter[124] = PFILE(203) * 1.0;
    //q_LN_nivolumab, mw805fd08e_f69b_4091_8a1a_945258e1c734, index: 125
    //Unit: second^(-1)
    _class_parameter[125] = PFILE(204) * 1.0;
    //q_LD_nivolumab, mw33e7929c_d964_4b08_a528_be6750fd1124, index: 126
    //Unit: second^(-1)
    _class_parameter[126] = PFILE(205) * 0.0166666666666667;
    //k_cl_nivolumab, mw8a7126ac_26db_44b0_951c_7900e9ca4a95, index: 127
    //Unit: second^(-1)
    _class_parameter[127] = PFILE(206) * 1.0;
    //gamma_C_nivolumab, mw6beecf82_9e28_450c_9bb5_d72a4f09cca2, index: 128
    //Unit: dimensionless^(1)
    _class_parameter[128] = PFILE(207) * 1.0;
    //gamma_P_nivolumab, mw1f5933e4_7b6b_4cdb_8185_7cf6a57ebc1d, index: 129
    //Unit: dimensionless^(1)
    _class_parameter[129] = PFILE(208) * 1.0;
    //gamma_T_nivolumab, mw05c69e80_7ac3_41a2_a526_4a151df7558a, index: 130
    //Unit: dimensionless^(1)
    _class_parameter[130] = PFILE(209) * 1.0;
    //gamma_LN_nivolumab, mw2560fed7_0d9d_4e14_b48d_ea32c4879030, index: 131
    //Unit: dimensionless^(1)
    _class_parameter[131] = PFILE(210) * 1.0;
    //T_PD1_total, mw3b38b76d_e81a_4a73_90f9_cfaa20268b11, index: 132
    //Unit: metre^(-2)mole^(1)
    _class_parameter[132] = PFILE(211) * 1.66053872801495e-12;
    //C_PDL1_total, mw75d601f7_ce54_44d8_8582_35fb83fb4da8, index: 133
    //Unit: metre^(-2)mole^(1)
    _class_parameter[133] = PFILE(212) * 1.66053872801495e-12;
    //C_PDL2_total, mw19c11cb2_327f_49f9_9659_0f5ce2e775d3, index: 134
    //Unit: metre^(-2)mole^(1)
    _class_parameter[134] = PFILE(213) * 1.66053872801495e-12;
    //A_CD80_total, mw0c872ea0_97b3_4ebd_b0a8_c4980a70ac9b, index: 135
    //Unit: metre^(-2)mole^(1)
    _class_parameter[135] = PFILE(214) * 1.66053872801495e-12;
    //A_CD86_total, mw5fb44543_1fd2_4ceb_9708_00d5f6494c0d, index: 136
    //Unit: metre^(-2)mole^(1)
    _class_parameter[136] = PFILE(215) * 1.66053872801495e-12;
    //N_CD28_total, mw93d4170c_3627_4a5d_8757_d71ccfdc2e8f, index: 137
    //Unit: metre^(-2)mole^(1)
    _class_parameter[137] = PFILE(216) * 1.66053872801495e-12;
    //N_PDL1_total, mw4809468e_6300_443b_a69b_9d8340830877, index: 138
    //Unit: metre^(-2)mole^(1)
    _class_parameter[138] = PFILE(217) * 1.66053872801495e-12;
    //N_CTLA4_total, mwc46af5da_bba0_432e_9908_e5ee38923c4d, index: 139
    //Unit: metre^(-2)mole^(1)
    _class_parameter[139] = PFILE(218) * 1.66053872801495e-12;
    //Treg_CTLA4_total, mw083eee70_8eb9_46ac_9b18_7002a6e5f323, index: 140
    //Unit: metre^(-2)mole^(1)
    _class_parameter[140] = PFILE(219) * 1.66053872801495e-12;
    //kon_PD1_PDL1, mwde9f074e_f7d0_4402_b9b2_fa3ffeec0622, index: 141
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[141] = PFILE(220) * 602214199000.0;
    //koff_PD1_PDL1, mwe1de3847_7d74_4d95_acda_9ee8d90a8ab4, index: 142
    //Unit: second^(-1)
    _class_parameter[142] = PFILE(221) * 1.0;
    //kon_PD1_PDL2, mwe419e5d7_011c_4f1c_ab81_98e4f81187a4, index: 143
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[143] = PFILE(222) * 602214199000.0;
    //koff_PD1_PDL2, mwc6d4ae27_3cbf_4cb5_a055_600fb2c5055d, index: 144
    //Unit: second^(-1)
    _class_parameter[144] = PFILE(223) * 1.0;
    //kon_CD28_CD80, mw00a5f90a_286a_4489_b93f_0a6a4a8cdee0, index: 145
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[145] = PFILE(224) * 602214199000.0;
    //koff_CD28_CD80, mw8460b8ef_e478_43f0_b914_dd0090d60c89, index: 146
    //Unit: second^(-1)
    _class_parameter[146] = PFILE(225) * 1.0;
    //kon_CD28_CD86, mw84512afc_20dd_4376_b96a_adaa3561d0ae, index: 147
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[147] = PFILE(226) * 602214199000.0;
    //koff_CD28_CD86, mw58b3980d_406c_49d1_b574_29096238e516, index: 148
    //Unit: second^(-1)
    _class_parameter[148] = PFILE(227) * 1.0;
    //kon_CTLA4_CD80, mwe02bbd54_35a8_458d_8189_c0e17f30c545, index: 149
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[149] = PFILE(228) * 602214199000.0;
    //koff_CTLA4_CD80, mw35381a46_d268_4d5d_867a_a5824dc41f0c, index: 150
    //Unit: second^(-1)
    _class_parameter[150] = PFILE(229) * 1.0;
    //kon_CTLA4_CD86, mw1a117675_55f4_4485_aa0f_5482d4bd6d6d, index: 151
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[151] = PFILE(230) * 602214199000.0;
    //koff_CTLA4_CD86, mw1d5c2a09_2906_455d_a942_1103ecba0101, index: 152
    //Unit: second^(-1)
    _class_parameter[152] = PFILE(231) * 1.0;
    //kon_PDL1_CD80, mw2d14a91a_59f0_4b44_a784_298156b9a399, index: 153
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[153] = PFILE(232) * 602214199000.0;
    //koff_PDL1_CD80, mwcf124048_87c8_46f7_bfe5_3c4635fd5b04, index: 154
    //Unit: second^(-1)
    _class_parameter[154] = PFILE(233) * 1.0;
    //kon_PD1_aPD1, mw48776c88_8886_46c0_b622_0be20e388672, index: 155
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[155] = PFILE(234) * 0.0010000000000000007;
    //koff_PD1_aPD1, mw1d11c17e_f8aa_4f05_be97_5bbea8a063e5, index: 156
    //Unit: second^(-1)
    _class_parameter[156] = PFILE(235) * 1.0;
    //Chi_PD1, mw304ffd24_9644_4054_af99_6ed2fcca3ce3, index: 157
    //Unit: metre^(-1)
    _class_parameter[157] = PFILE(236) * 999999999.9999999;
    //PD1_50, mw5b66577f_0bad_47e9_b436_ee165dd02194, index: 158
    //Unit: metre^(-2)mole^(1)
    _class_parameter[158] = PFILE(237) * 1.66053872801495e-12;
    //n_PD1, mwba153ce8_52a1_4c61_bafb_a89ee90a2e5e, index: 159
    //Unit: dimensionless^(1)
    _class_parameter[159] = PFILE(238) * 1.0;
    //CD28_50, mwb0be24a6_0dca_4763_a3ee_79b8b187b980, index: 160
    //Unit: metre^(-2)mole^(1)
    _class_parameter[160] = PFILE(239) * 1.66053872801495e-12;
    //n_CD28, mw8aaaad35_99f0_4eda_80fe_1dfef5ab2da9, index: 161
    //Unit: dimensionless^(1)
    _class_parameter[161] = PFILE(240) * 1.0;
    //N_CD28, mw837a6ff4_aa60_4715_aff1_7158308681bf, index: 162
    //Unit: dimensionless^(1)
    _class_parameter[162] = PFILE(241) * 1.0;
    //CTLA4_50, mwafdcfee7_38c7_41f8_8b9f_5fc83677c5fb, index: 163
    //Unit: metre^(-2)mole^(1)
    _class_parameter[163] = PFILE(242) * 1.66053872801495e-12;
    //n_ADCC, mwff6854f4_bc7a_4655_9a4f_dca17abf5b1d, index: 164
    //Unit: dimensionless^(1)
    _class_parameter[164] = PFILE(243) * 1.0;
    //k_ADCC, mweb73df96_2949_4662_bfb7_fa58e6f7ffe8, index: 165
    //Unit: second^(-1)
    _class_parameter[165] = PFILE(244) * 1.15740740740741e-05;
    //H_ADCC_P, mw1e482d45_5b24_4fda_8c34_bb2ce3f7aa1e, index: 166
    //Unit: dimensionless^(1)
    _class_parameter[166] = PFILE(245) * 1.0;
    //H_ADCC_T, mwc461a835_6702_488f_bc85_0b3821469ea3, index: 167
    //Unit: dimensionless^(1)
    _class_parameter[167] = PFILE(246) * 1.0;
    //k_a1_cabozantinib, mw59ed841c_001e_4cb2_9cbe_db4399a53688, index: 168
    //Unit: second^(-1)
    _class_parameter[168] = PFILE(247) * 0.000277777777777778;
    //k_a2_cabozantinib, mwf5cca1a1_2b31_481c_be71_1b798355cc9b, index: 169
    //Unit: second^(-1)
    _class_parameter[169] = PFILE(248) * 0.000277777777777778;
    //k_cln_cabozantinib, mw677116c7_d653_4e89_b698_2f7335fb2bb8, index: 170
    //Unit: metre^(-3)mole^(1)second^(-1)
    _class_parameter[170] = PFILE(249) * 0.277777777777778;
    //Kc_cabozantinib, mw3cd7b56e_fb10_4cbc_a85f_b092cfd11c58, index: 171
    //Unit: metre^(-3)mole^(1)
    _class_parameter[171] = PFILE(250) * 999.9999999999994;
    //lagP1_cabozantinib, mw5db48776_e402_4706_9320_37cb87ed95ef, index: 172
    //Unit: second^(1)
    _class_parameter[172] = PFILE(251) * 3600.0;
    //lagP2_cabozantinib, mw7d019554_1330_4001_b5e7_e7c1a3e0cb30, index: 173
    //Unit: second^(1)
    _class_parameter[173] = PFILE(252) * 3600.0;
    //F_cabozantinib, mw3b57c773_152b_4de7_9f3c_5fa27a42fe69, index: 174
    //Unit: dimensionless^(1)
    _class_parameter[174] = PFILE(253) * 1.0;
    //q_P_cabozantinib, mwdd72852e_cdd0_440e_83c2_7ce3af365a72, index: 175
    //Unit: second^(-1)
    _class_parameter[175] = PFILE(254) * 1.0;
    //q_T_cabozantinib, mwda278504_bf72_4dd5_af58_e991a990196b, index: 176
    //Unit: second^(-1)
    _class_parameter[176] = PFILE(255) * 1.0;
    //q_LN_cabozantinib, mw45fc2fac_ca3d_4277_bbf1_0c3112a1ccbb, index: 177
    //Unit: second^(-1)
    _class_parameter[177] = PFILE(256) * 1.0;
    //q_LD_cabozantinib, mw4e275c94_b9a6_403c_8fa3_1f84966252fc, index: 178
    //Unit: second^(-1)
    _class_parameter[178] = PFILE(257) * 0.0166666666666667;
    //gamma_C_cabozantinib, mwf29cc72b_370a_42d6_b305_a1f9e4fd1767, index: 179
    //Unit: dimensionless^(1)
    _class_parameter[179] = PFILE(258) * 1.0;
    //gamma_P_cabozantinib, mw55bc45c4_a014_4e21_b61a_8decab0a6507, index: 180
    //Unit: dimensionless^(1)
    _class_parameter[180] = PFILE(259) * 1.0;
    //gamma_T_cabozantinib, mwd441f95a_55f3_499c_860c_e276f418755c, index: 181
    //Unit: dimensionless^(1)
    _class_parameter[181] = PFILE(260) * 1.0;
    //gamma_LN_cabozantinib, mwa12ed763_a42c_4014_901c_966d051ecaf1, index: 182
    //Unit: dimensionless^(1)
    _class_parameter[182] = PFILE(261) * 1.0;
    //IC50_MET, mw514262c5_8b3d_4d42_835b_83e498a8875e, index: 183
    //Unit: metre^(-3)mole^(1)
    _class_parameter[183] = PFILE(262) * 1.0000000000000008e-06;
    //IC50_RET, mwf253a32d_f5fc_4e4a_b933_6ba92dd7b3eb, index: 184
    //Unit: metre^(-3)mole^(1)
    _class_parameter[184] = PFILE(263) * 1.0000000000000008e-06;
    //IC50_AXL, mwb9bb52b9_e309_49f6_aaca_6a64d8e80ed2, index: 185
    //Unit: metre^(-3)mole^(1)
    _class_parameter[185] = PFILE(264) * 1.0000000000000008e-06;
    //IC50_VEGFR2, mw29f5637c_8e11_4b57_b5cf_58bd76280f1d, index: 186
    //Unit: metre^(-3)mole^(1)
    _class_parameter[186] = PFILE(265) * 1.0000000000000008e-06;
    //cabo0, mw9ea32ef3_e919_434d_a61c_e4c9a54efe6f, index: 187
    //Unit: metre^(-3)mole^(1)
    _class_parameter[187] = PFILE(266) * 1.0000000000000008e-06;
    //kv_cabo, mw11cb8be5_a237_4212_aa45_1be1c7665188, index: 188
    //Unit: second^(-1)
    _class_parameter[188] = PFILE(267) * 1.15740740740741e-05;
    //kr_cabo, mwd845598f_f8b9_4583_baa0_b00e9ecda42f, index: 189
    //Unit: second^(-1)
    _class_parameter[189] = PFILE(268) * 1.15740740740741e-05;
    //lambda_C_cabo, mw400bcfd9_01c8_49ab_b0de_00e7c46728de, index: 190
    //Unit: dimensionless^(1)
    _class_parameter[190] = PFILE(269) * 1.0;
    //lambda_v_cabo, mw2cdf97ee_0771_4222_b4cf_08f94482531f, index: 191
    //Unit: second^(-1)
    _class_parameter[191] = PFILE(270) * 1.15740740740741e-05;
    //lambda_q_cabo, mwd938fbcd_624c_49ac_aff9_7311ec846e66, index: 192
    //Unit: dimensionless^(1)
    _class_parameter[192] = PFILE(271) * 1.0;
}

void ODE_system::setupVariables(void){

    _species_var = std::vector<realtype>(52, 0);
    _nonspecies_var = std::vector<realtype>(0, 0);
    //species not part of ode left-hand side
    _species_other =  std::vector<realtype>(20, 0);
    
    return;
}


void ODE_system::setup_instance_varaibles(Param& param){

    //V_C.nT1, mw9bb199c9_6324_46d7_a8fc_403834770ab1, index: 0
    //Unit: mole^(1)
    _species_var[0] = PFILE(16) * 1.66053872801495e-24;
    //V_C.T1, mw72ce8928_e928_45be_ba7f_50adaff6ec59, index: 1
    //Unit: mole^(1)
    _species_var[1] = PFILE(17) * 1.66053872801495e-24;
    //V_C.nT0, mw81afe5ba_b937_4d30_8ba3_68471e2d3f8e, index: 2
    //Unit: mole^(1)
    _species_var[2] = PFILE(18) * 1.66053872801495e-24;
    //V_C.T0, mw50086af5_e8e7_4a6e_a824_cfee76b4e3d3, index: 3
    //Unit: mole^(1)
    _species_var[3] = PFILE(19) * 1.66053872801495e-24;
    //V_C.nivolumab, mwd9d24e12_92e5_4231_9769_4b0acca89fa9, index: 4
    //Unit: mole^(1)metre^(-3)
    _species_var[4] = PFILE(20) * 1000.0;
    //V_C.cabozantinib, mwb57fedeb_cae2_43b4_9fd2_81f942925f49, index: 5
    //Unit: mole^(1)metre^(-3)
    _species_var[5] = PFILE(21) * 1000.0;
    //V_C.A_site1, mw9c8f5fb8_feb0_4d19_a408_576259d4144d, index: 6
    //Unit: mole^(1)metre^(-3)
    _species_var[6] = PFILE(22) * 1000.0;
    //V_C.A_site2, mw8c121cbe_831d_43eb_bbfc_37be87655074, index: 7
    //Unit: mole^(1)metre^(-3)
    _species_var[7] = PFILE(23) * 1000.0;
    //V_P.nT1, mwc5399f17_f08d_4fc6_af92_5925c7ad933c, index: 8
    //Unit: mole^(1)
    _species_var[8] = PFILE(24) * 1.66053872801495e-24;
    //V_P.T1, mwb0a70bc8_5c10_473b_a674_b21d1c3d0dd4, index: 9
    //Unit: mole^(1)
    _species_var[9] = PFILE(25) * 1.66053872801495e-24;
    //V_P.nT0, mwb859fce6_5efa_48be_b412_431bde044d3c, index: 10
    //Unit: mole^(1)
    _species_var[10] = PFILE(26) * 1.66053872801495e-24;
    //V_P.T0, mw347e60d1_f520_47fb_b09e_6937d7783599, index: 11
    //Unit: mole^(1)
    _species_var[11] = PFILE(27) * 1.66053872801495e-24;
    //V_P.nivolumab, mwe4917312_e3ba_4a69_a3ec_7c3e9aa5c172, index: 12
    //Unit: mole^(1)metre^(-3)
    _species_var[12] = PFILE(28) * 1000.0;
    //V_P.cabozantinib, mw5cadc79f_ae57_42a5_a88f_58806fe70cec, index: 13
    //Unit: mole^(1)metre^(-3)
    _species_var[13] = PFILE(29) * 1000.0;
    //V_T.C_x, mwd685144b_83b2_498b_8547_abf3d9428d2f, index: 14
    //Unit: mole^(1)
    _species_var[14] = PFILE(30) * 1.66053872801495e-24;
    //V_T.T_exh, mw2cc6bf89_0c93_4d3f_af0c_d4f2d928252c, index: 15
    //Unit: mole^(1)
    _species_var[15] = PFILE(31) * 1.66053872801495e-24;
    //V_T.C1, mw1d60b639_c7bd_450c_9e54_72821ff3349f, index: 16
    //Unit: mole^(1)
    _species_var[16] = PFILE(32) * 1.66053872801495e-24;
    //V_T.T1, mw906c79fd_b986_4ff0_af05_91a39eae8d0f, index: 17
    //Unit: mole^(1)
    _species_var[17] = PFILE(33) * 1.66053872801495e-24;
    //V_T.T0, mw87896245_9b3b_433a_a82a_c8eddee1c43e, index: 18
    //Unit: mole^(1)
    _species_var[18] = PFILE(34) * 1.66053872801495e-24;
    //V_T.APC, mw649182e7_be1a_445d_aef7_d001f409f451, index: 19
    //Unit: mole^(1)
    _species_var[19] = PFILE(35) * 1.66053872801495e-24;
    //V_T.mAPC, mwa8a7ab07_f1b7_4fe5_8fc8_17b2b92db889, index: 20
    //Unit: mole^(1)
    _species_var[20] = PFILE(36) * 1.66053872801495e-24;
    //V_T.c, mw8a621476_141c_42b5_99e0_25045997f46f, index: 21
    //Unit: mole^(1)metre^(-3)
    _species_var[21] = PFILE(37) * 999.9999999999999;
    //V_T.P0, mw97808c17_de36_4616_bb30_be5d2f5e3e1b, index: 22
    //Unit: mole^(1)metre^(-3)
    _species_var[22] = PFILE(38) * 999.9999999999999;
    //V_T.P1, mw911b473b_aac9_4e08_8d2b_01b2c5c00ae8, index: 23
    //Unit: mole^(1)metre^(-3)
    _species_var[23] = PFILE(39) * 999.9999999999999;
    //V_T.nivolumab, mw96a94d95_c06c_404f_894d_29fe61426273, index: 24
    //Unit: mole^(1)metre^(-3)
    _species_var[24] = PFILE(40) * 999.9999999999999;
    //V_T.cabozantinib, mw9c857630_844c_4b31_aa09_6d0ccfcaf650, index: 25
    //Unit: mole^(1)metre^(-3)
    _species_var[25] = PFILE(41) * 999.9999999999999;
    //V_T.v_cabo, mw646e551c_520d_4f81_8af4_78f1ca4c6032, index: 26
    //Unit: mole^(1)
    _species_var[26] = PFILE(42) * 1.66053872801495e-24;
    //V_LN.nT1, mw673efd58_3f22_49bc_b0e5_19885fcb8a11, index: 27
    //Unit: mole^(1)
    _species_var[27] = PFILE(43) * 1.66053872801495e-24;
    //V_LN.aT1, mw27fef97e_9fae_4550_8d4c_85ac90871439, index: 28
    //Unit: mole^(1)
    _species_var[28] = PFILE(44) * 1.66053872801495e-24;
    //V_LN.T1, mwbfaf7da3_3c95_4a1e_90ae_8177b328085e, index: 29
    //Unit: mole^(1)
    _species_var[29] = PFILE(45) * 1.66053872801495e-24;
    //V_LN.IL2, mw9e25f7c6_5c55_48af_b02d_540c574f5f04, index: 30
    //Unit: mole^(1)metre^(-3)
    _species_var[30] = PFILE(46) * 1000.0;
    //V_LN.nT0, mwc6a200b5_6f10_4357_b46b_fe65f94fee6a, index: 31
    //Unit: mole^(1)
    _species_var[31] = PFILE(47) * 1.66053872801495e-24;
    //V_LN.aT0, mw12bd819b_3ad5_4636_a0f4_4bec446de5e8, index: 32
    //Unit: mole^(1)
    _species_var[32] = PFILE(48) * 1.66053872801495e-24;
    //V_LN.T0, mw0f4c2873_69b4_4796_8c83_426ba8921f44, index: 33
    //Unit: mole^(1)
    _species_var[33] = PFILE(49) * 1.66053872801495e-24;
    //V_LN.APC, mw3080e883_efd4_4eaf_a3e4_0218479f2a4e, index: 34
    //Unit: mole^(1)
    _species_var[34] = PFILE(50) * 1.66053872801495e-24;
    //V_LN.mAPC, mwaf315482_94f3_4afc_9f22_b8a401bdfd15, index: 35
    //Unit: mole^(1)
    _species_var[35] = PFILE(51) * 1.66053872801495e-24;
    //V_LN.nivolumab, mwf9e0907b_3f65_4439_8238_304a0bedbb90, index: 36
    //Unit: mole^(1)metre^(-3)
    _species_var[36] = PFILE(52) * 1000.0;
    //V_LN.cabozantinib, mw0328e128_134e_433c_a781_6acf8a722422, index: 37
    //Unit: mole^(1)metre^(-3)
    _species_var[37] = PFILE(53) * 1000.0;
    //V_e.P0, mw15a4c371_a16a_49be_bfd7_7fc15cd186d3, index: 38
    //Unit: mole^(1)metre^(-3)
    _species_var[38] = PFILE(54) * 1000.0;
    //V_e.p0, mwde092df0_141c_4a9e_8dff_3cb80a2164db, index: 39
    //Unit: mole^(1)metre^(-3)
    _species_var[39] = PFILE(55) * 1000.0;
    //V_e.P1, mwf4dcc1ea_91dc_4bcd_9277_cb9a241643bb, index: 40
    //Unit: mole^(1)metre^(-3)
    _species_var[40] = PFILE(56) * 1000.0;
    //V_e.p1, mw92c103d6_76f4_41c0_9523_6cf1802d303b, index: 41
    //Unit: mole^(1)metre^(-3)
    _species_var[41] = PFILE(57) * 1000.0;
    //A_e.M1, mw6464f082_c5ad_4f59_b863_94ed2733f985, index: 42
    //Unit: mole^(1)metre^(-2)
    _species_var[42] = PFILE(58) * 1.66053872801495e-12;
    //A_e.M1p0, mw5ade50c9_9e3c_4c35_894b_9343d07234d3, index: 43
    //Unit: mole^(1)metre^(-2)
    _species_var[43] = PFILE(59) * 1.66053872801495e-12;
    //A_e.M1p1, mw1573bb71_02cb_4b5d_b0a3_54aae81df811, index: 44
    //Unit: mole^(1)metre^(-2)
    _species_var[44] = PFILE(60) * 1.66053872801495e-12;
    //A_s.M1, mw0fbde24a_5ddd_4fe0_bcbc_12764ad9595a, index: 45
    //Unit: mole^(1)metre^(-2)
    _species_var[45] = PFILE(61) * 1.66053872801495e-12;
    //A_s.M1p0, mw22bf7924_3742_4750_9820_c5a531c12b09, index: 46
    //Unit: mole^(1)metre^(-2)
    _species_var[46] = PFILE(62) * 1.66053872801495e-12;
    //A_s.M1p1, mw51e8e4a6_ae35_470b_9911_19da7ecaba55, index: 47
    //Unit: mole^(1)metre^(-2)
    _species_var[47] = PFILE(63) * 1.66053872801495e-12;
    //syn_CT.PD1_PDL1, mw8d946c5b_232e_47b8_b378_9cdeeb14fc60, index: 48
    //Unit: mole^(1)metre^(-2)
    _species_var[48] = PFILE(64) * 1.66053872801495e-12;
    //syn_CT.PD1_PDL2, mwc1c073ab_8af6_466d_9e5b_ce9c2b094249, index: 49
    //Unit: mole^(1)metre^(-2)
    _species_var[49] = PFILE(65) * 1.66053872801495e-12;
    //syn_CT.PD1_aPD1, mwd97c43c4_cbb0_43dc_b79f_516df8d61a67, index: 50
    //Unit: mole^(1)metre^(-2)
    _species_var[50] = PFILE(66) * 1.66053872801495e-12;
    //syn_CT.PD1_aPD1_PD1, mw014081f7_8bd8_4f18_aef2_c33344e79658, index: 51
    //Unit: mole^(1)metre^(-2)
    _species_var[51] = PFILE(67) * 1.66053872801495e-12;
    
    return;
}
    


void ODE_system::setup_instance_tolerance(Param& param){

    //Tolerance
    realtype reltol = PFILE(3);
    realtype abstol_base = PFILE(4);
    N_Vector abstol = N_VNew_Serial(_neq);

    for (size_t i = 0; i < 52; i++)
    {
        NV_DATA_S(abstol)[i] = abstol_base * get_unit_conversion_species(i);
    }
    int flag = CVodeSVtolerances(_cvode_mem, reltol, abstol);
    check_flag(&flag, "CVodeSVtolerances", 1);

    
    return;
}

void ODE_system::eval_init_assignment(void){
    //Assignment Rules required before IA
    //InitialAssignment
    _class_parameter[0] = _class_parameter[38];
    _class_parameter[1] = _class_parameter[39];
    _class_parameter[2] = _class_parameter[40];

    updateVar();
    
    return;
}
void ODE_system::setupEvents(void){

    _nevent = 7;
    _nroot = 7;

    _trigger_element_type = std::vector<EVENT_TRIGGER_ELEM_TYPE>(_nroot, TRIGGER_NON_INSTANT);
    _trigger_element_satisfied = std::vector<bool>(_nroot, false);
    _event_triggered = std::vector<bool>(_nevent, false);

    //V_T.C1 < (0.9 * cell)
    _trigger_element_type[0] = TRIGGER_NON_INSTANT;

    //V_T.T1 < (1e-10 * cell)
    _trigger_element_type[1] = TRIGGER_NON_INSTANT;

    //V_T.c < (1e-10 * molarity)
    _trigger_element_type[2] = TRIGGER_NON_INSTANT;

    //V_C.nivolumab < (1e-10 * molarity)
    _trigger_element_type[3] = TRIGGER_NON_INSTANT;

    //V_P.nivolumab < (1e-10 * molarity)
    _trigger_element_type[4] = TRIGGER_NON_INSTANT;

    //V_T.nivolumab < (1e-10 * molarity)
    _trigger_element_type[5] = TRIGGER_NON_INSTANT;

    //V_LN.nivolumab < (1e-10 * molarity)
    _trigger_element_type[6] = TRIGGER_NON_INSTANT;

    _event_triggered[0] = true;

    _event_triggered[1] = true;

    _event_triggered[2] = true;

    _event_triggered[3] = true;

    _event_triggered[4] = true;

    _event_triggered[5] = true;

    _event_triggered[6] = true;

    return;
}
int ODE_system::f(realtype t, N_Vector y, N_Vector ydot, void *user_data){

    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);

    //Assignment rules:

    realtype AUX_VAR_V_T = PARAM(37) + PARAM(35) * SPVAR(14) + PARAM(36) * SPVAR(15) + PARAM(35) * SPVAR(16) + PARAM(36) * SPVAR(17) + PARAM(36) * SPVAR(18);

    realtype AUX_VAR_C_total = 0.0 * PARAM(31) + SPVAR(16);

    realtype AUX_VAR_T_T = 0.0 * PARAM(31) + SPVAR(17);

    realtype AUX_VAR_H_PD1 = std::pow((SPVAR(48) + SPVAR(49)) / PARAM(158), PARAM(159)) / (std::pow((SPVAR(48) + SPVAR(49)) / PARAM(158), PARAM(159)) + 1.0);

    realtype AUX_VAR_T_C = 0.0 * PARAM(31) + SPVAR(1);

    realtype AUX_VAR_T_P = 0.0 * PARAM(31) + SPVAR(9);

    realtype AUX_VAR_T_LN = 0.0 * PARAM(31) + SPVAR(29);

    realtype AUX_VAR_mAPC_LN = SPVAR(35);

    realtype AUX_VAR_N_aT = PARAM(63) + PARAM(64) * SPVAR(30) / (PARAM(62) + SPVAR(30));

    realtype AUX_VAR_Tregs_C = SPVAR(3);

    realtype AUX_VAR_Tregs_P = SPVAR(11);

    realtype AUX_VAR_Tregs_T = SPVAR(18);

    realtype AUX_VAR_Tregs_LN = SPVAR(33);

    realtype AUX_VAR_APC_LN = SPVAR(34);

    realtype AUX_VAR_pTCR_p0_MHC_tot = 0.5 * (SPVAR(46) / PARAM(69) + PARAM(108) + PARAM(107) / PARAM(106) - PARAM(108) * std::pow(std::pow((SPVAR(46) / PARAM(69) + PARAM(108) + PARAM(107) / PARAM(106)) / PARAM(108), 2.0) - 4.0 * SPVAR(46) / PARAM(69) / PARAM(108), 0.5));

    realtype AUX_VAR_pTCR_p1_MHC_tot = PARAM(118) / (PARAM(118) + PARAM(120)) * std::pow(PARAM(119) / (PARAM(118) + PARAM(119)), PARAM(121)) * 0.5 * (SPVAR(47) / PARAM(46) + PARAM(122) + PARAM(118) / PARAM(117) - PARAM(122) * std::pow(std::pow((SPVAR(47) / PARAM(46) + PARAM(122) + PARAM(118) / PARAM(117)) / PARAM(122), 2.0) - 4.0 * SPVAR(47) / PARAM(46) / PARAM(122), 0.5));

    realtype AUX_VAR_R_cabo = SPVAR(25) / (SPVAR(25) + PARAM(185));

    realtype AUX_VAR_C_max = SPVAR(26);

    realtype AUX_VAR_frac = SPVAR(25) / (SPVAR(25) + PARAM(183));

    realtype AUX_VAR_R_Tcell = 0.0 * PARAM(31) / PARAM(32) + PARAM(66) * SPVAR(17) * SPVAR(16) / (AUX_VAR_C_total + AUX_VAR_T_T + PARAM(31)) * (1.0 - AUX_VAR_H_PD1);

    realtype AUX_VAR_H_act1 = PARAM(92) * AUX_VAR_mAPC_LN / (PARAM(92) * AUX_VAR_mAPC_LN + SPVAR(27) + PARAM(31));

    realtype AUX_VAR_H_act0 = PARAM(92) * AUX_VAR_APC_LN / (PARAM(92) * AUX_VAR_APC_LN + SPVAR(31) + PARAM(31));

    realtype AUX_VAR_H_P0 = AUX_VAR_pTCR_p0_MHC_tot / (AUX_VAR_pTCR_p0_MHC_tot + PARAM(101));

    realtype AUX_VAR_H_P1 = AUX_VAR_pTCR_p1_MHC_tot / (AUX_VAR_pTCR_p1_MHC_tot + PARAM(115));

    //Reaction fluxes:

    realtype ReactionFlux1 = PARAM(30) * SPVAR(14);

    realtype ReactionFlux2 = PARAM(30) * SPVAR(15);

    realtype ReactionFlux3 = PARAM(42) * SPVAR(16) * (1.0 - AUX_VAR_C_total / AUX_VAR_C_max) * (1.0 - PARAM(190) * SPVAR(25) / (SPVAR(25) + PARAM(183) + PARAM(184)) * AUX_VAR_R_cabo);

    realtype ReactionFlux4 = PARAM(43) * SPVAR(16);

    realtype ReactionFlux5 = PARAM(55) / PARAM(45);

    realtype ReactionFlux6 = PARAM(56) / PARAM(45) * SPVAR(8) / (PARAM(57) / PARAM(45) + SPVAR(8));

    realtype ReactionFlux7 = PARAM(56) / PARAM(45) * SPVAR(0) / (PARAM(57) / PARAM(45) + SPVAR(0));

    realtype ReactionFlux8 = PARAM(58) * SPVAR(8);

    realtype ReactionFlux9 = PARAM(58) * SPVAR(0);

    realtype ReactionFlux10 = PARAM(52) * PARAM(1) * SPVAR(0);

    realtype ReactionFlux11 = PARAM(53) * SPVAR(8);

    realtype ReactionFlux12 = PARAM(47) * PARAM(2) * SPVAR(0);

    realtype ReactionFlux13 = PARAM(48) * SPVAR(27);

    realtype ReactionFlux14 = PARAM(49) * AUX_VAR_H_act1 * AUX_VAR_H_P1 * SPVAR(27);

    realtype ReactionFlux15 = PARAM(49) * AUX_VAR_H_act1 * AUX_VAR_H_P1 * PARAM(46) * SPVAR(27);

    realtype ReactionFlux16 = PARAM(50) * SPVAR(28);

    realtype ReactionFlux17 = PARAM(50)  / AUX_VAR_N_aT * std::pow(2.0, AUX_VAR_N_aT) * SPVAR(28);

    realtype ReactionFlux18 = PARAM(51) * SPVAR(1);

    realtype ReactionFlux19 = PARAM(51) * SPVAR(9);

    realtype ReactionFlux20 = PARAM(51) * SPVAR(17);

    realtype ReactionFlux21 = PARAM(51) * SPVAR(29);

    realtype ReactionFlux22 = PARAM(67) * SPVAR(17) * AUX_VAR_Tregs_T / (AUX_VAR_C_total + AUX_VAR_T_T + PARAM(31));

    realtype ReactionFlux23 = PARAM(65) * SPVAR(17) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_T + PARAM(31)) * AUX_VAR_H_PD1;

    realtype ReactionFlux24 = PARAM(52) * PARAM(1) * SPVAR(1);

    realtype ReactionFlux25 = PARAM(53) * SPVAR(9);

    realtype ReactionFlux26 = PARAM(54) * PARAM(35) * AUX_VAR_C_total * SPVAR(1) * (1.0 + PARAM(192) * SPVAR(25) / (SPVAR(25) + PARAM(186))) * SPVAR(26) / AUX_VAR_C_total;

    realtype ReactionFlux27 = PARAM(48) * SPVAR(29);

    realtype ReactionFlux28 = PARAM(59) * SPVAR(30) * PARAM(2);

    realtype ReactionFlux29 = PARAM(60) * AUX_VAR_T_LN * SPVAR(30) / (PARAM(62) + SPVAR(30));

    realtype ReactionFlux30 = std::pow(2.0, AUX_VAR_N_aT) * PARAM(61) * SPVAR(28);

    realtype ReactionFlux31 = PARAM(66) * SPVAR(17) * SPVAR(16) / (AUX_VAR_C_total + AUX_VAR_T_T + PARAM(31)) * (1.0 - AUX_VAR_H_PD1);

    realtype ReactionFlux32 = PARAM(78) / PARAM(68);

    realtype ReactionFlux33 = PARAM(79) / PARAM(68) * SPVAR(10) / (PARAM(80) / PARAM(68) + SPVAR(10));

    realtype ReactionFlux34 = PARAM(79) / PARAM(68) * SPVAR(2) / (PARAM(80) / PARAM(68) + SPVAR(2));

    realtype ReactionFlux35 = PARAM(81) * SPVAR(10);

    realtype ReactionFlux36 = PARAM(81) * SPVAR(2);

    realtype ReactionFlux37 = PARAM(75) * PARAM(1) * SPVAR(2);

    realtype ReactionFlux38 = PARAM(76) * SPVAR(10);

    realtype ReactionFlux39 = PARAM(70) * PARAM(2) * SPVAR(2);

    realtype ReactionFlux40 = PARAM(71) * SPVAR(31);

    realtype ReactionFlux41 = PARAM(72) * AUX_VAR_H_act0 * AUX_VAR_H_P0 * SPVAR(31);

    realtype ReactionFlux42 = PARAM(72) * AUX_VAR_H_act0 * AUX_VAR_H_P0 * PARAM(69) * SPVAR(31);

    realtype ReactionFlux43 = PARAM(73) * SPVAR(32);

    realtype ReactionFlux44 = PARAM(73) / AUX_VAR_N_aT * std::pow(2.0, AUX_VAR_N_aT) * SPVAR(32);

    realtype ReactionFlux45 = PARAM(74) * SPVAR(3);

    realtype ReactionFlux46 = PARAM(74) * SPVAR(11);

    realtype ReactionFlux47 = PARAM(74) * SPVAR(18);

    realtype ReactionFlux48 = PARAM(74) * SPVAR(33);

    realtype ReactionFlux49 = PARAM(75) * PARAM(1) * SPVAR(3);

    realtype ReactionFlux50 = PARAM(76) * SPVAR(11);

    realtype ReactionFlux51 = PARAM(77) * PARAM(35) * AUX_VAR_C_total * SPVAR(3) * (1.0 + PARAM(192) * SPVAR(25) / (SPVAR(25) + PARAM(186))) * SPVAR(26) / AUX_VAR_C_total;

    realtype ReactionFlux52 = PARAM(71) * SPVAR(33);

    realtype ReactionFlux53 = PARAM(60) * SPVAR(33) * SPVAR(30) / (PARAM(62) + SPVAR(30));

    realtype ReactionFlux54 = PARAM(84) * (PARAM(86) * AUX_VAR_V_T - SPVAR(19));

    realtype ReactionFlux55 = PARAM(84) * (PARAM(87) * PARAM(2) - SPVAR(34));

    realtype ReactionFlux56 = PARAM(82) * SPVAR(21) / (SPVAR(21) + PARAM(90)) * SPVAR(19);

    realtype ReactionFlux57 = PARAM(83) * SPVAR(20);

    realtype ReactionFlux58 = PARAM(85) * SPVAR(20);

    realtype ReactionFlux59 = PARAM(85) * SPVAR(35);

    realtype ReactionFlux60 = PARAM(88) * (PARAM(89) - SPVAR(21)) * AUX_VAR_V_T;

    realtype ReactionFlux61 = AUX_VAR_R_Tcell * PARAM(91);

    realtype ReactionFlux62 = PARAM(94) * SPVAR(42) * PARAM(4) - PARAM(93) * SPVAR(45) * PARAM(5);

    realtype ReactionFlux63 = PARAM(69) * PARAM(102) * (PARAM(43) + PARAM(66) * SPVAR(17) / (AUX_VAR_C_total + AUX_VAR_T_T + PARAM(31)) * (1.0 - AUX_VAR_H_PD1)) * SPVAR(16) * AUX_VAR_V_T;

    realtype ReactionFlux64 = PARAM(96) * SPVAR(22) * AUX_VAR_V_T;

    realtype ReactionFlux65 = PARAM(95) * SPVAR(35) * SPVAR(22) * AUX_VAR_V_T;

    realtype ReactionFlux66 = PARAM(95) * PARAM(31) * SPVAR(22) * PARAM(3);

    realtype ReactionFlux67 = PARAM(97) * SPVAR(38) * PARAM(3);

    realtype ReactionFlux68 = PARAM(98) * SPVAR(39) * PARAM(3);

    realtype ReactionFlux69 = PARAM(99) * SPVAR(39) * SPVAR(42) * PARAM(4);

    realtype ReactionFlux70 = PARAM(100) * PARAM(99) * SPVAR(43) * PARAM(4);

    realtype ReactionFlux71 = PARAM(100) * PARAM(99) * SPVAR(46) * PARAM(5);

    realtype ReactionFlux72 = PARAM(94) * SPVAR(43) * PARAM(4);

    realtype ReactionFlux73 = PARAM(46) * PARAM(116) * (PARAM(43) + PARAM(66) * SPVAR(17) / (AUX_VAR_C_total + AUX_VAR_T_T + PARAM(31)) * (1.0 - AUX_VAR_H_PD1)) * SPVAR(16) * AUX_VAR_V_T;

    realtype ReactionFlux74 = PARAM(110) * SPVAR(23) * AUX_VAR_V_T;

    realtype ReactionFlux75 = PARAM(109) * SPVAR(35) * SPVAR(23) * AUX_VAR_V_T;

    realtype ReactionFlux76 = PARAM(109) * PARAM(31) * SPVAR(23) * PARAM(3);

    realtype ReactionFlux77 = PARAM(111) * SPVAR(40) * PARAM(3);

    realtype ReactionFlux78 = PARAM(112) * SPVAR(41) * PARAM(3);

    realtype ReactionFlux79 = PARAM(113) * SPVAR(41) * SPVAR(42) * PARAM(4);

    realtype ReactionFlux80 = PARAM(114) * PARAM(113) * SPVAR(44) * PARAM(4);

    realtype ReactionFlux81 = PARAM(114) * PARAM(113) * SPVAR(47) * PARAM(5);

    realtype ReactionFlux82 = PARAM(94) * SPVAR(44) * PARAM(4);

    realtype ReactionFlux83 = PARAM(123) * (SPVAR(4) / PARAM(128) - SPVAR(12) / PARAM(129)) * PARAM(0);

    realtype ReactionFlux84 = PARAM(124) * (SPVAR(4) / PARAM(128) - SPVAR(24) / PARAM(130)) * (1.0 + PARAM(192) * SPVAR(25) / (SPVAR(25) + PARAM(186))) * SPVAR(26) / AUX_VAR_C_total * PARAM(0);

    realtype ReactionFlux85 = PARAM(125) * (SPVAR(4) / PARAM(128) - SPVAR(36) / PARAM(131)) * PARAM(0);

    realtype ReactionFlux86 = PARAM(126) * SPVAR(24) / PARAM(130) * AUX_VAR_V_T;

    realtype ReactionFlux87 = PARAM(126) * SPVAR(36) / PARAM(131) * PARAM(2);

    realtype ReactionFlux88 = PARAM(127) * SPVAR(4) * PARAM(0);

    realtype ReactionFlux89 = PARAM(141) * (PARAM(132) - SPVAR(48) - SPVAR(49) - SPVAR(50) - 2.0 * SPVAR(51)) * (PARAM(133) - SPVAR(48) - PARAM(10) - 2.0 * PARAM(11)) * PARAM(6);

    realtype ReactionFlux90 = PARAM(142) * SPVAR(48) * PARAM(6);

    realtype ReactionFlux91 = PARAM(143) * (PARAM(132) - SPVAR(48) - SPVAR(49) - SPVAR(50) - 2.0 * SPVAR(51)) * (PARAM(134) - SPVAR(49)) * PARAM(6);

    realtype ReactionFlux92 = PARAM(144) * SPVAR(49) * PARAM(6);

    realtype ReactionFlux93 = PARAM(155) * (PARAM(132) - SPVAR(48) - SPVAR(49) - SPVAR(50) - 2.0 * SPVAR(51)) * SPVAR(24) * PARAM(6);

    realtype ReactionFlux94 = PARAM(156) * SPVAR(50) * PARAM(6);

    realtype ReactionFlux95 = PARAM(157) * PARAM(155) * (PARAM(132) - SPVAR(48) - SPVAR(49) - SPVAR(50) - 2.0 * SPVAR(51)) * SPVAR(50) * PARAM(6);

    realtype ReactionFlux96 = PARAM(156) * SPVAR(51) * PARAM(6);

    realtype ReactionFlux97 = PARAM(165) * PARAM(166) * SPVAR(11);

    realtype ReactionFlux98 = PARAM(165) * PARAM(167) * SPVAR(18);

    realtype ReactionFlux99 = PARAM(174) * PARAM(169) * SPVAR(7) * PARAM(0);

    realtype ReactionFlux100 = PARAM(174) * PARAM(168) * SPVAR(6) * PARAM(0);

    realtype ReactionFlux101 = PARAM(175) * (SPVAR(5) / PARAM(179) - SPVAR(13) / PARAM(180)) * PARAM(0);

    realtype ReactionFlux102 = PARAM(176) * (SPVAR(5) / PARAM(179) - SPVAR(25) / PARAM(181)) * (1.0 + PARAM(192) * SPVAR(25) / (SPVAR(25) + PARAM(186))) * SPVAR(26) / AUX_VAR_C_total * PARAM(0);

    realtype ReactionFlux103 = PARAM(177) * (SPVAR(5) / PARAM(179) - SPVAR(37) / PARAM(182)) * PARAM(0);

    realtype ReactionFlux104 = PARAM(178) * SPVAR(25) / PARAM(181) * AUX_VAR_V_T;

    realtype ReactionFlux105 = PARAM(178) * SPVAR(37) / PARAM(182) * PARAM(2);

    realtype ReactionFlux106 = PARAM(170) * SPVAR(5) / (SPVAR(5) + PARAM(171)) * PARAM(0);

    realtype ReactionFlux107 = PARAM(188) * PARAM(31) * std::pow(AUX_VAR_C_total / PARAM(31), 2.0 / 3.0) * std::pow(SPVAR(26) / PARAM(31), 1.0 / 3.0) - PARAM(191) * (SPVAR(25) / (SPVAR(25) + PARAM(186))) * AUX_VAR_R_cabo * SPVAR(26);

    //dydt:

    //d(V_C.nT1)/dt
    NV_DATA_S(ydot)[0] = ReactionFlux5 + ReactionFlux7 - ReactionFlux9 - ReactionFlux10 + ReactionFlux11 - ReactionFlux12 + ReactionFlux13;

    //d(V_C.T1)/dt
    NV_DATA_S(ydot)[1] =  - ReactionFlux18 - ReactionFlux24 + ReactionFlux25 - ReactionFlux26 + ReactionFlux27;

    //d(V_C.nT0)/dt
    NV_DATA_S(ydot)[2] = ReactionFlux32 + ReactionFlux34 - ReactionFlux36 - ReactionFlux37 + ReactionFlux38 - ReactionFlux39 + ReactionFlux40;

    //d(V_C.T0)/dt
    NV_DATA_S(ydot)[3] =  - ReactionFlux45 - ReactionFlux49 + ReactionFlux50 - ReactionFlux51 + ReactionFlux52;

    //d(V_C.nivolumab)/dt
    NV_DATA_S(ydot)[4] = 1/PARAM(0)*( - ReactionFlux83 - ReactionFlux84 - ReactionFlux85 + ReactionFlux87 - ReactionFlux88);

    //d(V_C.cabozantinib)/dt
    NV_DATA_S(ydot)[5] = 1/PARAM(0)*(ReactionFlux99 + ReactionFlux100 - ReactionFlux101 - ReactionFlux102 - ReactionFlux103 + ReactionFlux105 - ReactionFlux106);

    //d(V_C.A_site1)/dt
    NV_DATA_S(ydot)[6] = 1/PARAM(0)*( - ReactionFlux100);

    //d(V_C.A_site2)/dt
    NV_DATA_S(ydot)[7] = 1/PARAM(0)*( - ReactionFlux99);

    //d(V_P.nT1)/dt
    NV_DATA_S(ydot)[8] = ReactionFlux6 - ReactionFlux8 + ReactionFlux10 - ReactionFlux11;

    //d(V_P.T1)/dt
    NV_DATA_S(ydot)[9] =  - ReactionFlux19 + ReactionFlux24 - ReactionFlux25;

    //d(V_P.nT0)/dt
    NV_DATA_S(ydot)[10] = ReactionFlux33 - ReactionFlux35 + ReactionFlux37 - ReactionFlux38;

    //d(V_P.T0)/dt
    NV_DATA_S(ydot)[11] =  - ReactionFlux46 + ReactionFlux49 - ReactionFlux50 - ReactionFlux97;

    //d(V_P.nivolumab)/dt
    NV_DATA_S(ydot)[12] = 1/PARAM(1)*(ReactionFlux83);

    //d(V_P.cabozantinib)/dt
    NV_DATA_S(ydot)[13] = 1/PARAM(1)*(ReactionFlux101);

    //d(V_T.C_x)/dt
    NV_DATA_S(ydot)[14] =  - ReactionFlux1 + ReactionFlux4 + ReactionFlux31;

    //d(V_T.T_exh)/dt
    NV_DATA_S(ydot)[15] =  - ReactionFlux2 + ReactionFlux20 + ReactionFlux22 + ReactionFlux23 + ReactionFlux47;

    //d(V_T.C1)/dt
    NV_DATA_S(ydot)[16] = ReactionFlux3 - ReactionFlux4 - ReactionFlux31;

    //d(V_T.T1)/dt
    NV_DATA_S(ydot)[17] =  - ReactionFlux20 - ReactionFlux22 - ReactionFlux23 + ReactionFlux26;

    //d(V_T.T0)/dt
    NV_DATA_S(ydot)[18] =  - ReactionFlux47 + ReactionFlux51 - ReactionFlux98;

    //d(V_T.APC)/dt
    NV_DATA_S(ydot)[19] = ReactionFlux54 - ReactionFlux56;

    //d(V_T.mAPC)/dt
    NV_DATA_S(ydot)[20] = ReactionFlux56 - ReactionFlux57 - ReactionFlux58;

    //d(V_T.c)/dt
    NV_DATA_S(ydot)[21] = 1/AUX_VAR_V_T*(ReactionFlux60 + ReactionFlux61);

    //d(V_T.P0)/dt
    NV_DATA_S(ydot)[22] = 1/AUX_VAR_V_T*(ReactionFlux63 - ReactionFlux64 - ReactionFlux65);

    //d(V_T.P1)/dt
    NV_DATA_S(ydot)[23] = 1/AUX_VAR_V_T*(ReactionFlux73 - ReactionFlux74 - ReactionFlux75);

    //d(V_T.nivolumab)/dt
    NV_DATA_S(ydot)[24] = 1/AUX_VAR_V_T*(ReactionFlux84 - ReactionFlux86);

    //d(V_T.cabozantinib)/dt
    NV_DATA_S(ydot)[25] = 1/AUX_VAR_V_T*(ReactionFlux102 - ReactionFlux104);

    //d(V_T.v_cabo)/dt
    NV_DATA_S(ydot)[26] = ReactionFlux107;

    //d(V_LN.nT1)/dt
    NV_DATA_S(ydot)[27] = ReactionFlux12 - ReactionFlux13 - ReactionFlux14;

    //d(V_LN.aT1)/dt
    NV_DATA_S(ydot)[28] = ReactionFlux15 - ReactionFlux16;

    //d(V_LN.T1)/dt
    NV_DATA_S(ydot)[29] = ReactionFlux17 - ReactionFlux21 - ReactionFlux27;

    //d(V_LN.IL2)/dt
    NV_DATA_S(ydot)[30] = 1/PARAM(2)*( - ReactionFlux28 - ReactionFlux29 + ReactionFlux30 - ReactionFlux53);

    //d(V_LN.nT0)/dt
    NV_DATA_S(ydot)[31] = ReactionFlux39 - ReactionFlux40 - ReactionFlux41;

    //d(V_LN.aT0)/dt
    NV_DATA_S(ydot)[32] = ReactionFlux42 - ReactionFlux43;

    //d(V_LN.T0)/dt
    NV_DATA_S(ydot)[33] = ReactionFlux44 - ReactionFlux48 - ReactionFlux52;

    //d(V_LN.APC)/dt
    NV_DATA_S(ydot)[34] = ReactionFlux55;

    //d(V_LN.mAPC)/dt
    NV_DATA_S(ydot)[35] = ReactionFlux57 - ReactionFlux59;

    //d(V_LN.nivolumab)/dt
    NV_DATA_S(ydot)[36] = 1/PARAM(2)*(ReactionFlux85 + ReactionFlux86 - ReactionFlux87);

    //d(V_LN.cabozantinib)/dt
    NV_DATA_S(ydot)[37] = 1/PARAM(2)*(ReactionFlux103 + ReactionFlux104 - ReactionFlux105);

    //d(V_e.P0)/dt
    NV_DATA_S(ydot)[38] = 1/PARAM(3)*(ReactionFlux66 - ReactionFlux67);

    //d(V_e.p0)/dt
    NV_DATA_S(ydot)[39] = 1/PARAM(3)*(ReactionFlux67 - ReactionFlux68 - ReactionFlux69 + ReactionFlux70);

    //d(V_e.P1)/dt
    NV_DATA_S(ydot)[40] = 1/PARAM(3)*(ReactionFlux76 - ReactionFlux77);

    //d(V_e.p1)/dt
    NV_DATA_S(ydot)[41] = 1/PARAM(3)*(ReactionFlux77 - ReactionFlux78 - ReactionFlux79 + ReactionFlux80);

    //d(A_e.M1)/dt
    NV_DATA_S(ydot)[42] = 1/PARAM(4)*( - ReactionFlux62 - ReactionFlux69 + ReactionFlux70 - ReactionFlux79 + ReactionFlux80);

    //d(A_e.M1p0)/dt
    NV_DATA_S(ydot)[43] = 1/PARAM(4)*(ReactionFlux69 - ReactionFlux70 - ReactionFlux72);

    //d(A_e.M1p1)/dt
    NV_DATA_S(ydot)[44] = 1/PARAM(4)*(ReactionFlux79 - ReactionFlux80 - ReactionFlux82);

    //d(A_s.M1)/dt
    NV_DATA_S(ydot)[45] = 1/PARAM(5)*(ReactionFlux62 + ReactionFlux71 + ReactionFlux81);

    //d(A_s.M1p0)/dt
    NV_DATA_S(ydot)[46] = 1/PARAM(5)*( - ReactionFlux71 + ReactionFlux72);

    //d(A_s.M1p1)/dt
    NV_DATA_S(ydot)[47] = 1/PARAM(5)*( - ReactionFlux81 + ReactionFlux82);

    //d(syn_CT.PD1_PDL1)/dt
    NV_DATA_S(ydot)[48] = 1/PARAM(6)*(ReactionFlux89 - ReactionFlux90);

    //d(syn_CT.PD1_PDL2)/dt
    NV_DATA_S(ydot)[49] = 1/PARAM(6)*(ReactionFlux91 - ReactionFlux92);

    //d(syn_CT.PD1_aPD1)/dt
    NV_DATA_S(ydot)[50] = 1/PARAM(6)*(ReactionFlux93 - ReactionFlux94 - ReactionFlux95 + ReactionFlux96);

    //d(syn_CT.PD1_aPD1_PD1)/dt
    NV_DATA_S(ydot)[51] = 1/PARAM(6)*(ReactionFlux95 - ReactionFlux96);

    return(0);
}
int ODE_system::g(realtype t, N_Vector y, realtype *gout, void *user_data){

    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);

    //Assignment rules:

    //V_T.C1 < (0.9 * cell)
    gout[0] = 0.9 * PARAM(31) - (SPVAR(16));

    //V_T.T1 < (1e-10 * cell)
    gout[1] = 1e-10 * PARAM(31) - (SPVAR(17));

    //V_T.c < (1e-10 * molarity)
    gout[2] = 1e-10 * PARAM(34) - (SPVAR(21));

    //V_C.nivolumab < (1e-10 * molarity)
    gout[3] = 1e-10 * PARAM(34) - (SPVAR(4));

    //V_P.nivolumab < (1e-10 * molarity)
    gout[4] = 1e-10 * PARAM(34) - (SPVAR(12));

    //V_T.nivolumab < (1e-10 * molarity)
    gout[5] = 1e-10 * PARAM(34) - (SPVAR(24));

    //V_LN.nivolumab < (1e-10 * molarity)
    gout[6] = 1e-10 * PARAM(34) - (SPVAR(36));

    return(0);
}

bool ODE_system::triggerComponentEvaluate(int i, realtype t, bool curr) {

    bool discrete = false;
    realtype diff = 0;
    bool eval = false;
    //Assignment rules:

    switch(i)
    {
    case 0:
        //V_T.C1 < (0.9 * cell)
        diff = 0.9 * _class_parameter[31] - (NV_DATA_S(_y)[16]);
        break;
    case 1:
        //V_T.T1 < (1e-10 * cell)
        diff = 1e-10 * _class_parameter[31] - (NV_DATA_S(_y)[17]);
        break;
    case 2:
        //V_T.c < (1e-10 * molarity)
        diff = 1e-10 * _class_parameter[34] - (NV_DATA_S(_y)[21]);
        break;
    case 3:
        //V_C.nivolumab < (1e-10 * molarity)
        diff = 1e-10 * _class_parameter[34] - (NV_DATA_S(_y)[4]);
        break;
    case 4:
        //V_P.nivolumab < (1e-10 * molarity)
        diff = 1e-10 * _class_parameter[34] - (NV_DATA_S(_y)[12]);
        break;
    case 5:
        //V_T.nivolumab < (1e-10 * molarity)
        diff = 1e-10 * _class_parameter[34] - (NV_DATA_S(_y)[24]);
        break;
    case 6:
        //V_LN.nivolumab < (1e-10 * molarity)
        diff = 1e-10 * _class_parameter[34] - (NV_DATA_S(_y)[36]);
        break;
    default:
        break;
    }
    if (!discrete){
        eval = diff == 0 ? curr : (diff > 0);
    }
    return eval;
}

bool ODE_system::eventEvaluate(int i) {
    bool eval = false;
    switch(i)
    {
    case 0:
        eval = getSatisfied(0);
        break;
    case 1:
        eval = getSatisfied(1);
        break;
    case 2:
        eval = getSatisfied(2);
        break;
    case 3:
        eval = getSatisfied(3);
        break;
    case 4:
        eval = getSatisfied(4);
        break;
    case 5:
        eval = getSatisfied(5);
        break;
    case 6:
        eval = getSatisfied(6);
        break;
    default:
        break;
    }
    return eval;
}

bool ODE_system::eventExecution(int i, bool delayed, realtype& dt){

    bool setDelay = false;

    //Assignment rules:

    switch(i)
    {
    case 0:
        NV_DATA_S(_y)[16] = 0.0 * _class_parameter[31];
        break;
    case 1:
        NV_DATA_S(_y)[17] = 0.0 * 1e-08 * _class_parameter[31];
        break;
    case 2:
        NV_DATA_S(_y)[21] = 0.0 * 1e-08 * _class_parameter[34];
        break;
    case 3:
        NV_DATA_S(_y)[4] = 0.0 * 1e-08 * _class_parameter[34];
        break;
    case 4:
        NV_DATA_S(_y)[12] = 0.0 * 1e-08 * _class_parameter[34];
        break;
    case 5:
        NV_DATA_S(_y)[24] = 0.0 * 1e-08 * _class_parameter[34];
        break;
    case 6:
        NV_DATA_S(_y)[36] = 0.0 * 1e-08 * _class_parameter[34];
        break;
    default:
        break;
    }
    return setDelay;
}
void ODE_system::update_y_other(void){

    //syn_CT.PDL1_aPDL1
    _species_other[0] = _class_parameter[10];

    //syn_CT.PDL1_aPDL1_PDL1
    _species_other[1] = _class_parameter[11];

    //syn_DN.CD28_CD80
    _species_other[2] = _class_parameter[12];

    //syn_DN.CD28_CD80_CD28
    _species_other[3] = _class_parameter[13];

    //syn_DN.CD28_CD86
    _species_other[4] = _class_parameter[14];

    //syn_DN.PDL1_CD80
    _species_other[5] = _class_parameter[15];

    //syn_DN.CTLA4_CD80
    _species_other[6] = _class_parameter[16];

    //syn_DN.CTLA4_CD80_CTLA4
    _species_other[7] = _class_parameter[17];

    //syn_DN.CD80_CTLA4_CD80
    _species_other[8] = _class_parameter[18];

    //syn_DN.CTLA4_CD80_CTLA4_CD80
    _species_other[9] = _class_parameter[19];

    //syn_DN.CTLA4_CD86
    _species_other[10] = _class_parameter[20];

    //syn_DN.CD86_CTLA4_CD86
    _species_other[11] = _class_parameter[21];

    //syn_DN.CTLA4_aCTLA4
    _species_other[12] = _class_parameter[22];

    //syn_DN.CTLA4_aCTLA4_CTLA4
    _species_other[13] = _class_parameter[23];

    //syn_DN.PDL1_aPDL1
    _species_other[14] = _class_parameter[24];

    //syn_DN.PDL1_aPDL1_PDL1
    _species_other[15] = _class_parameter[25];

    //Treg_P.CTLA4_aCTLA4
    _species_other[16] = _class_parameter[26];

    //Treg_P.CTLA4_aCTLA4_CTLA4
    _species_other[17] = _class_parameter[27];

    //Treg_T.CTLA4_aCTLA4
    _species_other[18] = _class_parameter[28];

    //Treg_T.CTLA4_aCTLA4_CTLA4
    _species_other[19] = _class_parameter[29];

    return;
}
std::string ODE_system::getHeader(){

    std::string s = "";
    s += ",V_C.nT1";
    s += ",V_C.T1";
    s += ",V_C.nT0";
    s += ",V_C.T0";
    s += ",V_C.nivolumab";
    s += ",V_C.cabozantinib";
    s += ",V_C.A_site1";
    s += ",V_C.A_site2";
    s += ",V_P.nT1";
    s += ",V_P.T1";
    s += ",V_P.nT0";
    s += ",V_P.T0";
    s += ",V_P.nivolumab";
    s += ",V_P.cabozantinib";
    s += ",V_T.C_x";
    s += ",V_T.T_exh";
    s += ",V_T.C1";
    s += ",V_T.T1";
    s += ",V_T.T0";
    s += ",V_T.APC";
    s += ",V_T.mAPC";
    s += ",V_T.c";
    s += ",V_T.P0";
    s += ",V_T.P1";
    s += ",V_T.nivolumab";
    s += ",V_T.cabozantinib";
    s += ",V_T.v_cabo";
    s += ",V_LN.nT1";
    s += ",V_LN.aT1";
    s += ",V_LN.T1";
    s += ",V_LN.IL2";
    s += ",V_LN.nT0";
    s += ",V_LN.aT0";
    s += ",V_LN.T0";
    s += ",V_LN.APC";
    s += ",V_LN.mAPC";
    s += ",V_LN.nivolumab";
    s += ",V_LN.cabozantinib";
    s += ",V_e.P0";
    s += ",V_e.p0";
    s += ",V_e.P1";
    s += ",V_e.p1";
    s += ",A_e.M1";
    s += ",A_e.M1p0";
    s += ",A_e.M1p1";
    s += ",A_s.M1";
    s += ",A_s.M1p0";
    s += ",A_s.M1p1";
    s += ",syn_CT.PD1_PDL1";
    s += ",syn_CT.PD1_PDL2";
    s += ",syn_CT.PD1_aPD1";
    s += ",syn_CT.PD1_aPD1_PD1";
    s += ",syn_CT.PDL1_aPDL1";
    s += ",syn_CT.PDL1_aPDL1_PDL1";
    s += ",syn_DN.CD28_CD80";
    s += ",syn_DN.CD28_CD80_CD28";
    s += ",syn_DN.CD28_CD86";
    s += ",syn_DN.PDL1_CD80";
    s += ",syn_DN.CTLA4_CD80";
    s += ",syn_DN.CTLA4_CD80_CTLA4";
    s += ",syn_DN.CD80_CTLA4_CD80";
    s += ",syn_DN.CTLA4_CD80_CTLA4_CD80";
    s += ",syn_DN.CTLA4_CD86";
    s += ",syn_DN.CD86_CTLA4_CD86";
    s += ",syn_DN.CTLA4_aCTLA4";
    s += ",syn_DN.CTLA4_aCTLA4_CTLA4";
    s += ",syn_DN.PDL1_aPDL1";
    s += ",syn_DN.PDL1_aPDL1_PDL1";
    s += ",Treg_P.CTLA4_aCTLA4";
    s += ",Treg_P.CTLA4_aCTLA4_CTLA4";
    s += ",Treg_T.CTLA4_aCTLA4";
    s += ",Treg_T.CTLA4_aCTLA4_CTLA4";
    return s;
}
realtype ODE_system::get_unit_conversion_species(int i) const{

    static std::vector<realtype> scalor = {
        //sp_var
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000.0,
        1000.0,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000.0,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        //sp_other
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
    };
    return scalor[i];
}
realtype ODE_system::get_unit_conversion_nspvar(int i) const{

    static std::vector<realtype> scalor = {
    };
    return scalor[i];
}
};
