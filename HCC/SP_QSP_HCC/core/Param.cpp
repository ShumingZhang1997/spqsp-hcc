#include "Param.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <math.h>
#include "../ode/Param.h"
#include "../ode/ODE_system.h"

// parameters from QSP module. 
extern CancerVCT::Param qsp_params;

#define QP(x) CancerVCT::ODE_system::get_class_param(x)
#define AVOGADROS 6.022140857E23 
#define PI 3.1415926525897932384626

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

namespace pt = boost::property_tree;

static int SEC_PER_DAY = 86400;
static int HOUR_PER_DAY = 24;

// must match the order in the ParamInt and ParamFloat enums
#define PARAM_DESCRIPTION_FIELD_COUNT 3
const char* _description[][PARAM_DESCRIPTION_FIELD_COUNT] =
{
	//{"fullpath", "desc", "constraint"}

	//------------------------ float --------------------------//
	/* QSP */
	{ "Param.QSP.simulation.weight_qsp", "", "prob" },
	{ "Param.QSP.simulation.t_steadystate", "days", "pos" },
	{ "Param.QSP.simulation.t_resection", "days", "pos" },
	{ "Param.QSP.simulation.presimulation_diam_frac", "", "pos" },
	/* ABM */
	//environmental
	{ "Param.ABM.Environment.SecPerSlice", "", "pos" },
	{ "Param.ABM.Environment.recSiteFactor", "number of adhesion site per port voxel", "pos" },
	//Adhesion molecule density Reference
	//Increased ICAM-1 Expression Causes Endothelial Cell Leakiness, Cytoskeletal Reorganizationand Junctional Alterations
	// 657 molecules per field --> 1 field = 2.06um * 2.06um --> 657/(2.06 * 2.06) um
	{ "Param.ABM.Environment.adhSiteDens", "total adhesion site density on tumor vasculature, molecule/um^3", "pos" },
	//pharmacokinetics	
	{ "Param.ABM.Pharmacokinetics.nivoDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.nivoDose", "mole/m^3", "pos" },
	{ "Param.ABM.Pharmacokinetics.ipiDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.ipiDose", "mole/m^3", "pos" },
	{ "Param.ABM.Pharmacokinetics.caboDoseIntervalTime", "days", "pos" },
	{ "Param.ABM.Pharmacokinetics.caboDose", "mole/m^3", "pos" },
	//T cell
	{ "Param.ABM.TCell.lifespanSD", "days", "pos" },
	{ "Param.ABM.TCell.moveProb", "", "pr" },
	{ "Param.ABM.TCell.IL2_release_time", "amount of time to release IL2 after stimulation, sec", "pos" },
	{ "Param.ABM.TCell.IL2_prolif_th", "accumulative IL2 exposure to proliferate, sec*ng/mL", "pos" },
	{ "Param.ABM.TCell.TGFB_EC50", "Half-Maximal TGFb level for CD8 T cell inhibition, nanomolarity", "pos" },
	{ "Param.ABM.TCell.IFNg_recruit_Half", "Half-Maximal IFNg level for T cell recruitment, ng/mL", "pos" },
	// TCD4
	{ "Param.ABM.TCD4.moveProb", "", "pr" },
	{ "Param.ABM.TCD4.k_Th_Treg", "Th differentiation rate to Treg, 1/day", "pr" },
	{ "Param.ABM.TCD4.Th_frac", "Th ratio in central compartment, 1/day", "pr" },
	{ "Param.ABM.TCD4.TGFB_release_time", "amount of time to release TGFB after stimulation, sec", "pos"  },
	// MDSC
	{ "Param.ABM.MDSC.lifespanMean", "", "pos" },
	{ "Param.ABM.MDSC.lifespanSD", "", "pos" },
	{ "Param.ABM.MDSC.moveProb", "", "pr" },
	{ "Param.ABM.MDSC.k_rec_MDSC", "1/day", "pos" },
	{ "Param.ABM.MDSC.EC50_CCL2_rec", "molarity", "pos" },
	{ "Param.ABM.MDSC.IC50_ArgI_CTL", "mU", "pos" },
	{ "Param.ABM.MDSC.EC50_ArgI_Treg", "mU", "pos" },
	{ "Param.ABM.MDSC.IC50_NO_CTL", "molarity", "pos" },
	{ "Param.ABM.MDSC.Myeloid_scale", "", "pos" },
	//Mac
	{ "Param.ABM.Mac.lifespanMean", "", "pos" },
	{ "Param.ABM.Mac.lifespanSD", "", "pos" },
	{ "Param.ABM.Mac.k_rec_MAC", "cell/(milliliter*day)", "pos" },
	{ "Param.ABM.Mac.moveProb", "", "pr" },
	{ "Param.ABM.Mac.PD1_total", "molecules", "pos"},
	{ "Param.ABM.Mac.A_Mac", "surface area of macrophage, micrometer^2", "pos"},
	{ "Param.ABM.Mac.A_Mac_syn", "synapse area of macrophage, micrometer^2", "pos"},
	{ "Param.ABM.Mac.m2_pol", "Rate of M1 to M2 macrophage polarization, 1/day", "pos"},
	{ "Param.ABM.Mac.m1_pol", "Rate of M2 to M1 macrophage polarization, 1/day", "pos"},
	{ "Param.ABM.Mac.k_M1_phago", "Rate of M1 macrophage-mediated killing of cancer cell, 1/day", "pos"},
	{ "Param.ABM.Mac.TGFB_EC50", "Half-maximal TGFB level for M1 to M2 polarization / maintaining Treg function, nanomolarity", "pos"},
	{ "Param.ABM.Mac.IL10_EC50", "Half-maximal IL10 level for M1 to M2 polarization / maintaining Treg function, picomolarity", "pos"},
	{ "Param.ABM.Mac.IL10_half_phago", "Half-maximal IL-10 level for phagocytosis by macrophage , picomolarity", "pos"},
	{ "Param.ABM.Mac.IFNg_EC50", "molecules", "pos"},
	{ "Param.ABM.Mac.IL12_EC50", "surface area of macrophage, micrometer^2", "pos"},
	{ "Param.ABM.Mac.SIRPa_total", "SIRPa expression on macrophages, molecule/micrometer^2", "pos"},
	{ "Param.ABM.Mac.SIRPa_half", "SIRPa occupancy for half-maximal CD47/SIRPa inhibition on phagocytosis of cancer cells by macrophages, molecule/micrometer^2", "pos"},
	{ "Param.ABM.Mac.kon_SIRPa_CD47", "kon of CD47-SIRPa binding , 1/(micromolarity*minute*nanometer)", "pos"},
	{ "Param.ABM.Mac.koff_SIRPa_CD47", "koff of CD47-SIRPa binding , 1/minute", "pos"},
	{ "Param.ABM.Mac.n_SIRPa_CD47", "Hill coefficient for CD47/SIRPa inhibition on phagocytosis of cancer cells by macrophage, dimensionless", "pos"},
	//cancer cell
	{"Param.ABM.CancerCell.asymmetricDivProb", "", "pr"},
	{"Param.ABM.CancerCell.progGrowthRate", "per day", "pos"},
	{"Param.ABM.CancerCell.senescentDeathRate", "per day", "pos"},
	{"Param.ABM.CancerCell.moveProb", "", "pr"},
	{"Param.ABM.CancerCell.moveProb_csc", "", "pr"},
	{"Param.ABM.CancerCell.Tkill_scaler", "", "pos"},
	{"Param.ABM.CancerCell.mincc", "", "pos"},
	{"Param.ABM.CancerCell.C1_CD47", "CD47 expression on cancer cells, molecule/micrometer^2", "pos"},
	{"Param.ABM.CancerCell.IFNgUptake", "per sec", "pos"},
	{"Param.ABM.CancerCell.hypoxia_th", "ng/mL", "pos"},
	{"Param.ABM.CancerCell.density_csc", "cancer cell IC density in core", "prob"},
	//Vas
	{"Param.ABM.Vas.maxPerVoxel", "", "pos"},
	{"Param.ABM.Vas.vas_50", "half maximum of VEGFA of induction of angiogensis pg/ml", "pos"},
	{"Param.ABM.Vas.O2_conc", "", "pos"},
	{"Param.ABM.Vas.Rc", "radius of tumor capillary", "pos"},
	{"Param.ABM.Vas.sigma", "dimensionless ratio of intracapillary to extracapillary transport resistance from Sharan et al; D*alpha / k*Rc = 0.84", "pos"},
	//reference: Assessment of Blood Flow in Hepatocellular Carcinoma Correlations of Computed Tomography Perfusion Imagingan Circulating Angiogenic Factors
	//estimated: 20 / 0.62 mm^2
	{"Param.ABM.Vas.ref_vas_frac", "reference fraction of endotheilal cell per voxel", "pos"},
	//agent chemokine interaction
	{"Param.ABM.cell.PDL1_th", "percent of max PDL1_syn to be detectable", "prob"},
	{"Param.ABM.cell.IFNg_PDL1_half", "c of IFNg to induce PDL1 to half maximal level, ng/mL", "pos"},
	{"Param.ABM.cell.IFNg_PDL1_n", "hill coef for PDL1 expression", "pos"},
	/* molecular level */
	// diffusion grid
	{"Param.Molecular.biofvm.IFNg.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IFNg.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.IFNg.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.IL_2.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IL_2.release","ng/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.IL_2.uptake","1/s", "pos"},
	{"Param.Molecular.biofvm.IL_2.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.CCL2.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.CCL2.release","moles/day, one cell", "pos"},
	{"Param.Molecular.biofvm.CCL2.uptake","1/s", "pos"},
	{"Param.Molecular.biofvm.CCL2.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.CCL2.molecularWeight","kDa", "pos"},		
	{"Param.Molecular.biofvm.ArgI.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.ArgI.release","mU * Liter/sec, one cell", "pos"},
	{"Param.Molecular.biofvm.ArgI.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.ArgI.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.NO.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.NO.release","mg/sec", "pos"},
	{"Param.Molecular.biofvm.NO.decayRate","1/sec", "pos"},	
	{"Param.Molecular.biofvm.NO.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.TGFB.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.TGFB.release.Mac","moles/day, one cell", "pos"},
	{"Param.Molecular.biofvm.TGFB.release.Treg","moles/day, one cell", "pos"},
	{"Param.Molecular.biofvm.TGFB.release.CancerStem","mg/sec", "pos"},
	{"Param.Molecular.biofvm.TGFB.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.TGFB.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.IL10.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IL10.release.Treg","mg/sec", "pos"},
	{"Param.Molecular.biofvm.IL10.release.Mac","ng/sec", "pos"},
	{"Param.Molecular.biofvm.IL10.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.IL10.molecularWeight","kDa", "pos"},
	{"Param.Molecular.biofvm.IL12.diffusivity","cm^2/sec", "pos"},
	{"Param.Molecular.biofvm.IL12.release","ng/sec", "pos"},
	{"Param.Molecular.biofvm.IL12.decayRate","1/sec", "pos"},
	{"Param.Molecular.biofvm.IL12.molecularWeight","kDa", "pos"},
	{ "Param.Molecular.biofvm.VEGFA.diffusivity","cm^2/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.release.Mac","ng/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.release.CancerStem", "ng/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.release.CancerProgenitor", "ng/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.uptake", "1/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.decayRate","1/sec", "pos" },
	{ "Param.Molecular.biofvm.VEGFA.molecularWeight","kDa", "pos" },
	{ "Param.Molecular.biofvm.O2.diffusivity","cm^2/sec", "pos" },
	{ "Param.Molecular.biofvm.O2.uptake", "mg/sec", "pos" },
	{ "Param.Molecular.biofvm.O2.decayRate","1/sec", "pos" },
	{ "Param.Molecular.biofvm.O2.molecularWeight","kDa", "pos" },

	//------------------------ int ----------------------------//
	{"Param.ABM.Environment.Tumor.XSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.YSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.ZSize", "", "pos"},
	{"Param.ABM.Environment.Tumor.VoxelSize", "voxel resolution, microns", "pos"},
	{"Param.ABM.Environment.Tumor.nr_T_voxel", "", "pos"},
	{"Param.ABM.Environment.Tumor.nr_T_voxel_C", "", "pos"},
	{"Param.ABM.Environment.Tumor.stem_mode", "", "pos" },
	{"Param.ABM.Environment.ShuffleInterval", "", "pos"},
	{"Param.ABM.Environment.gridshiftInterval", "", "pos"},
	{"Param.ABM.TCell.div_interval", "", "pos"},
	{"Param.ABM.TCell.div_limit", "", "pos"},
	{"Param.ABM.TCD4.div_interval", "", "pos"},
	{"Param.ABM.TCD4.div_limit", "", "pos"},		
	{"Param.ABM.CancerCell.progenitorDivMax", "", "pos"},
	{"Param.Molecular.stepPerSlice","", "pos"},

	// ---------------------- bool -----------------------------//
	{"Param.QSP.simulation.use_resection", "", "" },
	{"Param.Molecular.allMolecularOff", "", ""},
	{"Param.Molecular.diffusionOff", "", ""},
	{"Param.Molecular.cellOdeOff", "", ""},
	{"Param.ABM.Pharmacokinetics.nivoOn", "", "" },
	{"Param.ABM.Pharmacokinetics.ipiOn", "", "" },
	{"Param.ABM.Pharmacokinetics.entOn", "", "" },
	{"Param.ABM.Pharmacokinetics.caboOn", "", "" },

};

Param::Param()
	:ParamBase()
{
	setupParam();
}

/*! Setup parameter storage
	instantiation of pure virtual member of the base class.
	setup description vector;
	initialize parameter value vectors with 0/false, 
	with size determined by enums, 
	so that other base class members can access vector sizes
*/
void Param::setupParam(){

	size_t nrExternalParam = PARAM_FLOAT_COUNT + PARAM_INT_COUNT + PARAM_BOOL_COUNT;
	for (size_t i = 0; i < nrExternalParam; i++)
	{
		_paramDesc.push_back(std::vector<std::string>(_description[i], 
			_description[i]+ PARAM_DESCRIPTION_FIELD_COUNT));
	}
	_paramFloat = std::vector<double>(PARAM_FLOAT_COUNT, 0);
	_paramInt= std::vector<int>(PARAM_INT_COUNT, 0);
	_paramBool= std::vector<bool>(PARAM_BOOL_COUNT, false);
	_paramFloatInternal = std::vector<double>(PARAM_FLOAT_INTERNAL_COUNT, 0);
	_paramIntInternal = std::vector<int>(PARAM_INT_INTERNAL_COUNT, 0);
	_paramBoolInternal = std::vector<bool>(PARAM_BOOL_INTERNAL_COUNT, false);
}

/*! Calculate internal parameters
*/
void Param::processInternalParams(){
	
	_paramFloatInternal[PARAM_AVOGADROS] = AVOGADROS;

	//micrometer to cm
	_paramFloatInternal[PARAM_VOXEL_SIZE_CM] = _paramInt[PARAM_VOXEL_SIZE] / 1e4;

	/*
	_paramFloatInternal[PARAM_T_CELL_LIFE_MEAN_SLICE] = _paramFloat[PARAM_T_CELL_LIFE_MEAN]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;
	*/
	_paramFloatInternal[PARAM_T_CELL_LIFE_SD_SLICE] = _paramFloat[PARAM_T_CELL_LIFE_SD]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;

	_paramBoolInternal[PARAM_MOLECULAR_MODULES_ON]
		= !getVal(PARAM_ALL_MOLECULAR_OFF);

	_paramBoolInternal[PARAM_DIFFUSION_ON] 
		= getVal(PARAM_MOLECULAR_MODULES_ON) && !getVal(PARAM_DIFFUSION_OFF);
	
	_paramBoolInternal[PARAM_CELL_ODE_ON] 
		= getVal(PARAM_MOLECULAR_MODULES_ON) && !getVal(PARAM_ALL_CELL_ODE_OFF);

}

//! update from QSP parameters
void Param::update_from_qsp(void){

	// cell (need to be 1,  not in mole)
	_paramFloatInternal[PARAM_CELL] = QP(31) * AVOGADROS;

	// minimum tumor volume, V_Tmin
	_paramFloatInternal[PARAM_VT_MIN] = QP(37);

	// volume of cancer cells, vol_cell
	_paramFloatInternal[PARAM_V_CELL] = QP(35);

	// volume of T cells, vol_Tcell
	_paramFloatInternal[PARAM_V_TCELL] = QP(36);

	// surface area of T cells, area_Tcell
	_paramFloatInternal[PARAM_A_TCELL] = QP(104);

	// initial tumor volume calculated by initial tumor diameter.
	_paramFloatInternal[PARAM_INIT_TUM_VOL] = (PI * std::pow(QP(44), 3)) / 6;

	_paramFloatInternal[PARAM_INIT_TUM_DIAM] = QP(44);

	std::cout << "cell: " << _paramFloatInternal[PARAM_CELL]
		<< ", minimum tumor size : " << _paramFloatInternal[PARAM_VT_MIN]
		<< ", cancer cell size: " << _paramFloatInternal[PARAM_V_CELL]
		<< ", T cell size: " << _paramFloatInternal[PARAM_V_TCELL] << std::endl;
//  ==========================================================
//|| Treg and MDSC does not have maximum number this version  ||
//  ===========================================================
	// maximum concentration of Tregs per volume in the tumor
	//_paramFloatInternal[PARAM_TREGMAX] = QP(178);		

	// maximum concentration of MDSC per volume
	//_paramFloatInternal[PARAM_MDSCMAX] = QP(177);


	//==============================================
	//  MDSC module are only in ABM  in this version
	//*================================================
	 
	// MDSC base recruitment
	//_paramFloatInternal[PARAM_K_BASE_REC] = QP(180);

	// MDSC recruitment by CCL2
	//_paramFloatInternal[PARAM_K_REC_BY_CCL2] = QP(162);		

	// half maximal effective concentration of CCL2 on recruitment of MDSC into the tumor (ng/ml)
	_paramFloat[PARAM_MDSC_EC50_CCL2_REC] =_paramFloat[PARAM_MDSC_EC50_CCL2_REC] * _paramFloat[PARAM_CCL2_MOLECULAR_WEIGHT] * 1e3 / 1e3;

	// half maximal inhibitory concentration of NO on inhibition of CD8+ T cell cytotoxic activity (ng/ml)
	_paramFloat[PARAM_MDSC_IC50_NO_CTL] = _paramFloat[PARAM_MDSC_IC50_NO_CTL] * _paramFloat[PARAM_NO_MOLECULAR_WEIGHT] * 1e3 / 1e3;

	// half maximal inhibitory concentration of Arg I on inhibition of CD8+ T cell cytotoxic activity (ng/ml) 
	//_paramFloat[PARAM_MDSC_IC50_ArgI_CTL] = _paramFloat[PARAM_MDSC_IC50_ArgI_CTL] * _paramFloat[PARAM_ARGI_MOLECULAR_WEIGHT] * 1e3 / 1e3;

	// half maximal effective concentration of arginase I on Treg expansion (ng/ml) (molecular weights in kDa
	//_paramFloat[PARAM_MDSC_EC50_ArgI_Treg] = _paramFloat[PARAM_MDSC_EC50_ArgI_Treg] * _paramFloat[PARAM_ARGI_MOLECULAR_WEIGHT] * 1e3 / 1e3;
	
	std::cout << "PARAM_MDSC_EC50_CCL2_REC: " << _paramFloat[PARAM_MDSC_EC50_CCL2_REC]
		<< ", PARAM_MDSC_IC50_NO_CTL: " << _paramFloat[PARAM_MDSC_IC50_NO_CTL]
		<< ", PARAM_MDSC_IC50_ArgI_CTL: " << _paramFloat[PARAM_MDSC_IC50_ArgI_CTL]
		<< ", PARAM_MDSC_EC50_ArgI_Treg: " << _paramFloat[PARAM_MDSC_EC50_ArgI_Treg] << std::endl;

	// half maximal inhibitory concentration of entinostat on anti-proliferation of tumor cell (mol/m^3)
	//_paramFloatInternal[PARAM_IC50_ENT_C] = QP(164);	

	// half maximal inhibitory concentration of entinostat on inhibition of CCL2 production (mol/m^3)
	//_paramFloatInternal[PARAM_IC50_ENT_CCL2] = QP(179);	

	// half maximal inhibitory concentration of entinostat on inhibition of NO production (mol/m^3)
	//_paramFloatInternal[PARAM_IC50_ENT_NO] = QP(171);

	// half maximal inhibitory concentration of entinostat on inhibition of Arg I production (mol/m^3)
	//_paramFloatInternal[PARAM_IC50_ENT_ARGI] = QP(181);	
	
	// Area of synapse
	_paramFloatInternal[PARAM_A_SYN] = QP(103);
	// number of PD1/PDL1 binding for half maximal inhibition, PD1_50
	// Change the PD1_50 value by 10 folds, since it's a fitted value in QSP model
	_paramFloatInternal[PARAM_PD1_PDL1_HALF] = QP(158) / 10;

	// total number of PD1 per synapse on T cell = T_PD1_total*A_syn; T_PD1_total in density (molecule per micrometer^2)
	_paramFloatInternal[PARAM_PD1_SYN] = QP(132) * QP(104);
	// total number of PD1 per synapse on Macrophages = T_PD1_total*A_syn; T_PD1_total in density (molecule per micrometer^2)
	// _paramFloat[PARAM_MAC_PD1_AREA] in molecule / micrometer^2, 1m^2 = 1e12 micrometer^2
	_paramFloatInternal[PARAM_MAC_PD1_SYN] = (_paramFloat[PARAM_MAC_PD1_TOTAL] / AVOGADROS) / (_paramFloat[PARAM_MAC_AREA] / 1e12) * (_paramFloat[PARAM_MAC_AREA_SYN] / 1e12);

	// _paramFloat[PARAM_C1_CD47_SYN] in molecule / micrometer^2, 1m^2 = 1e12 micrometer^2
	_paramFloatInternal[PARAM_C1_CD47_SYN] = (_paramFloat[PARAM_C1_CD47] / AVOGADROS * 1e12) ;

	// _paramFloat[PARAM_MAC_SIRPa] in molecule / micrometer^2, mole/m^2 1m^2 = 1e12 micrometer^2
	_paramFloatInternal[PARAM_MAC_SIRPa_SYN] = (_paramFloat[PARAM_MAC_SIRPa] / AVOGADROS * 1e12);

	// Number of PDL1 per synapse = C1_PDL1_total*A_syn; C1_PDL1_total in density (molecule per micrometer^2)
	_paramFloatInternal[PARAM_PDL1_SYN_MAX] = QP(133) * QP(104);
	
	// total number of PDL1 per cancer cell = C1_PDL1_total*A_Cancer_cell; C1_PDL1_total in density (molecule per micrometer^2)
	_paramFloatInternal[PARAM_PDL1_CELL] = QP(133) * QP(105);

	// k1 for PDL1-PD1 calculation =  kon_PD1_PDL1 / (koff_PD1_PDL1* A_syn)
	_paramFloatInternal[PARAM_PDL1_K1] = QP(155) / (QP(156) * QP(103));

	// k2 for PDL1-PD1 calculation = 2* kon_PD1_aPD1 / (koff_PD1_aPD1 * gamma_T_nivo)
	_paramFloatInternal[PARAM_PDL1_K2] = 2 * QP(155)  / (QP(156) * QP(130));

	// k3 for PDL1-PD1 calculation = (Chi_PD1 * kon_PD1_aPD1) / (2 * koff_PD1_aPD1)
	_paramFloatInternal[PARAM_PDL1_K3] = QP(157) * QP(155) / (2 * QP(156));
	// hill coefficient
	_paramFloatInternal[PARAM_N_PD1_PDL1] = QP(159);

	
	std::cout << "k1, k2, k3, T1, PDL1_tot" 
		<< ": " << _paramFloatInternal[PARAM_PDL1_K1]
		<< ", " << _paramFloatInternal[PARAM_PDL1_K2]
		<< ", " << _paramFloatInternal[PARAM_PDL1_K3]
		<< ", " << _paramFloatInternal[PARAM_PD1_SYN]
		<< ", " << _paramFloatInternal[PARAM_MAC_PD1_SYN]
		<< ", " << _paramFloatInternal[PARAM_PDL1_SYN_MAX]
		<< std::endl;
	std::cout << "k50: " << _paramFloatInternal[PARAM_PD1_PDL1_HALF]<< std::endl;
	
	// Binding rate between kon_CTLA4 and ipi
	//_paramFloatInternal[PARAM_KON_CTLA4_IPI] = QP(24);

	// Unbinding rate between CTLA4 and ipi
	//_paramFloatInternal[PARAM_KOFF_CTLA4_IPI] = QP(150);
	_paramFloatInternal[PARAM_KOFF_CTLA4_IPI] = 6.96e-06;

	// Volume fraction available to ipi in tumor compartment
	//_paramFloatInternal[PARAM_GAMMA_T_IPI] = QP(121);
	_paramFloatInternal[PARAM_GAMMA_T_IPI] = 0.718;

	// Antibody cross-arm binding efficiency  that also includes the conversion of kon from 3D to 2D (estimated)
	//_paramFloatInternal[PARAM_CHI_CTLA4_IPI] = QP(151);
	_paramFloatInternal[PARAM_CHI_CTLA4_IPI] = 3;

	// CTLA4 occupancy for half-maximal Treg inactivation by macrophages (estimated)  = CTLA_50 * A_treg
	_paramFloatInternal[PARAM_TREG_CTLA4_50] = QP(163) * QP(104);

	// total number of CTLA4 per synapse = Treg_CTLA4_total * A_treg; T_CTLA_total in density (molecule per micrometer^2)
	_paramFloatInternal[PARAM_CTLA4_TREG] = QP(140) * QP(104);

	// Anti-CTLA4 ADCC (antibody-dependent cellular cytotoxicity) rate of Treg (Richards 2008, PMID: 18723496)
	//_paramFloatInternal[PARAM_K_ADCC] = QP(159) * _paramFloat[PARAM_SEC_PER_TIME_SLICE];
	_paramFloatInternal[PARAM_K_ADCC] = 0.1 * _paramFloat[PARAM_SEC_PER_TIME_SLICE];
	/*
	std::cout << "kon_CTLA_ipi: " << _paramFloatInternal[PARAM_KON_CTLA4_IPI]
		<< ", koff_CTLA_ipi: " << _paramFloatInternal[PARAM_KOFF_CTLA4_IPI]
		<< ", gamma_T_ipi: " << _paramFloatInternal[PARAM_GAMMA_T_IPI]
		<< ", chi_CTLA_ipi: " << _paramFloatInternal[PARAM_CHI_CTLA4_IPI]
		<< ", CTLA4 total in Treg: " << _paramFloatInternal[PARAM_CTLA4_TREG] << std::endl;
	*/
	// Parameters calculated from QSP parameter values
	double t_step_sec = _paramFloat[PARAM_SEC_PER_TIME_SLICE];

	// Update cabozantinib module from QSP model (values are temperary)
	// IC50	for receptors inhibited by cabozantinib
	_paramFloatInternal[PARAM_IC50_AXL] = QP(185);
	_paramFloatInternal[PARAM_IC50_VEGFR2] = QP(186);
	_paramFloatInternal[PARAM_IC50_MET] = QP(183);
	_paramFloatInternal[PARAM_IC50_RET] = QP(184);

	// theraputic effects parameters elicited by cabozantinib
	_paramFloatInternal[PARAM_LAMBDA_C_CABO] = QP(190);
	_paramFloatInternal[PARAM_LAMBDA_V_CABO] = 0.025; //QP(191) * t_step_sec;
	_paramFloatInternal[PARAM_LAMBDA_Q_CABO] = QP(192);
	// T cell killing of Cancer cell
	// QP(48): k_C_death_by_T (day^-1, sec^-1 internal)
	// Becuase the ABM contain modules that QSP, which contains extra immunesuppresive modules
	// 5 is a arbitary coefficient to compensate immuno-suppresive module not present in QSP model
	_paramFloatInternal[PARAM_ESCAPE_BASE] = std::exp(-t_step_sec * QP(66));
	std::cout << "increased TK: " << _paramFloat[PARAM_ABM_KILL_SCALE] << std::endl;
	// Macrophage killing of Cancer cell
	//_paramFloat[PARAM_MAC_K_M1_PHAGO] k_C_death_by_Macrophage (day^-1, sec^-1 internal)
	// Becuase the ABM contain modules that QSP, which contains extra immunesuppresive modules
	// 5 is a arbitary coefficient to compensate immuno-suppresive module not present in QSP model
	_paramFloatInternal[PARAM_ESCAPE_MAC_BASE] = std::exp(-t_step_sec * _paramFloat[PARAM_MAC_K_M1_PHAGO] / SEC_PER_DAY);

	// T cell exhaustion from PDL1; 2 is arbitary
	_paramFloatInternal[PARAM_EXHUAST_BASE_PDL1] = std::exp(-t_step_sec * QP(65));

	// T cell exhaustion from Treg inhibition, k_Treg ; 50 is arbitary
	_paramFloatInternal[PARAM_EXHUAST_BASE_TREG] = std::exp(-t_step_sec * QP(67) / 50);

	// rate of Th to Treg transformation , units: 1/days -> 1/timestep
	_paramFloat[PARAM_K_TH_TREG] = _paramFloat[PARAM_K_TH_TREG] / SEC_PER_DAY * t_step_sec;
	// Macrophage module related parameter covnersion
	// Rate of M1 to M2 macrophage polarization, units: 1/days -> 1/timestep
	_paramFloat[PARAM_MAC_M2_POL] = _paramFloat[PARAM_MAC_M2_POL] / SEC_PER_DAY * t_step_sec;
	// Rate of M2 to M1 macrophage polarization, units: 1/days -> 1/timestep
	_paramFloat[PARAM_MAC_M1_POL] = _paramFloat[PARAM_MAC_M1_POL] / SEC_PER_DAY * t_step_sec;
	// TGFb_50 Half-Maximal TGFb level for Th-to-Treg differentiation / chemoresistance development / M1-to-M2 polarization
    //unit: nanomolarity (1e-9 mole/L) -> ng/ml
	_paramFloat[PARAM_MAC_TGFB_EC50] = _paramFloat[PARAM_MAC_TGFB_EC50] * _paramFloat[PARAM_TGFB_MOLECULAR_WEIGHT] * 1e3 / 1e3;
	//unit: nanomolarity (1e-9 mole/L) -> ng/ml
	_paramFloat[PARAM_TEFF_TGFB_EC50] = _paramFloat[PARAM_TEFF_TGFB_EC50] * _paramFloat[PARAM_TGFB_MOLECULAR_WEIGHT] * 1e3 / 1e3;
	//unit: picomolarity (1e-12 mole/L) -> ng/ml
	_paramFloat[PARAM_MAC_IL_10_EC50] = _paramFloat[PARAM_MAC_IL_10_EC50] * _paramFloat[PARAM_IL_10_MOLECULAR_WEIGHT] * 1e3 / 1e6;
	//unit: picomolarity (1e-12 mole/L) -> ng/ml
	_paramFloat[PARAM_MAC_IL_10_HALF_PHAGO] = _paramFloat[PARAM_MAC_IL_10_HALF_PHAGO] * _paramFloat[PARAM_IL_10_MOLECULAR_WEIGHT] * 1e3 / 1e6;
	//unit: picomolarity (1e-12 mole/L) -> ng/ml
	_paramFloat[PARAM_MAC_IFN_G_EC50] = _paramFloat[PARAM_MAC_IFN_G_EC50] * 55 * 1e3 / 1e3 / 1e3;
    //unit: picomolarity (1e-12 mole/L) -> ng/ml
	_paramFloat[PARAM_MAC_IL_12_EC50] = _paramFloat[PARAM_MAC_IL_12_EC50] * _paramFloat[PARAM_IL_12_MOLECULAR_WEIGHT] * 1e3 / 1e3 / 1e3;
	//unit: molecule/micrometer^2 -> moles/m^2
	_paramFloat[PARAM_MAC_SIRPa_HALF] = _paramFloat[PARAM_MAC_SIRPa_HALF] / AVOGADROS / 1e-12;
	// 1/(micromolarity*minute*nanometer) ->   0.001 (micromolarity to mole/m^3), 60 (minute to second), 1e-9 (nanometer to meter)
	_paramFloat[PARAM_KON_SIRPa_CD47] = _paramFloat[PARAM_KON_SIRPa_CD47] / (0.001 * 60 * 1e-9);
	// 1/minute
	_paramFloat[PARAM_KOFF_SIRPa_CD47] = _paramFloat[PARAM_KOFF_SIRPa_CD47] / 60;

	std::cout << "M1 to M2 Mac polarization: " << _paramFloat[PARAM_MAC_M2_POL]
		<< ", M2 to M1 Mac polarization: " << _paramFloat[PARAM_MAC_M1_POL]
		<< ", EC50 TGFB: " << _paramFloat[PARAM_MAC_TGFB_EC50] 
		<< ", PARAM_TEFF_TGFB_EC50" << _paramFloat[PARAM_TEFF_TGFB_EC50]
		<< ", PARAM_MAC_IL_10_EC50: " << _paramFloat[PARAM_MAC_IL_10_EC50]
		<< ", PARAM_MAC_IFN_G_EC50: " << _paramFloat[PARAM_MAC_IFN_G_EC50] << std::endl;
	// time for resection
	_paramFloatInternal[PARAM_RESECT_TIME_STEP] = _paramFloat[PARAM_QSP_T_RESECTION] * SEC_PER_DAY / t_step_sec;


	// Recruitment
	// T effector recruitment
	/* for each mole of adhesion site, the amount of T cell recruited (in unit of mole):
	dt * q_T1_T_in * Cent.T * V_T
	# Weighted QSP:
	# Units:
	All parameters come in SI units.
	Cent.T should use SI units --- in mole

	The result is mole recruited per mole site, or number per site, so no conversion needed.
	When calculating recruitment probability:
	p = 1 / (s * m^3) * cell  * (m^3 / cell) * dt
	*/

	//The number of adhesion site per voxel is:
	double site_per_voxel = _paramFloat[PARAM_ADH_SITE_DENSITY] * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]), 3) * _paramFloat[PARAM_VAS_REF_FRAC];
	//number of adhesion sites needed to recruit a single cell (becoming a recruitment port).
	double site_per_port = _paramFloat[PARAM_REC_SITE_FACTOR];
	//how many port per voxel have
	_paramFloatInternal[PARAM_REC_PORT] = site_per_voxel / site_per_port;
	std::cout << "PARAM_REC_PORT_PROB : " << _paramFloatInternal[PARAM_REC_PORT] << std::endl;
	/*When calculating recruitment probability:
	*/
	double  w = _paramFloat[PARAM_WEIGHT_QSP];
	// Teff -> k (1/mol) // p = k (1/mol) * Cent.T (mol), q_T1_T_in
	//_paramFloatInternal[PARAM_TEFF_RECRUIT_K] = QP(36) * site_per_port * t_step_sec  / w / _paramFloat[PARAM_ADH_SITE_DENSITY];
	// 2 is the arbitary number
	_paramFloatInternal[PARAM_TEFF_RECRUIT_K] = QP(54) *  t_step_sec * AVOGADROS * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]) / 1e6, 3) * _paramFloatInternal[PARAM_REC_PORT];

	// TCD4 -> k (1/mol) // p = k (1/mol) * Cent.T (mol), q_T0_T_in
	//_paramFloatInternal[PARAM_TCD4_RECRUIT_K] = QP(59) * site_per_port * t_step_sec  / w / _paramFloat[PARAM_ADH_SITE_DENSITY];
	_paramFloatInternal[PARAM_TCD4_RECRUIT_K] = QP(77) * t_step_sec * AVOGADROS * std::pow(double(_paramInt[PARAM_VOXEL_SIZE]) / 1e6, 3) * _paramFloatInternal[PARAM_REC_PORT];
	std::cout << "keff" << _paramFloatInternal[PARAM_TEFF_RECRUIT_K] << std::endl;
	std::cout << "kcd4" << _paramFloatInternal[PARAM_TCD4_RECRUIT_K] << std::endl;
	// The recruitment of MDSC does not depend on central compartment MDSC (central compartment does not have MDSC)
	// so the recruitment equation is different from Teff and Treg
	// 11000 * cells/(ml*day) -> cells / (m^3 * s)
	// 1 m^3 = 1e6 ml
	_paramFloatInternal[PARAM_MDSC_RECRUIT_BY_CCL2_K] = _paramFloat[PARAM_MDSC_K_REC] * t_step_sec / SEC_PER_DAY / 1e-6 * _paramFloat[PARAM_MDSC_MYELOID_SCALE];
	_paramFloat[PARAM_MAC_RECRUIT_BY_CCL2_K] = _paramFloat[PARAM_MAC_RECRUIT_BY_CCL2_K] * t_step_sec / SEC_PER_DAY / 1e-6 * _paramFloat[PARAM_MDSC_MYELOID_SCALE];

	// APC -> k (m^3/mol) // p = k (m^3/mol) * (APC0_T*V_T-V_T.APC)
	_paramFloatInternal[PARAM_APC_RECRUIT_K] = 0;
	
	// APC density in the tumour
	_paramFloatInternal[PARAM_APC0_T] = QP(86);

	// APC transmigration rate from tumor to lymph node
	_paramFloatInternal[PARAM_APC_TRANSMIG] = QP(83);

	// Number of T0 cell Clonality to lymph node
	_paramFloatInternal[PARAM_T0_CLONE] = QP(69);

	// Number of T1 cell Clonality to lymph node
	_paramFloatInternal[PARAM_T1_CLONE] = QP(46);

	// antigen concentration in cancer cell
	_paramFloatInternal[PARAM_ANTIGEN_PER_CELL] = QP(116);

	// antigen uptake rate by mAPC
	_paramFloatInternal[PARAM_ANTIGEN_UP] = QP(95);

	// antigen degrdation rate in the tumor
	_paramFloatInternal[PARAM_K_xP_DEG] = QP(110);
	
	// mean life of Tcell, unit: time step, 
	// divided by 100 since in ABM there is Tcell division in the tumor compartment, 
	// which is different from the QSP model.
	_paramFloatInternal[PARAM_T_CELL_LIFE_MEAN_SLICE] = 1 / QP(51) / t_step_sec;
	// mean life of TCD4, unit: time step, 10 is a arbitary number
	_paramFloatInternal[PARAM_TCD4_LIFE_MEAN] = 1 / QP(74) / t_step_sec / 5;
	std::cout << "T cell life mean: " << _paramFloatInternal[PARAM_T_CELL_LIFE_MEAN_SLICE]
		<< ", TCD4_LIFE_MEAN: " << _paramFloatInternal[PARAM_TCD4_LIFE_MEAN] << std::endl;
	// mean life of MDSC, unit: time step
	//_paramFloatInternal[PARAM_MDSC_LIFE_MEAN] = 1 / QP(163) / t_step_sec;
	// mean life of APC, unit: time step
	_paramFloatInternal[PARAM_APC_LIFE_MEAN] = 1 / QP(85) / t_step_sec;
	std::cout << "Internal param: " << QP(76)  <<  ", "  << getVal(PARAM_TCD4_LIFE_MEAN) << std::endl;
	
	// Cytokine concentration for half - maximal APC maturation
	_paramFloatInternal[PARAM_CKINE_50] = QP(90);
	// Maximum rate of APC maturation
	_paramFloatInternal[PARAM_K_APC_MAT] = QP(82);
	std::cout
	<< t_step_sec << ", " << QP(58) << ", " << QP(157) << "\n"
	<< "PARAM_ESCAPE_BASE, " << _paramFloatInternal[PARAM_ESCAPE_BASE] << "\n"
	<< "PARAM_EXHUAST_BASE_PDL1, " << _paramFloatInternal[PARAM_EXHUAST_BASE_PDL1] << "\n"
	<< "PARAM_EXHUAST_BASE_TREG, " << _paramFloatInternal[PARAM_EXHUAST_BASE_TREG] << "\n"
	<< std::endl;
	

	/*Cancer cell dynamics parameters*/

	// stem cell division rate is calculated from QSP parameter
	// unit: s^-1, k_C1_growth
	double rs = QP(42) / (1 - _paramFloat[PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB]);
	// unit: day^-1
	_paramFloatInternal[PARAM_CSC_GROWTH_RATE] = rs * SEC_PER_DAY;

	_paramFloatInternal[PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE]
		= std::log(2)/rs / getVal(PARAM_SEC_PER_TIME_SLICE);

	_paramFloatInternal[PARAM_CANCER_SENESCENT_MEAN_LIFE] =
		1 / _paramFloat[PARAM_CANCER_SENESCENT_DEATH_RATE]
		/ _paramFloat[PARAM_SEC_PER_TIME_SLICE] * SEC_PER_DAY;
		
	
	_paramFloatInternal[PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE]
		= std::log(2)/ getVal(PARAM_CANCER_PROG_GROWTH_RATE) 
		* SEC_PER_DAY / getVal(PARAM_SEC_PER_TIME_SLICE) + .5;
	
	//_paramFloatInternal[PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE]
		//= std::log(2) / rs / getVal(PARAM_SEC_PER_TIME_SLICE);
	
	std::cout << "PARAM_CSC_GROWTH_RATE: "<< getVal(PARAM_CSC_GROWTH_RATE) << ", k_C1_growth: " << QP(42)
		<< ", PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE: " << _paramFloatInternal[PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE] << std::endl;
	/*
	std::cout << getVal(PARAM_INT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE)
		<< ", "<< getVal(PARAM_CANCER_SENESCENT_MEAN_LIFE)
		<< ", "<< getVal(PARAM_INT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE)
		<< std::endl;
	*/

	return;
}

};
};