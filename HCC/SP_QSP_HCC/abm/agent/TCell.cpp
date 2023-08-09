//#include <boost/serialization/export.hpp>
#include "TCell.h"

//BOOST_CLASS_EXPORT_IMPLEMENT(TCell)
#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"


namespace SP_QSP_IO{
namespace SP_QSP_HCC{

//#define IFNG_ID 0

using std::string;
using std::stringstream;

static int TCellSize = 1;

TCell::TCell()
	: _count_neighbor_cancer(0)
	, _count_neighbor_Treg(0)
	, _count_neighbor_all(0)
	, _max_neighbor_PDL1(0)
{
}

TCell::TCell(SpatialCompartment* c)
	:Cell_Tumor(c)
	, _divide_flag(false)
	, _divide_cd(0)
	, _divide_limit(params.getVal(PARAM_T_DIV_LIMIT))
	, _IL2_exposure(0)
	, _IL2_release_remain(params.getVal(PARAM_IL_2_RELEASE_TIME))
	, _source_IFNg(NULL)
	, _source_IL_2(NULL)
	, _sink_IL_2(NULL)
	, _count_neighbor_cancer(0)
	, _count_neighbor_Treg(0)
	, _count_neighbor_all(0)
	, _max_neighbor_PDL1(0)
{
	_state = AgentStateEnum::T_CELL_EFF;
	_life = TCell::getTCellLife();
	//cout << getID() << ", " << getType() << endl;
}


TCell::TCell( const TCell & c)
	:Cell_Tumor(c)
	, _divide_flag(c._divide_flag)
	, _divide_cd(c._divide_cd)
	, _divide_limit(c._divide_limit)
	, _IL2_exposure(c._IL2_exposure)
	, _IL2_release_remain(params.getVal(PARAM_IL_2_RELEASE_TIME))
	, _source_IFNg(NULL)
	, _source_IL_2(NULL)
	, _sink_IL_2(NULL)
	, _count_neighbor_cancer(0)
	, _count_neighbor_Treg(0)
	, _count_neighbor_all(0)
	, _max_neighbor_PDL1(0)
{
	_state = c._state;
	//cout << getID() << ", " << getType() << endl;
	_life = TCell::getTCellLife();
}


TCell::~TCell()
{
}
string TCell::toString()const{
	stringstream ss;
	ss << CellAgent::toString();
	ss << "division flag: " << _divide_flag << ", division cool down: " 
		<< _divide_cd  << std::endl;
	return ss.str();
}

bool TCell::agent_movement_step(double t, double dt, Coord& c){
	bool move = false;
	if (rng.get_unif_01() < params.getVal(PARAM_T_CELL_MOVE_PROB))
	{
		// move
		int idx;
		const auto shape = getCellShape();
		if (_compartment->getOneOpenVoxel(shape->getMoveDestinationVoxels(), 
			shape->getMoveDirectionAnchor(), getCoord(), getType(), idx, rng))
		{
			move = true;
			c = getCellShape()->getMoveDirectionAnchor()[idx] + getCoord();
		}
	}
	return move;

}

/*! Scan neighborhood and count neighbor cell type
	number of Teff neighbor for Cancer cell; 
	number of Cancer cell neighbor for Teff; 
	number of Treg neighbor for Teff; 
*/
void TCell::agent_state_scan(void){

	const auto shape = getCellShape();
	if (_state == AgentStateEnum::T_CELL_CYT || _state == AgentStateEnum::T_CELL_EFF){
		/**/
		// scan cancer cells
		//std::cout << "T cell at: "<< getCoord() << std::endl;
		//int nr_cancer_neighbor = 0;
		_compartment->for_each_neighbor_ag(shape->getEnvironmentLocations(),
			getCoord(), [&](BaseAgent* ag){
			_count_neighbor_all++;
			auto pCell = dynamic_cast<Cell_Tumor*>(ag);
			update_neighbor_PDL1(pCell->get_PDL1());
			if (ag->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
				auto pCancer = dynamic_cast<CancerCell*>(ag);
				inc_neighbor_cancer();
				pCancer->inc_neighbor_Teff();
				//release TGFB when encounter T cells
				if (pCancer->getType() == AgentTypeEnum::CELL_TYPE_CANCER && pCancer->getState() == AgentStateEnum::CANCER_PROGENITOR) {
					double TGFB_pro_release = params.getVal(PARAM_CANCER_TGFB_RELEASE) / 100;
					pCancer->update_chem_source(pCancer->get_source_TGFB(), TGFB_pro_release);
				}
				//nr_cancer_neighbor += 1;
			}
			else if (ag->getType() == AgentTypeEnum::CELL_TYPE_TCD4 && ag->getState() == AgentStateEnum::TCD4_TREG){
				inc_neighbor_Treg();
			}
			return false;
		});
		//std::cout << "T cell neighbor PDL1: "<< _max_neighbor_PDL1 << std::endl;
		//std::cout << "T cell neighbor Cancer: "<< _count_neighbor_cancer << std::endl;
		//std::cout << "T cell neighbor Treg: " << _count_neighbor_Treg << std::endl;
	}
	return;
}

bool TCell::agent_state_step(double t, double dt, Coord& c){
	bool divide = false;
	if (!isDead())
	{
		_life--;
		if (_life <= 0)
		{
			setDead();
			// remove source when cell die
			return divide;
		}
	}

	const auto shape = getCellShape();

	Cell_Tumor::agent_state_step(t, dt, c);

	auto tumor = dynamic_cast<Tumor*>(_compartment);

	double IL2 = get_tumor().get_chem(getCoord(), CHEM_IL_2);
	_IL2_exposure += params.getVal(PARAM_SEC_PER_TIME_SLICE) * IL2;

	// effector cells to proliferate on IL2 exposure
	if (_IL2_exposure > params.getVal(PARAM_IL_2_PROLIF_TH))
	{
		//std::cout << "T cell divide due to IL_2 exposure. " << std::endl;
		_divide_flag = true;
		_IL2_exposure = 0;
	}

	// Look for tumor Ag
	if (_state == AgentStateEnum::T_CELL_EFF)
	{
		bool cancer_found = _compartment->hasTypeStateInTarget(shape->getEnvironmentLocations(), 
			getCoord(), AgentTypeEnum::CELL_TYPE_CANCER, AgentStateEnum::CANCER_PROGENITOR);

		//antigen presented by M1 macrophage
		//bool ag_found = _compartment->hasTypeStateInTarget(shape->getEnvironmentLocations(), getCoord(), AgentTypeEnum::CELL_TYPE_MAC, AgentStateEnum::MAC_M1);

		//std::cout << "Tcell detecting: " << _state << std::endl;
		if (cancer_found)
		{
			_state = AgentStateEnum::T_CELL_CYT;
			//std::cout << "Cancer found: "<<"Tcell status: " << _state << std::endl;
			_divide_flag = true;
			setup_chem_source(_source_IFNg, CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
			setup_chem_source(_source_IL_2, CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
			setup_chem_sink(_sink_IL_2, CHEM_IL_2, params.getVal(PARAM_IL_2_UPTAKE));
		}
	}

	if (_state == AgentStateEnum::T_CELL_CYT){
		/**/
		double nivo = tumor->get_Nivo();

		// kill one cancer cell
		// now handled from cancer cells

		// IL-2 release time limit
		if (_IL2_release_remain > 0)
		{
			_IL2_release_remain -= params.getVal(PARAM_SEC_PER_TIME_SLICE);
		}
		else{
			// set IL-2 source to 0
			update_chem_source(_source_IL_2, 0.0);
		}
		
		// exhaustion
		if (_count_neighbor_Treg > 0)// suppresion by Treg
		{
			double q = double(_count_neighbor_Treg) / (_count_neighbor_all + params.getVal(PARAM_CELL));
			double p_exhaust = get_exhaust_prob_Treg(q);

			/*
			std::cout << "RNG check: ID=" << getID() <<
				"(Teff exhaust by Treg): " << p_exhaust << std::endl;
			*/

			if (rng.get_unif_01() < p_exhaust)
			{

				set_suppressed();
			}
		}
		else if (_count_neighbor_all > 0)// suppresion by PDL1 
		{
			double bond = get_PD1_PDL1(_max_neighbor_PDL1, nivo);

			// PD1_50 is in molecule / um^2, therefore, the input to the hill function has to be in bond / area as well
			// multiplication is set get_PD1_supp
			double supp = get_PD1_supp(bond / params.getVal(PARAM_A_SYN), params.getVal(PARAM_N_PD1_PDL1));
			double p_exhaust = get_exhaust_prob_PDL1(supp, 1.0);
			
			/*
			std::cout << "RNG check: ID=" << getID() <<
				" bond : " << bond <<
				" supp :" << supp <<
				" (Teff exhaust by PDL1): " << p_exhaust << std::endl;
			*/
			
			if (rng.get_unif_01() < p_exhaust)
			{

				set_suppressed();
			}
		}
	}

	if (_divide_cd > 0)
	{
		_divide_cd--;
	}
	
	if (_divide_limit > 0 && _divide_flag && _divide_cd == 0) // && _state == AgentStateEnum::T_CELL_CYT)
	{
		int idx;
		if (_compartment->getOneOpenVoxel(shape->getProlifDestinationVoxels(), 
			shape->getProlifDestinationAnchor(), getCoord(), getType(), idx, rng))
		{
			divide = true;
			//cout << "idx: " << idx << ", " << getCellShape()->getProlif()[idx] << endl;
			c = getCellShape()->getProlifDestinationAnchor()[idx] + getCoord();

			_divide_flag = true;
			_divide_limit -= 1;
			_divide_cd = params.getVal(PARAM_T_DIV_INTERVAL);
		}

	}

	_count_neighbor_cancer = 0;
	_count_neighbor_Treg = 0;
	_count_neighbor_all = 0;
	_max_neighbor_PDL1 = 0;
	return divide;
}

/*! set to suppressed state
*/
void TCell::set_suppressed(void){
	//set the supressed T cell life to half when exhausted
	//_life = int(_life / 2 + 1);
	//std::cout << "suppressed T cell life: " << _life << std::endl;
	_divide_flag = false;
	_divide_limit = 0;
	_state = AgentStateEnum::T_CELL_SUPP;
	remove_source_sink(_source_IFNg);
	remove_source_sink(_source_IL_2);
	remove_source_sink(_sink_IL_2);
	return;
}

//! move sources (IFN and IL2)
void TCell::move_all_source_sink(void)const{
	//std::cout << "moving sources: " << getCoord() << std::endl;
	move_source_sink(_source_IFNg);
	move_source_sink(_source_IL_2);
	move_source_sink(_sink_IL_2);
	return;
}

//! remove sources (IFN and IL2)
void TCell::remove_all_source_sink(void){
	remove_source_sink(_source_IFNg);
	remove_source_sink(_source_IL_2);
	remove_source_sink(_sink_IL_2);
	return;
}

void TCell::update_neighbor_PDL1(double PDL1){
	if (PDL1 > _max_neighbor_PDL1)
	{
		_max_neighbor_PDL1 = PDL1;
	}
	return;
}

//! PD1_PDL1 bond in synapse
double TCell::get_PD1_PDL1(double PDL1, double Nivo){

	static double T1 = params.getVal(PARAM_PD1_SYN);
	double T2 = PDL1;
	static double k1 = params.getVal(PARAM_PDL1_K1);
	static double k2 = params.getVal(PARAM_PDL1_K2);
	static double k3 = params.getVal(PARAM_PDL1_K3);
	/*
	double c = T1 + PDL1 + (1 + k2 * Nivo) / k1;
	double discriminant = c * c - 4 * T1 * PDL1;
	double bond = .5*(c - std::sqrt(discriminant));
	*/

	double a = 1;
	double b = (Nivo*k2/k1*(2*k3/k1-1) - 2*T2 - T1 - 1/k1)/T2;
	double c = (Nivo*k2/k1 + 1/k1  +T2 + 2*T1 )/T2;
	double d = -T1/T2;

	/*
	std::cout << k1 << "," << k2 << "," << k3 << ","
		<< T1 << "," << T2 << std::endl;
	std::cout << a << "," << b << "," << c << "," << d << std::endl;
	*/

	//Newton_Raphson_root
	int max_iter = 20;
	double tol_rel = 1E-5;
	double root = 0;
	double res, root_new, f, f1;
	int i = 0;
	while (i < max_iter){
		f = a*std::pow(root, 3) + b*std::pow(root, 2)+ c*root + d;
		f1 = 3.0*a*std::pow(root, 2) + 2.0*b*root + c;
		root_new = root - f/f1;
		res = std::abs(root_new - root) / root_new;
		if (res > tol_rel){
			i++;
			root = root_new;
		}
		else{
			break;
		}
	}
	
	/*
	std::cout << "PDL1: " << PDL1 << ", nivo: " << Nivo 
		<< ", root: " << root << std::endl;
	*/

	double bond = T2*root;
	//std::cout << "bond : " << bond << std::endl;

	return bond;
}


//! get suppression from PD1_PDL1 bond
double TCell::get_PD1_supp(double concentration, double n) {
	// bond is in molecule (moles), therefore, the EC50 to the hill function has to be in bond as well
	// params.getVal(PARAM_PD1_PDL1_HALF) is in bonds/ um^2. 
    // PD1_bond is in bonds Need to divided by params.getVal(PARAM_A_SYN)
	double k50 = params.getVal(PARAM_PD1_PDL1_HALF);
	return get_Hill_equation(concentration, k50, n);
}



int TCell::getTCellLife(){
	double lifeMean = params.getVal(PARAM_T_CELL_LIFE_MEAN_SLICE);
	double lifeSd = params.getVal(PARAM_T_CELL_LIFE_SD_SLICE);
	double tLifeD = lifeMean + rng.get_norm_std() * lifeSd;

	int tLife = int(tLifeD + 0.5);
	tLife = tLife > 0 ? tLife : 0;
	//std::cout << "Base T cell life: " << lifeMean << std::endl;
	//std::cout << "Tcyt T cell life: " << tLife << std::endl;
	return tLife;
	
}

/*! probability of Cancer cell getting killed by Teff
	\param[in] supp: suppression
	p_kill = 1 - B_esc^((1-supp)*q),
	where B_esc = exp(-tstep * k_C_death_by_T)
*/
double TCell::get_kill_prob(double supp, double q){
	double T_cell_killing_rate = params.getVal(PARAM_ESCAPE_BASE);
	return 1 - std::pow(T_cell_killing_rate, q*(1-supp));
}

/*! probability of turning exhausted from PDL1-PD1 interaction
	p_exhaust_PDL1 = 1 - B^(supp*q)
*/
double TCell::get_exhaust_prob_PDL1(double supp, double q){
	return 1 - std::pow(params.getVal(PARAM_EXHUAST_BASE_PDL1), q*supp);
}

/*! probability of turning exhausted from Treg
	p_exhaust_Treg = 1 - B^q
*/
double TCell::get_exhaust_prob_Treg(double q){
	return 1 - std::pow(params.getVal(PARAM_EXHUAST_BASE_TREG), q);
}


};
};
