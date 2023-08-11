//#include <boost/serialization/export.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>
#include "Tumor.h"
#include "LymphCentral.h"
//BOOST_CLASS_EXPORT_IMPLEMENT(Tumor);

//#include "../agent/CancerCell.h"`
#include <algorithm>
#include <exception>
#include <numeric>
#include <queue>
#include <cmath>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include "../../core/GlobalUtilities.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//#include "../../MemoryUsage.h"


Tumor::Tumor(int x, int y, int z)
	: SpatialCompartment(x, y, z)
	, _stats()
	, _voxel_ic()
	, _chem()
	, _init_vas()
	, _t_source()
	, _mdsc_source()
	, _mac_source()
	, _tInitDummy(new TCell(this))
	, _cInitDummy(new CancerCell(this))
	, _macInitDummy(new Mac(this))
	, _fibInitDummy(new Fib(this))
	, _TCD4InitDummy(new TCD4(this))
	, _MDSCInitDummy(new MDSC(this))
	, _VasInitDummy(new Vas(this))
	//, _voxel_size(params.getVal(PARAM_VOXEL_SIZE))
	, _allow_shift_grid(false)
	, _center_target()
	, _agGrid_temp(AgentGrid(x, y, z, NULL))
	, _var_abm_to_qsp()
	, _concentration_cc(0)
	, _concentration_t_cyt(0)
	, _concentration_t_CD4(0)
	, _concentration_t_eff_tum(0)
	, _concentration_t_CD4_tum(0)
	, _concentration_mdsc(0)
	, _concentration_nivo(0)
	, _concentration_ipi(0)
	, _concentration_ent(0)
	, _concentration_cabo(0)
	, _concentration_cx(0)
	, _concentration_t_exh(0)
	, _concentration_mat_ckine(0)
	, _tumor_volume(0)
	, _f_tum_cap(0)
	, _R_cabo(0)
	, _cancer_counts(0)
{
	//CC total, CC death total, CC death Teff, Teff recruit, TCD4 recruit
	_var_abm_to_qsp = std::vector<double>(TUMEX_VAR_NUM, 0);

	// dummy cells do not act as source
	//_tInitDummy->get_source_IFNg()->set_remove();
	initAgentGrid();
}

/*!
  _agGrid uses pointer to AgentGridVoxel (or its derived classes)
  as its element. This has a higher computational cost than directly
  hosting the instantications of the derived class (base won't work).
  The benefit is that this way the voxel objects can be accesses from
  the SpatialCompartment abstract class.
*/
bool Tumor::initAgentGrid() {
	//std::cout << "init agGrid:" << std::endl;
	for (int i = 0; i < _sizeX; i++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int k = 0; k < _sizeZ; k++)
			{
				_agGrid(i, j, k) = new TumorGridVoxel(i, j, k);
				//std::cout << _agGrid(i, j, k) << std::endl;
			}
		}
	}

	return true;
}

Tumor::~Tumor()
{
	//std::cout << "destructor:" << std::endl;
	for (int i = 0; i < _sizeX; i++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int k = 0; k < _sizeZ; k++)
			{
				//std::cout << _agGrid(i, j, k) << std::endl;
				delete _agGrid(i, j, k);
			}
		}
	}

	delete _tInitDummy;
	delete _cInitDummy;
	delete _macInitDummy;
	delete _fibInitDummy;
	delete _TCD4InitDummy;
	delete _MDSCInitDummy;
	delete _VasInitDummy;
}




// add initial vasculature cell
void Tumor::add_init_vas(const Coord3D& c) {
	_init_vas.push_back(c);
}
// add vasculature entry points for all immune cells (each cell have one chance)
void Tumor::update_vas(){
	//remove vasculature cell from previous step
	CellVec::iterator last = stable_partition(_vecAgent.begin(), _vecAgent.end(), [&](CellAgent* p) {
		auto pC = dynamic_cast<Cell_Tumor*>(p);
		return (pC->getType() != AgentTypeEnum::CELL_TYPE_VAS);
		});// sortLiveCellAddToStats);

	for (CellVec::iterator it = last; it < _vecAgent.end(); it++)
	{
		auto p = dynamic_cast<Cell_Tumor*>(*it);
		if (p->getType() == AgentTypeEnum::CELL_TYPE_VAS)
		{
			auto pVas = dynamic_cast<Vas*>(p);
			//remove vas from previous step
			removeAgentFromGrid(pVas->getCoord(), pVas);
			_stats.incDropOut(pVas->getType(), pVas->getState());
			//nr_dropout++;
		}
		else {
			std::cout << "Wrong cell removed" << std::endl;
		}
		// Remove from grid when step is finished so they can still function during the state_change step
		delete p;
	}
	//remove vas from previous step
	_vecAgent.erase(last, _vecAgent.end());

	int neighborhood_radius = 3;
	int neighborhood_size = neighborhood_radius * 2 + 1;
	//Construct an agent grid with size X+2n-1, Y+2n-1, Z+2n-1
	//The grid is used to store the number of certain cell type in the moore neighborhood with radius n of coordinate C
	//access number of cell around Coord C by calling neighborhood_cancer_counts(c.x + neighborhood_radius, c.y + neighborhood_radius, c.z + neighborhood_radius)
	//Much faster compare to naive method
	Grid3D<int> neighborhood_cancer_counts(_sizeX + neighborhood_size - 1, _sizeY + neighborhood_size - 1, _sizeZ + neighborhood_size - 1, 0);
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		int count = 0;
		_agGrid(c)->countNumAgentByType(AgentTypeEnum::CELL_TYPE_CANCER, count, false);
		neighborhood_cancer_counts(c) = count;
	});
	neighborhood_cancer_counts.window_counts_inplace(neighborhood_size);
	//count all voxel around coord c
	Grid3D<int> neighborhood_voxel_counts(_sizeX + neighborhood_size - 1, _sizeY + neighborhood_size - 1, _sizeZ + neighborhood_size - 1, 0);
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		int count = 0;
		if (_agGrid.inGrid(c)) {
			neighborhood_voxel_counts(c) = 1;
		}
		});
	neighborhood_voxel_counts.window_counts_inplace(neighborhood_size);

	//update the cabozantinib resistance
	update_R_cabo();

	int large_p_vas = 0;

	std::vector<Coord> c_tumor;
	unsigned int nr_source_tumor;
	for_each_grid_coord(true, true, true, [&](Coord& c) {
		c_tumor.push_back(c);
		});
	nr_source_tumor = params.getVal(PARAM_TUMOR_X)* params.getVal(PARAM_TUMOR_Y)* params.getVal(PARAM_TUMOR_Z);
	
	rng.shuffle_first_k(c_tumor, nr_source_tumor);

	for (size_t i = 0; i < nr_source_tumor; i++)
	{
		int neighborhood_cancer_count = neighborhood_cancer_counts(c_tumor[i].x + neighborhood_radius, c_tumor[i].y + neighborhood_radius, c_tumor[i].z + neighborhood_radius);
		double max_num_cell = neighborhood_voxel_counts(c_tumor[i].x + neighborhood_radius, c_tumor[i].y + neighborhood_radius, c_tumor[i].z + neighborhood_radius);

		double local_cancer_ratio = neighborhood_cancer_count / max_num_cell;
		double VEGFA = get_chem(c_tumor[i], CHEM_VEGFA);
		double H_VEGFA = VEGFA / (VEGFA + params.getVal(PARAM_VAS_50));
		
		double vas_density = H_VEGFA * (1 - local_cancer_ratio) * (1 - params.getVal(PARAM_LAMBDA_V_CABO)*(_concentration_cabo /(_concentration_cabo + params.getVal(PARAM_IC50_VEGFR2))) *_R_cabo);
		//double vas_density = (raw_vas_density > 0.001 ? raw_vas_density : 0);
		/*
		std::cout << "cabo: " << _concentration_cabo << ", PARAM_IC50_VEGFR2: " << params.getVal(PARAM_IC50_VEGFR2) << ", PARAM_IC50_AXL: " << params.getVal(PARAM_IC50_AXL)
			<< ", PARAM_IC50_MET: " << params.getVal(PARAM_IC50_MET) << ", PARAM_LAMBDA_V_CABO: " << params.getVal(PARAM_LAMBDA_V_CABO) << ", R_cabo: " << _R_cabo
			<< ", vas_factor: " << (1 - params.getVal(PARAM_LAMBDA_V_CABO) * (_concentration_cabo / (_concentration_cabo + params.getVal(PARAM_IC50_VEGFR2))) * _R_cabo) << std::endl;
		*/
		auto p_agGrid = dynamic_cast<TumorGridVoxel*>(_agGrid(c_tumor[i]));
		p_agGrid->setVasDensity(vas_density);
		
		double O2 = get_chem(c_tumor[i], CHEM_O2);

		if (rng.get_unif_01() < vas_density) {
			//Adding O2 source to grid
			CellAgent* ptrCell = createOneInitCell(AgentTypeEnum::CELL_TYPE_VAS, AgentStateEnum::DEFAULT_STATE, c_tumor[i]);
			auto pVas = dynamic_cast<Vas*>(ptrCell);
			// Oxygen transport mechanism adopted from Sharan et al.
			//calculate oxygen transport rate based on vasulature density and oxygen level
			double pi = 3.1415926;
			// dimensionless ratio of intracapillary to extracapillary transport resistance from Sharan et al; D*alpha / k*Rc = 0.84
			double sigma = params.getVal(PARAM_VAS_SIGMA);
			// radius of capillary
			double Rc = params.getVal(PARAM_VAS_RC);
			double voxel_volume_cm = std::pow(params.getVal(PARAM_VOXEL_SIZE_CM), 3);
			//vasculature density (in cm/ml)
			double Lv = vas_density * voxel_volume_cm / (std::pow(Rc, 2) * pi);
			// radius of tissue in Krogh Cylinder 
			double Rt = 1 / std::pow(Lv * pi, 0.5);
			double w = Rc / Rt;
			double lambda = 1 - w * w;
			double Kv = 2 * pi * params.getVal(PARAM_O2_DIFFUSIVITY) * (lambda / (sigma * lambda - (2 + lambda) / 4 + (1 / lambda) * std::log(1 / w)));
			double O2_transport_raw = Kv * Lv * (params.getVal(PARAM_VAS_O2_CONC) - O2);
			double O2_transport = O2_transport_raw > 0 ? O2_transport_raw : 0;
			//std::cout << "Lv: " << Lv << ", Rt: " << Rt << ", w: " << w << ", lambda: " << lambda << ", Kv: " << Kv << ", Current O2 level: " << O2 << ", O2 transport rate: " << O2_transport << std::endl;
			pVas->setup_chem_source(pVas->get_source_O2(), CHEM_O2, O2_transport);
			pVas->setup_chem_sink(pVas->get_sink_VEGFA(), CHEM_VEGFA, params.getVal(PARAM_VEGFA_UPTAKE));
		}
		
		if (O2 > 0.03) {
			double ifng_conc = get_chem(c_tumor[i], CHEM_IFN);
			//double entry_density = (vas_density > 0.001 ? vas_density : 0);

			//double pRec = get_T_recruitment_prob(_concentration_t_cyt, params.getVal(PARAM_TEFF_RECRUIT_K) * (1 + recruitment_cabo_factor));
			//std::cout << "cabo: " << _concentration_cabo << ", PARAM_LAMBDA_Q_CABO: " << params.getVal(PARAM_LAMBDA_Q_CABO) << ", recruitment_cabo_factor: " << recruitment_cabo_factor << std::endl;
			
			double p_entry = vas_density;
		
			//double p_Tcell_source = p_entry;

			if (rng.get_unif_01() < p_entry) {
				_t_source.push_back(c_tumor[i]);
			}
			
			//Monocyte recruitment
			double ccl2_conc = get_chem(c_tumor[i], CHEM_CCL2);
			double H_CCL2_raw = ccl2_conc / (ccl2_conc + params.getVal(PARAM_MDSC_EC50_CCL2_REC));
			double H_CCL2 = (H_CCL2_raw > 0.01 ? H_CCL2_raw : 0);
			double p_MDSC_source = H_CCL2 * (1 - local_cancer_ratio);
			/*
			std::cout << "Coord: " << c << ", vas_density: " << vas_density << ", PARAM_REC_PORT: " << params.getVal(PARAM_REC_PORT) << ", vas_density: " << vas_density
				<< ", local_cancer_ratio: " << local_cancer_ratio << ", p_MDSC_source: " << p_MDSC_source << std::endl;
			*/

			if (rng.get_unif_01() < p_MDSC_source) {
				_mdsc_source.push_back(c_tumor[i]);
			}


			// The CCL2 EC50 for Mac and MDSC are the same, not a bug
			double p_Mac_source = H_CCL2 * (1 - local_cancer_ratio);
			if (rng.get_unif_01() < p_Mac_source) {
				_mac_source.push_back(c_tumor[i]);
			}
			//}
		}
			//std::cout << "prob T cell source: " << p_Tcell_source << ", prob MDSC source: " << p_MDSC_source << ", prob MAC source: " << p_Mac_source << std::endl;
		
	}
	std::cout << "t cell vector size: " << _t_source.size() << '\n';
	std::cout << "MDSC cell vector size: " << _mdsc_source.size() << '\n';
	std::cout << "MAC cell vector size: " << _mac_source.size() << '\n';
	return;
}

void Tumor::reset_immune_source() {
	_t_source.clear();
	_mdsc_source.clear();
	_mac_source.clear();
}

static bool sortNonTCell(CellAgent * ptrCell) { return ptrCell->getType() != AgentTypeEnum::CELL_TYPE_T; }

/*! simulate generaic tumor compartment for one time slice.
	Steps:
	-# process ABM step (cellular scale events)
		-# recruit cells
		-# for each cell,
			-# progress one step and determine consequences.
			-# process consequence of the step.
		-# go through cell list, remove dead cells, add live cells to ABM snapshot stats
		-# refresh MHC
	-# process molecular scale substeps
	-# shuffle cell list.
*/
void Tumor::timeSlice(unsigned long slice) {

	static double dt = params.getVal(PARAM_SEC_PER_TIME_SLICE);
	const double t = slice * dt;
	//static unsigned int nrMolSlice = params.getVal(PARAM_MOL_STEP_PER_SLICE);
	//static double dtMol = dt / nrMolSlice;
	// reset vasculature entry point based on

	update_vas();


	// reset time step stats
	_stats.resetStats();


	// Reset to zeros
	std::fill(_var_abm_to_qsp.begin(), _var_abm_to_qsp.end(), 0);

	//------------------------   Cellular scale events   -----------------------//
	//BOOST_TEST_MESSAGE("slice:" << slice);
	//cout << "slice:" << slice << endl;

	// Recruitment
	// model comparison: find invasive front to recruit T cells and MDSC

	time_slice_recruitment();

	//cout << "after recruitment:" << endl;

	//std::cout << "RNG check (" << slice << ") rec: " << rng.get_unif_01() << std::endl;

	time_slice_movement(t, dt);
	//cout << "after movement:" << endl;

	//std::cout << "RNG check (" << slice << ") move: " << rng.get_unif_01() << std::endl;
	time_slice_state_change(t, dt);

	//cout << "after state:" << endl;
	//std::cout << "RNG check (" << slice << ") state: " << rng.get_unif_01() << std::endl;

	if (_allow_shift_grid && slice % params.getVal(PARAM_SHIFTGRID_INTERVAL) == 0)
	{
		shift_adjust_center();
	}
	//cout << "after shift:" << endl;
	

	//std::cout << "RNG check (" << slice << ") shift: " << rng.get_unif_01() << std::endl;
	time_slice_final_scan();

	//cout << "after clear dead:" << endl;


	// ------------------------  molecular scale events     ---------------------//
	/*
	for (auto p : _vecAgent)
	{
		if (p->getState() == T_CELL_CYT)
		{
			auto pt = dynamic_cast<TCell*>(p);
			if (pt->get_source_IFNg()->get_index() != _chem.get_voxel_idx(pt->getCoord()))
				if (!pt->get_source_IFNg()->is_to_remove())
					std::cout << "location mismatch (IFNg):" << pt->getCoord() << ", "
						<< _chem.idx_to_coords(pt->get_source_IFNg()->get_index())
						<< std::endl;
		}
	}*/
	time_slice_molecular(t);

	//std::cout << "RNG check (" << slice << ") mol: " << rng.get_unif_01() << std::endl;

	//cout << "after molecular:" << endl;
	// shuffle cell vector
	if (slice % params.getVal(PARAM_SHUFFLE_CELL_VEC_INTERVAL) == 0)
	{
		std::shuffle(_vecAgent.begin(), _vecAgent.end(), rng.getRNG());
	}


	/*
	std::cout << "Num sources: " 
		<< _chem.get_num_source_sink() << std::endl;
	std::cout << "CELL INFO" << std::endl;
	}*/

	// remove all immune source 
	reset_immune_source();
}

/*! Adjust the center of the tumor compartment
*/
void Tumor::shift_adjust_center(void){
	Coord3D c = get_cam_shift();
	// only shift when distance is larger than c_tol
	int c_tol = 0;
	//std::cout << c << "," << c.length() << std::endl;
	//if (c.length()>c_tol)
	if (c.z != 0)
	{
		//shift_grid(c);
		//std::cout << "Shifting: " << c << "," << c.length() << std::endl;
		Coord3D c_zshift(0, 0, c.z);
		shift_grid(c_zshift);
		std::cout << "Shifting: " << c_zshift << std::endl;
	}
	return;
}

/*!
*/
Coord3D Tumor::get_cam_shift(void) {

//! get center of mass
	unsigned int nr_cancer = 0;
	Coord c(0, 0, 0);

	for (auto p : _vecAgent){
		if (p->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
			c = c + p->getCoord();
			nr_cancer += 1;
		}
	}
	//std::cout << "center: " << _center_target << "CoM: "<< c/nr_cancer << "," << c.length() << std::endl;
	if (nr_cancer) {
		return c / nr_cancer - _center_target;
	}
	else {
		return c;// (0,0,0)
	}
}
/*! Shift content of grid by crd
	In current version, grid size change is not allowed.
	When we need to shift center of grid by (x, y, z), 
	we move the content of the grid instead.
	The following contents need to be relocated:
	# Diffusion grid
		# concentration
		# source/sink
	# Agent grid
		# cells
		# structures
			# recruitment sources
	\param [in] crd: shifting distance (x, y, z) (of "camera")
*/
void Tumor::shift_grid(Coord3D& c_shift) {
	// k_step = 1 if k >0 else -1 (direction)
	bool x_pos = c_shift.x >= 0;
	bool y_pos = c_shift.y >= 0;
	bool z_pos = c_shift.z >= 0;

				
	int nr_shift = _sizeX*_sizeY*_sizeZ
		- (_sizeX - abs(c_shift.x))*(_sizeY - abs(c_shift.y))*(_sizeZ - abs(c_shift.z));
	int i_shift = 0;
	//cout << "shift: " << c_shift << ", nr: " << nr_shift << endl;
	CoordVec drop_out(nr_shift, Coord(0, 0, 0));
	/* shift agent grid:
		move camera by (x, y, z): 
		agents: (x0, y0, z0) -> (x0-x, y0-y, z-z0)
		*/
	// first round: save voxels to be moved out of sight to temp grid
	//	and shift remaining voxels to new location 

	for_each_grid_coord(x_pos, y_pos, z_pos, [&](Coord3D& c){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		//BOOST_TEST_MESSAGE(c);
		if (!_agGrid.inGrid(c_new)){
			//add to temp grid
			_agGrid_temp(c_target) = _agGrid(c);
			drop_out[i_shift] = c;
			i_shift++;
		}
		else{
			_agGrid(c_target) = _agGrid(c);
			for (auto pAg : _agGrid(c_target)->get_agents()){
				//BOOST_TEST_MESSAGE(c_target);
				auto pCell = dynamic_cast<Cell_Tumor*>(pAg);
				pCell->setCoord(c_target);
			}
		}
		return;
	});

	//cout << "added to drop_out: " << i_shift << endl;

	// populate new voxels (recycling out-of-site ones)
	/*
	for_each_voxel(x_pos, y_pos, z_pos, [&](Coord3D& c){
		return;
	});*/

	for (auto& c : drop_out){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		if (!_agGrid.inGrid(c_new)){
			// overwrite info and put back to grid
			_agGrid(c_target) = _agGrid_temp(c_target);
			// apply changes to voxels entering the grid (recycled)
			// remove old, 
			for (auto pAg : _agGrid(c_target)->get_agents()){
				auto pCell = dynamic_cast<Cell_Tumor*>(pAg);
				pCell->set_drop_out();
			}
			_agGrid(c_target)->remove_all_agents();
			//populate new
			populate_voxel_random(c_target);
		}
	}

	// remove dropped cells
	// now done together with dead cell

	// shift recruitment sources? 

	/* shift diffusion grid
		# shift values
		# shift point source/sink
		*/
	unsigned int num_sub = _chem.get_num_substrates();
	BioFVMGrid::chem_grid _chem_temp = _chem.get_concentrations();
	for_each_grid_coord(x_pos, y_pos, z_pos, [&](Coord3D& c){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		unsigned int idx = _chem.get_voxel_idx(c_target);
		for (size_t i = 0; i < num_sub; i++)
		{
			if (_agGrid.inGrid(c_new)){
				_chem_temp[idx][i] = _chem(c, i);
			}
			else{
				_chem_temp[idx][i] = 0;
			}
		}
		return;
	});
	_chem.setup_concentrations(_chem_temp);

	// shift source/sink: already taken care of.
	// Sources moving out of grid: set removed with agents
	// Other sources: moved together with agents.
	/* old code
	for (auto pS : _chem.get_sink_source()){
		auto c = _chem.idx_to_coords(pS->get_index());
		auto c_target = c - c_shift;
		_chem.move_source_sink(pS, c_target);
	}*/

	return;
}


void Tumor::time_slice_recruitment(){
	// enhanced T cell recruitment by vasculature normalization induced by cabozantinib
	double recruitment_cabo_factor = params.getVal(PARAM_LAMBDA_Q_CABO)* (_concentration_cabo / (_concentration_cabo + params.getVal(PARAM_IC50_VEGFR2)));
	std::cout << "cabo: " << _concentration_cabo << ", PARAM_LAMBDA_Q_CABO: " << params.getVal(PARAM_LAMBDA_Q_CABO) <<  ", recruitment_cabo_factor: " << recruitment_cabo_factor << std::endl;

	double max_cancer = _sizeX * _sizeY * _sizeZ;
	double tumor_scaler = max_cancer / (_concentration_cc * params.getVal(PARAM_AVOGADROS));
	//std::cout << "Tumor volume: " << _tumor_volume << ", _concentration_cc: " << _concentration_cc * params.getVal(PARAM_AVOGADROS)
	//	<< ", max_cancer: " << max_cancer << ", tumor_scaler: " << tumor_scaler << ", PARAM_REC_PORT: " << params.getVal(PARAM_REC_PORT) << std::endl;
	double T_cell_recruitment_rate = params.getVal(PARAM_TEFF_RECRUIT_K) * (1 + recruitment_cabo_factor) * params.getVal(PARAM_REC_PORT) * tumor_scaler;
	// CD8
	const auto dummy = _tInitDummy;
	double pRec = get_T_recruitment_prob(_concentration_t_cyt, T_cell_recruitment_rate);

	std::cout << "T cell concentration:" << _concentration_t_cyt << "\n"
		<< "recruitment rate: " << params.getVal(PARAM_TEFF_RECRUIT_K) << "\n"
		<< "rec prob (Teff): " << pRec << std::endl;
	for (auto& crd : _t_source){
		if (rng.get_unif_01() < pRec) {
			bool rec = recruitOneCellInMooreNeighborhood(dummy, crd, rng);
			if (rec)
			{
				_stats.incRecruit(dummy->getType(), dummy->getState());
				auto pT = dynamic_cast<TCell*>(_vecAgent.back());
				inc_abm_var_exchange(TUMEX_TEFF_REC);
			}
		}
	}
	// TCD4
	/**/
	const auto TCD4_dummy = _TCD4InitDummy;
	double pRec_TCD4 = get_T_recruitment_prob(_concentration_t_CD4, T_cell_recruitment_rate);
	std::cout << "T cell concentration:" << _concentration_t_CD4 << "\n"
		<< "recruitment rate: " << params.getVal(PARAM_TCD4_RECRUIT_K) << "\n"
		<< "rec prob (TCD4): " << pRec_TCD4 << std::endl;

	for (auto& crd : _t_source){
		if (rng.get_unif_01() < pRec_TCD4) {
			bool rec = recruitOneCellInMooreNeighborhood(TCD4_dummy, crd, rng);
			if (rec)
			{
				double ifng_conc = get_chem(crd, CHEM_IFN);
				double p_Th = params.getVal(PARAM_TH_FRAC);
				if (rng.get_unif_01() < p_Th) {
					//std::cout << "recruited a Th cell" << std::endl;
					_stats.incRecruit(TCD4_dummy->getType(), TCD4_dummy->getState());
					auto pT = dynamic_cast<TCD4*>(_vecAgent.back());
					pT->setup_chem_source(pT->get_source_IFNg(), CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
					pT->setup_chem_source(pT->get_source_IL_2(), CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
					pT->setup_chem_source(pT->get_source_IL_10(), CHEM_IL_10, 0);
					pT->setup_chem_source(pT->get_source_TGFB(), CHEM_TGFB, 0);
					inc_abm_var_exchange(TUMEX_TCD4_REC);
					
				}
				else 
				{
					//std::cout << "recruited a Treg cell" << std::endl;
					auto pT = dynamic_cast<TCD4*>(_vecAgent.back());
					pT->setup_chem_source(pT->get_source_IFNg(), CHEM_IFN, 0);
					pT->setup_chem_source(pT->get_source_IL_2(), CHEM_IL_2, 0);
					pT->setup_chem_source(pT->get_source_IL_10(), CHEM_IL_10, params.getVal(PARAM_TREG_IL_10_RELEASE));
					pT->setup_chem_source(pT->get_source_TGFB(), CHEM_TGFB, params.getVal(PARAM_TREG_TGFB_RELEASE));
					pT->setTreg();

					_stats.incRecruit(TCD4_dummy->getType(), TCD4_dummy->getState());
					inc_abm_var_exchange(TUMEX_TCD4_REC);
				}

			}
		}
	}


	//std::cout << "number of th recruited: " << cd4_recruited << std::endl;

	// MDSC
	/**/
	const auto MDSC_dummy = _MDSCInitDummy;
	//std::cout << "rec prob (MDSC): " << pRec_MDSC << std::endl;	

	const double voxel_volume = std::pow(double(params.getVal(PARAM_VOXEL_SIZE)) / 1e6, 3);

	for (auto& crd : _mdsc_source){	
		// for each coordinate, get the local concentration of CCL2 
		double ccl2_conc = get_chem(crd, CHEM_CCL2);

		// number of MDSC in the ABM
		// free spaces for MDSC 

		
		//double mdsc_free_fraction = ((params.getVal(PARAM_MDSC_MAX) * abm_volume) - _concentration_mdsc) / abm_volume;
		
		//MDSC Recruitment induced by CCL2
		double k_Rec_MDSC_ccl2 = params.getVal(PARAM_MDSC_RECRUIT_BY_CCL2_K);
		
		/*
		std::cout << "PARAM_MDSC_RECRUIT_BY_CCL2_K: " << params.getVal(PARAM_MDSC_RECRUIT_BY_CCL2_K)
			<< ",  abm_volume: " << abm_volume
			<< ", ccl2 : " << ccl2_conc << std::endl;
		*/
		double pRec_MDSC = get_Myeloid_recruitment_prob(k_Rec_MDSC_ccl2, voxel_volume);


		if (rng.get_unif_01() < pRec_MDSC) {
			bool rec = recruitOneCellInMooreNeighborhood(MDSC_dummy, crd, rng);		
			if (rec)
			{
				_stats.incRecruit(MDSC_dummy->getType(), MDSC_dummy->getState());
				auto pT = dynamic_cast<MDSC*>(_vecAgent.back());
				pT->setup_chem_source(pT->get_source_ArgI(), CHEM_ARGI, params.getVal(PARAM_ARGI_RELEASE));
				pT->setup_chem_source(pT->get_source_NO(), CHEM_NO, params.getVal(PARAM_NO_RELEASE));
				//sink
				pT->setup_chem_sink(pT->get_sink_CCL2(), CHEM_CCL2, params.getVal(PARAM_CCL2_UPTAKE));
			}
		}
	}	

	const auto mac_dummy = _macInitDummy;

	for (auto& crd : _mac_source) {
		// for each coordinate, get the local concentration of CCL2 
		double ccl2_conc = get_chem(crd, CHEM_CCL2);

		//double mdsc_free_fraction = ((params.getVal(PARAM_MDSC_MAX) * abm_volume) - _concentration_mdsc) / abm_volume;

		//MDSC Recruitment constant induced by CCL2
		double k_Rec_Mac_ccl2 = params.getVal(PARAM_MAC_RECRUIT_BY_CCL2_K);

		/*
		std::cout << "PARAM_MDSC_RECRUIT_BY_CCL2_K: " << params.getVal(PARAM_MDSC_RECRUIT_BY_CCL2_K)
			<< ",  abm_volume: " << abm_volume
			<< ", ccl2 : " << ccl2_conc << std::endl;
		*/
		double pRec_Mac = get_Myeloid_recruitment_prob(k_Rec_Mac_ccl2, voxel_volume);


		if (rng.get_unif_01() < pRec_Mac) {
			bool rec = recruitOneCellInMooreNeighborhood(mac_dummy, crd, rng);
			if (rec)
			{
				_stats.incRecruit(mac_dummy->getType(), mac_dummy->getState());
				auto pT = dynamic_cast<Mac*>(_vecAgent.back());
				if (pT->getState() == AgentStateEnum::MAC_M1)
				{
					//once M1 macrophage in contact with cancer cell, they release proinflammatory cytokine
					pT->setup_chem_source(pT->get_source_IFNg(), CHEM_IFN, 0);
					pT->setup_chem_source(pT->get_source_IL_12(), CHEM_IL_12, 0);
					pT->setup_chem_source(pT->get_source_TGFB(), CHEM_TGFB, 0);
					pT->setup_chem_source(pT->get_source_IL_10(), CHEM_IL_10, 0);
					pT->setup_chem_source(pT->get_source_VEGFA(), CHEM_VEGFA, 0);
					//sink
					pT->setup_chem_sink(pT->get_sink_CCL2(), CHEM_CCL2, params.getVal(PARAM_CCL2_UPTAKE));

					if (rng.get_unif_01() < 0.3) {
						pT->setM2();
					}

				}
				else
				{
					pT->setup_chem_source(pT->get_source_IFNg(), CHEM_IFN, 0);
					pT->setup_chem_source(pT->get_source_IL_12(), CHEM_IL_12, 0);
					pT->setup_chem_source(pT->get_source_TGFB(), CHEM_TGFB, params.getVal(PARAM_MAC_TGFB_RELEASE));
					pT->setup_chem_source(pT->get_source_IL_10(), CHEM_IL_10, params.getVal(PARAM_MAC_IL_10_RELEASE));
					pT->setup_chem_source(pT->get_source_VEGFA(), CHEM_VEGFA, params.getVal(PARAM_MAC_VEGFA_RELEASE));
					//sink
					pT->setup_chem_sink(pT->get_sink_CCL2(), CHEM_CCL2, params.getVal(PARAM_CCL2_UPTAKE));
				}
			}
		}
	}
}

void Tumor::time_slice_movement(double t, double dt){
	//BOOST_TEST_MESSAGE("total cells: " << _vecAgent.size());
	for (CellVec::size_type i = 0; i < _vecAgent.size(); i++)
	{
		auto ptP = dynamic_cast<Cell_Tumor*>(_vecAgent[i]);
		Coord c(0, 0, 0);
		//cout << "move: "  <<  i << ", type: " << ptP->getType() << endl;
		//cout.flush();

		//BOOST_TEST_MESSAGE(i);
		bool move = ptP->agent_movement_step(t, dt, c);
		if (move)
		{
			// need to do this first or problems will occur when relocating other cells to current location
			try{
				removeAgentFromGrid(ptP->getCoord(), ptP);
			}
			catch (...){
				cout << ptP->toString() << endl;
				cout << ptP->getCoord() << "->" << c << endl;
				std::rethrow_exception(std::current_exception());
			}

			//move current cell
			ptP->setCoord(c);
			addAgentToGrid(c, ptP);

			//printCellInfo(slice, ptP, "after moving cell");
			//ptP->move_all_source_sink();

			_stats.incMove(ptP->getType(), ptP->getState());
		}
		
	}
	return;

}

void Tumor::time_slice_state_change(double t, double dt){

	/*Scan neighbor*/
	for (auto const p: _vecAgent)
	{
		auto pCell = dynamic_cast<Cell_Tumor*>(p);
		pCell->agent_state_scan();
	}
	
	//std::cout << "RNG check "<< " scan: " << rng.get_unif_01() << std::endl;
	/*process state change*/
	CellVec::size_type originalCount = _vecAgent.size();
	for (CellVec::size_type i = 0; i < originalCount; i++)
	{
		CellAgent *ptP = _vecAgent[i];
		Coord c(0, 0, 0);

		//std::cout << "processing: " << ptP->getID() << ", " << ptP->getType() << ", " << ptP->getState() << std::endl;
		bool divide = ptP->agent_state_step(t, dt, c);
		//std::cout << "divide: " << divide << std::endl;
		if (divide)
		{
			//std::cout << "adding daughter cell" << std::endl;
			auto ptD = ptP->createCellCopy();
			//add daugther cell
			ptD->setCoord(c);
			addAgentToGrid(c, ptD);
			_vecAgent.push_back(ptD);

			if (ptP->getType() == AgentTypeEnum::CELL_TYPE_CANCER)
			{
				//std::cout << "cancer cell" << std::endl;
				auto ptPCancer = dynamic_cast<CancerCell*>(ptP);
				auto ptDCancer = dynamic_cast<CancerCell*>(ptD);

				ptDCancer->_stemID = ptPCancer->_stemID;
				if (ptP->getState() == AgentStateEnum::CANCER_STEM)
				{
					bool asymmetric = rng.get_unif_01() < params.getVal(PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB);

					//source
					ptDCancer->setup_chem_source(ptDCancer->get_source_CCL2(), CHEM_CCL2, params.getVal(PARAM_CCL2_RELEASE));
					ptDCancer->setup_chem_source(ptDCancer->get_source_TGFB(), CHEM_TGFB, params.getVal(PARAM_CANCER_TGFB_RELEASE));
					ptDCancer->setup_chem_source(ptDCancer->get_source_VEGFA(), CHEM_VEGFA, params.getVal(PARAM_CANCER_STEM_VEGFA_RELEASE));
					//sink
					ptDCancer->setup_chem_sink(ptDCancer->get_sink_O2(), CHEM_O2, params.getVal(PARAM_O2_UPTAKE));
					ptDCancer->setup_chem_sink(ptDCancer->get_sink_IFNg(), CHEM_IFN, params.getVal(PARAM_CANCER_IFN_G_UPTAKE));

					if (asymmetric)
					{
						ptDCancer->setProgenitor();
					}
					else{
						ptDCancer->_stemID = ptPCancer->getID();
					}
				} 
				else if (ptP->getState() == AgentStateEnum::CANCER_PROGENITOR)
				{
					//source
					ptDCancer->setup_chem_source(ptDCancer->get_source_CCL2(), CHEM_CCL2, params.getVal(PARAM_CCL2_RELEASE));
					ptDCancer->setup_chem_source(ptDCancer->get_source_TGFB(), CHEM_TGFB, 0);
					ptDCancer->setup_chem_source(ptDCancer->get_source_VEGFA(), CHEM_VEGFA, params.getVal(PARAM_CANCER_PRO_VEGFA_RELEASE));
					//sink
					ptDCancer->setup_chem_sink(ptDCancer->get_sink_O2(), CHEM_O2, params.getVal(PARAM_O2_UPTAKE));
					ptDCancer->setup_chem_sink(ptDCancer->get_sink_IFNg(), CHEM_IFN, params.getVal(PARAM_CANCER_IFN_G_UPTAKE));
				}
				else {
					std::cerr << "Wrong cancer cell Division" << std::endl;
				}
			}
			
			if (ptP->getType() == AgentTypeEnum::CELL_TYPE_T)
			{
				//std::cout << "cancer cell" << std::endl;
				auto ptPTCell = dynamic_cast<TCell*>(ptP);
				auto ptDTCell = dynamic_cast<TCell*>(ptD);

				if (ptD->getState() == AgentStateEnum::T_CELL_CYT)
				{
					ptDTCell->setup_chem_source(ptDTCell->get_source_IFNg(), CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
					ptDTCell->setup_chem_source(ptDTCell->get_source_IL_2(), CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
					ptDTCell->setup_chem_sink(ptDTCell->get_sink_IL_2(), CHEM_IL_2, params.getVal(PARAM_IL_2_UPTAKE));
				}
			}

			if (ptP->getType() == AgentTypeEnum::CELL_TYPE_TCD4)
			{
				//std::cout << "cancer cell" << std::endl;
				auto ptPTCD4 = dynamic_cast<TCD4*>(ptP);
				auto ptDTCD4 = dynamic_cast<TCD4*>(ptD);
				if (ptD->getState() == AgentStateEnum::TCD4_Th)
				{
					ptDTCD4->setup_chem_source(ptDTCD4->get_source_IFNg(), CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
					ptDTCD4->setup_chem_source(ptDTCD4->get_source_IL_2(), CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
					ptDTCD4->setup_chem_source(ptDTCD4->get_source_TGFB(), CHEM_TGFB, 0);
					ptDTCD4->setup_chem_source(ptDTCD4->get_source_IL_10(), CHEM_IL_10, 0);
				}

				if (ptD->getState() == AgentStateEnum::TCD4_TREG)
				{
					ptDTCD4->setup_chem_source(ptDTCD4->get_source_IFNg(), CHEM_IFN, 0);
					ptDTCD4->setup_chem_source(ptDTCD4->get_source_IL_2(), CHEM_IL_2, 0);
					ptDTCD4->setup_chem_source(ptDTCD4->get_source_TGFB(), CHEM_TGFB, params.getVal(PARAM_TREG_TGFB_RELEASE));
					ptDTCD4->setup_chem_source(ptDTCD4->get_source_IL_10(), CHEM_IL_10, params.getVal(PARAM_TREG_IL_10_RELEASE));
				}
			}


			//auto pD_tumor = dynamic_cast<Cell_Tumor*>(ptD);
			//pD_tumor->move_all_source_sink();

			_stats.incProlif(ptD->getType(), ptD->getState());
			//std::cout << "added" << std::endl;
		}
	}
	return;
}
/*! Final scan during a slice
	Process cell counts and slice statistics 
	Remove dead cells from cell vector
*/
void Tumor::time_slice_final_scan(void){

	CellVec::iterator last = stable_partition(_vecAgent.begin(), _vecAgent.end(), [&](CellAgent* p){
		auto pC = dynamic_cast<Cell_Tumor*>(p);
		return !(pC->isDead() || pC->is_drop_out());
	});// sortLiveCellAddToStats);

	// count live cells (and cancer cells)
	int nrPDL1 = 0;
	for (auto it = _vecAgent.begin(); it < last; it++)
	{
		Cell_Tumor* ptrCell = dynamic_cast<Cell_Tumor*>(*it);
		//statistics
		_stats.incCellState(ptrCell->getType(), ptrCell->getState());
		if (ptrCell->getType() == CELL_TYPE_CANCER)
		{
			inc_abm_var_exchange(TUMEX_CC);
		}
		if (ptrCell->is_PDL1_pos())
		{
			nrPDL1++;
		}
	}
	_cancer_counts = get_stats().getCancerCell();
	_stats.set_stats_misc(STATS_MISC_PDL1_POS, nrPDL1);
	int nr_death = 0;
	int nr_dropout = 0;
	// delete removed cells
	for (CellVec::iterator it = last; it < _vecAgent.end(); it++)
	{
		auto p = dynamic_cast<Cell_Tumor*>(*it);
		if (p->isDead())
		{

			// Cancer cell dies -> antigen presentation
			if (p->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
				inc_abm_var_exchange(TUMEX_CC_DEATH);
				//Release antigen after death
			}
			_stats.incDeath(p->getType(), p->getState());
			if (!p->is_drop_out())
			{
				// drop out cells are already removed from voxel
				removeAgentFromGrid(p->getCoord(), p);
			}
			nr_death++;
		}
		else if (p->is_drop_out() && !p->isDead()){
			// dead and dropout are only counted as dead

			removeAgentFromGrid(p->getCoord(), p);
			_stats.incDropOut(p->getType(), p->getState());
			nr_dropout++;
		}
		else{
			std::cout << "Error: cells neither dead nor drop out are removed" << std::endl;
		}
		// Remove from grid when step is finished so they can still function during the state_change step
		delete p;
	}
	_vecAgent.erase(last, _vecAgent.end());
	//std::cout << "Death: " << nr_death << ", dropout: " << nr_dropout << std::endl;
	return;
}

void Tumor::time_slice_molecular(double t){

	static double dt = params.getVal(PARAM_SEC_PER_TIME_SLICE);
	_chem.remove_dead_source_sink();
	if (params.getVal(PARAM_DIFFUSION_ON) && !_chem.grid_skippable())
	{
		for (auto const pS : _chem.get_sink_source()){
			auto c = _chem.idx_to_coords(pS->get_index());
		}
		static double dt_mol = dt / params.getVal(PARAM_MOLECULAR_STEP_PER_SLICE);
		double t_mol = t;
		double t_mol_end = t + dt;
		while (t_mol < t_mol_end)
		{
			// diffusion time step, including decay/point source/sink
			_chem.timestep(dt_mol);
			t_mol += dt_mol;
		}
	}
}




//! default header for extra remark column when writing cell grid to file
std::string Tumor::getExtraRemarkHeader() const {
	std::stringstream ss;
	ss << "extra";
	return ss.str();
}

/*! Print grid snapshot to a snapshot file.
	Info includes: Element type, IL2 concentration, and MHC concentration.
*/
std::string Tumor::printGridToFile() const {
//void Tumor::printGridToFile(unsigned long slice) const {
	std::stringstream ss;
	// header
	ss << "x,y,z,"<<_chem.get_substrate_names() << std::endl;
	// content
	ss << _chem;
	return ss.str();
}

/*! print gradient of cytokines along selected direction
* x=0; y=1; z=2 (axis)
*/
std::string Tumor::printGradientToFile(int axis) {
	//void Tumor::printGridToFile(unsigned long slice) const {
	std::stringstream ss;
	// header
	ss << "x,y,z," << _chem.get_substrate_names() << std::endl;
	// content
	for (size_t i = 0; i < _chem.get_num_voxel(); i++) {
		Coord3D c = _chem.idx_to_coords(i);
		std::vector< std::vector<double> > gradient = _chem.get_gradient(i);
		ss << c.x << ',' << c.y << ',' << c.z << ",";
		auto delim = "";
		for (size_t j = 0; j < _chem.get_num_substrates(); j++) {
			ss << delim << gradient[j][axis];
			delim = ",";
		}
		ss << std::endl;
	}

	return ss.str();
}

// Ask this question about const (current cannot add const to the function
std::string Tumor::printVasToFile() {
	//void Tumor::printGridToFile(unsigned long slice) const {
	std::stringstream ss;
	// header
	ss << "x,y,z,density" << std::endl;
	// content
	
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {

		TumorGridVoxel* p_agGrid = dynamic_cast<TumorGridVoxel*>(_agGrid(c));
		double density = p_agGrid->getVasDensity();

		//std::cout << c.x << "," << c.y << "," << c.z << "," << density << std::endl;
		ss << c.x << "," << c.y << "," << c.z << "," << density << std::endl;
		return;
	});
	return ss.str();
	
}

/*! Print state of ODE system of all T cells to ODE stats file,
	includeing a time stamp and cell ID.
	\param [in] slice: time slice
*/
void Tumor::printCellOdeToFile(unsigned long slice) const {
}

/*! print grid snapshot to screen
	\param [in] slice: time slice
*/
void Tumor::printGridToScreen(unsigned long slice)const {

	// header
	cout << "Time slice: "<< slice << endl;
	cout << "Grid size: " << getGridSize() << ", "<< _sizeX << " by "
		<< _sizeY << " by " << _sizeY <<  endl;
	//cout << getGridContent() << endl;

	cout << "x, y, z, nr_ag," << endl;
	for (int k = 0; k < _sizeZ; k++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int i = 0; i < _sizeX; i++)
			{
				Coord c(i, j, k);
				cout << c.x << "," << c.y << "," << c.z << ",";
				AgentGridVoxel* voxel = _agGrid.get(i, j, k);
				TumorGridVoxel * v2 = dynamic_cast<TumorGridVoxel*>(voxel);
				cout << voxel->getNumAgents();
				//cout << v2->getDistToOrigin();
				//cout << v2->getLoc();
				cout << endl;
			}
		}
	}
}

double Tumor::get_chem(const Coord3D&c, chem_ID i)const{
	if (!_agGrid.inGrid(c))
	{
		return 0;
	}
	else 
	{ 
		return _chem(c, i); 
	}
}
//-----------------------------------------  Protected --------------------------------------//

//-----------------------------------------  Private  --------------------------------------//

/*! Setup compartment environment. This includes: sources of T cells and MDSCs.
	create vasculature by mapping graph file to the grid.
	all locations mapped onto the graph are designated as T cell or MDSC sources.
*/
void Tumor::initEnvironment() {

	int nrSource = 0;

}

/*! Setup initial cells. Called when simulation starts.
	Create initial cells and set their coordinates;
	put initial cells on the grid;
	put initial cells into cell vector.
	After all intial cells are generated, go through the vector and record initial stats
	\param [in] initialCellFileName: initial arrangement of cells

	read initial cell configuration from a separate file
	structure of initial condition file:
	<initialCondition>
	  <cell>
		<x></x>
		<y></y>
		<z></z>
		<celltype></celltype>
		<cellstate></cellstate>
	  </cell>
	  <cell>
		...
	  </cell>
	</initialCondition>

	<cellCluster/> elements will have an additional subelement <count/>,
	indicating the number of cells in that cluster.
*/
void Tumor::initCell(string initialCellFileName){

	namespace pt = boost::property_tree;
	const std::string cellPath = "initialCondition";
	const std::string cellTag = "cell";
	const std::string cellClusterTag = "cellCluster";
	using std::pair;
	pt::ptree tree;

	// in case no initial condition file provided
	if (initialCellFileName.empty())
	{
		init_cell_fill_grid();

		//init_cell_single_center();
	}
	else {
		try {
			pt::read_xml(initialCellFileName, tree, pt::xml_parser::trim_whitespace);
			// get nodes
			BOOST_FOREACH(pt::ptree::value_type const& cell, tree.get_child(cellPath)) {
				if (cell.first == cellTag) {
					icProperty ic = cell.second;
					unsigned int e = ic.get<unsigned int>("celltype");
					unsigned int s = ic.get<unsigned int>("cellstate", 0);
					unsigned int x = ic.get<unsigned int>("x");
					unsigned int y = ic.get<unsigned int>("y");
					unsigned int z = ic.get<unsigned int>("z");
					AgentType type = static_cast<AgentType>(e);
					AgentState state = static_cast<AgentState>(s);
					createOneInitCell(type, state, Coord3D(x,y,z));
				}
				else if (cell.first == cellClusterTag)
				{
					icProperty ic = cell.second;
					createClusterInitCell(ic);
				}
			}
		}
		catch (std::exception & e) {
			std::cerr << "Error creating initial cells" << std::endl;
			//std::cerr << e.what() << std::endl;
			throw e;
		}
	}

	// check and update stats
	for (auto ptrCell : _vecAgent)
	{
		_stats.incCellState(ptrCell->getType(), ptrCell->getState());
		if (ptrCell->isDead())
		{
			cerr << "dead cell initiated" << endl;
			exit(1);
		}
	}
}

/*! 
*/
void Tumor::init_cell_single_center(void){
	auto crd = Coord3D(_sizeX, _sizeY, _sizeZ) / 2;
	AgentType type = AgentTypeEnum::CELL_TYPE_CANCER;
	AgentState state = AgentStateEnum::CANCER_STEM;
	createOneInitCell(type, state, crd);
	//_stats.incRecruit(type, state);
	return;
}

/*! initial cell: fill grid with random population
*/
void Tumor::init_cell_fill_grid(void){
	int nr_voxel = _agGrid.getSize();
	/*
	// For testing
	unsigned int k = 100;
	std::vector<int> voxel_id(nr_voxel);
	std::iota(std::begin(voxel_id), std::end(voxel_id), 0); 
	//std::shuffle(voxel_id.begin(), voxel_id.begin()+k, rng.getRNG());
	rng.shuffle_first_k<int>(voxel_id, k);
	for (size_t i = 0; i < k; i++)
	{
		Coord3D crd = _agGrid.get_coord(voxel_id[i]);
		//std::cout << "Voxel_id: " << voxel_id[i] << ", " << crd << endl;
		createOneInitCell(AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_CYT, crd);
	}*/ 

	//std::cout << "start population" << std::endl;
	Coord c_sum(0, 0, 0);
	unsigned int n = 0;

	for_each_grid_coord(true, true, true, [&](Coord3D& c){
		if (populate_voxel_random(c))
		{
			c_sum = c_sum + c;
			n += 1;
		}
		return;
	});

	_center_target = n ? c_sum / n : c_sum;
	//std::cout << _center_target << std::endl;
	return;
}

Cell_Tumor* Tumor::populate_voxel_random(const Coord3D& crd) {
	AgentType type = AgentTypeEnum::AGENT_DUMMY;
	AgentState state = AgentStateEnum::DEFAULT_STATE;
	//bool create_cancer_cell = false;
	Cell_Tumor* pTumorCell = NULL;
	int div = 0;
	//Adding cells based on initial condition
	if (_voxel_ic.get_type_state(crd, rng, type, state, div)) {
		CellAgent* ptrCell = createOneInitCell(type, state, crd);
		if (ptrCell && type == AgentTypeEnum::CELL_TYPE_CANCER)
		{
			//create_cancer_cell = true;
			CancerCell* pCancerCell = dynamic_cast<CancerCell*>(ptrCell);
			pTumorCell = pCancerCell;
			if (state == AgentStateEnum::CANCER_PROGENITOR)
			{
				pCancerCell->setDivCounter(div);
				pCancerCell->randomize_div_cd(int(params.getVal(PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) + .5));
			}
			else if (state == AgentStateEnum::CANCER_STEM)
			{
				pCancerCell->randomize_div_cd(int(params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE) + .5));
				
			}
		}
	}

	return pTumorCell;
}

/*! recruit a cluster of cells to the lattice
	\param [in] ic: initial condition

	-# create an empty queue; count = 0
	-# push indicated location to the queue
	-# recruit first cell to lattice, count++
	-# iterate until count matched or queue is empty
	  -# deque, get coordinate
	  -# iterate until no space found
		-# push found location to queue
		-# recruit, count++
*/
void Tumor::createClusterInitCell(icProperty &ic) {

	unsigned int nrCellToCreate = ic.get<unsigned int>("count");
	//ElementType e; not working, cannot directly map string in property_tree to enum
	unsigned int e = ic.get<unsigned int>("celltype");
	unsigned int s = ic.get<unsigned int>("cellstate", 0);
	unsigned int x = ic.get<unsigned int>("x");
	unsigned int y = ic.get<unsigned int>("y");
	unsigned int z = ic.get<unsigned int>("z");


	AgentType type = static_cast<AgentType>(e);
	AgentState state = static_cast<AgentState>(s);
	const ShapeBase *shape;
	switch (type)
	{
	case AgentTypeEnum::AGENT_DUMMY:
		break;
	case AgentTypeEnum::CELL_TYPE_CANCER:
		shape = _cInitDummy->getCellShape();
		break;
	case AgentTypeEnum::CELL_TYPE_T:
		shape = _tInitDummy->getCellShape();
		break;	
	case AgentTypeEnum::CELL_TYPE_TCD4:
		shape = _TCD4InitDummy->getCellShape();
		break;
	case AgentTypeEnum::CELL_TYPE_MDSC:
		shape = _MDSCInitDummy->getCellShape();
		break;	
	case AgentTypeEnum::CELL_TYPE_FIB:
		shape = _fibInitDummy->getCellShape();
		break;
	case AgentTypeEnum::CELL_TYPE_VAS:
		shape = _VasInitDummy->getCellShape();
		break;
	default:
		throw std::invalid_argument("unknown cell type in initial cells");
	}

	unsigned int count = 0;
	std::queue<Coord> nextCenter;
	auto crd = Coord(x, y, z);
	nextCenter.push(crd);
	createOneInitCell(type, state, crd);
	count++;

	while (count < nrCellToCreate || nextCenter.empty())
	{
		Coord c = nextCenter.front();
		nextCenter.pop();
		bool spaceFound = true;
		Coord cNew;
		while (count < nrCellToCreate)
		{
			// create one cell
			int idx;
			//spaceFound = getOneOpenDestinationByScan(shape->getProlifNewOccupy(),
			//	shape->getProlifSearchSeq(), shape->getProlifOccupyMap(),
			//	shape->getProlifRelocateMap(), c, ElementType(e), idx);
			spaceFound = getOneOpenVoxel(shape->getProlifDestinationVoxels(),
				shape->getProlifDestinationAnchor(), c, type, idx, rng);
			if (spaceFound)
			{
				cNew = shape->getProlifDestinationAnchor()[idx] + c;
				createOneInitCell(type, state, cNew);
				count++;
				nextCenter.push(cNew);
			}
			else {
				break;
			}
		}
	}
}

/*! recruit one cell to grid
	\param [in] e: cell type to create
	\param [in] state: cell state
	\param [in] crd: 3D coordinate, const
*/
CellAgent* Tumor::createOneInitCell(AgentType type, AgentState state, const Coord3D& crd) {

	CellAgent * ptrCell = NULL;

	//std::cout << "type: " << type << "; state: " << state << std::endl;
	// only add cell if the target voxel can take it
	if (_agGrid(crd)->isOpenToType(type))
	{
		switch (type)
		{
		case AgentTypeEnum::AGENT_DUMMY:
			throw std::invalid_argument("dummy type in initial cells");
			break;
		case AgentTypeEnum::CELL_TYPE_CANCER:
			ptrCell = _cInitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_T:
			ptrCell = _tInitDummy->createCellCopy();
			break;	
		case AgentTypeEnum::CELL_TYPE_TCD4:
			ptrCell = _TCD4InitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_MAC:
			ptrCell = _macInitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_FIB:
			ptrCell = _fibInitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_MDSC:
			ptrCell = _MDSCInitDummy->createCellCopy();
			break;	
		case AgentTypeEnum::CELL_TYPE_VAS:
			ptrCell = _VasInitDummy->createCellCopy();
			break;
		default:
			throw std::invalid_argument("unknown cell type in initial cells");
		}

		ptrCell->setCoord(crd);
		ptrCell->setAgentState(state);
		addAgentToGrid(crd, ptrCell);
		// reset state
		// now all starts from base state; otherwise report error
		/*
		if (state != AgentStateEnum::DEFAULT_STATE)
		{
		throw std::invalid_argument("invalide initial cell state");
		}
		*/
		// randomize life
		switch (type)
		{
		case AgentTypeEnum::AGENT_DUMMY:
			break;
		case AgentTypeEnum::CELL_TYPE_CANCER:{
			auto pCancerCell = dynamic_cast<CancerCell*>(ptrCell);
			
			//When cancer cells are created initally, its assumed to have Stem cell property (from the constructor)
			//set up CCL2, TGFB, VEGFA source
			pCancerCell->setup_chem_source(pCancerCell->get_source_CCL2(), CHEM_CCL2, params.getVal(PARAM_CCL2_RELEASE));
			pCancerCell->setup_chem_source(pCancerCell->get_source_TGFB(), CHEM_TGFB, params.getVal(PARAM_CANCER_TGFB_RELEASE));
			pCancerCell->setup_chem_source(pCancerCell->get_source_VEGFA(), CHEM_VEGFA, params.getVal(PARAM_CANCER_STEM_VEGFA_RELEASE));
			// O2, IFG_g sink
			pCancerCell->setup_chem_sink(pCancerCell->get_sink_O2(), CHEM_O2, params.getVal(PARAM_O2_UPTAKE));
			pCancerCell->setup_chem_sink(pCancerCell->get_sink_IFNg(), CHEM_IFN, params.getVal(PARAM_CANCER_IFN_G_UPTAKE));
			
			//Then, if they need to be set to progenitor cells, they are set by following code
			if (state == AgentStateEnum::CANCER_PROGENITOR)
			{
				pCancerCell->setProgenitor();
		
			}
			else if (state == AgentStateEnum::CANCER_SENESCENT)
			{
				pCancerCell->setSenescent();
			}
			break;
		}
		case AgentTypeEnum::CELL_TYPE_T:{
			auto pTCell = dynamic_cast<TCell*>(ptrCell);
			int life = (int)(rng.get_unif_01() * pTCell->getCellLife() + 0.5);
			pTCell->setCellLife(life);
			if (state == AgentStateEnum::T_CELL_CYT)
			{
				//pTCell->setup_source_IFNg();
				pTCell->setup_chem_source(pTCell->get_source_IFNg(), 
					CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
				pTCell->setup_chem_source(pTCell->get_source_IL_2(), 
					CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
			}
			break;
		}
		default:
			break;
		}
		// add to cell vector
		_vecAgent.push_back(ptrCell);
		_stats.incDropIn(type, state);
	}

	return ptrCell;
}

/*! return variables needed for QSP module 
	already updated during simulation.
*/
const std::vector<double>& Tumor::get_var_exchange(void){
	return _var_abm_to_qsp;
}
/*! update ABM module with variables from QSP 
*/
void Tumor::update_abm_with_qsp(const std::vector<double>& qsp_var){
	_concentration_cc = qsp_var[LymphCentral::QSPEX_TUM_C] / params.getVal(PARAM_AVOGADROS);
	_concentration_t_cyt = qsp_var[LymphCentral::QSPEX_CENT_TEFF];
	_concentration_t_CD4 = qsp_var[LymphCentral::QSPEX_CENT_TCD4];
	_concentration_t_eff_tum = qsp_var[LymphCentral::QSPEX_TUM_TEFF];
	_concentration_t_CD4_tum = qsp_var[LymphCentral::QSPEX_TUM_TCD4];
	_concentration_nivo = qsp_var[LymphCentral::QSPEX_TUM_NIVO];
	_concentration_ipi = qsp_var[LymphCentral::QSPEX_TUM_IPI];
	_concentration_cabo = qsp_var[LymphCentral::QSPEX_TUM_CABO];
	_concentration_cx = qsp_var[LymphCentral::QSPEX_TUM_CX];
	_concentration_t_exh = qsp_var[LymphCentral::QSPEX_TUM_TEXH];
	_concentration_mat_ckine = qsp_var[LymphCentral::QSPEX_TUM_MAT_CKINE];
	_concentration_mdsc = double(_stats.getMDSC()) / params.getVal(PARAM_AVOGADROS);
	//_concentration_ent = qsp_var[LymphCentral::QSPEX_TUM_ENT];
	_tumor_volume = params.getVal(PARAM_VT_MIN) + params.getVal(PARAM_V_CELL) * (_concentration_cc + _concentration_cx) + params.getVal(PARAM_V_TCELL) * (_concentration_t_eff_tum + _concentration_t_CD4_tum + _concentration_t_exh);
	_cmax = qsp_var[LymphCentral::QSPEX_TUM_CMAX];
	_f_tum_cap = _concentration_cc / _cmax;
	return;
}

/*! T recruitment probability
*/
double Tumor::get_T_recruitment_prob(double c, double k_rec) const{
	/*
	p = k (1/mol) * Cent.T (mol)
	*/	
	double num_rec = c * k_rec;
	//std::cout << "num_rec: " << num_rec  << std::endl;
	double p = (num_rec < 1 ? num_rec : 1);

	return p;
}

/*! MDSC recruitment probability
*/
// calculate recruitment probability for both macrophage and MDSCs (Myeloid cells)
double Tumor::get_Myeloid_recruitment_prob(double c, double k_rec) const{
	//p = k (m^3/mol) * Tum.Vol (m^3) * (Tum.MDSCmax * Tum.Vol - Tum.MDSC) (mol)

	//double num_rec = c * k_rec / _tumor_volume;
	double num_rec = c * k_rec;
	//std::cout << "Myeloid Recruitment rate: " << c << ", volume: " << k_rec << ", prob: " << num_rec << std::endl;
	double p;
	if (num_rec > 0){
		p = (num_rec < 1 ? num_rec : 1);
	}
	else{
		p = 0;
	}

	return p;
	
}

void Tumor::update_R_cabo() {
	_R_cabo = _concentration_cabo / (_concentration_cabo + params.getVal(PARAM_IC50_AXL));
	std::cout << "_R_cabo: " << _R_cabo << ", PARAM_IC50_AXL: " << params.getVal(PARAM_IC50_AXL) << std::endl;
}


};
};
