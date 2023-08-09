#include "TumorGridVoxel.h"
#include <math.h>

#include "../../core/GlobalUtilities.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

TumorGridVoxel::TumorGridVoxel(int x, int y, int z)
	: AgentGridVoxel()
	, _distToOrigin(sqrt((double)x*x+y*y+z*z))
	, _vas_density(0)
	//, _str_loc()
{
	//_str_loc = "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
}


TumorGridVoxel::~TumorGridVoxel()
{
}

/*!
  Check if this voxel is open to type
  \param[in] AgentType t: type to check against.
*/
bool TumorGridVoxel::isOpenToType(AgentType t)const {
	bool res = false;
	int count;
	if (t == AgentTypeEnum::AGENT_DUMMY)
	{
		res = AgentGridVoxel::isOpenToType(t);
	}
	else if (t == AgentTypeEnum::CELL_TYPE_CANCER)
	{
		res = !countNumAgentByType(AgentTypeEnum::CELL_TYPE_CANCER, count, true);
	}
	else if (t == AgentTypeEnum::CELL_TYPE_MDSC)
	{
		res = !countNumAgentByType(AgentTypeEnum::CELL_TYPE_MDSC, count, true);
	}	
	else if (t == AgentTypeEnum::CELL_TYPE_FIB)
	{
		res = !countNumAgentByType(AgentTypeEnum::CELL_TYPE_FIB, count, true);
	}
	else if (t == AgentTypeEnum::CELL_TYPE_APC)
	{
		res = !countNumAgentByType(AgentTypeEnum::CELL_TYPE_APC, count, true);
	}
	else if (t == AgentTypeEnum::CELL_TYPE_MAC)
	{
		res = !countNumAgentByType(AgentTypeEnum::CELL_TYPE_MAC, count, true);
	}
	else if (t == AgentTypeEnum::CELL_TYPE_T ||
		t == AgentTypeEnum::CELL_TYPE_TCD4)
	{
		int c = countNumAgentByType(AgentTypeEnum::CELL_TYPE_CANCER, count, true);
		int teff = countNumAgentByType(AgentTypeEnum::CELL_TYPE_T, count, false);
		//int tcd4 = countNumAgentByType(AgentTypeEnum::CELL_TYPE_TCD4, count, false);
		res = (c && count < params.getVal(PARAM_N_T_VOXEL_C)) ||
			(!c && count < params.getVal(PARAM_N_T_VOXEL));
	}
	else if (t == AgentTypeEnum::CELL_TYPE_VAS)
	{
		int endo = countNumAgentByType(AgentTypeEnum::CELL_TYPE_VAS, count, false);
		//std::cout << "number of Vas in the voxel: " << count << std::endl;
		res = (count < params.getVal(PARAM_VAS_MAX_PER_VOXEL));
	}
	return res;
}

/*!
  distance to origin (0,0,0).
  used to determine whether this voxel is in range of 
  invasive front when comparing simulation with model.
*/
double TumorGridVoxel::getDistToOrigin(void)const{
	return _distToOrigin;
}

/*!
  content of this voxel. currently a placeholder.
*/
int TumorGridVoxel::getVoxelContent() const{
	int sum = 0;
	for (auto && ag : _agents) {
		sum += ag->getTestValue();
	}
	return sum;
}

void TumorGridVoxel::setVasDensity(double new_density) {
	_vas_density = new_density;
}

};
};