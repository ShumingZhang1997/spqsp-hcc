#pragma once

#include "Cell_Tumor.h"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include "../../pde/DiffuseGrid.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

class MDSC:
	public Cell_Tumor
{
public:
	MDSC(){};
	MDSC(SpatialCompartment* c );
	MDSC(const MDSC& c);
	virtual ~MDSC();

	virtual CellAgent* createCellCopy() const { return new MDSC(*this); };

	//! print cancer cell information
	virtual std::string toString() const;

	static int getMDSCLife();
	//void setDead(void);

	//! step function for cancer cell
	bool agent_movement_step(double t, double dt, Coord& c);
	bool agent_state_step(double t, double dt, Coord& c);

	//! move sources (NO and ArgI)
	void move_all_source_sink(void)const;

	//! remove all sources
	void remove_all_source_sink(void);

	//! get cell agent type, in this case cancer cell.
	virtual AgentType getType() const { return AgentTypeEnum::CELL_TYPE_MDSC; };

	BioFVMSinkSource*& get_source_ArgI(void) { return _source_ArgI; };
	BioFVMSinkSource*& get_source_NO(void) { return _source_NO; };
	BioFVMSinkSource*& get_sink_CCL2(void) { return _sink_CCL2; };

private:
	friend class boost::serialization::access;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);

	BioFVMSinkSource* _source_ArgI;
	BioFVMSinkSource* _source_NO;
	BioFVMSinkSource* _sink_CCL2;
	//BioFVMSinkSource* _source_IL_10;
};

//BOOST_CLASS_EXPORT_KEY(MDSC)

template<class Archive>
inline void MDSC::serialize(Archive & ar, const unsigned int /* version */){
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Cell_Tumor);
	ar & BOOST_SERIALIZATION_NVP(_source_ArgI);
	ar & BOOST_SERIALIZATION_NVP(_source_NO);
	ar & BOOST_SERIALIZATION_NVP(_sink_CCL2);
	//ar & BOOST_SERIALIZATION_NVP(_source_IL_10);
}

};
};
