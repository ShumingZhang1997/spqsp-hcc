#pragma once

#include "Cell_Tumor.h"


#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

class CancerCell :
	public Cell_Tumor
{
public:
	static int getSenescentLife(void);

	//! default constructor for serialization
	CancerCell();
	CancerCell(SpatialCompartment* c );
	CancerCell(const CancerCell& c);
	virtual ~CancerCell();

	virtual CellAgent* createCellCopy() const { return new CancerCell(*this); };

	//! print cancer cell information
	virtual std::string toString() const;
	//! step function for cancer cell
	//virtual void agentStep(double t, double dt, AgentStep & as);
	virtual bool agent_movement_step(double t, double dt, Coord& c);
	virtual bool agent_state_step(double t, double dt, Coord& c);

	//! move sources (CCL2)
	void move_all_source_sink(void)const;

	//! remove all sources
	void remove_all_source_sink(void);

	//! molecular step; doing nothing currently
	virtual void molecularStep(double t, double dt){};
	//! ODE step overiding base class empty version
	virtual void odeStep(double t, double dt){};
	
	//! Cancer cell to stream
	friend std::ostream & operator<<(std::ostream &os, const CancerCell & c);

	//! set cancer cell state
	//! get cell agent type, in this case cancer cell.
	virtual AgentType getType() const { return AgentTypeEnum::CELL_TYPE_CANCER; };

	void inc_neighbor_Teff(void) { _count_neighbor_Teff++; };
	void inc_neighbor_MacM1(void){ _count_neighbor_MacM1++; };

	//! remarks to string
	virtual std::string getRemark() const;
	//! return count down of time steps to the next division
	int getDivideCD(void) const{ return _divideCD; } ;
	//! change division countdown
	void setDivideCD(int cd) { _divideCD = cd; };
	//! set cancer cell to progenitor
	void setProgenitor();
	//! set cancer cell to senescent
	void setSenescent();
	//! set remaining division for progenitors
	void setDivCounter(int div){ _divideCountRemaining = div; };
	//! randomize division cooldown 
	void randomize_div_cd(int mean);

	BioFVMSinkSource*& get_source_TGFB(void) { return _source_TGFB; };
	BioFVMSinkSource*& get_source_CCL2(void) { return _source_CCL2; };
	BioFVMSinkSource*& get_source_VEGFA(void) { return _source_VEGFA; };
	BioFVMSinkSource*& get_sink_O2(void) { return _sink_O2; };
	BioFVMSinkSource*& get_sink_IFNg(void) { return _sink_IFNg; };
	//const ShapeSingle* getCellShape()const{ return &CancerCell::_shapeC; };

	unsigned int _stemID;

private:
	friend class boost::serialization::access;
	//! boost serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);

	BioFVMSinkSource* _source_CCL2;
	BioFVMSinkSource* _source_TGFB;
	BioFVMSinkSource* _source_VEGFA;
	BioFVMSinkSource* _sink_O2;
	BioFVMSinkSource* _sink_IFNg;
	//! next division cool down timer
	int _divideCD;
	//! allowed to divide (currently always true for cancer cells)
	bool _divideFlag;
	//! remaining divisions, for progenitor cells.
	int _divideCountRemaining;

	//static ShapeSingle _shapeC;

	/* These variables do not need to be serialized,
		but should be set to 0 with default constructor
		during deserialization.
	*/
	//! number of Teff in neighborhood
	int _count_neighbor_Teff;
	//! number of M1 macrophage in neighborhood
	int _count_neighbor_MacM1;

};

//BOOST_CLASS_EXPORT_KEY(CancerCell)

inline std::ostream & operator<<(std::ostream &os, const CancerCell & c){
	os << c.getID();
	return os;
}

template<class Archive>
inline void CancerCell::serialize(Archive & ar, const unsigned int /* version */){
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Cell_Tumor);
	ar & BOOST_SERIALIZATION_NVP(_source_CCL2);
	ar & BOOST_SERIALIZATION_NVP(_source_TGFB);
	ar & BOOST_SERIALIZATION_NVP(_source_VEGFA);
	ar & BOOST_SERIALIZATION_NVP(_sink_O2);
	ar & BOOST_SERIALIZATION_NVP(_sink_IFNg);
	ar & BOOST_SERIALIZATION_NVP(_stemID);
	ar & BOOST_SERIALIZATION_NVP(_divideCD);
	ar & BOOST_SERIALIZATION_NVP(_divideFlag);
	ar & BOOST_SERIALIZATION_NVP(_divideCountRemaining);
}

};
};


