#pragma once

#include "HCC/SP_QSP_HCC/abm/compartment/Tumor.h"
#include "HCC/SP_QSP_HCC/abm/compartment/LymphCentral.h"
#include "FileOutputHub.h"

#include <boost/serialization/nvp.hpp>
//! Simulation central control
/*! Host all compartments (currently only one: tumor).
	Handle interaction between compartments; process global statistics output.
*/
class HCC_Core
{
public:
	HCC_Core();
	HCC_Core(std::string core, std::string margin);
	virtual ~HCC_Core();
	//! setup qsp module
	void setup_qsp(CancerVCT::Param& p);
	//! Initialize all compartments: from file
	void initializeSimulation(std::string, std::string);
	//! Initialize all compartments: random population 
	void initializeSimulation(void);
	//! Proceed on slice in time for all compartments
	virtual void timeSlice(const long slice);
	//! Write ABM stats header to stats file
	void write_stats_header(void) const;
	//! Write ABM stats of current time slice to stats file
	void write_stats_slice(unsigned long slice) const;
	//! write QSP content to file
	void write_QSP(unsigned long slice, bool header)const;
	//! Write ODE stats of current time slice to ode stats file
	virtual void writeOde(unsigned long slice);
	//! Write grid snapshots of current time slice to file 
	virtual void writeGrids(unsigned long slice, unsigned int option);
	//! Print brief summary of current slice
	virtual void briefStats(unsigned long slice);

private:
	friend class boost::serialization::access;
	//! boost serialization 
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);

	//! tumor compartment instance
	SP_QSP_IO::SP_QSP_HCC::Tumor _tumor_core;
	SP_QSP_IO::SP_QSP_HCC::Tumor _tumor_margin;
	SP_QSP_IO::SP_QSP_HCC::LymphCentral _lymph;
};

template<class Archive>
inline void HCC_Core::serialize(Archive & ar, const unsigned int version){
	ar & BOOST_SERIALIZATION_NVP(_tumor_core);
	ar & BOOST_SERIALIZATION_NVP(_tumor_margin);
	ar & BOOST_SERIALIZATION_NVP(_lymph);
	SP_QSP_IO::CellAgent::classSerialize(ar, version);
}