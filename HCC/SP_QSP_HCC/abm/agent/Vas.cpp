//#include <boost/serialization/export.hpp>
#include "Vas.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(Vas)
#include <iostream>
#include <sstream>

#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"
#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"
#include "../compartment/LymphCentral.h"

#include <iostream>
#include <sstream>

namespace SP_QSP_IO {
namespace SP_QSP_HCC {

	Vas::Vas(SpatialCompartment* c)
		:Cell_Tumor(c)
		, _source_O2(NULL)
		, _sink_VEGFA(NULL)
		, _vas_density(0)
	{

	}

	Vas::Vas(const Vas& c)
		: Cell_Tumor(c)
		, _source_O2(NULL)
		, _sink_VEGFA(NULL)
		, _vas_density(0)
	{

	}

	Vas::~Vas()
	{
	}

	std::string Vas::toString()const {
		std::stringstream ss;
		ss << Cell_Tumor::toString();
		return ss.str();
	}


	bool Vas::agent_movement_step(double t, double dt, Coord& c) {
		bool move = false;
		
		return move;
	}

	bool Vas::agent_state_step(double t, double dt, Coord& c) {
		bool divide = false;
		
		return divide;
	}

	void Vas::remove_all_source_sink(void) {

		remove_source_sink(_source_O2);
		remove_source_sink(_sink_VEGFA);

		return;
	}
};
};