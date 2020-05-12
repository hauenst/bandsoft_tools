#include "BBand.h"
#include "TVector3.h"
// ==============================================================
//Constructor for old hipo files, overloaded with extra input integer
BBand::BBand(hipo::schema __schema, int __old) : hipo::bank(__schema){

  bool oldfile = true;
  std::cout << "BBand constructor using old Band::hits structure files " << std::endl;
	id_order           = __schema.getEntryOrder("id"          );
	sector_order       = __schema.getEntryOrder("sector"      );
	layer_order        = __schema.getEntryOrder("layer"       );
	component_order    = __schema.getEntryOrder("component"   );
	x_order            = __schema.getEntryOrder("x"           );
	y_order            = __schema.getEntryOrder("y"           );
	z_order            = __schema.getEntryOrder("z"           );
	ux_order           = __schema.getEntryOrder("ux"          );
	uy_order           = __schema.getEntryOrder("uy"          );
	uz_order           = __schema.getEntryOrder("uz"          );
	meantimeTdc_order  = __schema.getEntryOrder("meantimeTdc" );
	meantimeFadc_order = __schema.getEntryOrder("meantimeFadc");
	difftimeTdc_order  = __schema.getEntryOrder("difftimeTdc" );
	difftimeFadc_order = __schema.getEntryOrder("difftimeFadc");
	adcLcorr_order     = __schema.getEntryOrder("adcLcorr"    );
	adcRcorr_order     = __schema.getEntryOrder("adcRcorr"    );
	tFadcLcorr_order   = __schema.getEntryOrder("tFadcLcorr"  );
	tFadcRcorr_order   = __schema.getEntryOrder("tFadcRcorr"  );
	tTdcLcorr_order    = __schema.getEntryOrder("tTdcLcorr"   );
	tTdcRcorr_order    = __schema.getEntryOrder("tTdcRcorr"   );
	energy_order       = -1;
	indexLpmt_order	   = -1;
	indexRpmt_order    = -1;
	status_order			 = -1;

}
//STANDARD constructor for new files
BBand::BBand(hipo::schema __schema) : hipo::bank(__schema){

  bool oldfile = false;
  std::cout << "BBand constructor using new Band::hits structure files " << std::endl;
	id_order           = __schema.getEntryOrder("id"          );
	sector_order       = __schema.getEntryOrder("sector"      );
	layer_order        = __schema.getEntryOrder("layer"       );
	component_order    = __schema.getEntryOrder("component"   );
	x_order            = __schema.getEntryOrder("x"           );
	y_order            = __schema.getEntryOrder("y"           );
	z_order            = __schema.getEntryOrder("z"           );
	ux_order           = __schema.getEntryOrder("ex"          );
	uy_order           = __schema.getEntryOrder("ey"          );
	uz_order           = __schema.getEntryOrder("ez"          );
	meantimeTdc_order  = __schema.getEntryOrder("time"        );
	meantimeFadc_order = __schema.getEntryOrder("timeFadc"    );
	difftimeTdc_order  = __schema.getEntryOrder("difftime"    );
	difftimeFadc_order = __schema.getEntryOrder("difftimeFadc");
	energy_order       = __schema.getEntryOrder("energy"      );
	indexLpmt_order	   = __schema.getEntryOrder("indexLpmt"   );
	indexRpmt_order    = __schema.getEntryOrder("indexRpmt"   );
	status_order			 = __schema.getEntryOrder("status"      );
	adcLcorr_order     = -1;
	adcRcorr_order     = -1;
	tFadcLcorr_order   = -1;
	tFadcRcorr_order   = -1;
	tTdcLcorr_order    = -1;
	tTdcRcorr_order    = -1;

}

// ==============================================================
BBand::~BBand(){}
// ==============================================================
int BBand::getBarKey(int index)
{
	int s = getSector   (index);
	int l = getLayer    (index);
	int c = getComponent(index);
	return s*100+l*10+c;
}
