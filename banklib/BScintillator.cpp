#include "BScintillator.h"
#include "TVector3.h"
// ==============================================================
BScintillator::BScintillator(hipo::schema __schema) : hipo::bank(__schema){

	index_order     = __schema.getEntryOrder("index"    );
	pindex_order    = __schema.getEntryOrder("pindex"   );
	detector_order  = __schema.getEntryOrder("detector" );
	sector_order    = __schema.getEntryOrder("sector"   );
	layer_order     = __schema.getEntryOrder("layer"    );
	component_order = __schema.getEntryOrder("component");
	energy_order    = __schema.getEntryOrder("energy"   );
	time_order      = __schema.getEntryOrder("time"     );
	path_order      = __schema.getEntryOrder("path"     );
	chi2_order      = __schema.getEntryOrder("chi2"     );
	x_order         = __schema.getEntryOrder("x"        );
	y_order         = __schema.getEntryOrder("y"        );
	z_order         = __schema.getEntryOrder("z"        );
	hx_order        = __schema.getEntryOrder("hx"       );
	hy_order        = __schema.getEntryOrder("hy"       );
	hz_order        = __schema.getEntryOrder("hz"       );
	status_order    = __schema.getEntryOrder("status"   );

}
// ==============================================================
BScintillator::~BScintillator(){}
// ==============================================================
