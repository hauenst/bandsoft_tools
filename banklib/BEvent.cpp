#include "BEvent.h"
#include "TVector3.h"
// ==============================================================
BEvent::BEvent(hipo::schema __schema) : hipo::bank(__schema){

	Category_order  = __schema.getEntryOrder("category"  );
	Topo_order      = __schema.getEntryOrder("topology"   );
	BCG_order       = __schema.getEntryOrder("beamCharge");
	LiveTime_order  = __schema.getEntryOrder("liveTime"  );
	StartTime_order = __schema.getEntryOrder("startTime" );
	RFTime_order    = __schema.getEntryOrder("RFTime" );
	Helic_order     = __schema.getEntryOrder("helicity"  );
	HelicRaw_order  = __schema.getEntryOrder("helicityRaw"  );
	ProcTime_order  = __schema.getEntryOrder("procTime"  );

}
// ==============================================================
BEvent::~BEvent(){}
// ==============================================================
