#include "BConfig.h"
#include "TVector3.h"
// ==============================================================
BConfig::BConfig(hipo::schema __schema) : hipo::bank(__schema){

	Run_order  	= __schema.getEntryOrder("run"  );
	Event_order     = __schema.getEntryOrder("event"   );
	UTime_order     = __schema.getEntryOrder("unixtime");
	Trigger_order   = __schema.getEntryOrder("trigger"  );
	TimeStamp_order = __schema.getEntryOrder("timestamp" );
	Type_order      = __schema.getEntryOrder("type" );
	Mode_order      = __schema.getEntryOrder("mode"  );
	Torus_order     = __schema.getEntryOrder("torus"  );
	Solenoid_order  = __schema.getEntryOrder("solenoid"  );
}
// ==============================================================
BConfig::~BConfig(){}
// ==============================================================
