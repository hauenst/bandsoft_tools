#include "BScaler.h"
#include "TVector3.h"
// ==============================================================
BScaler::BScaler(hipo::schema __schema) : hipo::bank(__schema){

	FCupGated_order  = __schema.getEntryOrder("fcupgated"  );
	FCup_order       = __schema.getEntryOrder("fcup"   );
	LiveTime_order   = __schema.getEntryOrder("livetime");
}
// ==============================================================
BScaler::~BScaler(){}
// ==============================================================
