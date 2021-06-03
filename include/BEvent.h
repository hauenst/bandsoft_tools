#ifndef BEVENT_H
#define BEVENT_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BEvent : public hipo::bank {

	private:

		int Category_order  ;
		int Topo_order      ;
		int BCG_order       ;
		int LiveTime_order  ;
		int StartTime_order ;
		int RFTime_order    ;
		int Helic_order     ;
		int HelicRaw_order  ;
		int ProcTime_order  ;


	public:

		BEvent(){};

		BEvent(hipo::schema __schema);

		~BEvent();

		int   getCategory(int row) { return getInt   ( Category_order   ,row);}
		int   getTopo    (int row) { return getInt   ( Topo_order       ,row);}
		float getBCG     (int row) { return getFloat ( BCG_order        ,row);}
		float getLT      (int row) { return getFloat ( LiveTime_order   ,row);}
		float getSTTime  (int row) { return getFloat ( StartTime_order  ,row);}
		float getRFTime  (int row) { return getFloat ( RFTime_order     ,row);}
		int   getHelic   (int row) { return getInt   ( Helic_order      ,row);}
		int   getHelicRaw(int row) { return getInt   ( HelicRaw_order   ,row);}
		float getPTime   (int row) { return getFloat ( ProcTime_order   ,row);}
};

#endif
