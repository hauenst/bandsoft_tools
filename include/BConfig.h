#ifndef BCONFIG_H
#define BCONFIG_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BConfig : public hipo::bank {

	private:

		int Run_order  		 ;  //Run Number
		int Event_order    ;  //Event Number
		int UTime_order    ;  //Unix time
		int Trigger_order  ;  //Trigger Bits
		int TimeStamp_order;  //Time Stap from Trigger Interface board
		int Type_order     ;  //Type of the Run
		int Mode_order     ;  //Run mode
		int Torus_order    ;  //Torus setting value (-1 to 1)
		int Solenoid_order ;  //Solenoid setting value (-1 to 1)

	public:

		BConfig(){};

		BConfig(hipo::schema __schema);

		~BConfig();

		int   getRunNumber (int row) { return getInt ( Run_order        ,row);}
		int   getEvent     (int row) { return getInt ( Event_order      ,row);}
		int   getUnixTime  (int row) { return getInt ( UTime_order      ,row);}
		int   getTrigger   (int row) { return getInt ( Trigger_order    ,row);}
		int   getTimeStamp (int row) { return getInt ( TimeStamp_order  ,row);}
		int   getType      (int row) { return getInt ( Type_order       ,row);}
		int   getMode      (int row) { return getInt ( Mode_order       ,row);}
		float getTorus     (int row) { return getFloat ( Torus_order    ,row);}
		float getSolenoid  (int row) { return getFloat ( Solenoid_order ,row);}
};

#endif
