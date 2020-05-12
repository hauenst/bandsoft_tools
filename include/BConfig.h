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

		int   getRunNumber (int index) { return getInt ( Run_order        ,index);}
		int   getEvent     (int index) { return getInt ( Event_order      ,index);}
		int   getUnixTime  (int index) { return getInt ( UTime_order      ,index);}
		int   getTrigger   (int index) { return getInt ( Trigger_order    ,index);}
		int   getTimeStamp (int index) { return getInt ( TimeStamp_order  ,index);}
		int   getType      (int index) { return getInt ( Type_order       ,index);}
		int   getMode      (int index) { return getInt ( Mode_order       ,index);}
		float getTorus     (int index) { return getFloat ( Torus_order    ,index);}
		float getSolenoid  (int index) { return getFloat ( Solenoid_order ,index);}
};

#endif
