#ifndef BSCALER_H
#define BSCALER_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BScaler : public hipo::bank {

	private:

		int FCupGated_order  ;  //Beam charge, integrated from beginning of run, gated in nano-Coulombs
		int FCup_order   		 ;  //Beam charge, integrated from beginning of run, ungated in nano-Coulombs
		int LiveTime_order   ;  //Livetime during one scaler period

	public:

		BScaler(){};

		BScaler(hipo::schema __schema);

		~BScaler();


		float getFCupGated (int index) { return getFloat ( FCupGated_order,index);}
		float getFCup      (int index) { return getFloat ( FCup_order     ,index);}
		float getLiveTime  (int index) { return getFloat ( LiveTime_order ,index);}
};

#endif
