#ifndef BCALORIMETER_H
#define BCALORIMETER_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BCalorimeter : public hipo::bank {

	private:
		int pindex_order   ;
		int detector_order ;
		int sector_order   ;
		int layer_order    ;
		int energy_order   ;
		int time_order     ;
		int path_order     ;
		int x_order        ;
		int y_order        ;
		int z_order        ;
		int lu_order       ;
		int lv_order       ;
		int lw_order       ;

	public:

		BCalorimeter(){};

		BCalorimeter(hipo::schema __schema);

		~BCalorimeter();

		int   getIndex   (int index) { return getInt  (pindex_order    ,index);}
		int   getDetector(int index) { return getInt  (detector_order  ,index);}
		int   getSector  (int index) { return getInt  (sector_order    ,index);}
		int   getLayer   (int index) { return getInt  (layer_order     ,index);}
		float getEnergy  (int index) { return getFloat(energy_order    ,index);}
		float getLU      (int index) { return getFloat(lu_order        ,index);}
		float getLV      (int index) { return getFloat(lv_order        ,index);}
		float getLW      (int index) { return getFloat(lw_order        ,index);}

		float getPcalE (int pindex);
		float getECinE (int pindex);
		float getECoutE(int pindex);
		float getTotE  (int pindex);
};

#endif
