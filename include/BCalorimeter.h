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
		int index_order    ;
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

		int   getIndex	 (int row) { return getInt  (index_order     ,row);}
		int   getPindex  (int row) { return getInt  (pindex_order    ,row);}
		int   getDetector(int row) { return getInt  (detector_order  ,row);}
		int   getSector  (int row) { return getInt  (sector_order    ,row);}
		int   getLayer   (int row) { return getInt  (layer_order     ,row);}
		float getEnergy  (int row) { return getFloat(energy_order    ,row);}
		float getLU      (int row) { return getFloat(lu_order        ,row);}
		float getLV      (int row) { return getFloat(lv_order        ,row);}
		float getLW      (int row) { return getFloat(lw_order        ,row);}


		int   getElectronSector  (int pindex);
		int   getPcalRow( int pindex);
		float getPcalE (int pindex);
		float getECinE (int pindex);
		float getECoutE(int pindex);
		float getTotE  (int pindex);
};

#endif
