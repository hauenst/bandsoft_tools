#ifndef BBAND_H
#define BBAND_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BBand : public hipo::bank {

	private:

    bool oldfile ;
		//for both BAND::hits (old and new structure)
		int id_order           ;
		int sector_order       ;
		int layer_order        ;
		int component_order    ;
		int x_order            ;
		int y_order            ;
		int z_order            ;
		int ux_order           ; //corresponds to ex in new BAND::hits
		int uy_order           ; //corresponds to ey in new BAND::hits
		int uz_order           ; //corresponds to ez in new BAND::hits
		int meantimeTdc_order  ; //corresponds to time in new BAND::hits
		int meantimeFadc_order ; //corresponds to timeFadc in new BAND::hits
		int difftimeTdc_order  ; //corresponds to difftime in new BAND::hits
		int difftimeFadc_order ; //corresponds to difftimeFadc in new BAND::hits
		//only new BAND::hits structure
		int energy_order       ;
		int indexLpmt_order    ;
		int indexRpmt_order    ;
		int status_order       ;
		//only old BAND::hits structure
		int adcLcorr_order     ;
		int adcRcorr_order     ;
		int tFadcLcorr_order   ;
		int tFadcRcorr_order   ;
		int tTdcLcorr_order    ;
		int tTdcRcorr_order    ;


	public:

		BBand(){};

		BBand(hipo::schema __schema,int __old);
		BBand(hipo::schema __schema);

		~BBand();

		bool  isOldfile       () 					{ return oldfile;}
		int   getId           (int row) { return getInt   (id_order           ,row);}
		int   getSector       (int row) { return getInt   (sector_order       ,row);}
		int   getLayer        (int row) { return getInt   (layer_order        ,row);}
		int   getComponent    (int row) { return getInt   (component_order    ,row);}
		float getMeantimeTdc  (int row) { return getFloat (meantimeTdc_order  ,row);}
		float getMeantimeFadc (int row) { return getFloat (meantimeFadc_order ,row);}
		float getDifftimeTdc  (int row) { return getFloat (difftimeTdc_order  ,row);}
		float getDifftimeFadc (int row) { return getFloat (difftimeFadc_order ,row);}
		float getX            (int row) { return getFloat (x_order            ,row);}
		float getY            (int row) { return getFloat (y_order            ,row);}
		float getZ            (int row) { return getFloat (z_order            ,row);}
		float getUx           (int row) { return getFloat (ux_order           ,row);}
		float getUy           (int row) { return getFloat (uy_order           ,row);}
		float getUz           (int row) { return getFloat (uz_order           ,row);}
		//next 6 functions are exclusively for old BAND::hits structure
		float getAdcLcorr     (int row) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getAdcLcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (adcLcorr_order    ,row);
		}
		float getAdcRcorr     (int row) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getAdcRcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (adcRcorr_order     ,row);
		}
		float getTfadcLcorr   (int row) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getTfadcLcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (tFadcLcorr_order   ,row);
		}
		float getTfadcRcorr   (int row) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getTfadcRcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (tFadcRcorr_order   ,row);
		}
		float getTtdcLcorr    (int row) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getTtdcLcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (tTdcLcorr_order    ,row);
		}
		float getTtdcRcorr    (int row) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getTtdcRcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (tTdcRcorr_order    ,row);
		}
		//next functions are for BAND::hits new structure. some of them are mapped to the same class members than previous get-functions
		float getEx           (int row) { return getFloat (ux_order           ,row);}
		float getEy           (int row) { return getFloat (uy_order           ,row);}
		float getEz           (int row) { return getFloat (uz_order           ,row);}
		float getTime         (int row) { return getFloat (meantimeTdc_order  ,row);}
		float getTimeFadc     (int row) { return getFloat (meantimeFadc_order ,row);}
		//exclusively for new BAND::hits structure
		float getEnergy       (int row) {
			if (oldfile == true) {
				std::cout << "Warning from BBand class: Usage of getEnergy with old file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (energy_order       ,row);
		}
		int   getLpmtindex    (int row) {
			if (oldfile == true) {
				std::cout << "Warning from BBand class: Usage of getLpmtindex with old file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getInt   (indexLpmt_order    ,row);
		}
		int   getRpmtindex    (int row) {
			if (oldfile == true) {
				std::cout << "Warning from BBand class: Usage of getRpmtindex with old file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getInt   (indexRpmt_order    ,row);
		}
		int   getStatus       (int row) {
			if (oldfile == true) {
				std::cout << "Warning from BBand class: Usage of getStatus with old file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getInt   (status_order       ,row);
		}

		int   getBarKey(int row);
};

#endif
