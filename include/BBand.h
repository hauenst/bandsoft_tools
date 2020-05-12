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
		int   getId           (int index) { return getInt   (id_order           ,index);}
		int   getSector       (int index) { return getInt   (sector_order       ,index);}
		int   getLayer        (int index) { return getInt   (layer_order        ,index);}
		int   getComponent    (int index) { return getInt   (component_order    ,index);}
		float getMeantimeTdc  (int index) { return getFloat (meantimeTdc_order  ,index);}
		float getMeantimeFadc (int index) { return getFloat (meantimeFadc_order ,index);}
		float getDifftimeTdc  (int index) { return getFloat (difftimeTdc_order  ,index);}
		float getDifftimeFadc (int index) { return getFloat (difftimeFadc_order ,index);}
		float getX            (int index) { return getFloat (x_order            ,index);}
		float getY            (int index) { return getFloat (y_order            ,index);}
		float getZ            (int index) { return getFloat (z_order            ,index);}
		float getUx           (int index) { return getFloat (ux_order           ,index);}
		float getUy           (int index) { return getFloat (uy_order           ,index);}
		float getUz           (int index) { return getFloat (uz_order           ,index);}
		//next 6 functions are exclusively for old BAND::hits structure
		float getAdcLcorr     (int index) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getAdcLcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (adcLcorr_order    ,index);
		}
		float getAdcRcorr     (int index) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getAdcRcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (adcRcorr_order     ,index);
		}
		float getTfadcLcorr   (int index) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getTfadcLcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (tFadcLcorr_order   ,index);
		}
		float getTfadcRcorr   (int index) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getTfadcRcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (tFadcRcorr_order   ,index);
		}
		float getTtdcLcorr    (int index) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getTtdcLcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (tTdcLcorr_order    ,index);
		}
		float getTtdcRcorr    (int index) {
			if (oldfile == false) {
				std::cout << "Warning from BBand class: Usage of getTtdcRcorr with new file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (tTdcRcorr_order    ,index);
		}
		//next functions are for BAND::hits new structure. some of them are mapped to the same class members than previous get-functions
		float getEx           (int index) { return getFloat (ux_order           ,index);}
		float getEy           (int index) { return getFloat (uy_order           ,index);}
		float getEz           (int index) { return getFloat (uz_order           ,index);}
		float getTime         (int index) { return getFloat (meantimeTdc_order  ,index);}
		float getTimeFadc     (int index) { return getFloat (meantimeFadc_order ,index);}
		//exclusively for new BAND::hits structure
		float getEnergy       (int index) {
			if (oldfile == true) {
				std::cout << "Warning from BBand class: Usage of getEnergy with old file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getFloat (energy_order       ,index);
		}
		int   getLpmtindex    (int index) {
			if (oldfile == true) {
				std::cout << "Warning from BBand class: Usage of getLpmtindex with old file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getInt   (indexLpmt_order    ,index);
		}
		int   getRpmtindex    (int index) {
			if (oldfile == true) {
				std::cout << "Warning from BBand class: Usage of getRpmtindex with old file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getInt   (indexRpmt_order    ,index);
		}
		int   getStatus       (int index) {
			if (oldfile == true) {
				std::cout << "Warning from BBand class: Usage of getStatus with old file or constructor was used wrong. Return 0 " << std::endl;
				return 0;
			}
			return getInt   (status_order       ,index);
		}

		int   getBarKey(int index);
};

#endif
