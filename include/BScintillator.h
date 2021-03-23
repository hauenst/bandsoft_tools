#ifndef BSCINTILLATOR_H
#define BSCINTILLATOR_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BScintillator : public hipo::bank {

	private:

		int index_order     ;
		int pindex_order    ;
		int detector_order  ;
		int sector_order    ;
		int layer_order     ;
		int component_order ;
		int energy_order    ;
		int time_order      ;
		int path_order      ;
		int chi2_order      ;
		int x_order         ;
		int y_order         ;
		int z_order         ;
		int hx_order        ;
		int hy_order        ;
		int hz_order        ;
		int status_order    ;

	public:

		BScintillator(){};

		BScintillator(hipo::schema __schema);

		~BScintillator();

		int   getIndex	  (int row) { return getInt   ( index_order     ,row);}
		int   getPindex   (int row) { return getInt   ( pindex_order    ,row);}
		int   getDetector (int row) { return getInt   ( detector_order  ,row);}
		int   getSector   (int row) { return getInt   ( sector_order    ,row);}
		int   getLayer    (int row) { return getInt   ( layer_order     ,row);}
		int   getComponent(int row) { return getInt   ( component_order ,row);}
		float getEnergy   (int row) { return getFloat ( energy_order    ,row);}
		float getTime     (int row) { return getFloat ( time_order      ,row);}
		float getPath     (int row) { return getFloat ( path_order      ,row);}
		float getChi2     (int row) { return getFloat ( chi2_order      ,row);}
		float getX        (int row) { return getFloat ( x_order         ,row);}
		float getY        (int row) { return getFloat ( y_order         ,row);}
		float getZ        (int row) { return getFloat ( z_order         ,row);}
		float getHx       (int row) { return getFloat ( hx_order        ,row);}
		float getHy       (int row) { return getFloat ( hy_order        ,row);}
		float getHz       (int row) { return getFloat ( hz_order        ,row);}
		int   getStatus   (int row) { return getInt   ( status_order    ,row);}

};

#endif
