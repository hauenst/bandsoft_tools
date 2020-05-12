#include "BCalorimeter.h"
#include "TVector3.h"
// ==============================================================
BCalorimeter::BCalorimeter(hipo::schema __schema) : hipo::bank(__schema){

	pindex_order   = __schema.getEntryOrder("pindex"   );
	detector_order = __schema.getEntryOrder("detector" );
	sector_order   = __schema.getEntryOrder("sector"   );
	layer_order    = __schema.getEntryOrder("layer"    );
	energy_order   = __schema.getEntryOrder("energy"   );
	time_order     = __schema.getEntryOrder("time"     );
	path_order     = __schema.getEntryOrder("path"     );
	x_order        = __schema.getEntryOrder("x"        );
	y_order        = __schema.getEntryOrder("y"        );
	z_order        = __schema.getEntryOrder("z"        );
	lu_order       = __schema.getEntryOrder("lu"       );
	lv_order       = __schema.getEntryOrder("lv"       );
	lw_order       = __schema.getEntryOrder("lw"       );
}

// ==============================================================
BCalorimeter::~BCalorimeter(){}
// ==============================================================
float BCalorimeter::getPcalE (int index)
{
	int nCal = getSize();
	for(int i = 0 ; i < nCal ; i++){
		if(BCalorimeter::getIndex(i)==index&&BCalorimeter::getLayer(i)==1){
			return BCalorimeter::getEnergy(i);
		}
	}
	return 0;
}
// ==============================================================
float BCalorimeter::getECinE (int index)
{
	int nCal = getSize();
        for(int i = 0 ; i < nCal ; i++){
                if(BCalorimeter::getIndex(i)==index&&BCalorimeter::getLayer(i)==4){
                        return BCalorimeter::getEnergy(i);
                }
        }
        return 0;
}
// ==============================================================
float BCalorimeter::getECoutE(int index)
{
	int nCal = getSize();
        for(int i = 0 ; i < nCal ; i++){
                if(BCalorimeter::getIndex(i)==index&&BCalorimeter::getLayer(i)==7){
                        return BCalorimeter::getEnergy(i);
                }
        }
        return 0;
}
// ==============================================================
float BCalorimeter::getTotE  (int index)
{
	return BCalorimeter::getPcalE(index) + BCalorimeter::getECinE(index) + BCalorimeter::getECoutE(index);
}
// ==============================================================
