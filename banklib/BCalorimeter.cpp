#include "BCalorimeter.h"
#include "TVector3.h"
// ==============================================================
BCalorimeter::BCalorimeter(hipo::schema __schema) : hipo::bank(__schema){

	index_order    = __schema.getEntryOrder("index"    );
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
int BCalorimeter::getElectronSector (int pindex)
{
	//Sector of PCAL hit layer 1
	int nCal = getRows();
	for(int i = 0 ; i < nCal ; i++){
		if(BCalorimeter::getPindex(i)==pindex&&BCalorimeter::getLayer(i)==1 && BCalorimeter::getDetector(i)==7){
			return BCalorimeter::getSector(i);
		}
	}
	return -1;
}
// ==============================================================
int BCalorimeter::getPcalRow (int pindex)
{
	//Find row of BANK with PCAL for Particle Index
	int nCal = getRows();
	for(int i = 0 ; i < nCal ; i++){
		if(BCalorimeter::getPindex(i)==pindex&&BCalorimeter::getLayer(i)==1 && BCalorimeter::getDetector(i)==7){
			return i;
		}
	}
	return -1;
}
// ==============================================================
float BCalorimeter::getPcalE (int pindex)
{
	int nCal = getRows();
	for(int i = 0 ; i < nCal ; i++){
		if(BCalorimeter::getPindex(i)==pindex&&BCalorimeter::getLayer(i)==1 && BCalorimeter::getDetector(i)==7){
			return BCalorimeter::getEnergy(i);
		}
	}
	return 0;
}
// ==============================================================
float BCalorimeter::getECinE (int pindex)
{
	int nCal = getRows();
	for(int i = 0 ; i < nCal ; i++){
		if(BCalorimeter::getPindex(i)==pindex&&BCalorimeter::getLayer(i)==4 && BCalorimeter::getDetector(i)==7){
			return BCalorimeter::getEnergy(i);
		}
	}
	return 0;
}
// ==============================================================
float BCalorimeter::getECoutE(int pindex)
{
	int nCal = getRows();
	for(int i = 0 ; i < nCal ; i++){
		if(BCalorimeter::getPindex(i)==pindex&&BCalorimeter::getLayer(i)==7 && BCalorimeter::getDetector(i)==7){
			return BCalorimeter::getEnergy(i);
		}
	}
	return 0;
}
// ==============================================================
float BCalorimeter::getTotE  (int pindex)
{
	return BCalorimeter::getPcalE(pindex) + BCalorimeter::getECinE(pindex) + BCalorimeter::getECoutE(pindex);
}
// ==============================================================
