#include "e_pid.h"

#include "TF1.h"
#include "TString.h"
#include "TMath.h"

#include "clashit.h"

#include <fstream>
#include <unistd.h>

#define PID_DIR _VAR

using namespace std;

e_pid::e_pid(){

  interval = 3.5;
  color = 1;
}

e_pid::~e_pid(){
}

void e_pid::setParamsRGB(double Ebeam){

  //Determine correct file name
  char fileName[100];
  if(std::abs(Ebeam-10.6) < 0.001){
    sprintf(fileName,"%s/SFCutParams_106_RGB.dat",std::string(PID_DIR).c_str()); 
  }
  else if(std::abs(Ebeam-10.2) < 0.001){
    sprintf(fileName,"%s/SFCutParams_102_RGB.dat",std::string(PID_DIR).c_str()); 
  }  
  else{
    std::cout<<"Attempting to set a beam energy\n"
	     <<"without a defined RGB fiducal cut\n"
	     <<"\t Using E=10.6 GeV as Default. \n\n\n";
    sprintf(fileName,"%s/SFCutParams_106_RGB.dat",std::string(PID_DIR).c_str()); 
  }      

  //Load file data
  ifstream paramFile(fileName);
  if(!paramFile.is_open()){
    cout<<"Fiducial cut parameter files failed to load.\n"
	<<"Aborting...\n\n\n";
    exit(-2);
  }

  for(int i = 0; i < 6; i++){
    for(int j = 0; j < 8; j++){
      params[i][j] = 0;
      paramFile >> params[i][j];
    }
  }

  paramFile.close();

}

void e_pid::setParamsRGA(){

  //Determine correct file name
  char fileName[100];
  sprintf(fileName,"%s/SFCutParams_RGA.dat",std::string(PID_DIR).c_str()); 

  //Load file data
  ifstream paramFile(fileName);
  if(!paramFile.is_open()){
    cout<<"Fiducial cut parameter files failed to load.\n"
	<<"Aborting...\n\n\n";
    exit(-2);
  }

  for(int i = 0; i < 6; i++){
    for(int j = 0; j < 8; j++){
      params[i][j] = 0;
      paramFile >> params[i][j];
    }
  }

  paramFile.close();

}

void e_pid::setInterval(double newInterval){
  interval = newInterval;
}

bool e_pid::isElectron(clashit * eHit){
  
    bool passSFtot = SF_Epcal_Cut(eHit->getSector()-1,eHit->getEoP(),eHit->getEpcal());
    bool passSFpi = SFpcal_SFecin_Cut(eHit->getMomentum(),eHit->getEpcal(),eHit->getEecin());
    if(eHit->getPID() != 11){ return false; }
    if(eHit->getCharge() != -1){ return false; }
    if(!passSFtot){ return false; }
    if(!passSFpi){ return false; }
    return true;

}

bool e_pid::SF_Epcal_Cut(int sector, double SF, double Epcal){
  if(SF < minSF(sector,Epcal)){
    return false;
  }  
  if(SF > maxSF(sector,Epcal)){
    return false;
  }  
  return true;
}

bool e_pid::SFpcal_SFecin_Cut(double p, double Epcal, double Eecin){
  
  //This cut is only for high momentum particles
  if(p < 4.5){ return true; }

  double SFpi = (Epcal + Eecin) / p;
  if(SFpi < 0.2){ return false; }
  return true;

}

double e_pid::meanSF(int sector, double Epcal){
  return FF(Epcal,params[sector][0],params[sector][1],params[sector][2],params[sector][3]);
}

double e_pid::sigmaSF(int sector, double Epcal){
  return FF(Epcal,params[sector][4],params[sector][5],params[sector][6],params[sector][7]);
}

double e_pid::maxSF(int sector, double Epcal){
  return meanSF(sector,Epcal) + interval * sigmaSF(sector,Epcal);
}

double e_pid::minSF(int sector, double Epcal){
  return meanSF(sector,Epcal) - interval * sigmaSF(sector,Epcal);
}

double e_pid::FF(double E, double sf1, double sf2, double sf3, double sf4){
  return sf1 * (sf2 + (sf3/E) + (sf4/(E*E)) ); 
}
