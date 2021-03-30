// I want to write this code
// To read hipo file and get the
// charge information

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TVector3.h"

#include "bank.h"
#include "clas12fiducial.h"
#include "reader.h"

#include "BCalorimeter.h"
#include "BEvent.h"
#include "BParticle.h"
#include "BScintillator.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"

using namespace std;

int main(int argc, char **argv) {
  // check number of arguments
  if (argc < 4) {
    cerr << "Incorrect number of arugments. Instead use:\n\t./code "
            "[outputFile] [MC/DATA] [inputFile] \n\n";
    cerr << "\t\t[outputFile] = ____.root\n";
    cerr << "\t\t[<MC,DATA> = <0, 1> \n";
    cerr << "\t\t[inputFile] = ____.hipo ____.hipo ____.hipo ...\n\n";
    return -1;
  }

  int MC_DATA_OPT = atoi(argv[2]);

  // Create output tree
  TFile *outFile = new TFile(argv[1], "RECREATE");

  TTree *outTree = new TTree("Q", "charge information");

  //	Event info:
  int Runno = 0;
  double Ebeam = 0;
  double gated_charge = 0;
  double livetime = 0;
  double starttime = 0;
  double current = 0;
  int eventnumber = 0;
  double fcupgated = 0;
  double fcup = 0;
  int row =0;

  // event branches:

  outTree->Branch("eventnumber", &eventnumber);
  outTree->Branch("gated_charge", &gated_charge);
  outTree->Branch("fcupgated", &fcupgated);
  outTree->Branch("fcup", &fcup);
  outTree->Branch("row", &row );

  // Connect to the RCDB
  rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

  // loading input file

  for (int i = 3; i < argc; i++) {
    TString inputFile = argv[i];

    cout << "working on file: " << inputFile << endl;
    hipo::reader reader;
    //Set tag event
    reader.setTags(0);
    reader.open(inputFile);
    hipo::dictionary factory;
    hipo::schema schema;
    reader.readDictionary(factory);

    // adding bank here
    hipo::bank run_config(factory.getSchema("RUN::config"));
    hipo::bank rec_event(factory.getSchema("REC::Event"));
    hipo::bank scaler(factory.getSchema("RUN::scaler"));

    hipo::event readevent;

    // loop over all even in file
    int event_counter = 0;
    while (reader.next() == true) {

      // clear all branches
      eventnumber = 0;
      gated_charge = 0;
      fcupgated = 0;
      fcup = 0;

      // Count events
      if (event_counter % 1000000 == 0)
        cout << "event: " << event_counter << endl;
      // if (event_counter > 100000)
      //  break;
      event_counter++;

      // loading data structure for this event
      reader.read(readevent);
      readevent.getStructure(run_config);
      readevent.getStructure(rec_event);
      readevent.getStructure(scaler);

      // Get Event number from RUN::config
      eventnumber = run_config.getInt(1, 0);
      gated_charge = rec_event.getFloat(2, 0);
      fcupgated = scaler.getFloat(0, 0);
      fcup = scaler.getFloat(1, 0);
      row = rec_event.getRows();

      outTree->Fill();
    } // end loop over event
  }   // end loop over file

  outFile->cd();
  outTree->Write();
  outFile->Close();

  return 0;
}
