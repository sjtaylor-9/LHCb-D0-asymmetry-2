/* 
 * Author: Suzanne Klaver (suzanne.klaver@cern.ch)
 * September 2016
 * This program simulates ccbar production and fragmentation.
 * It is used to simulate the production asymmetry of Ds mesons.
 *
 * Interesting things to test:
 * - influence of beam crossing angle (should be negligible)
 *   see default LHCb settings 
 * - change ColourReconnection:range (previously BeamRemnants:reconnectRange)
 * - try various ColourReconnection:mode (in particular 3 and 4)
 * - for the Fragmentation there are many models and parameters to understand
*/

// Headers and Namespaces.
#include "Pythia8/Pythia.h" // Include Pythia headers.
#include <TFile.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace Pythia8; // Let Pythia8:: be implicit.

// Begin main program.
int main(int argc, char* argv[]) {

  if (argc < 2){
    std::cout << "Please provide a seed, i.e. 57" << std::endl;
    return 1;
    }

  // Set up generation.
  Pythia pythia; // Declare Pythia object

  // Set seed
  pythia.readString("Random:setSeed = on");
  pythia.readString((std::string("Random:seed = ") + std::string(argv[1])));

  // Processes used to generate ccbar events
  pythia.readString("HardQCD:gg2ccbar = on"); // Switch on process.
  pythia.readString("HardQCD:qqbar2ccbar = on"); // add another process
  pythia.readString("PhaseSpace:pTHatMin = 1"); // 1 GeV
  pythia.readString("PhaseSpace:pTHatMax = 10"); // 10 GeV

  // Beam parameters, assuming no crossing angle
  pythia.readString("Beams:eCM = 13000."); // 8 TeV CM energy; default = 14
  pythia.readString("Beams:AllowMomentumSpread = on"); // same as in LHCb

  // Designed to correct the wrong meson/baryon ratio in LHCb
  // Advised to use mode = 1
  // See slides: https://urldefense.com/v3/__http://indico.cern.ch/event/305160/contributions/__;!!PDiH4ENfjr2_Jw!CYP-wkdYVjGn8O9hwXS00OMoO1YYXxcWWo82qUuMp_YptcVcenxvD1wRYU-CMrGAqtS6sdVkORj2GXG14Hm7pOOVfpluGfgDYmCN5vrCbof-CA$ [indico[.]cern[.]ch]
  // 1673812/attachments/579436/797864/presentation.pdf
  // Also: https://urldefense.com/v3/__https://arxiv.org/pdf/1505.01681v1.pdf__;!!PDiH4ENfjr2_Jw!CYP-wkdYVjGn8O9hwXS00OMoO1YYXxcWWo82qUuMp_YptcVcenxvD1wRYU-CMrGAqtS6sdVkORj2GXG14Hm7pOOVfpluGfgDYmCN5vohfD2g_Q$ [arxiv[.]org]
  // pythia.readString("MultiPartonInteractions:pT0Ref = 2.15"); // try?
  pythia.readString("ColourReconnection:mode = 1"); 
  pythia.readString("ColourReconnection:range = 1.8"); // default

  // Fragmentation: values are very dependent on pTmin!
  // See slides: https://urldefense.com/v3/__https://indico.cern.ch/event/226053/contributions/__;!!PDiH4ENfjr2_Jw!CYP-wkdYVjGn8O9hwXS00OMoO1YYXxcWWo82qUuMp_YptcVcenxvD1wRYU-CMrGAqtS6sdVkORj2GXG14Hm7pOOVfpluGfgDYmCN5vo4wV-yFg$ [indico[.]cern[.]ch]
  // 1532535/attachments/371720/517255/PPTS_LHCb040313.pdf 
  pythia.readString("StringZ:rFactC = 0.67"); // default value for b quarks  
 
  // Spacelike shower: advised to look at effect by Marco Adinolfi
  // Parton shower cut-off mass for pure QED branchings. 
  // Assumed smaller than (or equal to) pTminChgQ (default = 0.5)
  pythia.readString("SpaceShower:pTminChgL = 0.0005"); // default = 0.0005

  // Set number events you want to see in output; default = 1
  pythia.readString("Next:numberShowEvent = 1"); 


  pythia.init(); // Initialize; incoming pp beams is default.
  int nEvents = 2000000; // number of events

  // Create file to write events to
  std::ofstream outputfile;
  outputfile.open((std::string("/afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/pythia_hadronisation/simulated_data")+std::string(argv[1])+std::string(".csv")));
  outputfile << "Event , PID , PT , Y \n";

  // Set rapidity range limits
  double rapidityMin = 2.0; // Minimum rapidity 
  double rapidityMax = 5.0;  // Maximum rapidity

  // Generate events in event loop.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent){
    pythia.next(); // Generate an(other) event. Fill event record.
    // particle loop; properties of each pythia particle event
    int iDs = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
        {
        if (abs( pythia.event[i].id()) == 421){
            iDs = i;
            if (abs( pythia.event[iDs+1].id() ) != 421) {
                // cout << "i = " << iDs << ", pt = " << pythia.event[iDs].pT() << endl;
                // cout << "i = " << iDs << ", y  = " << pythia.event[iDs].y()  << endl;
                if (pythia.event[iDs].y() >= rapidityMin && pythia.event[iDs].y() <= rapidityMax) { // CHatGPT???
                    outputfile << iEvent << "," << pythia.event[iDs].id() << "," << pythia.event[iDs].pT() << "," << pythia.event[iDs].y() << "\n";
                }
            }
        }
    }
  }
  pythia.stat(); // show statistics, cross-sections, problems etc.
  outputfile.close();
 
  return 0;
  // End main program with error-free return.
}
