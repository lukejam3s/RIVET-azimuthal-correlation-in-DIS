// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Math/Constants.hh"

namespace Rivet {


  /// @Azimuthal correlation angle between scattered lepton and leading jet at ZEUS
  class ZEUS_2024_I2794054 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_2024_I2794054);
    /// @name Analysis methods

    /// @author Luke Jones


    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
	declare(FinalState(), "FS");
	declare(DISFinalState(DISFrame::LAB), "DISFS");
	declare(FastJets(DISFinalState(DISFrame::LAB), fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::Et_scheme, 1.0), "DISFSJets");
        const DISLepton dl;
	declare(dl, "Lepton");
	declare(DISKinematics(), "Kinematics");

	// Booking all d sigma / d phi graphs
        book(_h_dsigdphi[0], 1, 1, 1);
	book(_h_dsigdphi[1], 2, 1, 1);
	book(_h_dsigdphi[2], 3, 1, 1);
	book(_h_dsigdphi[3], 4, 1, 1);
        book(_h_dsigdphi[4], 5, 1, 1);
        book(_h_dsigdphi[5], 6, 1, 1);
        book(_h_dsigdphi[6], 7, 1, 1);
        book(_h_dsigdphi[7], 8, 1, 1);
        book(_h_dsigdphi[8], 9, 1, 1);
        book(_h_dsigdphi[9], 10, 1, 1);
        book(_h_dsigdphi[10], 11, 1, 1);
        book(_h_dsigdphi[11], 12, 1, 1);
        book(_h_dsigdphi[12], 13, 1, 1);
        book(_h_dsigdphi[13], 14, 1, 1);
        book(_h_dsigdphi[14], 15, 1, 1);
        book(_h_dsigdphi[15], 16, 1, 1);
        book(_h_dsigdphi[16], 17, 1, 1);
        book(_h_dsigdphi[17], 18, 1, 1);
        book(_h_dsigdphi[18], 19, 1, 1);

	}

	void analyze(const Event& event) {



 // Get the DIS Kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      if ( dk.failed() ) vetoEvent;
      const int orientation = dk.orientation();
      double y = dk.y();
      double Q2  = dk.Q2();


	//  Momentum of scattered lepton
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");
      if ( dl.failed() )  vetoEvent;
      const FourMomentum leptonMom = dl.out();
      const double enel = leptonMom.E();
      const double thel = leptonMom.angle(dk.beamHadron().mom())/degree;


        // Jet Pseudorapidity
        double etamin = -1.5;
        double etamax = 1.8;
        if (orientation  < 0) {
            etamin = -1.8;
            etamax = 1.5;

	}

        // Basic kinematic cuts
        if (!inRange(Q2, 10, 350)) vetoEvent;
        if (!inRange(y, 0.04, 0.7)) vetoEvent;
        if (enel < 10) vetoEvent;
        if (thel < 140.0 || thel > 180.0) vetoEvent;  // Cut on lepton polar angle

      const Jets jets = apply<FastJets>(event, "DISFSJets").jets(Cuts::Et > 2.5*GeV && Cuts::Et < 30*GeV && Cuts::etaIn(etamin, etamax), cmpMomByEt);

	double dPhi2 = 0.0;
        double dPhi3 = 0.0;
	if (jets.size() < 1) vetoEvent;

	const Jet& firstJet = jets[0];
	const FourMomentum jetMom = firstJet.mom();

	if (jets.size()  >= 2) {
        dPhi2 = deltaPhi(leptonMom, jetMom);

	}

        if (jets.size()  >= 3) {
        dPhi3 = deltaPhi(leptonMom, jetMom);

	}


	//Azimuthal angle calculations
	double dPhi = deltaPhi(leptonMom, jetMom);


	// Fill Histograms
	_h_dsigdphi[0]->fill(dPhi);

  if (firstJet.pT() >  2.5*GeV && firstJet.pT() < 7*GeV)
        _h_dsigdphi[1]->fill(dPhi);

  if (firstJet.pT() >  2.5*GeV && firstJet.pT() < 7*GeV && jets.size() >= 2)
        _h_dsigdphi[2]->fill(dPhi2);

  if (firstJet.pT() >  2.5*GeV && firstJet.pT() < 7*GeV && jets.size() >= 3)
        _h_dsigdphi[3]->fill(dPhi3);

  if (firstJet.pT() >  7*GeV && firstJet.pT() < 12*GeV)
        _h_dsigdphi[4]->fill(dPhi);

  if (firstJet.pT() >  7*GeV && firstJet.pT() < 12*GeV && jets.size() >= 2)
        _h_dsigdphi[5]->fill(dPhi2);

  if (firstJet.pT() >  7*GeV && firstJet.pT() < 12*GeV && jets.size() >= 3)
        _h_dsigdphi[6]->fill(dPhi3);

  if (firstJet.pT() >  12*GeV && firstJet.pT() < 30*GeV)
        _h_dsigdphi[7]->fill(dPhi);

  if (firstJet.pT() >  12*GeV && firstJet.pT() < 30*GeV && jets.size() >= 2)
        _h_dsigdphi[8]->fill(dPhi2);

  if (firstJet.pT() >  12*GeV && firstJet.pT() < 30*GeV && jets.size() >= 3)
        _h_dsigdphi[9]->fill(dPhi3);

	// Table 5.1 start here
if (Q2  > 10 && Q2 < 50)
        _h_dsigdphi[10]->fill(dPhi);

if (Q2  > 10 && Q2 < 50 && jets.size() >= 2)
        _h_dsigdphi[11]->fill(dPhi2);

if (Q2  > 10 && Q2 < 50 && jets.size() >= 3)
        _h_dsigdphi[12]->fill(dPhi3);

	// Table 6.1
if (Q2  > 50 && Q2 < 100)
        _h_dsigdphi[13]->fill(dPhi);

if (Q2  > 50 && Q2 < 100 && jets.size() >= 2)
        _h_dsigdphi[14]->fill(dPhi2);

if (Q2  > 50 && Q2 < 100 && jets.size() >= 3)
        _h_dsigdphi[15]->fill(dPhi3);

	//Table 7.1
if (Q2  > 100 && Q2 < 350)
        _h_dsigdphi[16]->fill(dPhi);

if (Q2  > 100 && Q2 < 350 && jets.size() >= 2)
        _h_dsigdphi[17]->fill(dPhi2);

if (Q2  > 100 && Q2 < 350 && jets.size() >= 3)
        _h_dsigdphi[18]->fill(dPhi3);


	}



    /// Normalise histograms etc., after the run
    void finalize() {


	const double sf = crossSection() / picobarn / sumOfWeights();
	for (auto& h : _h_dsigdphi)      scale(h, sf);


    }

    /// @}


	Histo1DPtr _h_dsigdphi[19];

  };


  RIVET_DECLARE_PLUGIN(ZEUS_2024_I2794054);

}
