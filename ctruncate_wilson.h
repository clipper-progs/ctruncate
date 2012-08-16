//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_WILSON_H
#define __CTRUNCATE_WILSON_H

#include "clipper/clipper.h"
#include "clipper/clipper-minimol.h"
#include "clipper/clipper-ccp4.h"
#include "ctruncate_utils.h"

namespace ctruncate {
	
	class Scattering {
	public:
		enum MODE { PROTEIN, NUCLEIC, COMPLEX };
		Scattering( MODE mode=PROTEIN) : _mode(mode) { init(); }
		~Scattering() {}
		clipper::ftype operator()(const int nres, clipper::ftype power=2.0 );
		clipper::ftype operator()(const clipper::MPolymerSequence& poly, clipper::ftype power=2.0);
		clipper::ftype operator()(const clipper::Cell& cell, const clipper::Spacegroup& spg, clipper::ftype solvent=0.5, clipper::ftype power=2.0 );
		clipper::ftype operator()(clipper::MMDBManager& mmdb, const clipper::Spacegroup& spg, clipper::ftype power=2.0);
		clipper::ftype operator()(std::vector<clipper::Atom>& mmdb, const clipper::Spacegroup& spg, clipper::ftype power=2.0 );
		
		static const clipper::ftype C, O, N, H, S, P;
		
	private:
		void init();
		static const std::string ProteinAtomNames[5];  //atom names
		static const char ProteinResidueNames[20];  //residue names
		static const int ProteinCatoms[20];         // number of C per residue
		static const int ProteinHatoms[20];         // number of H per residue
		static const int ProteinNatoms[20];         // number of N per residue
		static const int ProteinOatoms[20];         // number of O per residue
		static const int ProteinSatoms[20];         // number of S per residue
		static const clipper::ftype ProteinComp[5];          // average composition
		
		static const std::string NucleicAtomNames[5];
		static const char NucleicResidueNames[5];
		static const int NucleicCatoms[5];
		static const int NucleicHatoms[5];
		static const int NucleicNatoms[5];
		static const int NucleicOatoms[5];
		static const int NucleicPatoms[5];
		static const clipper::ftype NucleicComp[5]; 
		
        MODE _mode;                                   // mode of operation
		std::vector<std::string> atomtype;
		std::vector<clipper::ftype> scattering;
		
	};
	
	class WilsonB {
	public:
		enum MODE {STRAIGHT, BEST };
		
		WilsonB ( WilsonB::MODE _mode = BEST ) : mode(_mode), nresidue_supplied(false), sequence_supplied(false), _a(-1.0), 
		_bscale(1.0f), _boff(0.0f)
		{ numatoms.resize(5); }
		
		float operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings* ice = NULL);
		
		float operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, int nresidues, ctruncate::Rings* ice = NULL);
		
		float operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::MPolymerSequence& poly, ctruncate::Rings* ice = NULL);
		
		float a() {return _a;}
		float b() {return _b;}
		float B() {return 2.0*_a; }
		float intercept () { return exp(_b); }
		float siga() {return _siga;}
		float sibb() {return _sigb;}
		
		void summary();
		void plot ();
		void plot(clipper::HKL_data<clipper::data32::I_sigI>& ref, clipper::String& name);
		
	private:			
		MODE mode;
        clipper::ftype _a;  // slope
        clipper::ftype _b;  // intercept
        clipper::ftype _siga; // uncertainty in intercept
        clipper::ftype _sigb; //uncertainty in b-value
        clipper::ftype _bscale; // scale best to scattering
        clipper::ftype _boff; //offset on best to scattering
		clipper::HKL_data<clipper::data32::I_sigI>* intensity;
		
		bool nresidue_supplied, sequence_supplied; //book keeping
		int nresidues; // number of residues
		std::vector<int> numatoms;  // number of atoms
		int _totalscat;  //number of scattering electrons in cell (scale for BEST)
		float maxres;
		
		static const std::string AtomNames[5];  //atom names
		static const char ResidueNames[21];  //residue names
		static const int Catoms[21];         // number of C per residue
		static const int Hatoms[21];         // number of H per residue
		static const int Natoms[21];         // number of N per residue
		static const int Oatoms[21];         // number of O per residue
		static const int Satoms[21];         // number of S per residue
		
		void wilson_straight(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings& ice);
		void wilson_best(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings& ice);
	};
		
}
#endif
