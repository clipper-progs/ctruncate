#ifndef CTRUNCATE_UTILS_H
#define CTRUNCATE_UTILS_H

#include "clipper/clipper.h"
#include "clipper/clipper-contrib.h"
#include "clipper/clipper-ccp4.h"
#include "clipper/clipper-minimol.h"

#include <vector>
#include <cmath>

void tricart(clipper::Cell& cell, clipper::Mat33<clipper::ftype>& transf);
void MatrixToString( clipper::Mat33<int>& twinoper, clipper::String &s );
void straight_line_fit(std::vector<clipper::ftype>& x, std::vector<clipper::ftype>& y, std::vector<clipper::ftype>& w, int n, clipper::ftype &a, clipper::ftype &b, clipper::ftype &siga, clipper::ftype &sigb);


namespace ctruncate {
	
	// Close(a,b[,tol])  true if a == b within tolerance
	template<class T1, class T2> inline static bool
	Close(const T1& a, const T2& b,
		  const T1& tol=1.0e-6)
	{ return std::abs(a-b)<=tol; }
	
	class IceRing
	{
	public:
		// Resolution in A, width in 1/d^2 units
		IceRing(const double& Resolution, const double& width);
		
		// Resolution of ring: return centre of ring as d* = 1/d
		double Dstar() const {return std::sqrt(ring_invressqr);}
		
		// If in ring, returns true
		bool InRing(const double& invresolsq) const;
		// Clear intensity sums
		void ClearSums();
		// Add in IsigI, use clipper version for now rather than pointless version
		void AddObs(const clipper::datatypes::I_sigI<float>& I_sigI, const double& invresolsq);
		
		void SetReject() {reject=true;}
		void SetReject(const bool& Rej) {reject=Rej;}
		
		// Results
		double MeanI() const;
		double MeanSigI() const;
		double MeanSSqr() const;
		int N() const {return nI;}
		
		// Reject flag, true means reject reflection in this range
		bool Reject() const {return reject;}
		
	private:
		double ring_invressqr;  // centre of ring in 1/d^2
		double halfwidth_invressqr;  // halfwidth of ring in 1/d^2
		double sum_I;
		double sum_sigI;
		double sum_sSqr;
		int nI;
		bool reject;
	};
	//--------------------------------------------------------------
	class Rings
	{
	public:
		Rings() : nrings(0) {}
		
		// Resolution in A, width in 1/d^2 units
		void AddRing(const double& Resolution, const double& width);
		
		// Set up default ice rings: 3.90, 3.67, 3.44A
		void DefaultIceRings();
		
		// Clear list
		void Clear();
		
		int Nrings() const {return nrings;}
		
		// Resolution of iring'th ring as d*
		double Dstar(const int& iring) const
		{return rings.at(iring).Dstar();}
		
		// Copy rejected rings only
		void CopyRejRings(const Rings& other);
		
		// If in ring, returns ring number (0,n-1), else = -1 
		int InRing(const double& invresolsq) const;
		// Clear intensity sums
		void ClearSums();
		// Add in IsigI
		void AddObs(const int& Iring, const clipper::datatypes::I_sigI<float>& I_sigI, const double& invresolsq);
		
		void SetReject(const int& Iring);
		void SetReject(const int& Iring, const bool& Rej);
		
		// Results
		double MeanI(const int& Iring) const;
		double MeanSigI(const int& Iring) const;
		double MeanSSqr(const int& Iring) const;
		// Reject flag, true means reject reflection in this range
		bool Reject(const int& Iring) const; 
		int N(const int& Iring) const;
		
	private:
		int nrings;
		std::vector<IceRing> rings;
		
		void CheckRing(const int& Iring) const;
	};
	
	
}
#endif
