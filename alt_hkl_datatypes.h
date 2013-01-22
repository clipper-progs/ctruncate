/*! \file lib/hkl_datatypes.h
    Fundamental data types for the clipper libraries
*/
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#ifndef CLIPPER_ALT_HKL_DATATYPES
#define CLIPPER_ALT_HKL_DATATYPES

#include <complex>
#include "ctruncate_utils.h"
#include "clipper/core/hkl_data.h"


namespace clipper
{

  // Now define some actual datatypes


  namespace datatypes
  {

    // modified version of I_sigI_ano without covariance
    //! Reflection data type: I(+) I(+) sigI(+) sigI(-) 
    /*! Note that I_sigI_ano also has methods for returning I(),
      sigI(), so you can use this type in any template type where you
      would use I_sigI. */
    template<class dtype> class J_sigJ_ano : private Datatype_base
    { public:
      J_sigJ_ano() { set_null(); }
      void set_null() { Util::set_null(I_pl_); Util::set_null(I_mi_); Util::set_null(sigI_pl_); Util::set_null(sigI_mi_);  }
      static String type() { return "J_sigJ_ano"; }
      void friedel() { dtype I=I_pl_; I_pl_=I_mi_; I_mi_=I;
                       I=sigI_pl_; sigI_pl_=sigI_mi_; sigI_mi_=I; }
      void shift_phase(const ftype& dphi) {}
      bool missing() const { return (Util::is_nan(I_pl_) && Util::is_nan(I_mi_)); }
      static int data_size() { return 4; }
      static String data_names() { return "I+ sigI+ I- sigI-"; }
      void data_export( xtype a[] ) const { a[0] = I_pl(); a[1] = sigI_pl(); a[2] = I_mi(); a[3] = sigI_mi(); }
      void data_import( const xtype a[] ) { I_pl() = a[0]; sigI_pl() = a[1]; I_mi() = a[2]; sigI_mi() = a[3]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { I_pl_ *= (s*s); sigI_pl_ *= (s*s); I_mi_ *= (s*s); sigI_mi_ *= (s*s); }
      // accessors
      const dtype& I_pl() const { return I_pl_; }  //<! read access
      const dtype& sigI_pl() const { return sigI_pl_; }  //<! read access
      const dtype& I_mi() const { return I_mi_; }  //<! read access
      const dtype& sigI_mi() const { return sigI_mi_; }  //<! read access
		dtype d() const { if ( Util::is_nan(I_pl_) || Util::is_nan(I_mi_) ) return 0.0; else return I_pl_ - I_mi_; } //<! read access
		dtype sigd() const { if ( Util::is_nan(sigI_pl_) ) return sigI_mi_; if ( Util::is_nan(sigI_mi_) ) return sigI_pl_;
			else return std::sqrt(std::pow(sigI_pl_,dtype(2.0))+std::pow(sigI_mi_,dtype(2.0)) ); } //<! read access
      dtype& I_pl() { return I_pl_; }  //<! write access
      dtype& sigI_pl() { return sigI_pl_; }  //<! write access
      dtype& I_mi() { return I_mi_; }  //<! write access
      dtype& sigI_mi() { return sigI_mi_; }  //<! write access
      // nonanomalous-type accessors
      //dtype I() const { return Util::mean(I_pl_,I_mi_); } //use 0.5a+0.5b
      //dtype sigI() const { return Util::sig_mean(sigI_pl_,sigI_mi_,dtype(0.0)); } // use 0.5s1+0.5s2 
      dtype I() const { return ctruncate::Utils::mean(I_pl_,I_mi_,sigI_pl_,sigI_mi_); }  //<! read access as simple sigma weighted
      dtype sigI() const { return ctruncate::Utils::sig_mean(I_pl_,I_mi_,sigI_pl_,sigI_mi_,dtype(0.0)); }  //<! read access as simple
    private:
      dtype I_pl_, I_mi_, sigI_pl_, sigI_mi_ ;
    };

	// modified version of F_sigF_ano without covariance
    //! Reflection data type: F(+) F(+) sigF(+) sigF(-) 
    /*! Note that F_sigF_ano also has methods for returning f(),
      sigf(), so you can use this type in any template type where you
      would use F_sigF. */
    template<class dtype> class G_sigG_ano : private Datatype_base
    { public:
      G_sigG_ano() { set_null(); }
      void set_null() { Util::set_null(f_pl_); Util::set_null(f_mi_); Util::set_null(sigf_pl_); Util::set_null(sigf_mi_); }
      static String type() { return "G_sigG_ano"; }
      void friedel() { dtype f=f_pl_; f_pl_=f_mi_; f_mi_=f;
                       f=sigf_pl_; sigf_pl_=sigf_mi_; sigf_mi_=f; }
      void shift_phase(const ftype& dphi) {}
      bool missing() const { return (Util::is_nan(f_pl_) && Util::is_nan(f_mi_)); }
      static int data_size() { return 4; }
      static String data_names() { return "F+ sigF+ F- sigF-"; }
      void data_export( xtype a[] ) const { a[0] = f_pl(); a[1] = sigf_pl(); a[2] = f_mi(); a[3] = sigf_mi(); }
      void data_import( const xtype a[] ) { f_pl() = a[0]; sigf_pl() = a[1]; f_mi() = a[2]; sigf_mi() = a[3]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { f_pl_ *= s; sigf_pl_ *= s; f_mi_ *= s; sigf_mi_ *= s; }
      // accessors
      const dtype& f_pl() const { return f_pl_; }  //<! read access
      const dtype& sigf_pl() const { return sigf_pl_; }  //<! read access
      const dtype& f_mi() const { return f_mi_; }  //<! read access
      const dtype& sigf_mi() const { return sigf_mi_; }  //<! read access
		dtype d() const { if ( Util::is_nan(f_pl_) || Util::is_nan(f_mi_) ) return 0.0; else return f_pl_ - f_mi_; } //<! read access
		dtype sigd() const { if ( Util::is_nan(sigf_pl_) ) return sigf_mi_; if ( Util::is_nan(sigf_mi_) ) return sigf_pl_;
			else return std::sqrt(std::pow(sigf_pl_,dtype(2.0))+std::pow(sigf_mi_,dtype(2.0)) ); } //<! read access
      dtype& f_pl() { return f_pl_; }  //<! write access
      dtype& sigf_pl() { return sigf_pl_; }  //<! write access
      dtype& f_mi() { return f_mi_; }  //<! write access
      dtype& sigf_mi() { return sigf_mi_; }  //<! write access
      // nonanomalous-type accessors
      //dtype f() const { return Util::mean(f_pl_,f_mi_); }
      //dtype sigf() const { return Util::sig_mean(sigf_pl_,sigf_mi_,dtype(0.0)); } 
      dtype f() const { return ctruncate::Utils::mean(f_pl_,f_mi_,sigf_pl_,sigf_mi_); }  //<! read access as simple
      dtype sigf() const { return ctruncate::Utils::sig_mean(f_pl_,f_mi_,sigf_pl_,sigf_mi_,dtype(0.0)); }  //<! read access as simple
    private:
      dtype f_pl_, f_mi_, sigf_pl_, sigf_mi_ ;
    };


	class ISym : private Datatype_base
	{
	public:
		ISym() { bias_ = -1;}
		explicit ISym( const int& bias ) : bias_(bias) {}
		void set_null() { bias_ = -1; }
		static String type() { return "ISym"; }
		void friedel() {}
		void shift_phase(const ftype& dphi) {}
		bool missing() const { return (bias_ == -1); }
		static int data_size() { return 1; }
		static String data_names() { return "ISym"; }
		void data_export( xtype array[] ) const
		{ array[0] = xtype( bias_ ); }
		void data_import( const xtype array[] )
		{ bias_ = int(array[0]) ; }
		// accessors
		const int& isym() const { return bias_; }  //<! read access
		int& isym() { return bias_; }  //<! write access
	private:
		int bias_; // ISYM=0 normally but for those HKLs where all data
		//    refer to I+ or to I-  ISYM is set to 1 or 2 indicating that
		//    FMEAN value is biased.
	};
		  
}
  namespace data32
  {
    typedef clipper::datatypes::J_sigJ_ano<ftype32> J_sigJ_ano;  //!< datatype
	typedef clipper::datatypes::G_sigG_ano<ftype32> G_sigG_ano;  //!< datatype
	typedef clipper::datatypes::ISym ISym;						 //!< datatype
  }

  namespace data64
  {
    typedef clipper::datatypes::J_sigJ_ano<ftype64> J_sigJ_ano;  //!< datatype
	typedef clipper::datatypes::G_sigG_ano<ftype64> G_sigG_ano;  //!< datatype
	typedef clipper::datatypes::ISym ISym;						 //!< datatype
  }



} // namespace clipper

#endif
