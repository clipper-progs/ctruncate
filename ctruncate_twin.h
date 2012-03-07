//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_TWIN
#define __CTRUNCATE_TWIN

#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"

namespace ctruncate {
	
    clipper::ftype Htest_driver_fp(clipper::HKL_data<clipper::data32::I_sigI>& isig, bool debug);
    clipper::ftype Htest_driver_table(clipper::HKL_data<clipper::data32::I_sigI>& isig, bool debug);
	
    clipper::ftype Htest( clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::Mat33<int>& twinop, int scalefac, 
			   clipper::String s, bool debug );
	
    clipper::ftype Ltest_driver(clipper::HKL_data<clipper::data32::I_sigI>& isig, bool debug);
    
    void twin_summary(clipper::ftype, clipper::ftype);
	
}

#endif

