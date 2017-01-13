/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// A wrapper around BOINC. This header should not include any BOINC headers!

#include "util/BoincWrapper.h"
#include <iostream>
#include <string>
#include <fstream>
#include "util/errors.h"

#ifdef WITHOUT_BOINC
#include <stdlib.h>	// exit()
#include <cmath>		// ceil
#include <boost/format.hpp>
#else	// With BOINC
#ifdef _WIN32
#include "boinc_win.h"
#endif
#include "boinc_api.h"
#include "diagnostics.h"
#include "md5.h"
#include "filesys.h"    // lib/filesys.h if not installed
#include <sstream>
#endif

using namespace std;

namespace OM { namespace util {

#ifdef WITHOUT_BOINC
namespace BoincWrapper {
  void init () {
    cout << "BoincWrapper: not using BOINC" << endl;
  }
  void finish(int err) {
      cout << '\r' << flush;	// clean last line of progress-output
    exit(err);	// doesn't return
  }
  
  string resolveFile (const string& inName) {
    return inName;
  }
  
  bool fileExists (const char* path){
    ifstream test (path);
    return test.is_open();
  }
  
  int lastPercent = -1;	// last _integer_ percentage value
  void reportProgress (int now, int duration) {
      int percent = (now * 100) / duration;
      if( percent != lastPercent ){	// avoid huge amounts of output for performance/log-file size reasons
	  lastPercent = percent;
	// \r cleans line. Then we print progress as a percentage.
	cout << (boost::format("\r[%|3i|%%]\t") %percent) << flush;
      }
  }
  int timeToCheckpoint() {
    return 0;
  }
  void checkpointCompleted() {}
  
  void beginCriticalSection() {}
  void endCriticalSection() {}
}
Checksum Checksum::generate (istream& fileStream) {
    // Return a dummy checksum; making sure it is always the same.
    // Note: checkpoints from BOINC and non-BOINC versions are thus incompatible.
    Checksum ret;
    for(int i = 0; i < 16; ++i)
	ret.data[i] = 0;
    return ret;
}
// In non-BOINC mode we don't need checksums, so don't write one:
void Checksum::writeToFile (string filename) {}

#else	// With BOINC

namespace BoincWrapper {
  void init () {
    boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED|BOINC_DIAG_REDIRECTSTDERR);
    int err = boinc_init();
    if (err) {
      std::cerr << "APP. boinc_init() failed with code: "<<err<<std::endl;
      exit (err);
    }
    // Probably we don't need to know this. Disable to compress stderr.txt:
//     std::cout << "BoincWrapper: BOINC initialized" << std::endl;
  }
  void finish(int err) {
    boinc_finish(err);	// doesn't return
  }
  
  std::string resolveFile (const std::string& inName) {
    std::string ret;
    int err = boinc_resolve_filename_s(inName.c_str(),ret);
    if (err) {
      stringstream t;
      t << "APP. boinc_resolve_filename_s failed with code: "<<err;
      throw TRACED_EXCEPTION(t.str(),util::Error::FileIO);	// can't call finish/exit here; need to free memory
    }
    return ret;
  }
  
  bool fileExists (const char* path){
    return boinc_file_exists(path);
  }
  
    void reportProgress (int now, int duration) {
        double progress = static_cast<double>(sim::now().raw()) /
            static_cast<double>(totalSimDuration.raw());
        boinc_fraction_done (progress);
    }
  
  int timeToCheckpoint() {
    return boinc_time_to_checkpoint();
  }
  void checkpointCompleted() {
    boinc_checkpoint_completed();
  }
  
  void beginCriticalSection() {
      boinc_begin_critical_section();
  }
  void endCriticalSection() {
      boinc_end_critical_section();
  }
}
Checksum Checksum::generate (istream& fileStream) {
    fileStream.clear ();        // clear EOF flag (must happen before tellg/seekg)
    streampos firstLen = fileStream.tellg ();
    fileStream.seekg (0);
    // now reset to beginning; we can read it again without reopening it
    
    // Code copied from boinc/lib/md5_file.cpp: md5_file()
    // and modified to work with an input stream (so we only open the input file once).
    
    char buf[4096];
    Checksum output;
    md5_state_t state;
    streampos bytes_read = 0;
    
    md5_init(&state);
    while (1) {
	fileStream.read (buf, 4096);
	
#       ifndef NO_CHECKSUM_PERTURB
#       error "For BOINC builds: insert checksum perturbation (compile with -DNO_CHECKSUM_PERTURB to ignore)."
#       endif
	
	streamsize n = fileStream.gcount ();
	if (n<=0) break;
	bytes_read += n;
	md5_append(&state, (unsigned char*) buf, n);
	if (!fileStream.good ()) break;
    }
    md5_finish(&state, output.data);
    
    if (firstLen != bytes_read) {	// fileStream.tellg () returns -1 now, not what I'd expect
	ostringstream msg;
	msg << "Initialisation read error:\tfirst: "<<firstLen<<"\tsecond:"<<fileStream.tellg()<<"\tread:  "<<bytes_read;
	throw TRACED_EXCEPTION(msg.str(),Error::Checksum);
    }
    
    return output;
}
void Checksum::writeToFile (string filename) {
    ifstream test (filename.c_str());
    if (test.is_open())
	throw util::base_exception("File scenario.sum exists!",Error::Checksum);
    
    // Use C file commands, since these have clearer behaviour with binary data:
    FILE *f = fopen( filename.c_str(), "wb" );
    size_t written=0;
    if( f != NULL )
	written=fwrite( data, 1, 16, f );
    fclose( f );
    if( written != 16 )
	throw util::base_exception("Error writing scenario.sum",Error::Checksum);
}
#endif	// Without/with BOINC

} }
