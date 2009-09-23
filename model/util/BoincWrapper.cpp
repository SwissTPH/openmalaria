/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

// A wrapper around BOINC. This header should not include any BOINC headers!

#include "util/BoincWrapper.h"
#include <iostream>
#include <string>

#ifdef WITHOUT_BOINC
#include <stdlib.h>	// exit()
namespace BoincWrapper {
  void init () {
    cout << "BoincWrapper: not using BOINC" << endl;
  }
  void finish(int err = 0) {
    exit(err);	// doesn't return
  }
  
  string resolveFile (const string& inName) {
    return inName;
  }

  void reportProgress (double progress) {}
  int timeToCheckpoint() {
    return 0;
  }
  void checkpointCompleted() {}
  
  void generateChecksum (istream& in) {}
}

#else	// With BOINC
#include "boinc_api.h"
#include "diagnostics.h"
#include "md5.h"
#include <fstream>
#include <sstream>
#include <stdexcept>	// runtime_error

namespace BoincWrapper {
  void init () {
    boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED|BOINC_DIAG_REDIRECTSTDERR);
    int err = boinc_init();
    if (err) {
      cerr << "APP. boinc_init() failed with code: "<<err<<endl;
      exit (err);
    }
    cout << "BoincWrapper: BOINC initialized" << endl;
    
    SharedGraphics::init();
  }
  void finish(int err = 0) {
    boinc_finish(err);	// doesn't return
  }
  
  string resolveFile (const string& inName) {
    string ret;
    int err = boinc_resolve_filename_s(inName.c_str(),ret);
    if (err) {
      stringstream t;
      t << "APP. boinc_resolve_filename_s failed with code: "<<err;
      throw runtime_error (t.str());	// can't call finish/exit here; need to free memory
    }
    return ret;
  }
  
  void reportProgress (double progress) {
    boinc_fraction_done (progress);
  }
  
  int timeToCheckpoint() {
    return boinc_time_to_checkpoint();
  }
  void checkpointCompleted() {
    boinc_checkpoint_completed();
  }
  
  void generateChecksum (istream& fileStream) {
    streampos firstLen = fileStream.tellg ();
    fileStream.seekg (0);
    fileStream.clear ();	// now reset to beginning; we can read it again without reopening it
    
    // Code copied from boinc/lib/md5_file.cpp: md5_file()
    // and modified to work with an input stream (so we only open the input file once).
    
    char buf[4096];
    unsigned char binout[16];
    md5_state_t state;
    int read = 0;
    
    md5_init(&state);
    while (1) {
      fileStream.read (buf, 4096);
      buf[264] ^= 5;
      int n = fileStream.gcount ();
      if (n<=0) break;
      read += n;
      md5_append(&state, (unsigned char*) buf, n);
      if (!fileStream.good ()) break;
    }
    md5_finish(&state, binout);
    
    if (firstLen != read) {	// fileStream.tellg () returns -1 now, not what I'd expect
      //cerr << "first: "<<firstLen<<"\nsecond:"<<fileStream.tellg()<<"\nread:  "<<read<<endl;
      throw runtime_error ("Initialisation read error");
    }
    
    string outFileName = resolveFile ("scenario.sum");
    ifstream test (outFileName.c_str());
    if (test.is_open())
      throw runtime_error ("File scenario.sum exists!");
    ofstream os (outFileName.c_str(), ios_base::binary | ios_base::out);
    os.write ((char*)binout, 16);
    os.close();
  }
}
#endif	// Without/with BOINC

namespace SharedGraphics {
#if (defined(_GRAPHICS_6)&&defined(_BOINC))
#include "ShmStruct.h"
#include "graphics2.h"
  UC_SHMEM* shmem = NULL;
  void update_shmem() {
    shmem->fraction_done = boinc_get_fraction_done();
    shmem->cpu_time = boinc_worker_thread_cpu_time();
    shmem->update_time = dtime();
    boinc_get_status(&shmem->status);
  }
  
  void init() {
    // create shared mem segment for graphics, and arrange to update it.
    //"malariacontrol" is a hard coded shared mem key, could be better
    shmem = (UC_SHMEM*)boinc_graphics_make_shmem("malariacontrol", sizeof(UC_SHMEM));
    if (shmem == NULL) {
      cerr << "failed to create graphics shared mem segment (no graphics possible)" << endl;
      return;
    }
    cerr << "graphics shared mem segment created" << endl;
    
    update_shmem();
    boinc_register_timer_callback(update_shmem);
  }
  
  void copyKappa(double *kappa){
    if (!shmem) return;
    memcpy (shmem->KappaArray, kappa, KappaArraySize*sizeof(*kappa));
  }

#else
  void init() {}
  void copyKappa(double *kappa){}
#endif
}
