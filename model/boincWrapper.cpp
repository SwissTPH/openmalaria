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

#include "boincWrapper.h"
#include "boinc_bridge.h"
#include "ShmStruct.h"

#include <iostream>

namespace BoincWrapper {
  void init () {
    boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED|BOINC_DIAG_REDIRECTSTDERR);
    int err = boinc_init();
    if (err) {
      cerr << "APP. boinc_init() failed with code: "<<err<<endl;
      exit (err);
    }
    cout << "boinc initialized" << endl;
    
    SharedGraphics::init();
  }
  void finish(int err = 0) {
    boinc_finish(err);	// doesn't return
  }
  
  string resolveFile (const char* inName) {
    string ret;
    int err = boinc_resolve_filename_s(inName,ret);
    if (err) {
      cerr << "APP. boinc_resolve_filename failed with code: "<<err<<endl;
      finish(err);
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
}

namespace SharedGraphics {
#if (defined(_GRAPHICS_6)&&defined(_BOINC))
  UC_SHMEM* shmem;
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
    if (!shmem) {
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
