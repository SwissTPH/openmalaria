/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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


#include <iostream>
using namespace std;
#include <fstream>

//#include "graphics_api.h"
//#include "diagnostics.h"

//Include the graphicals libraries


#ifdef _WIN32
#ifndef _NO_GRAPHICS
#include <GL/glut.h>
#include "GraphicsBridge.h"
#include <windows.h>
#include <gl/gl.h>           // Header File For The OpenGL32 Library
#include <gl/glu.h>          // Header File For The GLu32 Library
#include <gl/glaux.h>        // Header File For The Glaux Library
#include <boinc_win.h>
#endif
#else
#define _GRAPHICS_6
#define _NO_GRAPHICS
#include <pthread.h>
#endif

#include "boinc_bridge.h"

#if (defined(_GRAPHICS_6))
#include "ShmStruct.h"
UC_SHMEM* shmem;
#endif

#include "inputData.h" //Include parser for the input
#include "GSLWrapper.h" //Include wrapper for GSL library

#ifdef _WIN32
//extern struct APP_INIT_DATA dataBOINC;
struct APP_INIT_DATA dataBOINC;
#endif

#ifdef __cplusplus
#include "simulation.h"
#endif

#ifndef __cplusplus
extern "C" {
  void simulate();
  void add_kappa(double*);
  void start_cp_timer();
  void stop_cp_timer();
}
#endif

#if (defined(_GRAPHICS_6))
void update_shmem() {
	if (!shmem) return;
	shmem->fraction_done = boinc_get_fraction_done();
	shmem->cpu_time = boinc_worker_thread_cpu_time();
	shmem->update_time = dtime();
	boinc_get_status(&shmem->status);
}

void add_kappa(double *kappa){
	if (!shmem) return;
	memcpy (shmem->KappaArray, kappa, KappaArraySize*sizeof(*kappa));
}

#else
void add_kappa(double *kappa){
#ifdef _WIN32
	//add_and_copy_data(kappa);
#endif
}
#endif

bool finishedCP=false;

#ifdef _WIN32

	HANDLE timer_threadCP;
    DWORD timer_threadId;
    DWORD thread_result;

DWORD WINAPI write_cp_timer(PVOID arg) {
	int counter=0;
	while (counter<10800){
	Sleep(1000);
	if(finishedCP){
		return 0;
	}
	counter++;
	}
	//Checkpoint write timed out
	cerr<<"cpw_to\n";
	exit(-6);
	//return 0;
}

void start_cp_timer(){
	finishedCP=false;
	//setup windows thread
    timer_threadCP = CreateThread(NULL, 0, write_cp_timer, NULL,0, &timer_threadId);
}

void stop_cp_timer(){
	finishedCP=true;
    if (WaitForSingleObject(timer_threadCP, INFINITE) != WAIT_OBJECT_0) {
        perror("Thread join failed");
        exit(EXIT_FAILURE);
    }
    // Retrieve the code returned by the thread.
    GetExitCodeThread(timer_threadCP, &thread_result);
}

#else
//Pthread version
int res;
pthread_t timer_thread;
void *thread_result;

void *write_cp_timer(void *arg) {
    int counter=0;
	while (counter<10800){
		sleep(1);
		if(finishedCP){
			return 0;
		}
		counter++;
	}
	//Checkpoint write timed out
	cerr<<"cpw_to\n";
	exit(-6);
	//return 0;
}

void start_cp_timer(){
	finishedCP=false;
	//setup pthread
	res = pthread_create(&timer_thread, NULL, write_cp_timer, NULL);
}

void stop_cp_timer(){
	finishedCP=true;
    res = pthread_join(timer_thread, &thread_result);
    if (res != 0) {
        perror("Thread join failed");
        exit(EXIT_FAILURE);
    }
}
#endif


int main2(int argc, char* argv[]){
  try {
    //Function to test output given by boinc
    int retval =0;	
    
    const int MAX_LENGTH = 1000;
    char scenario [MAX_LENGTH];

//Call the initialisation of boinc
#ifdef _WIN32	
    boinc_get_init_data(dataBOINC);
#endif
    
    string scenario_name = "scenario.xml";
    
    if (argc == 2){
      scenario_name = argv[1];

      fstream scenario_file(scenario_name.c_str(),fstream::in);
      if (scenario_file.is_open()){
        scenario_file.close();
      }
      else {
        cout << "Error: " << scenario_name << " file does not exist" << endl;
        exit(-1);
      }
    }
  
    //Resolve the scenario filename
    retval = boinc_resolve_filename(scenario_name.c_str(),scenario,MAX_LENGTH);
    if (retval){
      cerr << "APP. boinc_resolve_filename failed \n";
      boinc_finish(retval);
    }

	//Change it and read it with boinc
	bool succeeded;
	succeeded = createDocument(scenario);
	if (!succeeded){
		cerr << "APP. createDocument failed \n";
		boinc_finish(-1);
		exit(-1);
	}

#if (defined(_GRAPHICS_6)&&defined(_BOINC))
	// create shared mem segment for graphics, and arrange to update it.
	//"malariacontrol" is a hard coded shared mem key, could be better
	shmem = (UC_SHMEM*)boinc_graphics_make_shmem("malariacontrol", sizeof(UC_SHMEM));
	if (!shmem) {
		fprintf(stderr, "failed to create graphics shared mem segment\n");
		}
	else{
		fprintf(stderr, "graphics shared mem segment created\n");
	}
	update_shmem();
	boinc_register_timer_callback(update_shmem);
#endif

	GSL_SETUP();	
  Simulation* simulation = new Simulation();
	simulation->start();
  delete simulation;
	GSL_TEARDOWN();
	
  
	//Now we can exit with boinc. The function does not return

	cleanDocument();
	retval = boinc_finish(0);	
	if (retval){
		cerr << "APP. boinc_finish() failed \n";
		exit(retval);
	}
  } catch (...) {
    cerr << "Error occurred." << endl;
    return -1;
  }

  return 0;
}
#if (defined(_WIN32) && defined(_NO_GRAPHICS)) || !defined(_WIN32) || defined(_GRAPHICS_6)
/** main() - initializes BOINC and calls main2. */
int main(int argc, char* argv[]){
	int retval;
        retval = boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED|BOINC_DIAG_REDIRECTSTDERR);
	//Call the initialisation of boinc
	retval = boinc_init();
	if (retval){
		cerr << "APP. boinc_init() failed \n";
		exit(retval);
	}		
	cout << "boinc initialized\n"; 
	main2(0,NULL);
}
#endif

#if (!(defined(_GRAPHICS_6)))
void worker(){
	main2(0, (char **) NULL);
	return;
	}
#endif

//Graphic stuff
extern "C" {
	void graphic_progress_(int* num,int* numtotal){	}
}

