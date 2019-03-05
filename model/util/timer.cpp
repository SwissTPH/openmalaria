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

// Checkpoint timers

#ifdef _WIN32
#include <windows.h>
#else
#include <cstdio>	// perror
#include <pthread.h>
#include <unistd.h>     // sleep
#endif

#include "util/timer.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>

using namespace std;

namespace OM { namespace util {

bool finishedCP=false;

//NOTE: is calling exit from another thread even valid?

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
  cerr<<"cpw_to"<<endl;
  exit(-6);
  //return 0;
}

void timer::startCheckpoint (){
  finishedCP=false;
  //setup windows thread
  timer_threadCP = CreateThread(nullptr, 0, write_cp_timer, nullptr,0, &timer_threadId);
}

void timer::stopCheckpoint (){
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
  cerr<<"cpw_to"<<endl;
  exit(-6);
  //return 0;
}

void timer::startCheckpoint (){
  finishedCP=false;
  //setup pthread
  res = pthread_create(&timer_thread, nullptr, write_cp_timer, nullptr);
}

void timer::stopCheckpoint (){
  finishedCP=true;
  res = pthread_join(timer_thread, &thread_result);
  if (res != 0) {
    perror("Thread join failed");
    exit(EXIT_FAILURE);
  }
}

#endif
} }
