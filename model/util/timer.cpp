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

// Checkpoint timers

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#include "util/timer.h"
#include "boinc_bridge.h"
#include <iostream>
using namespace std;

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

void timer::startCheckpoint (){
  finishedCP=false;
  //setup windows thread
  timer_threadCP = CreateThread(NULL, 0, write_cp_timer, NULL,0, &timer_threadId);
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
#include <stdio.h>	// perror

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

void timer::startCheckpoint (){
  finishedCP=false;
  //setup pthread
  res = pthread_create(&timer_thread, NULL, write_cp_timer, NULL);
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
