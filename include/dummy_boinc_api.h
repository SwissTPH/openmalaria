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

#ifndef _DUMMY_BOINC_API_
#define _DUMMY_BOINC_API_
#include <string>

typedef void (*FUNC_PTR)();

struct APP_INIT_DATA;
struct BOINC_OPTIONS;
struct BOINC_STATUS;

int boinc_resolve_filename_s(const char*, std::string&) {return 0;};
int boinc_get_init_data(APP_INIT_DATA&) {return 0;};
int boinc_wu_cpu_time(double&) {return 0;};
int boinc_upload_file(std::string& name) {return 0;};
int boinc_upload_status(std::string& name) {return 0;};
int boinc_write_init_data_file(APP_INIT_DATA&) {return 0;};
int suspend_activities() {return 0;};   
int resume_activities() {return 0;};    
int restore_activities() {return 0;};  

int boinc_init(void) {return 0;};
int boinc_finish(int status) {return 0;};
int boinc_resolve_filename(const char*, char*, int len) {return 0;};
int boinc_get_init_data_p(struct APP_INIT_DATA*) {return 0;};
int boinc_parse_init_data_file(void) {return 0;};
int boinc_send_trickle_up(char* variety, char* text) {return 0;};
int boinc_checkpoint_completed(void) {return 0;};
int boinc_fraction_done(double) {return 0;};
int boinc_suspend_other_activities(void) {return 0;};
int boinc_resume_other_activities(void) {return 0;};
int boinc_report_app_status(double, double, double) {return 0;};
int boinc_time_to_checkpoint() {return 0;};
void boinc_begin_critical_section() {};
int boinc_try_critical_section() {return 0;};
void boinc_end_critical_section() {};
void boinc_need_network() {};
int boinc_network_poll() {return 0;};
void boinc_network_done() {};
int boinc_is_standalone(void) {return 0;};
void boinc_ops_per_cpu_sec(double fp, double integer) {};
void boinc_ops_cumulative(double fp, double integer) {};
int boinc_receive_trickle_down(char* buf, int len) {return 0;};
int boinc_init_options(BOINC_OPTIONS*) {return 0;};
int boinc_get_status(BOINC_STATUS*) {return 0;};
double boinc_get_fraction_done() {return 0;};
void boinc_register_timer_callback(FUNC_PTR) {};
double boinc_worker_thread_cpu_time() {return 0;};
void boinc_exit(int) {};    // deprecated

int setMacPList(void) {return 0;};
int setMacIcon(char *filename, char *iconData, long iconSize) {return 0;};

#endif
