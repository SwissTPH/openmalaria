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

/* This header either includes BOINC headers or defines dummy functions.
 * It's no longer used by the malaria simulation, but is still used by the
 * graphics application. */

#ifndef _BOINC_BRIDGE_
#define _BOINC_BRIDGE_


#ifndef WITHOUT_BOINC
#include "boinc_api.h"
#include "diagnostics.h"
#if (defined(_GRAPHICS_6))
#include "util.h"
#include "graphics2.h"
#endif
#else	// #ifndef WITHOUT_BOINC

#include <string>
#include <string.h>
#include <stdlib.h>

// flags for boinc_init_diagnostics()
//
#define BOINC_DIAG_DUMPCALLSTACKENABLED     0x00000001L
#define BOINC_DIAG_HEAPCHECKENABLED         0x00000002L
#define BOINC_DIAG_MEMORYLEAKCHECKENABLED   0x00000004L
#define BOINC_DIAG_ARCHIVESTDERR            0x00000008L
#define BOINC_DIAG_ARCHIVESTDOUT            0x00000010L
#define BOINC_DIAG_REDIRECTSTDERR           0x00000020L
#define BOINC_DIAG_REDIRECTSTDOUT           0x00000040L
#define BOINC_DIAG_REDIRECTSTDERROVERWRITE  0x00000080L
#define BOINC_DIAG_REDIRECTSTDOUTOVERWRITE  0x00000100L
#define BOINC_DIAG_TRACETOSTDERR            0x00000200L
#define BOINC_DIAG_TRACETOSTDOUT            0x00000400L
#define BOINC_DIAG_HEAPCHECKEVERYALLOC      0x00000800L
#define BOINC_DIAG_BOINCAPPLICATION         0x00001000L
#define BOINC_DIAG_DEFAULTS \
    BOINC_DIAG_DUMPCALLSTACKENABLED | \
    BOINC_DIAG_HEAPCHECKENABLED | \
    BOINC_DIAG_MEMORYLEAKCHECKENABLED | \
    BOINC_DIAG_REDIRECTSTDERR | \
    BOINC_DIAG_TRACETOSTDERR

typedef void (*FUNC_PTR)();

struct APP_INIT_DATA;
struct BOINC_OPTIONS;
struct BOINC_STATUS {
  int dummy;
};

inline int boinc_resolve_filename_s(const char *virtual_name, std::string& physical_name) {
    physical_name = virtual_name;
    return 0;
}
inline int boinc_get_init_data(APP_INIT_DATA&) {return 0;};
inline int boinc_wu_cpu_time(double&) {return 0;};
inline int boinc_upload_file(std::string& name) {return 0;};
inline int boinc_upload_status(std::string& name) {return 0;};
inline int boinc_write_init_data_file(APP_INIT_DATA&) {return 0;};
inline int suspend_activities() {return 0;};   
inline int resume_activities() {return 0;};    
inline int restore_activities() {return 0;};  

inline int boinc_init(void) {return 0;};
inline int boinc_finish(int status) {
  exit (status);
}
inline int boinc_resolve_filename(const char* filenamein, char* filenameout, int len) {
    strcpy(filenameout, filenamein); return 0;
}
inline int boinc_get_init_data_p(struct APP_INIT_DATA*) {return 0;};;
inline int boinc_parse_init_data_file(void) {return 0;};
inline int boinc_send_trickle_up(char* variety, char* text) {return 0;};
inline int boinc_checkpoint_completed(void) {return 0;};
inline int boinc_fraction_done(double) {return 0;};
inline int boinc_suspend_other_activities(void) {return 0;};
inline int boinc_resume_other_activities(void) {return 0;};
inline int boinc_report_app_status(double, double, double) {return 0;};
inline int boinc_time_to_checkpoint() {return 0;};
inline void boinc_begin_critical_section() {};
inline int boinc_try_critical_section() {return 0;};
inline void boinc_end_critical_section() {};
inline void boinc_need_network() {};
inline int boinc_network_poll() {return 0;};
inline void boinc_network_done() {};
inline int boinc_is_standalone(void) {return 0;};
inline void boinc_ops_per_cpu_sec(double fp, double integer) {};
inline void boinc_ops_cumulative(double fp, double integer) {};
inline int boinc_receive_trickle_down(char* buf, int len) {return 0;};
inline int boinc_init_options(BOINC_OPTIONS*) {return 0;};
inline int boinc_get_status(BOINC_STATUS*) {return 0;};
inline double boinc_get_fraction_done() {return 0;};
inline void boinc_register_timer_callback(FUNC_PTR) {};
inline double boinc_worker_thread_cpu_time() {return 0;};
inline void boinc_exit(int) {};    // deprecated

inline int setMacPList(void) {return 0;};
inline int setMacIcon(char *filename, char *iconData, long iconSize) {return 0;};

inline void app_graphics_render(int xs, int ys, double time_of_day) {};
inline void app_graphics_init(void) {};
inline void app_graphics_resize(int width, int height) {};
inline void boinc_app_mouse_button(int x, int y, int which, int is_down) {};
inline void boinc_app_mouse_move(int x, int y, int left, int middle, int right) {};
inline void boinc_app_key_press(int, int) {};
inline void boinc_app_key_release(int, int) {};

// Functions that the app can call
//
inline void boinc_graphics_loop(int, char**) {};
inline void* boinc_graphics_make_shmem(char*, int) {return 0;};
inline void* boinc_graphics_get_shmem(char*) {return 0;};
inline void boinc_set_windows_icon(const char* icon16,const char* icon48) {};
inline void boinc_close_window_and_quit(const char*) {};

inline int boinc_init_diagnostics( int flags ) {return 0;};
inline int boinc_init_graphics_diagnostics( int flags ) {return 0;};
inline int boinc_install_signal_handlers() {return 0;};
inline int boinc_finish_diag() {return 0;};
inline double dtime() {return 0;};

#endif	// #ifndef WITHOUT_BOINC #else

#endif	// #ifndef _BOINC_BRIDGE_
