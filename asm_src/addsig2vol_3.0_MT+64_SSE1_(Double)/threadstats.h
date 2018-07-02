#ifndef THREADSTATS_H
#define THREADSTATS_H

#define THREADSTATS_ON


//------------------------------
//Switches
#undef  THREADSTATS_ON              //Exclude code from compiling
//------------------------------

#include <pthread.h>
#include "as2v_array.h"
#include "timestats.h"

#ifdef BUILDMEX
    #include "mex.h"
#endif

void threadstats_markStartTask(unsigned int thread, unsigned int currentTask, unsigned int run);
void threadstats_markEndTask(unsigned int thread, unsigned int currentTask, unsigned int run);


void threadstats_init(unsigned int _nThreads, unsigned int _nTasks, unsigned int _nReruns);
void threadstats_free();

void threadstats_init_mexed(unsigned int _nThreads, unsigned int _nTasks, unsigned int _nReruns, double* _dataTimeStartTask, double* _dataTimeEndTask, double* _dataMoveToTask, double* _dataThreadnumber);
void threadstats_free_mexed();

#endif // THREADSTATS_H
