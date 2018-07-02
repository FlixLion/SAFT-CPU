#include "threadstats.h"

static pthread_mutex_t mutexStartTask;
static pthread_mutex_t mutexEndTask;

static unsigned int nThreads = 0;
static unsigned int nTasks = 0;
static unsigned int nReruns = 0;

double* dataTimeStartTask = NULL;
double* dataTimeEndTask = NULL;
double* dataMoveToTask = NULL;
double* dataThreadnumber = NULL;

cArrayDouble arrayTimeStartTask;
cArrayDouble arrayTimeEndTask;
cArrayDouble arrayMoveToTask;
cArrayDouble arrayThreadnumber;

static unsigned int* timeStartTracker = NULL;
static unsigned int* timeEndTracker = NULL;

static char isInitalized = 0;
static char isMexed = 0;

void threadstats_init(unsigned int _nThreads, unsigned int _nTasks, unsigned int _nReruns){
    #ifdef THREADSTATS_ON
        if(!isInitalized){
            nThreads = _nThreads;
            nTasks = _nTasks;
            nReruns = _nReruns;

            dataTimeStartTask = (double*) malloc(nTasks*nReruns*sizeof(double));
            dataTimeEndTask = (double*) malloc(nTasks*nReruns*sizeof(double));
            dataMoveToTask = (double*) malloc(nTasks*nReruns*sizeof(double));
            dataThreadnumber = (double*) malloc(nTasks*nReruns*sizeof(double));

            timeStartTracker = (unsigned int*) calloc(nReruns, sizeof(unsigned int));
            timeEndTracker = (unsigned int*) calloc(nReruns, sizeof(unsigned int));
            isInitalized = 0xff;
            isMexed = 0;
            //printf("threadstats| nThreads %i, nTasks %i, nReruns %i, elements %i\n", nThreads, nTasks, nReruns, nTasks*nReruns);
        }
    #endif
}

void threadstats_init_mexed(unsigned int _nThreads, unsigned int _nTasks, unsigned int _nReruns, double* _dataTimeStartTask, double* _dataTimeEndTask, double* _dataMoveToTask, double* _dataThreadnumber){
    #ifdef THREADSTATS_ON
        if(!isInitalized){
            nThreads = _nThreads;
            nTasks = _nTasks;
            nReruns = _nReruns;

            dataTimeStartTask = _dataTimeStartTask;
            dataTimeEndTask = _dataTimeEndTask;
            dataMoveToTask = _dataMoveToTask;
            dataThreadnumber = _dataThreadnumber;

            timeStartTracker = (unsigned int*) calloc(nReruns, sizeof(unsigned int));
            timeEndTracker = (unsigned int*) calloc(nReruns, sizeof(unsigned int));
            isInitalized = 0xff;
            isMexed = 0xff;
            //printf("threadstats| nThreads %i, nTasks %i, nReruns %i, elements %i\n", nThreads, nTasks, nReruns, nTasks*nReruns);
        }
    #endif
}


void threadstats_free(){
    #ifdef THREADSTATS_ON
        if(!isMexed && isInitalized){
            free(dataTimeStartTask);
            free(dataTimeEndTask);
            free(dataMoveToTask);
            free(dataThreadnumber);
            dataTimeStartTask = NULL;
            dataTimeEndTask = NULL;
            dataMoveToTask = NULL;
            dataThreadnumber = NULL;

            free(timeStartTracker);
            free(timeEndTracker);
            timeStartTracker = NULL;
            timeEndTracker = NULL;

            isInitalized = 0;
        }
    #endif
}

void threadstats_free_mexed(){
    #ifdef THREADSTATS_ON
        if(isMexed && isInitalized){
            free(timeStartTracker);
            free(timeEndTracker);
            timeStartTracker = NULL;
            timeEndTracker = NULL;
            #ifdef BUILDMEX
                mxFree(dataTimeStartTask);
                mxFree(dataTimeEndTask);
                mxFree(dataMoveToTask);
                mxFree(dataThreadnumber);
                dataTimeStartTask = NULL;
                dataTimeEndTask = NULL;
                dataMoveToTask = NULL;
                dataThreadnumber = NULL;
            #endif
            isMexed = 0;
            isInitalized = 0;
        }
    #endif
}

void threadstats_markStartTask(unsigned int thread, unsigned int currentTask, unsigned int run){
    #ifdef THREADSTATS_ON
        if((run < nReruns) && (currentTask < nTasks)){
            pthread_mutex_lock(&mutexStartTask);
            unsigned int pointer = run*nTasks+timeStartTracker[run];

            dataTimeStartTask[pointer] = getSumTime(thread, TS_MILI);
            dataMoveToTask[pointer] = (double) currentTask;
            dataThreadnumber[pointer] = (double) thread;
            //printf("threadstats_markStartTask| currentTask %i, thread %i, run %i => pointer %i\n", currentTask, thread, run, pointer);
            //printf("threadstats_markStartTask| currentTask %i, timestampStartTask %f, currentTask %f, thread %f\n", currentTask, *(dataTimeStartTask+pointer), *(dataMoveToTask+pointer), *(dataThreadnumber+pointer));
            timeStartTracker[run]++;
            pthread_mutex_unlock(&mutexStartTask);
        }
    #endif
}

void threadstats_markEndTask(unsigned int thread, unsigned int currentTask, unsigned int run){
    #ifdef THREADSTATS_ON
        if((run < nReruns) && (currentTask < nTasks)){
            pthread_mutex_lock(&mutexEndTask);
            unsigned int pointer = run*nTasks+timeEndTracker[run];

            *(dataTimeEndTask+pointer) = getSumTime(thread, TS_MILI);
            //printf("threadstats_markEndTask| currentTask %i, thread %i, run %i => pointer %i\n", currentTask, thread, run, pointer);
            timeEndTracker[run]++;
            pthread_mutex_unlock(&mutexEndTask);
        }
    #endif
}

cArrayDouble* threadstats_get_timeStartTask(){
    #ifdef THREADSTATS_ON
        return &arrayTimeStartTask;
    #endif
}

cArrayDouble* threadstats_get_timeEndTask(){
    #ifdef THREADSTATS_ON
        return &arrayTimeEndTask;
    #endif
}

cArrayDouble* threadstats_get_moveToTask(){
    #ifdef THREADSTATS_ON
        return &arrayMoveToTask;
    #endif
}

cArrayDouble* threadstats_get_threadnumber(){
    #ifdef THREADSTATS_ON
        return &arrayThreadnumber;
    #endif
}
