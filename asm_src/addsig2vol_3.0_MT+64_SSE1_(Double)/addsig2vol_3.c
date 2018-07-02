    #ifndef __WIN32__
    #ifdef WIN32
    #define __WIN32__
    #endif
    #ifdef _WIN32
    #define __WIN32__
    #endif
    #ifdef __WIN32
    #define __WIN32__
    #endif
    #ifdef _WIN64
    #define __WIN32__
    #endif
    #ifdef WIN64
    #define __WIN32__
    #endif
    #ifdef _WINDOWS
    #define __WIN32__
    #endif
#endif

#ifdef __WIN32__
    #include <windows.h>
    #include <float.h>
#else
    #define _POSIX_C_SOURCE 199309L
    #define _POSIX_TIMERS 1
#endif

#ifdef __linux__
    #include <time.h>
    #include <unistd.h>
#endif

#include <math.h>
#include <stdio.h>

#ifndef _PORTABLEFPU_H_INCLUDED
#include "portableFPU.h"
#endif

#ifndef _PSTDINT_H_INCLUDED
#include "pstdint_new.h"
#endif

#include "addsig2vol_3.h"


//extern __declspec(dllimport) int _imp_pthread_join();
//extern __declspec(dllimport) int _imp_pthread_create();
/*#ifdef __cplusplus
extern "C" {
#endif
 __declspec(dllimport) int __imp_pthread_join();

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
__declspec(dllimport) __imp_pthread_create();

#ifdef __cplusplus
}
#endif*/


//enable/disable debugging output
#define addsig2vol_debug
#undef addsig2vol_debug

// enable/disable pthreads
#define p_threads

//enable/disable C-CODE version (disabled is Asm-code)
#undef C_CODE

#ifdef p_threads
#include "pthread.h"
#define NUMCORES 4 //32 //potential maximum value for init structs etc.
#else
#define NUMCORES 1
#endif

//define interpol-ratio of AScan (for asm fix)
#define INTERP_RATIO 5    // resizes ascan from 3000x1 -> 15000x1 (FIXED SIZE lin. interp.)

//define min-Voxel Number for X (parallel voxels in X_pipe; for 64bit code = 4voxel, for 32 bit = 2voxel)
#define MIN_VOXEL 4


//global variable
unsigned int addsig2vol_mode=0; //mode = single ascan, constant speed, 1=blocked, constant speed, 2= unblocked with soundmap


//for MT handling
static int64_t latency[NUMCORES]={0}; //uint64_t nSec
static int64_t throughput[NUMCORES]={0}; //uint64 MVoxel/s
static uint32_t nCores_bench = -1;
static uint32_t nCores = NUMCORES; //used value might be reduced by imagesize, benchmark etc

//TEST-CASE 0
//tic; x=100; image_1=addsig2vol_2(rand([3000 1]),single(ones(3,1)),single(ones(3,1)),single(ones(3,1)),single(ones(1,1)),single(ones(1,1)),single(ones(3,1)),uint32([x,x,x]),zeros([x,x,x]),2); toc, j=reshape(image_1,[x x x]); imagesc(j(:,:,1));
//inter_p test
//ascan=[repmat(1000,[1 100]) repmat(0,[1 2900])]';
//tic; x=100; image_1=addsig2vol_2(ascan,single(ones(3,1)),single(ones(3,1)),single(ones(3,1)),single(ones(1,1)),single(ones(1,1)),single(ones(3,1)),uint32([x,x,x]),zeros([x,x,x]),2); toc, j=reshape(image_1,[x x x]); imagesc(j(:,:,1));

//TEST-CASE 1
//x=128; addsig2vol_3(4), image_1=zeros([x,x,x]); rand('seed',0); x=128; for i=1:2 ascan=rand([1 3000])'; [image_1,kkk]=addsig2vol_3(ascan,rand([3 1],'single'),10.*rand([3 1],'single'),400.*rand([3 1],'single'),rand([1 1],'single'),rand([1 1],'single'),rand([1 1],'single'),uint32([x,x,x]),image_1); end, j=reshape(image_1,[x x x]); imagesc(j(:,:,1));
//j=reshape(image_1-image_2,[x x x]); figure; imagesc(j(:,:,128));

//// TEST--CASE Blocked
//count=2; senderPos = 0.01.*rand(3,count); receiverPos = 0.01.*rand(3,count); IMAGE_STARTPOINT = [0,0,0]; IMAGE_RESOLUTION= 0.001; Speed=1500; TimeInterval=1e-7; DataLength=3000; Data=zeros(3000,count); Data(floor(DataLength.*rand(count,1)),1:count)=1;
//x=100; bild=addsig2vol_3(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),zeros([x,x,x]));

// 4x times same parameter to be have compatible win64&linux64 calling convention
 typedef   struct   /*struct reordered to circumvent alignment problems, 64pointer then 32bit values*/
        {
        double*outz;
        double*AScanz;
        double*out_complexz;
        double*bufferz;
        float*pix_vectz;
		double*buffer_complexz;
        float*rec_posz;
        float*send_posz;
        float*speedz;
        float*resz;
        float*timeintz;
        double*AScan_complexz;
		    double *IMAGE_SUMz;
        double *IMAGE_SUM_complexz;
        unsigned int n_Yz;
        unsigned int n_Zz;
        unsigned int n_AScanz;
        unsigned int n_Xz;
        unsigned int qwb0;
        unsigned int qwb1;
        unsigned int qwb2;
        unsigned int qwb3;
		    } Addsig2vol_param;

  /* typedef   struct   //struct reorded to circumvent alignment problems, 64pointer then 32bit values
        {
        double*outz;
        double*out_complexz;
        double*AScanz;
        double*AScan_complexz;
        double*bufferz;
        double*buffer_complexz;
        double *IMAGE_SUMz;
        double *IMAGE_SUM_complexz;
        float*pix_vectz;
        float*rec_posz;
        float*send_posz;
        float*speedz;
        float*resz;
        float*timeintz;
        unsigned int n_AScanz;
		    unsigned int n_Xz;
        unsigned int n_Yz;
        unsigned int n_Zz;
        } Addsig2vol_param;*/



static double* out_real = NULL;
static double* out_complex = NULL;
static double* buffer_real = NULL;
static double* buffer_complex = NULL;

typedef   struct   /*struct reordered to circumvent alignment problems, 64pointer then 32bit values*/
       {
       unsigned int x;
       unsigned int y;
       unsigned int z;
       unsigned int currentJob;
       unsigned int currentTotalJob;
   } coordinate;

enum axis { XAXIS, YAXIS, ZAXIS };

static unsigned int L3CACHE_SIZE = 3072000;
static unsigned int L3CACHE_LINESIZE = 64;
static unsigned int jobs = 0;
static unsigned int nextJobWaiting = 0;
static unsigned int currentAScanThreads = 0;
static coordinate threadInfo[NUMCORES];
pthread_mutex_t threadMutex;
static float globalPixPointer[3];
static int segmentedAxis = ZAXIS;

//char *myArray[100] = { "str1", "str2", ... "str100" };

// For interlaced threads
static unsigned int  stepZ = 0;
static unsigned int  stepY = 0;
static unsigned int  stepX = 0;
static unsigned int  posZ = 0;
static unsigned int  posY = 0;
static unsigned int  posX = 0;

static unsigned int imgX = 0;
static unsigned int imgY = 0;
static unsigned int imgZ = 0;
static unsigned int halfJobs = 0;
static unsigned int fullJobs = 0;
static unsigned int halfStepZ = 0;
static unsigned int halfStepY = 0;
static unsigned int halfStepX = 0;
int algostarts = 0;

static unsigned int totalAscan = 0;
static unsigned int NtotalAscan = 0;
static unsigned int currentAscan = 0;
static double* global_outz = NULL;
static double* global_imagez = NULL;
static double* flip_imagez = NULL;

static unsigned int global_segmentElements = 0;


//CPUcount
uint64_t CPUCount(void);
uint64_t TimeCounter(void);

//fpu
void fpu_check(void);

//ellipsoide backprojections
void as2v_complex(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void as2v_complex_sm(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void as2v_c(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void as2v_MT(double*outz, double*AScanz, unsigned int n_AScanz, double*bufferz, float*pix_vectz,
		    unsigned int n_Xz, float*rec_posz, float*send_posz, float*speedz, float*resz,
		    float*timeintz, double*AScan_complexz,
		    double*buffer_complexz, double*out_complexz, unsigned int n_Yz, unsigned int n_Zz,
		    double* IMAGE_SUMz, double* IMAGE_SUM_complexz);

//Xsum and interpol
void xsum_complex(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void xsum_c(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);

//thread function
void *thread_function(void *arg);
void *thread_function_Z(void *arg);

//thread benchmark
void as2v_bench( uint64_t lat[],  uint64_t through[]);

///////////////////End declarations



double* as2v_addsig2vol_3(cArrayDouble* AScan_realz, cArrayDouble* AScan_complexz,
            cArrayFloat* pix_vectz, cArrayFloat* rec_posz, cArrayFloat* send_posz, cArrayFloat* speedz, float* resz, float* timeintz,
		    cArrayDouble* IMAGE_SUM_realz, cArrayDouble* IMAGE_SUM_complexz, cArrayDouble* out_image, cArrayDouble* out_buffer)
            {

                tssettimer();

                double* out_real = NULL;
                double* buffer_real = NULL;


        double* AScan_realz_ptr = AScan_realz->data;
        double* IMAGE_SUM_realz_ptr = IMAGE_SUM_realz->data;
        float* pix_vectz_ptr = pix_vectz->data;
        float *rec_vec_ptr=rec_posz->data;
        float *send_vec_ptr=send_posz->data;
        float *speedz_vec_ptr=speedz->data;

        unsigned int i;
        unsigned int n_AScan = AScan_realz->x;
        unsigned int n_AScan_block = AScan_realz->y;
        unsigned int n_Z = IMAGE_SUM_realz->z;
        unsigned int n_Y = IMAGE_SUM_realz->y;
        unsigned int n_X = IMAGE_SUM_realz->x;
        unsigned int n_IMAGE = IMAGE_SUM_realz->len;
        unsigned int n_rec_vec_block = rec_posz->y;
        unsigned int n_send_vec_block= send_posz->y;
        unsigned int n_speed_vec_block = speedz->y;

    	//check if X_Dim >= MIN_VOXEL
		if (n_X < MIN_VOXEL)
		{
		//out     = mxCreateDoubleMatrix(0, 0, mxCOMPLEX); //out  == empty
	  printf("Error: X-Dim has to be at least %d (X is %u, Y is %u, Z %u). Use Y-dim for such small sizes.\n",MIN_VOXEL,n_X,n_Y,n_Z);
	  return;
	  }
    //number of voxel
		n_IMAGE   = n_X * n_Y * n_Z;




	 ////benchmark of performance, selecting number of threads
    if (nCores_bench==-1) as2v_bench(&throughput[0], &latency[0]);

    //select if use or not use multithreading
     #ifdef addsig2vol_debug
      printf("selectedCore perf: %f , %f\n",((double)n_IMAGE/((double)throughput[nCores_bench-1] * 1000000))+((double)latency[nCores_bench-1]/1000000000),(((double)n_IMAGE/((double)throughput[0]*1000000))+((double)latency[0]/1000000000)));
      printf("nimage %i, throughput: %f , latency %f, %f, %f\n",n_IMAGE, (double)throughput[nCores_bench-1],(double)latency[nCores_bench-1],(double)throughput[0],(double)latency[0]);
     #endif

    if ( ( ((double)n_IMAGE/((double)throughput[nCores_bench-1] *1000000))+((double)latency[nCores_bench-1]/1000000000)) <= (((double)n_IMAGE/((double)throughput[0] *1000000))+((double)latency[0]/1000000000)))
         {  nCores = nCores_bench;}
    else {
      #ifdef addsig2vol_debug
      printf("Overhead to big, switch to single thread.\n");
      #endif
      nCores = 1;}



      addsig2vol_mode = 0;


        if(!out_image){
            out_real =  malloc(n_IMAGE*sizeof(double));
        }
        else{
            out_real = out_image->data;
        }
        if(!out_buffer){
            buffer_real = malloc(n_AScan*INTERP_RATIO*n_AScan_block*sizeof(double));
        }
        else{
            buffer_real = out_buffer->data;
        }
         double* AScan_complexz_ptr= NULL;
         double* IMAGE_SUM_complexz_ptr =NULL;
         double* out_complex = NULL;
          double* buffer_complex = NULL;

          NtotalAscan = n_AScan;
          totalAscan = n_AScan_block;
          currentAscan = 0;
          currentAScanThreads = 0;

        ///first Ascan
        // combined REAL & COMPLEX VERSION
        tsclock(80);

        as2v_MT(out_real, AScan_realz_ptr, n_AScan, buffer_real,
            pix_vectz_ptr, n_X, rec_vec_ptr, send_vec_ptr,
            speedz_vec_ptr, resz, timeintz, AScan_complexz_ptr, buffer_complex, out_complex, n_Y, n_Z, IMAGE_SUM_realz_ptr, IMAGE_SUM_complexz_ptr);
            tsclock(80);
            tsprint(80, TS_SEC);

        // //loop over ascans > 1
        // for (i = 2; i <= n_AScan_block; i++) {
        //     #ifdef addsig2vol_debug
        //     printf("as2v_MT call %i:\n", i);
        //     #endif
        //     currentAscan++;
        //
        //     //check for complex ascan only increase if available because NULL-Pointer +something -> not anymore a nullpointer!
        //
        // //    if (AScan_complexz != NULL)	AScan_complexz = AScan_complexz + (n_AScan * (i - 1)); //set to next value
        //
        //     if (1<n_rec_vec_block)  rec_vec_ptr  = rec_posz->data  + (3 * (i - 1)); else rec_vec_ptr = rec_posz->data;
        //     if (1<n_send_vec_block) send_vec_ptr = send_posz->data + (3 * (i - 1)); else send_vec_ptr = send_posz->data;
        //     if (1<n_speed_vec_block) speedz_vec_ptr = speedz->data + (1  * (i - 1)); else speedz_vec_ptr = speedz->data;
        //
        //
        //     // NOTE: pass out_real/out_complex as IMAGE_SUM to sum up the single scans
        //     as2v_MT(out_real, AScan_realz_ptr + (n_AScan * (i - 1)), n_AScan, buffer_real,
        //     pix_vectz_ptr, n_X, rec_vec_ptr, send_vec_ptr,
        //     speedz_vec_ptr, resz, timeintz, AScan_complexz_ptr, buffer_complex, out_complex, n_Y, n_Z, out_real, out_complex);
        //
        // }
        tsclock(83);

        if(!out_buffer){
            free(buffer_real);
        }
        tsclock(83);
        tsprint(83, TS_SEC);
	return out_real;

}

void setSegmentation(unsigned int elements){
    printf("Given segmentation: %i elements per task\n", elements);
    global_segmentElements = elements;
}


#ifdef BUILDMEX

double*  as2v_addsig2vol_3_mex(mxArray* AScan_realz, mxArray* AScan_complexz,
            mxArray* pix_vectz, mxArray* rec_posz, mxArray* send_posz, mxArray* speedz, mxArray* resz, mxArray* timeintz,
		    mxArray* IMAGE_SUM_realz, mxArray* IMAGE_SUM_complexz, mxArray* out_image, mxArray* out_buffer)
            {
            tsclock(80);
        double* AScan_realz_ptr = mxGetPr(AScan_realz);
        double* IMAGE_SUM_realz_ptr = mxGetPr(IMAGE_SUM_realz);
        float* pix_vectz_ptr = mxGetPr(pix_vectz);
        float *rec_vec_ptr= mxGetPr(rec_posz);
        float *send_vec_ptr= mxGetPr(send_posz);
        float *speedz_vec_ptr= mxGetPr(speedz);

        mwSize* scanDim      = (mwSize*) mxGetDimensions(AScan_realz);
        unsigned int n_AScan         = scanDim[0];        //gesamtanzahl elemente IN EINEM ASCAN!!!
        unsigned int n_AScan_block   = scanDim[1];    //2 dim %number of parallel


        unsigned int i;

        unsigned int n_IMAGE;
        mwSize* dim = (mwSize*) mxGetDimensions(IMAGE_SUM_realz);

        //IMAGE_ptr       = *(uint32_t*) mxGetPr(IMAGE_XYZ);
		unsigned int  n_X              = dim[0];//*IMAGE_ptr;
		unsigned int  n_Y              = dim[1];//*(IMAGE_ptr + 1);
        unsigned int  n_Z;

		//% workaround for matlab behaviour for size(IMAGE(1 x1x1))->reduced to 1 x1 (some kind of squeeze)
		if (mxGetNumberOfDimensions(IMAGE_SUM_realz) < 3)
		     n_Z = 1;
		else
			n_Z  = dim[2];//*(IMAGE_ptr + 2);
        n_IMAGE = n_X*n_Y*n_Z;

        #ifdef addsig2vol_debug
       //fpu_status
        fpu_check();
   	 #endif

	 ////benchmark of performance, selecting number of threads
    if (nCores_bench==-1) as2v_bench(&throughput[0], &latency[0]);

    //select if use or not use multithreading
     #ifdef addsig2vol_debug
      printf("selectedCore perf: %f , %f\n",((double)n_IMAGE/((double)throughput[nCores_bench-1] * 1000000))+((double)latency[nCores_bench-1]/1000000000),(((double)n_IMAGE/((double)throughput[0]*1000000))+((double)latency[0]/1000000000)));
      printf("nimage %i, throughput: %f , latency %f, %f, %f\n",n_IMAGE, (double)throughput[nCores_bench-1],(double)latency[nCores_bench-1],(double)throughput[0],(double)latency[0]);
     #endif

    if ( ( ((double)n_IMAGE/((double)throughput[nCores_bench-1] *1000000))+((double)latency[nCores_bench-1]/1000000000)) <= (((double)n_IMAGE/((double)throughput[0] *1000000))+((double)latency[0]/1000000000)))
         {  nCores = nCores_bench;}
    else {
      #ifdef addsig2vol_debug
      printf("Overhead to big, switch to single thread.\n");
      #endif
      nCores = 1;}



      addsig2vol_mode = 0;
       double* out_real;
       double* buffer_real ;

        if(!out_image){
            out_real =  malloc(n_IMAGE*sizeof(double));
        }
        else{
            out_real = mxGetPr(out_image);
        }
        if(!out_buffer){
            // TODO this malloc ist bigger => try to make this not scale
            buffer_real = malloc(n_AScan*INTERP_RATIO*n_AScan_block*sizeof(double));
        }
        else{
            buffer_real = mxGetPr(out_buffer);
        }


         double* AScan_complexz_ptr= NULL;
         double* IMAGE_SUM_complexz_ptr =NULL;
         double* out_complex = NULL;
          double* buffer_complex = NULL;

          NtotalAscan = n_AScan;
          totalAscan = n_AScan_block;
          currentAscan = 0;
          currentAScanThreads = 0;
          tsclock(80);
          tsprint(80, TS_MILI);
        ///first Ascan
        // combined REAL & COMPLEX VERSION
        tsclock(81);

        as2v_MT(out_real, AScan_realz_ptr, n_AScan, buffer_real,
            pix_vectz_ptr, n_X, rec_vec_ptr, send_vec_ptr,
            speedz_vec_ptr, resz, timeintz, AScan_complexz_ptr, buffer_complex, out_complex, n_Y, n_Z, IMAGE_SUM_realz_ptr, IMAGE_SUM_complexz_ptr);


        // //loop over ascans > 1
        // for (i = 2; i <= n_AScan_block; i++) {
        //     #ifdef addsig2vol_debug
        //     printf("as2v_MT call %i:\n", i);
        //     #endif
        //     currentAscan++;
        //
        //     //check for complex ascan only increase if available because NULL-Pointer +something -> not anymore a nullpointer!
        //
        // //    if (AScan_complexz != NULL)	AScan_complexz = AScan_complexz + (n_AScan * (i - 1)); //set to next value
        //
        //     if (1<n_rec_vec_block)  rec_vec_ptr  = rec_posz->data  + (3 * (i - 1)); else rec_vec_ptr = rec_posz->data;
        //     if (1<n_send_vec_block) send_vec_ptr = send_posz->data + (3 * (i - 1)); else send_vec_ptr = send_posz->data;
        //     if (1<n_speed_vec_block) speedz_vec_ptr = speedz->data + (1  * (i - 1)); else speedz_vec_ptr = speedz->data;
        //
        //
        //     // NOTE: pass out_real/out_complex as IMAGE_SUM to sum up the single scans
        //     as2v_MT(out_real, AScan_realz_ptr + (n_AScan * (i - 1)), n_AScan, buffer_real,
        //     pix_vectz_ptr, n_X, rec_vec_ptr, send_vec_ptr,
        //     speedz_vec_ptr, resz, timeintz, AScan_complexz_ptr, buffer_complex, out_complex, n_Y, n_Z, out_real, out_complex);
        //
        // }
        tsclock(81);
        tsprint(81, TS_MILI);

        tsclock(82);
                if(!out_buffer){
            free(buffer_real);
                }
        tsclock(82);
        tsprint(82, TS_MILI);
	return out_real;
}

#endif



//////////////////////////////////////////////////////////////////////////////////

unsigned int as2v_determineSegmentation(unsigned int imgX, unsigned int imgY, unsigned int imgZ){
    unsigned int elementsPerPackage = global_segmentElements;
    if(elementsPerPackage == 0)
    {
        // Definition nach Cache
        float cacheFracture = 2;
        unsigned int imagePrecision = 8;
        //Anzahl doubles per thread
        //Anzahl der jobs mindestens, um L3 Größe einzuhalten (Anzahl der Stücke, die mind. geschnitten werden müssen)
        float workPackages = (float) imgX*imgY*imgZ*cacheFracture*imagePrecision*nCores/L3CACHE_SIZE;

        // Check that there is at least one package for each thread
        if(workPackages < nCores){
             workPackages = nCores;
         }
        //berechne Elemente für jedes package
        elementsPerPackage = (unsigned int) imgX*imgY*imgZ/workPackages;
        // TODO Hier kann man Größe auch forcen für tests
        //elementsPerPackage = 4;
    }
    // Falls es zuwenig tasks für die Cores gibt: schneide mindestens nCores-Teile
    if (elementsPerPackage < 4) elementsPerPackage = 4;
    if(((float) imgX*imgY*imgZ/nCores) < elementsPerPackage )
    {
        elementsPerPackage = (unsigned int) ceil((float) imgX*imgY*imgZ/nCores);
        printf("WARNING: Tasks waren zu groß für %i threads, schneide kleinere tasks\n", nCores);
        if (elementsPerPackage < 4){
            printf("ERROR Zuviele Cores für dieses Bild!\n");
        }
    }

    // force layout: Für Interlacing muss hier einfach weiter geteilt werden, s.d. fullZ und fullY 0 werden, fullX = 0;
    posZ = ceil((float)elementsPerPackage/(imgX*imgY));
    posY = ceil((float)(elementsPerPackage-imgX*imgY*posZ)/imgX);
    posX = elementsPerPackage-posZ*imgY*imgZ-posY*imgX;
    // TODO posX needs to multiple of 4;
    //posX -= posX%4;

    if (posZ > 0){
        segmentedAxis = ZAXIS;
        //printf("ZAXIS\n");
        posX = 0;
        posY = 0;
        jobs = ceil((float)imgZ/posZ);
        stepZ = posZ;
        stepX = imgX;
        stepY = imgY;
        fullJobs = floor((float)imgZ/posZ);
        halfStepZ = imgZ - fullJobs*stepZ;
        if (halfStepZ == 0) halfStepZ = stepZ;
        else halfJobs = 1;
        jobs = ceil((float)imgZ/posZ);
    }
    if (posY > 0){
        segmentedAxis = YAXIS;
        //printf("YAXIS\n");

        posX = 0;
        posZ = 0;
        stepY=posY;
        stepX = imgX;
        stepZ = 1;
        fullJobs = floor((float)imgY/posY);
        halfStepY = imgY - fullJobs*stepY;
        if (halfStepY == 0) halfStepY = stepY;
        else halfJobs = 1;
        jobs = ceil((float)imgY/posY)*imgZ;
    }
    if (posX > 0){
        segmentedAxis = XAXIS;
        //printf("XAXIS\n");

        posY = 0;
        posZ = 0;
        stepY = 1;
        stepZ = 1;
        stepX = posX;
        fullJobs = floor((float)imgX/posX);
        halfStepX = imgX - fullJobs*stepX;
        if (halfStepX == 0) halfStepX = stepX;
        else halfJobs = 1;
        jobs = ceil((float)imgX/posX)*imgY*imgZ;
    }

    //printf("elementsPerPackage %i (%fpercent), jobs %i | posX, posY, posZ: %i %i %i | x, y, z: %i %i %i \n", elementsPerPackage, ((double)elementsPerPackage)/(imgX*imgY*imgZ)*100., jobs, posX, posY, posZ, imgX, imgY, imgZ);
    return jobs;
}

void as2v_MT(double*outz, double*AScanz, unsigned int n_AScanz, double*bufferz, float*pix_vectz,
		    unsigned int n_Xz, float*rec_posz, float*send_posz, float*speedz, float*resz,
		    float*timeintz, double*AScan_complexz,
		    double*buffer_complexz, double*out_complexz, unsigned int n_Yz, unsigned int n_Zz,
		    double *IMAGE_SUMz, double *IMAGE_SUM_complexz)

{
        tsclock(83);

        //pthread variables
        #ifdef p_threads
            pthread_t mythread[NUMCORES]; //numCPU -1
            Addsig2vol_param threadArg[NUMCORES]; //numCPU -1
            int rc = 0; //return-value from thread functions
        #endif
        //NOTE: This must be set on 0, otherwise results doesnt match
        unsigned int i,j;

        // Hier muss pix pointer für aktuellen AScan rein
        globalPixPointer[0] = *pix_vectz;
        globalPixPointer[1] = *(pix_vectz+1);
        globalPixPointer[2] = *(pix_vectz+2);

        imgX = n_Xz;
        imgY = n_Yz;
        imgZ = n_Zz;

        as2v_determineSegmentation(n_Xz, n_Yz, n_Zz);


        #ifdef addsig2vol_debug
        printf("Z-Dim multithreading\n");
        #endif

        //Generate parameter structs for Z-multithreading

        global_outz = outz;
        global_imagez = IMAGE_SUMz;

        for (j = 0; j<nCores; j++)
        {
            //fill parameter struct
            threadArg[j].outz=outz;//
            threadArg[j].AScanz=AScanz;
            threadArg[j].n_AScanz=n_AScanz;
            threadArg[j].bufferz=bufferz;
            //threadArg[j].pix_vectz=&(pix_vecz_buffer[j][0]) //will be done by each thread individually
            threadArg[j].n_Xz=stepX;
            threadArg[j].rec_posz=rec_posz;
            threadArg[j].send_posz=send_posz;
            threadArg[j].speedz=speedz;
            threadArg[j].resz=resz;
            threadArg[j].timeintz=timeintz;
            threadArg[j].AScan_complexz=AScan_complexz;
            threadArg[j].buffer_complexz=buffer_complexz;
            threadArg[j].out_complexz=out_complexz;//
            threadArg[j].n_Yz=stepY;
            threadArg[j].n_Zz=stepZ;/*n_Zz*/
            threadArg[j].IMAGE_SUMz=IMAGE_SUMz;//
            threadArg[j].IMAGE_SUM_complexz=IMAGE_SUM_complexz; //
            threadArg[j].qwb0 = j; // Thread ID, this is a hack
        }

        threadstats_init(nCores, jobs, totalAscan);
        tsclock(83);
        tsprint(83, TS_MILI);
        tsclock(84);

        //    interpol & X-SUM (in the case of NUMCORE=1 only call)
        // Für gesamten Scan aufeinmal! So dass threadArg[0} am Schluss auf 0 start
        for (int k = totalAscan-1; k > -1; k--)
        {
            threadArg[0].AScanz = AScanz + k*NtotalAscan;
            threadArg[0].bufferz = bufferz + k*NtotalAscan*INTERP_RATIO;
            #ifdef C_CODE
            xsum_c(&threadArg[0],&threadArg[0],&threadArg[0],&threadArg[0]);
            #else
            xsum_complex(&threadArg[0],&threadArg[0],&threadArg[0],&threadArg[0]);
             #endif

         }

        nextJobWaiting = 0;
        tsclock(84);
        tsprint(84, TS_MILI);

        tsclock(85);

        for (i=0;i<nCores;i++)
        {
            if (segmentedAxis == ZAXIS){
                rc = pthread_create( &(mythread[i]), NULL, thread_function_Z, &(threadArg[i]));
                if (rc) { printf("ERROR: return code from pthread_create() is %d\n", rc); return;}
            }else{
            #ifdef p_threads
            rc = pthread_create( &(mythread[i]), NULL, thread_function, &(threadArg[i]));
            if (rc) { printf("ERROR: return code from pthread_create() is %d\n", rc); return;}
            #endif
            }
        }
        // AT THIS POINT IMAGE_SUM_1 not OUT_0 anymore (from 0 t0 8?)
        //-> Because of threading! Not all finish the same time

        //catches threads again
        for (i=0;i<nCores;i++)
        {
            #ifdef p_threads
            rc = pthread_join ( mythread[i], NULL );
            if (rc) { printf("ERROR: return code from pthread_join() is %d\n", rc); return;}
            #endif
        }

        threadstats_free();

        tsclock(85);
        tsprint(85, TS_MILI);

    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void as2v_c(Addsig2vol_param* tt, Addsig2vol_param* t1, Addsig2vol_param* t2, Addsig2vol_param* t3)
{

      ///TODO: COMPLEX part!!!!

float dist_sv[3] = {0,0,0};
float dist_rv[3] = {0,0,0};
float factor = 0;

unsigned int index, image_index = 0;
unsigned int z,y,x,sampl = 0;

//decompose variables from struct
  	    double*outz = tt->outz;
        double*AScanz = tt->AScanz;
        double*out_complexz=tt->out_complexz;
        double*bufferz = tt->bufferz;
        float*pix_vectz = tt->pix_vectz;
		    double*buffer_complexz=tt->buffer_complexz;
        float*rec_posz = tt->rec_posz;
        float*send_posz = tt->send_posz;
        float*speedz = tt->speedz;
        float*resz = tt->resz;
        float*timeintz = tt->timeintz;
        double*AScan_complexz = tt->AScan_complexz;
		    double *IMAGE_SUMz = tt->IMAGE_SUMz;
        double *IMAGE_SUM_complexz = tt->IMAGE_SUM_complexz;
        unsigned int n_Yz = tt->n_Yz;
        unsigned int n_Zz = tt->n_Zz;
        unsigned int n_AScanz = tt->n_AScanz;
        unsigned int n_Xz = tt->n_Xz;

///bildgebung

factor =  INTERP_RATIO / (*speedz * *timeintz);

for (z=1; z<=n_Zz; z++)
	{
		//dist_sv[2] = pow(send_posz[2] - (((float)z* *resz)+pix_vectz[2]) ,2);
        dist_sv[2]  = send_posz[2] - (((float)z* *resz)+pix_vectz[2]);
        dist_sv[2] *= dist_sv[2];
		//dist_rv[2] = pow(rec_posz[2] - (((float)z* *resz)+pix_vectz[2]) ,2);
        dist_rv[2]  = rec_posz[2] - (((float)z* *resz)+pix_vectz[2]);
        dist_rv[2] *= dist_rv[2];

	for (y=1; y<=n_Yz; y++)
		{
		//dist_sv[1] = pow(send_posz[1] - ((y* *resz)+pix_vectz[1]) ,2);
        dist_sv[1]  = send_posz[1] - (((float)y* *resz)+pix_vectz[1]);
        dist_sv[1] *= dist_sv[1];
		//dist_rv[1] =  pow(rec_posz[1] - ((y* *resz)+pix_vectz[1]) ,2);
        dist_rv[1]  = rec_posz[1] - (((float)y* *resz)+pix_vectz[1]);
        dist_rv[1] *= dist_rv[1];

		for (x=1; x<=n_Xz; x++)
			{
				//dist_sv[0] = pow(send_posz[0] - (x* *resz),2);
				//dist_rv[0] =  pow(rec_posz[0] - (x* *resz),2);
				dist_sv[0] = send_posz[0] - ((x* *resz)+pix_vectz[0]);
				dist_rv[0] =  rec_posz[0] - ((x* *resz)+pix_vectz[0]);

				    //dist = (sqrt(dist_sv[1]+ dist_sv[2] + (dist_sv[0]*dist_sv[0])) + sqrt(dist_rv[1]+ dist_rv[2]+ (dist_rv[0]*dist_rv[0])) );
				index = (unsigned int) floor( ( sqrt(dist_sv[1]+ dist_sv[2] + (dist_sv[0]*dist_sv[0])) + sqrt(dist_rv[1]+ dist_rv[2]+ (dist_rv[0]*dist_rv[0])) ) * factor);
				    //printf("index:  %i\n\n", index);
				    //outz[image_index] = (double)index;
				if ((index >= n_AScanz*INTERP_RATIO) | (index < 0))
					outz[image_index] = IMAGE_SUMz[image_index]; //nix addiert
				else
					outz[image_index] = IMAGE_SUMz[image_index] + bufferz[index];//AScanz[index];

				image_index++;
			}
		}
	}

}

///////////////////////////////////////////////

void xsum_c(Addsig2vol_param* tt, Addsig2vol_param* t1, Addsig2vol_param* t2, Addsig2vol_param* t3)
{  ///TODO: COMPLEX part!!!!

float dist_sv[3] = {0,0,0};
float dist_rv[3] = {0,0,0};
float factor = 0;

unsigned int image_index = 0;
unsigned int i,j,sampl = 0;

double i_buffer=0.0;

double* sec_buffer;

//decompose variables from struct
  /* as2v_c(tt->outz, tt->AScanz, tt->n_AScanz, tt->bufferz, tt->pix_vectz,
		    tt->n_Xz, tt->rec_posz, tt->send_posz, tt->speedz, tt->resz,
		    tt->timeintz, tt->AScan_complexz,
		    tt->buffer_complexz, tt->out_complexz, tt->n_Yz, tt->n_Zz,
		    tt->IMAGE_SUMz, tt->IMAGE_SUM_complexz);*/

		double*outz = tt->outz;
        double*AScanz = tt->AScanz;
        double*out_complexz=tt->out_complexz;
        double*bufferz = tt->bufferz;
        float*pix_vectz = tt->pix_vectz;
		    double*buffer_complexz=tt->buffer_complexz;
        float*rec_posz = tt->rec_posz;
        float*send_posz = tt->send_posz;
        float*speedz = tt->speedz;
        float*resz = tt->resz;
        float*timeintz = tt->timeintz;
        double*AScan_complexz = tt->AScan_complexz;
		    double *IMAGE_SUMz = tt->IMAGE_SUMz;
        double *IMAGE_SUM_complexz = tt->IMAGE_SUM_complexz;
        unsigned int n_Yz = tt->n_Yz;
        unsigned int n_Zz = tt->n_Zz;
        unsigned int n_AScanz = tt->n_AScanz;
        unsigned int n_Xz = tt->n_Xz;

#ifdef addsig2vol_debug
for (i=0;i<n_AScanz*INTERP_RATIO;i++)
{	bufferz[i] = i; //set marking for NON-set or initalized values
}
#endif

//C-Version needs second buffer
sec_buffer = malloc(INTERP_RATIO * n_AScanz * sizeof(double));


//interp for first Samples
for (i=0;i<(unsigned int) floor(INTERP_RATIO/2);i++)
{	sec_buffer[i] = AScanz[0] * (unsigned int) ((i+1)/floor(INTERP_RATIO/2));
}

//almost all samples
for (i=0;i<n_AScanz-1;i++)
{    	for (j=0;j<INTERP_RATIO;j++)
	{	sec_buffer[(unsigned int) floor(INTERP_RATIO/2)+(i* INTERP_RATIO)+j] = AScanz[i+1] * j/INTERP_RATIO + AScanz[i] * (INTERP_RATIO-j)/INTERP_RATIO; //
	}
}

//interp for last Samples
for (i=0;i<(unsigned int) floor(INTERP_RATIO/2)+1;i++)
{	sec_buffer[n_AScanz*INTERP_RATIO-(unsigned int) floor(INTERP_RATIO/2)-1+i] = AScanz[n_AScanz-1] * ((floor(INTERP_RATIO/2)+1-i) / (floor(INTERP_RATIO/2)+1));
}
///end interp



/////xsum
sampl = (unsigned int)(ceil((float)1.7*(( *resz / *speedz)/ (*timeintz/INTERP_RATIO)) /2)); //halbe breite

i_buffer = 0;
for (i=0;i<sampl;i++)
{	i_buffer = i_buffer + sec_buffer[i]/(2*sampl);
}

for (i=0;i<sampl;i++)
{	if (i+sampl<n_AScanz*INTERP_RATIO){
    i_buffer = i_buffer +sec_buffer[i+sampl]/(2*sampl);
	bufferz[i] = i_buffer;}
}

for (i=sampl;i<(n_AScanz*INTERP_RATIO)-sampl;i++)
{ if (i+sampl<n_AScanz*INTERP_RATIO){
    i_buffer = i_buffer + sec_buffer[i+sampl]/(2*sampl) - sec_buffer[i-sampl]/(2*sampl);
	bufferz[i] = i_buffer;}
}

for (i=n_AScanz*INTERP_RATIO-sampl;i<n_AScanz*INTERP_RATIO;i++)
{	if (i-sampl>=0){
    i_buffer = i_buffer - sec_buffer[i-sampl]/(2*sampl);
	bufferz[i] = i_buffer / (sampl-(n_AScanz*INTERP_RATIO)-i);}
}
/////end xsum

//free sec_buffer (only buffer needed now)
free(sec_buffer);

}

///////////////////////////////////////////////

void *thread_function_Z(void *argument)
{
    // Place variables on each thread stack
    volatile Addsig2vol_param* arg = (Addsig2vol_param*) argument;
    volatile unsigned int id = arg->qwb0; // TODO change volatile depending on benchmarking
    tsclock(id);
    volatile unsigned int x = id+nCores; // TODO change volatile depending on benchmarking
    volatile unsigned int currentTotalJob = 0; //threadInfo[id].currentTotalJob;

    volatile float pixPtr[3] = { globalPixPointer[0],  globalPixPointer[1],  globalPixPointer[2] };
    arg->pix_vectz = &(pixPtr[0]);
    volatile float resz = *(arg->resz);
    volatile float pixZ = globalPixPointer[2];

    volatile unsigned int currentZn;
    volatile unsigned int totalElementJumps;
    volatile unsigned int nextStepZ;
    // TODO fix this AScan-incremental error
    volatile unsigned int localCurrentAscan = 0;
    volatile unsigned int deltaAscan = 0;
    volatile char switchvar = 0;

    // local variables
    volatile unsigned int _stepZ = stepZ;
    volatile unsigned int _diskXY = imgX*imgY;
    volatile unsigned int _taskCheck = fullJobs+halfJobs-1;
    volatile unsigned int _halfStepZ = halfStepZ;
    volatile unsigned int _totalScans = NtotalAscan;

    //printf("jobs %i | posX, posY, posZ: %i %i %i | x, y, z: %i %i %i \n", jobs, posX, posY, posZ, imgX, imgY, imgZ);

    tsclock(id);
    tsprint(id, TS_MILI);
    while(1){

        //tstimer(id);



        //Look up next task number (and next AScan number)
        pthread_mutex_lock(&threadMutex);
        // printf("T%i | calculated task was %i/%i, Ascan %i/%i, currentZn %i, nextStepZ %i totalElementJumps %i\n", id, currentTotalJob, jobs, localCurrentAscan, totalAscan, currentZn, nextStepZ, totalElementJumps);

        //threadstats_markStartTask(id, nextJobWaiting, currentAScanThreads); // get updated numbers
        currentTotalJob = nextJobWaiting;
        nextJobWaiting++;
        pthread_mutex_unlock(&threadMutex);

        // if(currentTotalJob >=jobs){
        //     return;
        // }

        // // If all AScans are done, return to main thread
        if(localCurrentAscan >= totalAscan){
            //printf("T%i | out (task %i/%i)\n", id, currentTotalJob, jobs);
            // tstimer(x);
            // threadstats_markEndTask(x, currentTotalJob, localCurrentAscan);

            return NULL;
        }

        currentZn = currentTotalJob*_stepZ;
        totalElementJumps = currentZn*_diskXY;
        if(currentTotalJob == _taskCheck) nextStepZ = _halfStepZ;
        else nextStepZ = _stepZ;




    //     // printf("T%i| .........setting up task  %i/%i \n", id, currentTotalJob, jobs);
    //     // printf("T%i| .........setting up Ascan %i/%i, delta %i \n", id, internAscanNew, totalAscan, deltaAscan);
        pixPtr[2] = pixZ + currentZn * resz;
        arg->n_Zz = nextStepZ;/*n_Zz*/
        arg->outz = global_outz + totalElementJumps;//
        arg->IMAGE_SUMz = global_imagez + totalElementJumps;//


        #ifdef C_CODE
            as2v_c(arg,arg,arg,arg); //compatible win64 & linxu64 function-call
        #else
            if (addsig2vol_mode==0) as2v_complex(arg,arg,arg,arg);
            if (addsig2vol_mode==2) as2v_complex_sm(arg,arg,arg,arg);
        #endif


    // //    arg->IMAGE_SUM_complexz+=totalElementJumps; //
//    printf("T%i| Next task: %i currentZ %i nextStepZ %i imgX %i imgY %i posZ %i totalElementJumps %i \n", id, currentTotalJob, currentZn, arg->n_Zz, imgX, imgY, posZ, totalElementJumps);

        //tstimer(x);
        //threadstats_markEndTask(x, currentTotalJob, localCurrentAscan);

    }
    return NULL;
}


void *thread_function(void *argument)
{
    Addsig2vol_param* arg = (Addsig2vol_param*) argument;
    unsigned int id = arg->qwb0;
    unsigned int currentTotalJob = threadInfo[id].currentTotalJob;
    tstimer(id);
    threadstats_markStartTask(id, currentTotalJob, currentAscan);
    unsigned int firstJobDone = 0;
    float pix_vecz_buffer[3];
    float resz = *(arg->resz);
    unsigned int currentXn = threadInfo[id].x;
    unsigned int currentYn = threadInfo[id].y;
    unsigned int currentZn = threadInfo[id].z;
    unsigned int currentJob = threadInfo[id].currentJob;

    unsigned int internJob;
    unsigned int totalElementJumps;
    unsigned int nextStepX;
    unsigned int nextStepY;
    unsigned int nextStepZ;
    unsigned int oldCurrentX;
    unsigned int oldCurrentY;
    unsigned int oldCurrentZ;
    unsigned int jumps = 0;
    unsigned int jumpsZ;
    unsigned int delta;
    int x = id+nCores;



    while(1){
        if(firstJobDone > 0)
        {
            tstimer(id);
            threadstats_markStartTask(id, currentTotalJob, currentAscan);
        }
        #ifdef C_CODE
            as2v_c(arg,arg,arg,arg); //compatible win64 & linxu64 function-call
        #else
            if (addsig2vol_mode==0) as2v_complex(arg,arg,arg,arg);
            if (addsig2vol_mode==2) as2v_complex_sm(arg,arg,arg,arg);
        #endif
        firstJobDone= 1;

        pthread_mutex_lock(&threadMutex);
        internJob = nextJobWaiting;
        nextJobWaiting++;
        pthread_mutex_unlock(&threadMutex);

        // Anzahl der Jobsprünge (global)
        delta = (internJob-currentTotalJob);
        //printf("T%i| internJob %i currentTotalJob %i delta %i\n", id, internJob, currentTotalJob, delta);

        // nächste Jobnummer dieses threads
        currentTotalJob += delta;
        currentJob += delta;
        jumps = currentJob/(fullJobs+halfJobs);
        currentJob = currentJob % (fullJobs+halfJobs);

        // printf("T%i: delta %i, nextinternJob %i, modulo %i, jumps %i  \n", id, delta, currentTotalJob, currentJob, jumps);
        // printf("T%i, lastCurrentZN %i, lastNextStep: %i \n", id, currentZn, nextStepZ);

        switch(segmentedAxis){
            case ZAXIS:;

            currentZn = currentTotalJob*posZ;

            if(currentZn >= imgZ){
                tstimer(x);
                threadstats_markEndTask(x, currentTotalJob, currentAscan);
                return NULL;
            } // finished, last case
            nextStepZ = posZ;
            totalElementJumps = delta*posZ*imgX*imgY;
            if(currentTotalJob == fullJobs+halfJobs-1){
                nextStepZ = halfStepZ;
            }
            nextStepX = imgX;
            nextStepY = imgY;
            //wait(0.000000001);
            break;

                case YAXIS:;
                //print("YAXIS\n");
                oldCurrentY = currentYn;
                currentYn = currentJob*posY;
                currentZn += jumps;
                if(currentZn >= imgZ){
                    tstimer(x);
                    threadstats_markEndTask(x, currentTotalJob, currentAscan);
                    return NULL;} // finished, last case

                    nextStepY = stepY;
                    if(jumps == 0) totalElementJumps = delta*stepY*imgX;
                    else totalElementJumps = (imgY-oldCurrentY)*imgX + currentYn*imgX +  (jumps-1)*imgY*imgX;
                    if(currentJob ==  fullJobs+halfJobs-1){
                        nextStepY = halfStepY;

                    }
                    nextStepX = imgX;
                    nextStepZ = 1;
                    break;
                    case XAXIS:;
                    //print("XAXIS\n");
                    oldCurrentX = currentXn;
                    currentXn = currentJob*stepX;
                    // Insgesamte Anzahl an jumps (overflow jumps)
                    jumpsZ = floor(currentYn+jumps/imgY);
                    if(jumpsZ == 0) currentYn = currentYn + jumps;
                    else currentYn = jumps%imgY;
                    currentZn += jumpsZ;
                    if(currentZn >= imgZ){
                        tstimer(x);
                        threadstats_markEndTask(x, currentTotalJob, currentAscan);
                        return NULL;} // finished, last case

                        nextStepX = stepX;
                        if(jumps == 0) totalElementJumps = nCores*stepX;
                        // TODO ausrechnen!!
                        else totalElementJumps = (imgX-oldCurrentX) + currentXn + jumps*imgX;

                        if(currentJob== fullJobs+halfJobs-1){
                            nextStepX = halfStepX;
                        }
                        nextStepY = 1;
                        nextStepZ = 1;
                        break;

                    }

                    pix_vecz_buffer[0]=globalPixPointer[0]+currentXn* (resz);
                    pix_vecz_buffer[1]=globalPixPointer[1]+currentYn* (resz);
                    pix_vecz_buffer[2]=globalPixPointer[2]+currentZn* (resz);

                    // print("T%i: current x,y,z: %i %i %i\n", id, currentXn, currentYn, currentZn);
                    // print("T%i: nextsteps x,y,z: %i %i %i, totalElementJumps: %i\n", id, nextStepX, nextStepY, nextStepZ, totalElementJumps);
                    arg->outz+=totalElementJumps;//
                    arg->pix_vectz=&(pix_vecz_buffer[0]);//
                    arg->n_Xz=nextStepX;
                    arg->out_complexz+=totalElementJumps;//
                    arg->n_Yz=nextStepY;
                    arg->n_Zz=nextStepZ;/*n_Zz*/
                    arg->IMAGE_SUMz+=totalElementJumps;//
                    arg->IMAGE_SUM_complexz+=totalElementJumps; //
                    arg->n_Xz=nextStepX;

                    //printf("T%i| currentTotalJob %i moduloJob %i jumps %i delta %i\n", id, currentTotalJob, currentJob, jumps, delta);
                    //printf("T%i| totalElementJumps %i currentZn %i nextStepZ %i\n", id, totalElementJumps,currentZn, nextStepZ);
                    //printf("T%i| currentZn %i nextStepZ %i\n", id, currentZn, nextStepZ);

                    tstimer(x); // Octave crashes, if variable not from stack(?)
                    threadstats_markEndTask(x, currentTotalJob, currentAscan);

                }
                return NULL;
            }

////////////////////////////////////////////////////////////////////
void as2v_benchLocal(){
    as2v_bench(&throughput[0], &latency[0]);
}

void as2v_bench(uint64_t throughput[], uint64_t latency[])
{

    #pragma fenv_access (on)
    uint32_t i,j,l,k;
	uint32_t n_AScan;
	uint32_t n_AScan_block;
	int n_X;
	int n_Y;
	int n_Z;
	int n_IMAGE;
	// mwSize setImageDim[3];
	// mwSize setBufferDim[2];
	// mwSize setAScanDim[2];
	float pix_vec_bench[3]  = {2,1,77};
	//float pix_vec_bench[3]  = {2,1,77};
	float send_vec_bench[3] = {1,2,2};
	float rec_vec_bench[3]  = {5,6,9};
    float float_bench = 1;

	double* pr;
   /* double*pr1;
    double* pr2;
    double*pr3;
    double*pr4;*/

    #define MAX_AVERAGE 128
    uint64_t average_buffer[MAX_AVERAGE]={0};
    uint64_t throughput_sort_buffer[NUMCORES]={0};
    uint64_t counter, counter2, stdabw, mittelw;
    uint64_t bench_ref_time=0; //ns

    uint64_t minBenchTime = 100000000; //ns

    uint32_t minAverage = 7;
    uint32_t average = 8;
    uint32_t nVoxel_throughput = 128;

    double *AScan_bench;
	double *buffer_bench;
	double *out_bench;
	double *image_sum_bench;
	double *time_bench;

    //autodetect free memory
   /*  do { nVoxel_throughput=(uint32_t) nVoxel_throughput/2;
       pr = mxMalloc(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput * sizeof(double));
     if (pr!=NULL)  mxFree(pr);
       nVoxel_throughput=(uint32_t) nVoxel_throughput/2; //nocnmal halbieren f�r 2Bilder
   } while (pr==NULL);*/

  printf("Benchmarking System:\nLatency with %d Byte Image\nThroughput with %d kByte Image.\n",(int)(MIN_VOXEL*NUMCORES*sizeof(double)),(int)(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput*sizeof(double))/1024);
  printf("Threads | Throughput in MVoxel/s                     | Latency in nsec                              | Median-Size\n");
  printf("        | Median | Mean   | Min    | Max    | Std    | Median    | Min       | Max       | Std      | \n");
  printf("-------------------------------------------------------------------------------------------------------------------\n");

  //check for minimum amount
  if (nVoxel_throughput<NUMCORES*MIN_VOXEL) nVoxel_throughput=NUMCORES*MIN_VOXEL;
  /////kill eventually running tic-toc!!! (not needed)
   //mexCallMATLAB(0, NULL, 0, NULL, "toc");

	////fix variables
    time_bench = malloc(sizeof(double));
	n_AScan        = 3000;        //gesamtanzahl elemente IN EINEM ASCAN!!!
	n_AScan_block  = 100;    //2 dim %number of parallel

    //buffer

    printf("buffer size %i\n", (n_AScan*INTERP_RATIO));
    buffer_bench = malloc(n_AScan*INTERP_RATIO *sizeof(double));

    //pr1=pr;

	AScan_bench = malloc(n_AScan *n_AScan_block*sizeof(double));




    ///////////////throughput calc (assumption for Latency = Zero because working in saturation)

	n_X             = nVoxel_throughput;
	n_Y             = nVoxel_throughput;
    n_Z             = nVoxel_throughput;

    //number of voxel
		n_IMAGE         = n_X * n_Y * n_Z;

	out_bench = malloc(n_IMAGE*sizeof(double));
    image_sum_bench = malloc(n_IMAGE*sizeof(double));

        printf("Voxels in x, y, z: [%i, %i, %i]\n", n_X, n_Y, n_Z );
        printf("Voxels total: %i\n", n_IMAGE);
        printf("Buffer total: %i\n", INTERP_RATIO * n_AScan);

       //first call to get REAL memory from system (WARMUP later on faster)
      //set nCore!!!
	  nCores = 1;
      do{
      counter = TimeCounter();
      as2v_MT((out_bench), (AScan_bench), n_AScan, (buffer_bench),
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL, NULL, NULL, n_Y, n_Z, (image_sum_bench), NULL);
       counter2 = TimeCounter();
       } while (counter2<counter); //retry on error like used wrong core

        //get througput time time for setup average
        bench_ref_time= counter2-counter;
        //print("benchreftime %llu c2:%llu c:%llu\n",bench_ref_time, counter, counter2);
        average = (uint32_t) ((minBenchTime/bench_ref_time)+1); //integer ceil

        if (average < minAverage) average = minAverage;
        if (average > MAX_AVERAGE) average = MAX_AVERAGE;


	for (i=1;i<=NUMCORES;i++)
	{
	  //set nCore!!!
		nCores = i;

	    //benchmark throughput
        //mexCallMATLAB(0, NULL, 0, NULL, "tic");

		for (j=0;j<average;j++)
         {
         do{ counter=TimeCounter();
          //no sizeof(double) needed because compilers assumes already double as datatype for pointer!!!
		  as2v_MT((out_bench), (AScan_bench), n_AScan, (buffer_bench),
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL, NULL, NULL, n_Y, n_Z, (image_sum_bench), NULL);
          /*as2v_MT(pr3, pr2, n_AScan, pr1,
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL,NULL, NULL, n_Y, n_Z, pr4,NULL);*/
		  //mexCallMATLAB(1, &time_bench, 0, NULL, "toc");
          //print("%llu ms\n",((TimeCounter()-counter)/1000000));
          counter2=TimeCounter();
         } while(counter2<counter); //retry on error like used wrong core
          average_buffer[j]= counter2-counter;
        }
        //tsprint(1,TS_MILI);
        //tsclear(1);

        //bubblesort (small time top) for median
		  for (k=average-1;k>0;k--)
          {  for (l=average-1;l>0;l--){
                 if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;}
             }
          }

         mittelw=0;
         for  (k=0;k<average;k++)  {mittelw=average_buffer[k]+mittelw;}
         mittelw=(uint64_t) mittelw/average;

         stdabw=0;
         for  (k=0;k<average;k++)
         {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);}
         stdabw = (uint64_t) sqrt( (int64_t)stdabw/(int64_t)average);

         #ifdef addsig2vol_debug
          //print("minTime: %f, timeCounter: %llu counter: %llu\n",(float)(uint64_t)bench_ref_time,counter2,counter);
          // for (k=0;k<average;k++){ print("%llu average\n",average_buffer[k]);}
         #endif

        //minimum selection
        bench_ref_time = (uint64_t) average_buffer[0];
        //median selection
        bench_ref_time = (uint64_t) average_buffer[(uint32_t) ceil(average/2)];

        //throughput[i-1] = (uint64_t) ((double)n_IMAGE/ ((double) ((uint32_t) bench_ref_time)/1000)); //gekuerzt ns und MVoxel ->Mvoxel/s // index 0 = core1; index1 = core 2 etc...(!)
        throughput[i-1] =  ((1000*(uint64_t)n_IMAGE)/ bench_ref_time); //gekuerzt ns und MVoxel ->Mvoxel/s // index 0 = core1; index1 = core 2 etc...(!)

        printf("%7i |%7llu |%7llu |%7llu |%7llu |%7llu", i, (int64_t)throughput[i-1],((1000*(int64_t)n_IMAGE)/(int64_t)mittelw), (1000*(int64_t)n_IMAGE)/((int64_t)average_buffer[average-1]), (1000*(int64_t)n_IMAGE)/((int64_t)average_buffer[0]), (1000*(int64_t)n_IMAGE)/((int64_t)stdabw));


       //get throughput time time for setup average
       average = (uint32_t) ((minBenchTime/bench_ref_time)+1); //in ns; +1=integer-ceil
       if (average < minAverage) average = minAverage;
       if (average > MAX_AVERAGE) average = MAX_AVERAGE;

		 // printf("%e MVoxel/s for %i threads, %e sec through\n",(double)(4*n_IMAGE)/(*mxGetPr(time_bench)),n_IMAGE,*mxGetPr(time_bench));
      //  	  printf("%f",(float)(((4*2000000)/(*(double*)mxGetPr(time_bench)))));

//         //first call to get REAL memory form system (later on faster)
//       as2v_MT(mxGetPr(out_bench), mxGetPr(AScan_bench), n_AScan, mxGetPr(buffer_bench),
// 			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
// 			     &float_bench,
// 			     &float_bench, &float_bench, mxGetPi(AScan_bench), mxGetPi(
// 				     buffer_bench), mxGetPi(
// 				     out_bench), n_Y, n_Z, mxGetPr(image_sum_bench), mxGetPi(
// 				     image_sum_bench));

	   ///////////////////////////////////////////////////////benchmark latency
	//	 mexCallMATLAB(0, NULL, 0, NULL, "tic");

		for (j=0;j<minAverage;j++)
         {
          do{
          counter=TimeCounter();
			//no sizeof(double) needed because compilers assumes already double as datatype for pointer!!!
		  as2v_MT((out_bench), (AScan_bench), n_AScan, (buffer_bench),
			     &pix_vec_bench[0],(unsigned int) MIN_VOXEL, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL, NULL, NULL, (unsigned int) 1,(unsigned int) NUMCORES, (image_sum_bench), NULL);

		//  mexCallMATLAB(1, &time_bench, 0, NULL, "toc");
        //latency[i] =(*mxGetPr(time_bench)/minAverage)/1;	 //assumption for 1 PIXEL!
		 counter2 = TimeCounter(); } while(counter2<counter); //retry on error like used wrong core
         average_buffer[j]= counter2-counter;
        }
        //bubblesort (small time top)
		 for (k=minAverage-1;k>0;k--)
         {  for (l=minAverage-1;l>0;l--){
               if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;}
            }
         }
         #ifdef addsig2vol_debug
            //for (k=0;k<minAverage;k++){ print("%llu average\n",average_buffer[k]);  }
        #endif

         mittelw=0;
         for  (k=0;k<minAverage;k++)  {mittelw=average_buffer[k]+mittelw;}
         mittelw=(uint64_t) mittelw/minAverage;

         stdabw=0;
         for  (k=0;k<minAverage;k++)
         {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);}
         stdabw = (uint64_t) sqrt((int64_t) stdabw/minAverage);

        //minimum selection
        latency[i-1] = (uint64_t) average_buffer[0];	 //assumption for 1 PIXEL! // index 0 = core1; index1 = core 2 etc...(!)
        //median selection
        latency[i-1] = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
	    printf("|%10llu |%10llu |%10llu |%10llu |%8i\n", latency[i-1],average_buffer[0],average_buffer[minAverage-1],stdabw,average);
  		}

		for (i=NUMCORES;i>0;i--)
        { throughput_sort_buffer[i-1]=throughput[i-1];
		}

      //bubblesort biggest on top (bottom up)
      for (i=NUMCORES-1;i>=1;i--)
      { if (throughput_sort_buffer[i]>(throughput_sort_buffer[i-1])) throughput_sort_buffer[i-1] = throughput_sort_buffer[i]; }

      //find core-number according to perf.-value
      for (i=0;i<NUMCORES;i++)
      { if (throughput[i] == throughput_sort_buffer[0]) break; } //break on first core-number which fits

      //set up used Cores
      nCores = i+1;
      //backup in nCores_bench
      nCores_bench = i+1;

      switch (i+1)
      {
      case 1:
          printf("Detected Single-core System, 1 thread prefered\n");
          break;
      case 2:
          printf("Detected Dual-core or Hyperthreading System, 2 threads prefered\n");
          break;
      case 3:
          printf("Detected Triple-core or Hyperthreading System, 3 threads prefered\n");
          break;
      case 4:
          printf("Detected Quad-Core system, 4 threads prefered\n");
          break;
      case 8:
          printf("Detected Octa-Core system (or Quadcore with HT), 8 threads prefered\n");
          break;
      case 16:
          printf("Detected 16 Core system, 16 threads prefered\n");
          break;
      default:
          printf("Detected %i-core system, %i threads prefered (?)\n",i+1,i+1);
      }


      //benchmark perf-per size
      printf("\nPerformance for various imagesize in Voxel (with potentially %i Cores)\n",nCores);
      printf("     Voxel | Throughput in kVoxel/s         | Time in mikros  | Malloc time (mikro-sec)\n");
      printf("           | Median   | Mean     | Std      |           | Median | mean  | Std  | min  | max  \n");
      printf("------------------------------------------------------------------------------------------------\n");
      printf("voxel throughput maximal: %i\n", nVoxel_throughput*nVoxel_throughput*nVoxel_throughput);
      for (i=MIN_VOXEL;i<=(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput);i=i*2)
      {
          //print("i: %i, Nx: %i, Nz: %i \n", i, MIN_VOXEL, (uint32_t) floor(i/MIN_VOXEL) );


          for (j=0;j<minAverage;j++)
          {

              do {
                  counter=TimeCounter();
                  //no sizeof(double) needed because compilers assumes already double as datatype for pointer!!!
                  as2v_MT((out_bench), (AScan_bench), n_AScan, (buffer_bench),
                  &pix_vec_bench[0], (uint32_t) MIN_VOXEL, &rec_vec_bench[0], &send_vec_bench[0],
                  &float_bench,
                  &float_bench, &float_bench, NULL, NULL, NULL, (uint32_t) 1, (uint32_t) (i/MIN_VOXEL), (image_sum_bench), NULL);
                  counter2 = TimeCounter();  } while(counter2<counter); //retry on error like used wrong core
              average_buffer[j]= counter2-counter;

          }

          //bubblesort (small time top)
          for (k=minAverage-1;k>0;k--)
          {  for (l=minAverage-1;l>0;l--){
                     if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;}
                 }
          }
         mittelw=0;
         for  (k=0;k<minAverage;k++)  {mittelw=average_buffer[k]+mittelw;}
         mittelw=(uint64_t) mittelw/minAverage;

         stdabw=0;
         for  (k=0;k<minAverage;k++)
         {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);}
         stdabw = (uint64_t) sqrt((int64_t) stdabw/minAverage);

          #ifdef addsig2vol_debug
           //for (k=0;k<minAverage;k++){ print("%llu values-range \n",average_buffer[k]);}
          #endif

          //minimum selection
          counter2 = (uint64_t) average_buffer[0];	 //assumption for 1 PIXEL! // index 0 = core1; index1 = core 2 etc...(!)
          //median selection
          counter2 = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
          // print("%10luu , %10llu, %10llu \n",TimeCounter, counter, TimeCounter());
          //counter=pow(2,64)-counter;

          //mexCallMATLAB(1, &time_bench, 0, NULL, "toc");
          //print("%llu ms\n",((TimeCounter()-counter)/1000000));
          printf("%10i | %8llu | %8llu | %8llu |%8llu ",i,(1000000*(uint64_t)i)/(counter2),(1000000*(uint64_t)i)/(mittelw),((uint64_t)i*1000000)/(stdabw),(uint64_t)(counter2/1000) );

          //fix  (UGLY!!!!)
          nCores=nCores_bench;

          //benchmark mem-alloc
          minAverage=100;

              for (j=0;j<minAverage;j++)
              {
                  do {
                      free(image_sum_bench);

                      counter=TimeCounter();

                      image_sum_bench = malloc(n_IMAGE*sizeof(double));

                     counter2 = TimeCounter(); } while(counter2<counter); //retry on error like used wrong core
                     average_buffer[j]= counter2-counter;
              }
              //bubblesort (small time top)
              for (k=minAverage-1;k>0;k--)
              {  for (l=minAverage-1;l>0;l--){
                  if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;}
              }
              }

              mittelw=0;
             for  (k=0;k<minAverage;k++)  {mittelw=average_buffer[k]+mittelw;}
             mittelw=(uint64_t) mittelw/minAverage;

             stdabw=0;
             for  (k=0;k<minAverage;k++)
             {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);}
             stdabw = (uint64_t) sqrt((int64_t) stdabw/minAverage);

             //minimum selection
             counter2 = (uint64_t) average_buffer[0];	 //assumption for 1 PIXEL! // index 0 = core1; index1 = core 2 etc...(!)
             //median selection
             counter2 = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
             printf("| %8llu | %8llu | %8llu |%8llu |%8llu\n", (uint64_t)counter2/1000,(uint64_t)mittelw/1000,(uint64_t)stdabw/1000,(uint64_t) average_buffer[0]/1000,(uint64_t) average_buffer[minAverage-1]/1000);

          }

        //check status
        fpu_check();


        //free memory
		  free(out_bench);
		  free(image_sum_bench);
        //free memory
		  free(buffer_bench);
		  free(AScan_bench);
		  free(time_bench);
}


uint64_t TimeCounter(void) {
uint64_t counter;
    #ifndef __WIN32__
        #ifdef WIN32
        #define __WIN32__
        #endif
        #ifdef _WIN32
        #define __WIN32__
        #endif
        #ifdef __WIN32
        #define __WIN32__
        #endif
        #ifdef _WIN64
        #define __WIN32__
        #endif
        #ifdef WIN64
        #define __WIN32__
        #endif
        #ifdef _WINDOWS
        #define __WIN32__
        #endif
    #endif

    #ifdef __WIN32__ //fitting windows64 AND windows32
    #include <windows.h>
    //register uint64_t temp;
    uint64_t iFreq, iCount;
    QueryPerformanceFrequency((LARGE_INTEGER*)&iFreq);
    QueryPerformanceCounter((LARGE_INTEGER*)&iCount);
    //counter = (uint64_t) (1000000000*((double)iCount/(double)iFreq)); //f�r nSekunden (balancing der multiplikation)
    counter = (uint64_t) ((1000000000*iCount)/iFreq); //f�r nSekunden (balancing der multiplikation)

//#ifdef addsig2vol_debug
    //printff("iFreq:%llu, iCount:%llu, Counter:%llu\n",iFreq,iCount,counter);
    //#endif
    #endif

    #ifdef __linux__

    struct timespec time_str;
    struct timespec time_res;
    clock_gettime(CLOCK_REALTIME, &time_str);
    counter = (uint64_t) time_str.tv_nsec + (uint64_t) (time_str.tv_sec*1000000000);

    /*
    clock_gettime(CLOCK_REALTIME, &time_str); clock_getres(CLOCK_REALTIME, &time_res);
    printf("CLOCK_REALTIME clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec);
    clock_gettime(CLOCK_MONOTONIC, &time_str); clock_getres(CLOCK_MONOTONIC, &time_res);
    printf("CLOCK_MONOTONIC clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_str); clock_getres(CLOCK_PROCESS_CPUTIME_ID, &time_res);
    printf("CLOCK_PROCESS_CPUTIME_ID clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec);
    */
      //int clock_getres(CLOCK_REALTIME, struct timespec *res);
    //int clock_gettime(clockid_t clk_id, struct timespec *tp);
    #endif
    return counter;
 }

void fpu_check()
{
// http://www.christian-seiler.de/projekte/fpmath/
    #ifdef __WIN32__
    //#include <float.h>
    uint32_t control_wordfp=0, status_wordfp=0;

    //read out
    control_wordfp = _controlfp(0, 0);//,&control_wordfp, 0);
    control_wordfp = getFPUStateX86();
    //http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz/
    //DAZ and FTZ in MXCRS
    //FTZ: The FTZ bit (bit 15)& The underflow exception (bit 11)
    //DAZ: DAZ bit (bit 6)
    //control_wordfp = _controlfp( _CW_DEFAULT, 0xfffff);
    //control_wordfp =_controlfp(_DN_FLUSH|32 ,_MCW_DN|32 ); //_controlfp(control_wordfp | 34816, 4294967295);
    //control_wordfp =_controlfp(_PC_64, _MCW_PC); //_PC_24, _PC_53 _PC_64
    //control_wordfp =_controlfp(_RC_NEAR, _MCW_RC); //_RC_UP _RC_CHOP _RC_DOWN _RC_NEAR

    status_wordfp = _statusfp();
    _clearfp();
    #else
    #include <fpu_control.h>
    fpu_control_t control_wordfp=0, fpu_cw=0;
    _FPU_GETCW(control_wordfp);

    //Status WORD under linux????

    #endif


    printf("FPU/SSE Control-register(0x%.4x): ", control_wordfp );
    //RC round-control
     switch ((control_wordfp & __FPU_CW_ROUND_MASK__)>>0)//& 3072)
    {case __FPU_CW_ROUND_NEAR__:
         printf("nearest rounding");
         break;
     case __FPU_CW_ROUND_UP__:
          printf("ceil-rounding");
         break;
     case __FPU_CW_ROUND_DOWN__:
          printf("floor-rounding");
         break;
     case __FPU_CW_ROUND_CHOP__:
          printf("truncation rounding");
         break;
    }
    // printf("%x %x %x",control_wordfp & _MCW_PC, _MCW_PC,_MCW_RC);

    //PC Precision-control
    switch ((control_wordfp & __FPU_CW_PREC_MASK__ )>>0 )//& 768)
    {case __FPU_CW_PREC_SINGLE__:
         printf(", internal precision float (32bit)\n");
         break;
     case __FPU_CW_PREC_DOUBLE__:
         printf(", internal precision double (64bit)\n");
         break;
     case __FPU_CW_PREC_EXTENDED__:
         printf(", internal precision extended double (80bit)\n");
         break;
     default :
         printf(", internal precision invalid \n");
    }

#ifdef __WIN32__
    //FPU/SSE status
    printf("FPU/SSE status (0x%.4x): ", status_wordfp);
    if ((status_wordfp & 32)>>5)  printf("inexact ");
    if ((status_wordfp & 16)>>4)  printf("underflow ");
    if ((status_wordfp &  8)>>3)  printf("overflow ");
    if ((status_wordfp &  4)>>2)  printf("division-by-zero ");
    if ((status_wordfp &  2)>>1)  printf("denormal ");
    if ((status_wordfp &  1)>>0)  printf("invalid operation mask ");
    if ( status_wordfp == 0 )  printf("OK");
    printf("\n");
    #endif
}



// void as2v_addsig2vol_3(cArrayDouble* AScan_realz, cArrayDouble* AScan_complexz,
//             cArrayFloat* pix_vectz, cArrayFloat* rec_posz, cArrayFloat* send_posz, cArrayFloat* speedz, float* resz, float* timeintz,
// 		    cArrayDouble* IMAGE_SUM_realz, cArrayDouble* IMAGE_SUM_complexz, cArrayDouble* outputImage, cArrayDouble* buffer)
// {
//     //Check if minimal data is available...
//     if (AScan_realz && IMAGE_SUM_realz && pix_vectz && rec_posz && send_posz && speedz && resz && timeintz){
//
//         double* AScan_complexz_ptr = NULL;
//         double* IMAGE_SUM_complexz_ptr = NULL;
//
//         if ( (AScan_complexz==!NULL) & (IMAGE_SUM_complexz==!NULL) ){
//             AScan_complexz_ptr = AScan_complexz->data;
//             IMAGE_SUM_complexz_ptr = IMAGE_SUM_complexz->data;
//         } else if( (AScan_complexz==NULL) ^ (IMAGE_SUM_complexz==NULL) ){
//             printf("Error: Mismatch complex/real AScan vs complex/real sum image, break\n");
//         }
//
//         // Save pointers locally: Do NOT do pointer arithmetic in the structs themselves (used to reference/free data)
//         double* AScan_realz_ptr = AScan_realz->data;
//         double* IMAGE_SUM_realz_ptr = IMAGE_SUM_realz->data;
//         float* pix_vectz_ptr = pix_vectz->data;
//         float *rec_vec_ptr=rec_posz->data;
//         float *send_vec_ptr=send_posz->data;
//         float *speedz_vec_ptr=speedz->data;
//
//         unsigned int i;
//         unsigned int n_AScan = AScan_realz->x;
//         unsigned int n_AScan_block = AScan_realz->y;
//         unsigned int n_Z = IMAGE_SUM_realz->z;
//         unsigned int n_Y = IMAGE_SUM_realz->y;
//         unsigned int n_X = IMAGE_SUM_realz->x;
//         unsigned int n_IMAGE = IMAGE_SUM_realz->len;
//         unsigned int n_rec_vec_block = rec_posz->y;
//         unsigned int n_send_vec_block= send_posz->y;
//         unsigned int n_speed_vec_block = speedz->y;
//
//         #ifdef addsig2vol_debug
//         printf("AScan_real:\t"); caMiniprint(AScan_realz);
//         printf("IMAGE_SUM_real:\t"); caMiniprint(IMAGE_SUM_realz);
//         printf("send_pos:\t"); caMiniprint(send_posz);
//         printf("rec_pos:\t"); caMiniprint(rec_posz);
//         printf("speed:\t\t"); caMiniprint(speedz);
//         printf("n_AScan: %i\n", n_AScan);
//         printf("n_AScan_block: %i\n", n_AScan_block);
//         printf("n_Z: %i\n", n_Z);
//         printf("n_Y: %i\n", n_Y);
//         printf("n_X: %i\n", n_X);
//         printf("n_IMAGE: %i\n", n_IMAGE);
//         printf("res: %f\n", *resz);
//         printf("timeint: %f\n", *timeintz);
//         printf("INTERP_RATIO: %i\n", INTERP_RATIO);
//         printf("NUMCORES: %i\n", NUMCORES);
//         printf("MIN_VOXEL: %i\n", MIN_VOXEL);
//         fpu_check();
//         #endif
//
//         //check if X_Dim >= MIN_VOXEL
//         if (n_X < MIN_VOXEL)
//         {
//             printf("Error: X-Dim has to be at least %d (X is %u, Y is %u, Z %u). Use Y-dim for such small sizes.\n", MIN_VOXEL, n_X, n_Y, n_Z);
//         }
//
//         ////benchmark of performance, selecting number of threads
//         if (nCores_bench == -1) as2v_benchLocal();
//
//         //select if use or not use multithreading
//         #ifdef addsig2vol_debug
//         printf("selectedCore perf: %f, %f\n",((double)n_IMAGE/((double)throughput[nCores_bench-1] * 1000000))+((double)latency[nCores_bench-1]/1000000000),(((double)n_IMAGE/((double)throughput[0]*1000000))+((double)latency[0]/1000000000)));
//         printf("n_IMAGE: %i, throughput: %f, latency: %f, %f, %f\n",n_IMAGE, (double)throughput[nCores_bench-1],(double)latency[nCores_bench-1],(double)throughput[0],(double)latency[0]);
//         #endif
//
//         if ( ( ((double)n_IMAGE/((double)throughput[nCores_bench-1] *1000000))+((double)latency[nCores_bench-1]/1000000000)) <= (((double)n_IMAGE/((double)throughput[0] *1000000))+((double)latency[0]/1000000000)))
//         {
//             nCores = nCores_bench;
//         }
//         else
//         {
//             #ifdef addsig2vol_debug
//             printf("Overhead to big, switch to single thread.\n");
//             #endif
//             nCores = 1;
//         }
//
//         #ifdef addsig2vol_debug
//         printf("selectedNumCores: %i\n", nCores);
//         printf("savedNumCORE: %i\n", nCores_bench);
//         printf("perf_MT: %e\n", ( (throughput[nCores_bench] * n_IMAGE)+latency[nCores_bench]));
//         printf("perf_single: %e\n", ((throughput[1]*n_IMAGE)+latency[1]));
//         #endif
//
//
//
//         //soundmap version ?
//         if (((speedz->x == n_X) & (speedz->y == n_Y) ) & ((speedz->z == n_Z) | (speedz->z == 1)))
//         {
//             printf("Info: Soundmap version\n");
//             addsig2vol_mode = 2;
//         }
//         else
//         {
//             //not soundmapversion; blocked version?
//             if ( (rec_posz->x != 3) |  (send_posz->x != 3) ) //check first dimension
//             {
//                 printf("send_pos:\t"); caMiniprint(send_posz);
//                 printf("rec_pos:\t"); caMiniprint(rec_posz);
//                 printf("ascan_block: %d, blocksize rec_vec: %d, blocksize send_vec: %d",n_AScan_block,n_rec_vec_block,n_send_vec_block);
//
//                 printf("Error: 3-d vectors needed for emitter & receiver positions or transposed blocked pos (1x3 instead of 3x1), break\n");
//             }
//
//             if (!( (((n_AScan_block == n_rec_vec_block) & (n_AScan_block == n_send_vec_block)) | ((1 == n_rec_vec_block) & (n_AScan_block == n_send_vec_block)) | ((n_AScan_block == n_rec_vec_block) & (1 == n_send_vec_block))) & ((n_AScan_block == n_speed_vec_block) | (1 == n_speed_vec_block)) ))
//             {
//                 printf("Error: Blocked sizes parameter mismatch. Size(AScan,2) has to be size(rec_vec,2) and/or size(send_vec,2), n_rec_vec_block:%i n_send_vec_block:%i n_speed_vec_block:%i break\n",n_rec_vec_block, n_send_vec_block, n_speed_vec_block);
//             }
//             if (n_AScan_block == 1)
//             {
//                 //if here blocked mode successful initalized!
//                 addsig2vol_mode = 1;
//             }
//             else{
//                 //normal version
//                 addsig2vol_mode = 0;
//             }
//         }
//
//         // malloc space for results
//         if(outputImage){
//             out_real = outputImage->data;
//             //cAout_real = outputImage;
//         }else{
//             out_real = malloc(sizeof(double)*n_IMAGE);
//         }
//
//         if(buffer){
//             buffer_real = buffer->data;
//             //cAbuffer_real = buffer;
//         }else{
//             buffer_real = malloc(sizeof(double)*n_AScan*INTERP_RATIO);
//         }
//
//         if (AScan_complexz) {
//             //we know that both AScan_complexz and IMAGE_SUM_complexz exist
//             // TODO same as in REAL case
//             // out_complex = malloc(sizeof(double)*n_IMAGE);
//             // buffer_complex = malloc(sizeof(double)*n_AScan*INTERP_RATIO);
//             // cAout_complex = caNewArrayDoubleFromData(out_complex, n_X, n_Y, n_Z);
//             // cAbuffer_complex = caNewArrayDoubleFromData(buffer_complex, n_AScan*INTERP_RATIO, 1, 1);
//         }
//
//         #ifdef addsig2vol_debug
//         printf("Set addsig2vol_mode to %i\n", addsig2vol_mode);
//         printf("Elements in out_real: %i\n", n_IMAGE);
//         printf("Elements in out_complex: %i\n", n_IMAGE);
//         printf("Elements in buffer_real: %i\n", n_AScan*INTERP_RATIO);
//         printf("Elements in buffer_complex: %i\n", n_AScan*INTERP_RATIO);
//         #endif
//
//
//         totalAscan = n_AScan_block;
//         currentAscan = 0;
//
//         threadstats_clearAll();
//         tssettimer();
//
//         ////first Ascan
//         // combined REAL & COMPLEX VERSION
//         as2v_MT(out_real, AScan_realz_ptr, n_AScan, buffer_real,
//             pix_vectz_ptr, n_X, rec_vec_ptr, send_vec_ptr,
//             speedz_vec_ptr, resz, timeintz, AScan_complexz_ptr, buffer_complex, out_complex, n_Y, n_Z, IMAGE_SUM_realz_ptr, IMAGE_SUM_complexz_ptr);
//
//
//         //loop over ascans > 1
//         for (i = 2; i <= n_AScan_block; i++) {
//             #ifdef addsig2vol_debug
//             printf("as2v_MT call %i:\n", i);
//             #endif
//             currentAscan++;
//
//             //check for complex ascan only increase if available because NULL-Pointer +something -> not anymore a nullpointer!
//             if (AScan_complexz != NULL)	AScan_complexz = AScan_complexz + (n_AScan * (i - 1)); //set to next value
//             if (1<n_rec_vec_block)  rec_vec_ptr  = rec_posz->data  + (3 * (i - 1)); else rec_vec_ptr = rec_posz->data;
//             if (1<n_send_vec_block) send_vec_ptr = send_posz->data + (3 * (i - 1)); else send_vec_ptr = send_posz->data;
//             if (1<n_speed_vec_block) speedz_vec_ptr = speedz->data + (1  * (i - 1)); else speedz_vec_ptr = speedz->data;
//
//
//             // NOTE: pass out_real/out_complex as IMAGE_SUM to sum up the single scans
//             as2v_MT(out_real, AScan_realz_ptr + (n_AScan * (i - 1)), n_AScan, buffer_real,
//             pix_vectz_ptr, n_X, rec_vec_ptr, send_vec_ptr,
//             speedz_vec_ptr, resz, timeintz, AScan_complexz_ptr, buffer_complex, out_complex, n_Y, n_Z, out_real, out_complex);
//
//             //             mxSetPr(ascan_output_buffer,mxGetPr(AScan)+(n_AScan*(i-1)));
//             //             mxSetPr(recpos_output_buffer,(char*)mxGetPr(recpos_org)+(3*sizeof(float)*(i-1)));
//             //
//             //             mexCallMATLAB(0, NULL, 0, NULL, "figure");
//             //             mexCallMATLAB(0, NULL, 1, &ascan_output_buffer, "plot");
//             //             mexCallMATLAB(0, NULL, 1, &recpos_output_buffer, "disp");
//
//         }
//         //
//         // if (segmentedAxis == ZAXIS){
//         //     printf("ZAXIS\n");
//         // }
//         // if (segmentedAxis == YAXIS){
//         //     printf("YAXIS\n");
//         // }
//         // if (segmentedAxis == XAXIS){
//         //     printf("XAXIS\n");
//         // }
//         // printf("img size x,y,z: %i %i %i\n", imgX, imgY, imgZ);
//         // printf("pos x,y,z: %i %i %i\n", posX, posY, posZ);
//         // printf("halfstep x,y,z: %i %i %i\n", halfStepX, halfStepY, halfStepZ);
//
//         totalAscan = 0;
//         currentAscan = 0;
//             tsclearAll();
//
//     }
// ;
//     //Daten fehlen
// }

void printIntro(void){
    printf("\naddSig2Vol_2 SSE1 Assembler Optimized 64bit LINUX&Windows v3.1 (Multiple Rec-AScan Vers.)\n\n Calculate the ellip. backprojection.\nUses SSE. Features: Win&Linux, 32&64bit version, SSE1-Assembler&C-Implementation, Multithreaded by PosixThreads (under windows pthreadsVC2.dll might be needed)\n\t %s. M.Zapf KIT-IPE\n\n",__DATE__);
    // TODO print som flags for debugging
    #ifdef addsig2vol_debug
        printf("addsig2vol_debug: 1\n");
        #ifdef C_CODE
            printf("C_CODE: 1, use c implementation\n");
        #else
            printf("C_CODE: 0, use asm implementation\n");
        #endif

    #endif
}

void as2v_overwriteBenchresultToThreadcount_n(uint32_t n){
    if ((n <= NUMCORES) & (n >= 1)){
        printf("Clear benchmark results, manual set number of threads to %i\n", n);
        nCores_bench = n;
        //fill
        for (int i=NUMCORES; i>0; i--){
            throughput[i-1]=1;
            latency[i-1]=10000000;
        }
        throughput[nCores_bench-1]=299;
        latency[nCores_bench-1]=1;
    } else {
        printf("Number of threads %i\n is invalid, benchmark results remain unchanged\n", n);
    }
}

unsigned int as2v_getThreadsNumber(){
    return nCores;
}
