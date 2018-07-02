#include "addsig2vol_3.h"

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
    #include <stdlib.h>
#endif


#include <math.h>

#ifndef _PORTABLEFPU_H_INCLUDED
#include "portableFPU.h"
#endif

#ifndef _PSTDINT_H_INCLUDED
#include "pstdint_new.h"
#endif

#ifndef PRINT(string, ...)
    #define PRINT(string, ...) printf(string, ##__VA_ARGS__)
#endif



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



#define addsig2vol_debug
#undef addsig2vol_debug


#define p_threads






#ifdef p_threads
#include "pthread.h"
#define NUMCORES 16 
#else
#define NUMCORES 1
#endif


#define INTERP_RATIO 5    


#define MIN_VOXEL 4



unsigned int addsig2vol_mode=0; 



static int64_t latency[NUMCORES]={0}; 
static int64_t throughput[NUMCORES]={0}; 
static uint32_t nCores_bench = -1;
static uint32_t nCores = NUMCORES; 

static double* out_real = NULL;
static double* out_complex = NULL;
static double* buffer_real = NULL;
static double* buffer_complex = NULL;


















typedef   struct   /*struct reordered to circumvent alignment problems, 64pointer then 32bit values*/
{
    double *outz;
    double *AScanz;
    double *out_complexz;
    double *bufferz;
    float  *pix_vectz;
    double *buffer_complexz;
    float  *rec_posz;
    float  *send_posz;
    float  *speedz;
    float  *resz;
    float  *timeintz;
    double *AScan_complexz;
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




  /* typedef   struct   
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


uint64_t CPUCount(void);
uint64_t TimeCounter(void);


void fpu_check(void);


void as2v_complex(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void as2v_complex_sm(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void as2v_c(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void as2v_MT(double*outz, double*AScanz, unsigned int n_AScanz, double*bufferz, float*pix_vectz,
		    unsigned int n_Xz, float*rec_posz, float*send_posz, float*speedz, float*resz,
		    float*timeintz, double*AScan_complexz,
		    double*buffer_complexz, double*out_complexz, unsigned int n_Yz, unsigned int n_Zz,
		    double* IMAGE_SUMz, double* IMAGE_SUM_complexz);


void xsum_complex(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void xsum_c(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);


void *thread_function(void *arg);


void as2v_bench( uint64_t lat[],  uint64_t through[]);






void as2v_MT(double*outz, double*AScanz, unsigned int n_AScanz, double*bufferz, float*pix_vectz,
		    unsigned int n_Xz, float*rec_posz, float*send_posz, float*speedz, float*resz,
		    float*timeintz, double*AScan_complexz,
		    double*buffer_complexz, double*out_complexz, unsigned int n_Yz, unsigned int n_Zz,
		    double *IMAGE_SUMz, double *IMAGE_SUM_complexz)

{

  
#ifdef p_threads
    pthread_t mythread[NUMCORES]; 
    int rc = 0; 
#endif
    Addsig2vol_param threadArg[NUMCORES]; 
    float pix_vecz_buffer[NUMCORES][3]= {0,0,0};
    unsigned int n_Zz_start = 0;
    unsigned int n_Zz_num = 0;
    
    unsigned int i = 0;


    
    if (n_Zz>1)
      {
      
      if (n_Zz < nCores)  nCores = n_Zz;

      #ifdef addsig2vol_debug
      
      #endif
      int z = 0;

      
		  for (i=0;i<=nCores-1;i++) 
        {
        n_Zz_num  = (unsigned int)(floor(n_Zz/nCores));
        n_Zz_start = n_Zz_num*i*n_Xz*n_Yz;

         
        pix_vecz_buffer[i][0]= *pix_vectz;
        pix_vecz_buffer[i][1]= *(pix_vectz+1);
        pix_vecz_buffer[i][2]= (*(pix_vectz+2))+n_Zz_num*i* *resz;

        if (i==nCores-1)
        {      
        n_Zz_num = (unsigned int)(n_Zz-floor(n_Zz/nCores)*(nCores-1));}

        /* PRINT("delta_address :  %p\n\n", (pix_vectz));
        PRINT("z :  %f\n\n", *(pix_vectz+2));
        PRINT("z_num :  %i\n\n", n_Zz_num);
        PRINT("z_start :  %i\n\n", n_Zz_start);
        PRINT("x :  %f\n\n", pix_vecz_buffer[i][0]);
        PRINT("y :  %f\n\n", pix_vecz_buffer[i][1]);
        PRINT("z :  %f\n\n", pix_vecz_buffer[i][2]);
        PRINT("p_x :  %p\n\n", &(pix_vecz_buffer[i][0]));*/

        z+= n_Zz_num;
        
        threadArg[i].outz=outz+n_Zz_start;
        threadArg[i].AScanz=AScanz;
        threadArg[i].n_AScanz=n_AScanz;
        threadArg[i].bufferz=bufferz;
        threadArg[i].pix_vectz=&(pix_vecz_buffer[i][0]);/*mxGetPr(pix_vect)*/
        threadArg[i].n_Xz=n_Xz;
        threadArg[i].rec_posz=rec_posz;
        threadArg[i].send_posz=send_posz;
        threadArg[i].speedz=speedz;
        threadArg[i].resz=resz;
        threadArg[i].timeintz=timeintz;
        threadArg[i].AScan_complexz=AScan_complexz;
        threadArg[i].buffer_complexz=buffer_complexz;
        threadArg[i].out_complexz=out_complexz+n_Zz_start;
        threadArg[i].n_Yz=n_Yz;
        threadArg[i].n_Zz=n_Zz_num;/*n_Zz*/
        threadArg[i].IMAGE_SUMz=IMAGE_SUMz+n_Zz_start;
        threadArg[i].IMAGE_SUM_complexz=IMAGE_SUM_complexz+n_Zz_start; 
        }

        
        
        
        
        
        
        


      }
      else  {

     
      if (n_Yz>1)
      {
      
      if (n_Yz < nCores)  nCores = n_Yz;

      #ifdef addsig2vol_debug
      
      #endif

      
	  	for (i=0;i<=nCores-1;i++) 
        {
        n_Zz_num  = (unsigned int)(floor(n_Yz/nCores));
        n_Zz_start = n_Zz_num*i*n_Xz;

         
        pix_vecz_buffer[i][0]= *pix_vectz;
        pix_vecz_buffer[i][1]= (*(pix_vectz+1))+n_Zz_num*i* *resz;
        pix_vecz_buffer[i][2]= *(pix_vectz+2);

        if (i==nCores-1)
        {      
        n_Zz_num = (unsigned int)(n_Yz-floor(n_Yz/nCores)*(nCores-1));}

        
        threadArg[i].outz=outz+n_Zz_start;
        threadArg[i].AScanz=AScanz;
        threadArg[i].n_AScanz=n_AScanz;
        threadArg[i].bufferz=bufferz;
        threadArg[i].pix_vectz=&(pix_vecz_buffer[i][0]);/*mxGetPr(pix_vect)*/
        threadArg[i].n_Xz=n_Xz;
        threadArg[i].rec_posz=rec_posz;
        threadArg[i].send_posz=send_posz;
        threadArg[i].speedz=speedz;
        threadArg[i].resz=resz;
        threadArg[i].timeintz=timeintz;
        threadArg[i].AScan_complexz=AScan_complexz;
        threadArg[i].buffer_complexz=buffer_complexz;
        threadArg[i].out_complexz=out_complexz+n_Zz_start;
        threadArg[i].n_Yz=n_Zz_num; /*n_Yz*/
        threadArg[i].n_Zz=n_Zz;
        threadArg[i].IMAGE_SUMz=IMAGE_SUMz+n_Zz_start;
        threadArg[i].IMAGE_SUM_complexz=IMAGE_SUM_complexz+n_Zz_start; 
        }

      }else
        { 
        nCores = 1;

        #ifdef addsig2vol_debug
        
        #endif

        
        threadArg[nCores-1].outz=outz;
        threadArg[nCores-1].AScanz=AScanz;
        threadArg[nCores-1].n_AScanz=n_AScanz;
        threadArg[nCores-1].bufferz=bufferz;
        threadArg[nCores-1].pix_vectz=pix_vectz;
        threadArg[nCores-1].n_Xz=n_Xz;
        threadArg[nCores-1].rec_posz=rec_posz;
        threadArg[nCores-1].send_posz=send_posz;
        threadArg[nCores-1].speedz=speedz;
        threadArg[nCores-1].resz=resz;
        threadArg[nCores-1].timeintz=timeintz;
        threadArg[nCores-1].AScan_complexz=AScan_complexz;
        threadArg[nCores-1].buffer_complexz=buffer_complexz;
        threadArg[nCores-1].out_complexz=out_complexz;
        threadArg[nCores-1].n_Yz=n_Yz;
        threadArg[nCores-1].n_Zz=n_Zz;
        threadArg[nCores-1].IMAGE_SUMz=IMAGE_SUMz;
        threadArg[nCores-1].IMAGE_SUM_complexz=IMAGE_SUM_complexz; 
        }
      }

      

      
      #ifdef C_CODE
      xsum_c(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);
      #else
      xsum_complex(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);
      #endif

      
      for (i=0;i<nCores-1;i++)
	  {
      #ifdef p_threads
      rc = pthread_create( &(mythread[i]), NULL, thread_function, (void*)&(threadArg[i]));
      if (rc) { PRINT("ERROR: return code from pthread_create() is %d\n", rc); return;}
      #endif
      }

	
      #ifdef C_CODE
      as2v_c(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);
      #else
      if (addsig2vol_mode==0) as2v_complex(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);
      if (addsig2vol_mode==2) as2v_complex_sm(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);

      #endif
			
			for (i=0;i<nCores-1;i++)
		  {
              #ifdef p_threads
              rc = pthread_join ( mythread[i], NULL );
              if (rc) { PRINT("ERROR: return code from pthread_join() is %d\n", rc); return;}
              #endif
        }

      
      


}



void as2v_c(Addsig2vol_param* tt, Addsig2vol_param* t1, Addsig2vol_param* t2, Addsig2vol_param* t3)
{

      

float dist_sv[3] = {0,0,0};
float dist_rv[3] = {0,0,0};
float factor = 0;

unsigned int index, image_index = 0;
unsigned int z,y,x,sampl = 0;


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



factor =  INTERP_RATIO / (*speedz * *timeintz);
















for (z=1; z<=n_Zz; z++)
	{
		
        dist_sv[2]  = send_posz[2] - (((float)z* *resz)+pix_vectz[2]);
        dist_sv[2] *= dist_sv[2];
		
        dist_rv[2]  = rec_posz[2] - (((float)z* *resz)+pix_vectz[2]);
        dist_rv[2] *= dist_rv[2];

	for (y=1; y<=n_Yz; y++)
		{
		
        dist_sv[1]  = send_posz[1] - (((float)y* *resz)+pix_vectz[1]);
        dist_sv[1] *= dist_sv[1];
		
        dist_rv[1]  = rec_posz[1] - (((float)y* *resz)+pix_vectz[1]);
        dist_rv[1] *= dist_rv[1];

		for (x=1; x<=n_Xz; x++)
			{
				
				
				dist_sv[0] = send_posz[0] - ((x* *resz)+pix_vectz[0]);
				dist_rv[0] =  rec_posz[0] - ((x* *resz)+pix_vectz[0]);

				    
				index = (unsigned int) floor( ( sqrt(dist_sv[1]+ dist_sv[2] + (dist_sv[0]*dist_sv[0])) + sqrt(dist_rv[1]+ dist_rv[2]+ (dist_rv[0]*dist_rv[0])) ) * factor);
				    
				    
                    

				if ((index >= n_AScanz*INTERP_RATIO) | (index < 0)){
					outz[image_index] = IMAGE_SUMz[image_index]; 
                }
				else
					outz[image_index] = IMAGE_SUMz[image_index] + bufferz[index];

				image_index++; 
			}
		}
	}

}



void xsum_c(Addsig2vol_param* tt, Addsig2vol_param* t1, Addsig2vol_param* t2, Addsig2vol_param* t3)
{  

float dist_sv[3] = {0,0,0};
float dist_rv[3] = {0,0,0};
float factor = 0;

unsigned int image_index = 0;
unsigned int i,j,sampl = 0;

double i_buffer=0.0;

double* sec_buffer;


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
{	bufferz[i] = i; 
}
#endif


sec_buffer = malloc(INTERP_RATIO * n_AScanz * sizeof(double));



for (i=0;i<(unsigned int) floor(INTERP_RATIO/2);i++)
{	sec_buffer[i] = AScanz[0] * (unsigned int) ((i+1)/floor(INTERP_RATIO/2));
}


for (i=0;i<n_AScanz-1;i++)
{    	for (j=0;j<INTERP_RATIO;j++)
	{	sec_buffer[(unsigned int) floor(INTERP_RATIO/2)+(i* INTERP_RATIO)+j] = AScanz[i+1] * j/INTERP_RATIO + AScanz[i] * (INTERP_RATIO-j)/INTERP_RATIO; 
	}
}


for (i=0;i<(unsigned int) floor(INTERP_RATIO/2)+1;i++)
{	sec_buffer[n_AScanz*INTERP_RATIO-(unsigned int) floor(INTERP_RATIO/2)-1+i] = AScanz[n_AScanz-1] * ((floor(INTERP_RATIO/2)+1-i) / (floor(INTERP_RATIO/2)+1));
}





sampl = (unsigned int)(ceil((float)1.7*(( *resz / *speedz)/ (*timeintz/INTERP_RATIO)) /2)); 

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



free(sec_buffer);

}



 void *thread_function(void *arg)
{
    Addsig2vol_param* tt = (Addsig2vol_param*) arg;

    
    #ifdef C_CODE
    as2v_c(arg,arg,arg,arg); 
    #else
    if (addsig2vol_mode==0) as2v_complex(arg,arg,arg,arg);
    if (addsig2vol_mode==2) as2v_complex_sm(arg,arg,arg,arg);

    #endif

   
   /* as2v_c(tt->outz, tt->AScanz, tt->n_AScanz, tt->bufferz, tt->pix_vectz,
		    tt->n_Xz, tt->rec_posz, tt->send_posz, tt->speedz, tt->resz,
		    tt->timeintz, tt->AScan_complexz,
		    tt->buffer_complexz, tt->out_complexz, tt->n_Yz, tt->n_Zz,
		    tt->IMAGE_SUMz, tt->IMAGE_SUM_complexz);*/

     return NULL;
}



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
	
	
	
	float pix_vec_bench[3]  = {2,1,77};
	
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
    uint64_t bench_ref_time=0; 

    uint64_t minBenchTime = 100000000; 

    uint32_t minAverage = 7;
    uint32_t average = 8;
    uint32_t nVoxel_throughput = 128;

    double *AScan_bench;
	double *buffer_bench;
	double *out_bench;
	double *image_sum_bench;
	double *time_bench;

    
   /*  do { nVoxel_throughput=(uint32_t) nVoxel_throughput/2;
       pr = mxMalloc(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput * sizeof(double));
     if (pr!=NULL)  mxFree(pr);
       nVoxel_throughput=(uint32_t) nVoxel_throughput/2; 
   } while (pr==NULL);*/

  PRINT("Benchmarking System:\nLatency with %d Byte Image\nThroughput with %d kByte Image.\n",(int)(MIN_VOXEL*NUMCORES*sizeof(double)),(int)(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput*sizeof(double))/1024);
  PRINT("Threads | Throughput in MVoxel/s                     | Latency in nsec                              | Median-Size\n");
  PRINT("        | Median | Mean   | Min    | Max    | Std    | Median    | Min       | Max       | Std      | \n");
  PRINT("-----------------------------------------------------------------------------------------------------------\n");

  
  if (nVoxel_throughput<NUMCORES*MIN_VOXEL) nVoxel_throughput=NUMCORES*MIN_VOXEL;
  
   

	
    time_bench = malloc(sizeof(double));
	n_AScan        = 3000;        
	n_AScan_block  = 1;    

    

    PRINT("buffer size %i\n", (n_AScan*INTERP_RATIO));
    buffer_bench = malloc(n_AScan*INTERP_RATIO *sizeof(double));

    

	AScan_bench = malloc(n_AScan *n_AScan_block*sizeof(double));




    

	n_X             = nVoxel_throughput;
	n_Y             = nVoxel_throughput;
    n_Z             = nVoxel_throughput;

    
		n_IMAGE         = n_X * n_Y * n_Z;

	out_bench = malloc(n_IMAGE*sizeof(double));
    image_sum_bench = malloc(n_IMAGE*sizeof(double));

        PRINT("Voxels in x, y, z: [%i, %i, %i]\n", n_X, n_Y, n_Z );
        PRINT("Voxels total: %i\n", n_IMAGE);
        PRINT("Buffer total: %i\n", INTERP_RATIO * n_AScan);

       
      
	  nCores = 1;
      do{
      counter = TimeCounter();
      as2v_MT((out_bench), (AScan_bench), n_AScan, (buffer_bench),
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL, NULL, NULL, n_Y, n_Z, (image_sum_bench), NULL);
       counter2 = TimeCounter();
       } while (counter2<counter); 

        
        bench_ref_time= counter2-counter;
        
        average = (uint32_t) ((minBenchTime/bench_ref_time)+1); 

        if (average < minAverage) average = minAverage;
        if (average > MAX_AVERAGE) average = MAX_AVERAGE;


	for (i=1;i<=NUMCORES;i++)
	{
	  
		nCores = i;

	    
        

		for (j=0;j<average;j++)
         {
         do{ counter=TimeCounter();
          
		  as2v_MT((out_bench), (AScan_bench), n_AScan, (buffer_bench),
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL, NULL, NULL, n_Y, n_Z, (image_sum_bench), NULL);
          /*as2v_MT(pr3, pr2, n_AScan, pr1,
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL,NULL, NULL, n_Y, n_Z, pr4,NULL);*/
		  
          
          counter2=TimeCounter();
         } while(counter2<counter); 
          average_buffer[j]= counter2-counter;
        }

        
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
          
          
         #endif

        
        bench_ref_time = (uint64_t) average_buffer[0];
        
        bench_ref_time = (uint64_t) average_buffer[(uint32_t) ceil(average/2)];

        
        throughput[i-1] =  ((1000*(uint64_t)n_IMAGE)/ bench_ref_time); 

        PRINT("%7i |%7llu |%7llu |%7llu |%7llu |%7llu", i, (int64_t)throughput[i-1],((1000*(int64_t)n_IMAGE)/(int64_t)mittelw), (1000*(int64_t)n_IMAGE)/((int64_t)average_buffer[average-1]), (1000*(int64_t)n_IMAGE)/((int64_t)average_buffer[0]), (1000*(int64_t)n_IMAGE)/((int64_t)stdabw));


       
       average = (uint32_t) ((minBenchTime/bench_ref_time)+1); 
       if (average < minAverage) average = minAverage;
       if (average > MAX_AVERAGE) average = MAX_AVERAGE;

		 
      










	   
	

		for (j=0;j<minAverage;j++)
         {
          do{
          counter=TimeCounter();
			
		  as2v_MT((out_bench), (AScan_bench), n_AScan, (buffer_bench),
			     &pix_vec_bench[0],(unsigned int) MIN_VOXEL, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL, NULL, NULL, (unsigned int) 1,(unsigned int) NUMCORES, (image_sum_bench), NULL);

		
        
		 counter2 = TimeCounter(); } while(counter2<counter); 
         average_buffer[j]= counter2-counter;
        }
        
		 for (k=minAverage-1;k>0;k--)
         {  for (l=minAverage-1;l>0;l--){
               if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;}
            }
         }
         #ifdef addsig2vol_debug
            
        #endif

         mittelw=0;
         for  (k=0;k<minAverage;k++)  {mittelw=average_buffer[k]+mittelw;}
         mittelw=(uint64_t) mittelw/minAverage;

         stdabw=0;
         for  (k=0;k<minAverage;k++)
         {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);}
         stdabw = (uint64_t) sqrt((int64_t) stdabw/minAverage);

        
        latency[i-1] = (uint64_t) average_buffer[0];	 
        
        latency[i-1] = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
	    PRINT("|%10llu |%10llu |%10llu |%10llu |%8i\n", latency[i-1],average_buffer[0],average_buffer[minAverage-1],stdabw,average);
  		}

		for (i=NUMCORES;i>0;i--)
        { throughput_sort_buffer[i-1]=throughput[i-1];
		}

      
      for (i=NUMCORES-1;i>=1;i--)
      { if (throughput_sort_buffer[i]>(throughput_sort_buffer[i-1])) throughput_sort_buffer[i-1] = throughput_sort_buffer[i]; }

      
      for (i=0;i<NUMCORES;i++)
      { if (throughput[i] == throughput_sort_buffer[0]) break; } 

      
      nCores = i+1;
      
      nCores_bench = i+1;

      switch (i+1)
      {
      case 1:
          PRINT("Detected Single-core System, 1 thread prefered\n");
          break;
      case 2:
          PRINT("Detected Dual-core or Hyperthreading System, 2 threads prefered\n");
          break;
      case 3:
          PRINT("Detected Triple-core or Hyperthreading System, 3 threads prefered\n");
          break;
      case 4:
          PRINT("Detected Quad-Core system, 4 threads prefered\n");
          break;
      case 8:
          PRINT("Detected Octa-Core system (or Quadcore with HT), 8 threads prefered\n");
          break;
      case 16:
          PRINT("Detected 16 Core system, 16 threads prefered\n");
          break;
      default:
          PRINT("Detected %i-core system, %i threads prefered (?)\n",i+1,i+1);
      }


      
      PRINT("\nPerformance for various imagesize in Voxel (with potentially %i Cores)\n",nCores);
      PRINT("     Voxel | Throughput in kVoxel/s         | Time in mikros  | Malloc time (mikro-sec)\n");
      PRINT("           | Median   | Mean     | Std      |           | Median | mean  | Std  | min  | max  \n");
      PRINT("------------------------------------------------------------------------------------------------\n");
      PRINT("voxel throughput maximal: %i\n", nVoxel_throughput*nVoxel_throughput*nVoxel_throughput);
      for (i=MIN_VOXEL;i<=(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput);i=i*2)
      {
          

          for (j=0;j<minAverage;j++)
          {

              do {
                  counter=TimeCounter();
                  
                  as2v_MT((out_bench), (AScan_bench), n_AScan, (buffer_bench),
                  &pix_vec_bench[0], (uint32_t) MIN_VOXEL, &rec_vec_bench[0], &send_vec_bench[0],
                  &float_bench,
                  &float_bench, &float_bench, NULL, NULL, NULL, (uint32_t) 1, (uint32_t) (i/MIN_VOXEL), (image_sum_bench), NULL);
                  counter2 = TimeCounter(); } while(counter2<counter); 
              average_buffer[j]= counter2-counter;
          }
          
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
           
          #endif

          
          counter2 = (uint64_t) average_buffer[0];	 
          
          counter2 = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
          
          

          
          
          PRINT("%10i | %8llu | %8llu | %8llu |%8llu ",i,(1000000*(uint64_t)i)/(counter2),(1000000*(uint64_t)i)/(mittelw),((uint64_t)i*1000000)/(stdabw),(uint64_t)(counter2/1000) );

          
          nCores=nCores_bench;

          
          minAverage=100;

              for (j=0;j<minAverage;j++)
              {
                  do {
                      free(image_sum_bench);

                      counter=TimeCounter();

                      image_sum_bench = malloc(n_IMAGE*sizeof(double));

                     counter2 = TimeCounter(); } while(counter2<counter); 
                     average_buffer[j]= counter2-counter;
              }
              
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

             
             counter2 = (uint64_t) average_buffer[0];	 
             
             counter2 = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
             PRINT("| %8llu | %8llu | %8llu |%8llu |%8llu\n", (uint64_t)counter2/1000,(uint64_t)mittelw/1000,(uint64_t)stdabw/1000,(uint64_t) average_buffer[0]/1000,(uint64_t) average_buffer[minAverage-1]/1000);

          }

        
        fpu_check();


        
		  free(out_bench);
		  free(image_sum_bench);
        
		  free(buffer_bench);
		  free(AScan_bench);
		  free(time_bench);
}


uint64_t TimeCounter(void) {
uint64_t counter;


    #ifdef __WIN32__ 
    #include <windows.h>
    
    uint64_t iFreq, iCount;
    QueryPerformanceFrequency((LARGE_INTEGER*)&iFreq);
    QueryPerformanceCounter((LARGE_INTEGER*)&iCount);
    
    counter = (uint64_t) ((1000000000*iCount)/iFreq); 


    
    
    #endif

    #ifdef __linux__

    struct timespec time_str;
    struct timespec time_res;
    clock_gettime(CLOCK_REALTIME, &time_str);
    counter = (uint64_t) time_str.tv_nsec + (uint64_t) (time_str.tv_sec*1000000000);

    /*
    clock_gettime(CLOCK_REALTIME, &time_str); clock_getres(CLOCK_REALTIME, &time_res);
    PRINT("CLOCK_REALTIME clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec);
    clock_gettime(CLOCK_MONOTONIC, &time_str); clock_getres(CLOCK_MONOTONIC, &time_res);
    PRINT("CLOCK_MONOTONIC clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_str); clock_getres(CLOCK_PROCESS_CPUTIME_ID, &time_res);
    PRINT("CLOCK_PROCESS_CPUTIME_ID clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec);
    */
      
    
    #endif
    return counter;
 }

void fpu_check()
{

    #ifdef __WIN32__
    
    uint32_t control_wordfp=0, status_wordfp=0;

    
    control_wordfp = _controlfp(0, 0);
    control_wordfp = getFPUStateX86();
    
    
    
    
    
    
    
    

    status_wordfp = _statusfp();
    _clearfp();
    #else
    #include <fpu_control.h>
    fpu_control_t control_wordfp=0, fpu_cw=0;
    _FPU_GETCW(control_wordfp);

    

    #endif


    PRINT("FPU/SSE Control-register(0x%.4x): ", control_wordfp );
    
     switch ((control_wordfp & __FPU_CW_ROUND_MASK__)>>0)
    {case __FPU_CW_ROUND_NEAR__:
         PRINT("nearest rounding");
         break;
     case __FPU_CW_ROUND_UP__:
          PRINT("ceil-rounding");
         break;
     case __FPU_CW_ROUND_DOWN__:
          PRINT("floor-rounding");
         break;
     case __FPU_CW_ROUND_CHOP__:
          PRINT("truncation rounding");
         break;
    }
    

    
    switch ((control_wordfp & __FPU_CW_PREC_MASK__ )>>0 )
    {case __FPU_CW_PREC_SINGLE__:
         PRINT(", internal precision float (32bit)\n");
         break;
     case __FPU_CW_PREC_DOUBLE__:
         PRINT(", internal precision double (64bit)\n");
         break;
     case __FPU_CW_PREC_EXTENDED__:
         PRINT(", internal precision extended double (80bit)\n");
         break;
     default :
         PRINT(", internal precision invalid \n");
    }

#ifdef __WIN32__
    
    PRINT("FPU/SSE status (0x%.4x): ", status_wordfp);
    if ((status_wordfp & 32)>>5)  PRINT("inexact ");
    if ((status_wordfp & 16)>>4)  PRINT("underflow ");
    if ((status_wordfp &  8)>>3)  PRINT("overflow ");
    if ((status_wordfp &  4)>>2)  PRINT("division-by-zero ");
    if ((status_wordfp &  2)>>1)  PRINT("denormal ");
    if ((status_wordfp &  1)>>0)  PRINT("invalid operation mask ");
    if ( status_wordfp == 0 )  PRINT("OK");
    PRINT("\n");
    #endif
}



void as2v_overwriteBenchresultToThreadcount_n(uint32_t n){
    if ((n <= NUMCORES) & (n >= 1)){
        PRINT("Clear benchmark results, manual set number of threads to %i\n", n);
        nCores_bench = n;
        
        for (int i=NUMCORES; i>0; i--){
            throughput[i-1]=1;
            latency[i-1]=10000000;
        }
        throughput[nCores_bench-1]=299;
        latency[nCores_bench-1]=1;
    } else {
        PRINT("Number of threads %i\n is invalid, benchmark results remain unchanged\n", n);
    }
}

as2v_results as2v_addsig2vol_3(as2v_doubleArray* AScan_realz, as2v_doubleArray* AScan_complexz,
            as2v_floatArray* pix_vectz, as2v_floatArray* rec_posz, as2v_floatArray* send_posz, as2v_floatArray* speedz, float* resz, float* timeintz,
		    as2v_doubleArray* IMAGE_SUM_realz, as2v_doubleArray* IMAGE_SUM_complexz)
{
    as2v_results results;

    
    if (AScan_realz && IMAGE_SUM_realz && pix_vectz && rec_posz && send_posz && speedz && resz && timeintz){

        double* AScan_complexz_ptr = NULL;
        double* IMAGE_SUM_complexz_ptr = NULL;

        if ( (AScan_complexz==!NULL) & (IMAGE_SUM_complexz==!NULL) ){
            AScan_complexz_ptr = AScan_complexz->data;
            IMAGE_SUM_complexz_ptr = IMAGE_SUM_complexz->data;
      } else if( (AScan_complexz==NULL) ^ (IMAGE_SUM_complexz==NULL) ){
          PRINT("Error: Mismatch complex/real AScan vs complex/real sum image, break\n");
        return;
      }

      
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

    
    upstreamDoubleArray(&AScan_realz);


    #ifdef addsig2vol_debug
    PRINT("Elements in AScan_real:  %i\n", AScan_real.len);
    PRINT("n_AScan:  %i\n\n", AScan_real.x );
    PRINT("n_AScan_block :  %i\n\n", AScan_real.y );
    PRINT("n_Z:  %i\n\n", IMAGE_SUM_real.z);
    PRINT("n_Y:  %i\n\n", IMAGE_SUM_real.y);
    PRINT("n_X:  %i\n\n", IMAGE_SUM_real.x);
    PRINT("n_IMAGE:  %i\n\n", IMAGE_SUM_real.len);
    PRINT("rec_posz:  %i\n\n", rec_pos.y);
    PRINT("send_posz:  %i\n\n", send_pos.y);
    PRINT("speedz:  %i, %i, %i\n\n", speedz.x, speedz.y, speedz.z);
    PRINT("res:  %i\n\n", *res);
    PRINT("timeint:  %i\n\n", *timeint);
    fpu_check();
    #endif

	
	if (n_X < MIN_VOXEL)
	{
    	PRINT("Error: X-Dim has to be at least %d (X is %u, Y is %u, Z %u). Use Y-dim for such small sizes.\n", MIN_VOXEL, n_X, n_Y, n_Z);
        return;
	}

    	 
        if (nCores_bench == -1) as2v_benchLocal();

        
         #ifdef addsig2vol_debug
          PRINT("selectedCore perf: %f , %f\n",((double)n_IMAGE/((double)throughput[nCores_bench-1] * 1000000))+((double)latency[nCores_bench-1]/1000000000),(((double)n_IMAGE/((double)throughput[0]*1000000))+((double)latency[0]/1000000000)));
          PRINT("nimage %i, throughput: %f , latency %f, %f, %f\n",n_IMAGE, (double)throughput[nCores_bench-1],(double)latency[nCores_bench-1],(double)throughput[0],(double)latency[0]);
         #endif

        if ( ( ((double)n_IMAGE/((double)throughput[nCores_bench-1] *1000000))+((double)latency[nCores_bench-1]/1000000000)) <= (((double)n_IMAGE/((double)throughput[0] *1000000))+((double)latency[0]/1000000000)))
             {  nCores = nCores_bench;}
        else {
          #ifdef addsig2vol_debug
          PRINT("Overhead to big, switch to single thread.\n");
          #endif
          nCores = 1;}

        #ifdef addsig2vol_debug
    		PRINT("selectedNumCores:  %i\n", nCores);
    		PRINT("savedNumCORE:  %i\n", nCores_bench);
    		PRINT("perf_MT:  %e\n", ( (throughput[nCores_bench] * n_IMAGE)+latency[nCores_bench]));
    		PRINT("perf_single:  %e\n", ((throughput[1]*n_IMAGE)+latency[1]));
            #endif

        
        if (((speedz->x == n_X) & (speedz->y == n_Y) ) & ((speedz->z == n_Z) | (speedz->z == 1)))
        {
            PRINT("Info: Soundmap version\n");
            addsig2vol_mode = 2;
        }
        else{
        

            if ( (rec_posz->x != 3) |  (send_posz->x != 3) ) 
            {
            #ifdef addsig2vol_debug
             PRINT("Dim1 send_vec: %d, Dim1 rec_vec: %d",*(int*)rec_vec_ptr,*(int*)send_vec_ptr);
            PRINT("ascan_block: %d, blocksize rec_vec: %d, blocksize send_vec: %d",n_AScan_block,n_rec_vec_block,n_send_vec_block);
             #endif
            PRINT("Error: 3-d vectors needed for emitter & receiver positions or transposed blocked pos (1x3 instead of 3x1), break\n");
            return;}

             #ifdef addsig2vol_debug
                PRINT("nascan: %d,ascan_block: %d, blocksize rec_vec: %d, blocksize send_vec: %d",n_AScanz,n_AScan_block,n_rec_vec_block,n_send_vec_block);
             #endif

             if (!( (((n_AScan_block == n_rec_vec_block) & (n_AScan_block == n_send_vec_block)) | ((1 == n_rec_vec_block) & (n_AScan_block == n_send_vec_block)) | ((n_AScan_block == n_rec_vec_block) & (1 == n_send_vec_block))) & ((n_AScan_block == n_speed_vec_block) | (1 == n_speed_vec_block)) ))
              {
              PRINT("Error: Blocked sizes parameter mismatch. Size(AScan,2) has to be size(rec_vec,2) and/or size(send_vec,2), n_rec_vec_block:%i n_send_vec_block:%i n_speed_vec_block:%i break\n",n_rec_vec_block, n_send_vec_block, n_speed_vec_block);
              return;}

                
             if (!(((AScan_complexz==NULL) & (IMAGE_SUM_complexz==NULL)) | ((AScan_complexz!=NULL) & (IMAGE_SUM_complexz!=NULL))))
             {
                 
              
              addsig2vol_mode = 1;
            }
             else
                 
                  addsig2vol_mode = 0;
        }

        PRINT("Set addsig2vol_mode to %i\n", addsig2vol_mode);

        out_real = malloc(sizeof(double)*n_IMAGE);
        buffer_real = malloc(sizeof(double)*n_AScan*INTERP_RATIO);

    		if (AScan_complexz) {
                
                out_complex = malloc(sizeof(double)*n_IMAGE);
                buffer_complex = malloc(sizeof(double)*n_AScan*INTERP_RATIO);
    		}

            #ifdef addsig2vol_debug
            PRINT("Set addsig2vol_mode to %i\n", addsig2vol_mode);

            PRINT("Elements in out_real:  %i\n", ARRAYLENGTH(out_real));
    		PRINT("Elements in out_complex:  %i\n", ARRAYLENGTH(out_complex));
            PRINT("Elements in buffer_real:  %i\n", ARRAYLENGTH(buffer_real));
            PRINT("Elements in buffer_complex:  %i\n", ARRAYLENGTH(buffer_complex));
            #endif

    		
    		

    		
    		as2v_MT(out_real, AScan_realz_ptr, n_AScan, buffer_real,
    			     pix_vectz_ptr, n_X, rec_vec_ptr, send_vec_ptr,
    			     speedz_vec_ptr, resz, timeintz, AScan_complexz_ptr, buffer_complex, out_complex, n_Y, n_Z, IMAGE_SUM_realz_ptr, IMAGE_SUM_complexz_ptr);


    		
    		for (i = 2; i <= n_AScan_block; i++) {
    			
    			if (AScan_complexz != NULL)	AScan_complexz = AScan_complexz + (n_AScan * (i - 1)); 
                if (1<n_rec_vec_block)  rec_vec_ptr  = rec_posz->data  + (3 * sizeof(float) * (i - 1)); else rec_vec_ptr = rec_posz->data;
                if (1<n_send_vec_block) send_vec_ptr = send_posz->data + (3 * sizeof(float) * (i - 1)); else send_vec_ptr = send_posz->data;
                if (1<n_speed_vec_block) speedz_vec_ptr = speedz->data + (1 * sizeof(float) * (i - 1)); else speedz_vec_ptr = speedz->data;

    
    
    
    
    
    
    

    			
    			
    			

             as2v_MT(out_real, AScan_realz_ptr + (n_AScan * (i - 1)), n_AScan, buffer_real,
     			     pix_vectz_ptr, n_X, rec_vec_ptr, send_vec_ptr,
     			     speedz_vec_ptr, resz, timeintz, AScan_complexz_ptr, buffer_complex, out_complex, n_Y, n_Z, IMAGE_SUM_realz_ptr, IMAGE_SUM_complexz_ptr);

    		/*	as2v_MT(mxGetPr(out), mxGetPr(
    					     AScan) + (n_AScan * (i - 1)), n_AScan, mxGetPr(
    					     buffer), mxGetPr(
    					     pix_vect), n_X, (char*)mxGetPr(
    					     rec_pos) + (3 * sizeof(float) * (i - 1)),
    				     mxGetPr(send_pos), mxGetPr(
    					     speed), mxGetPr(res), mxGetPr(
    					     timeint), AScan_pi, mxGetPi(buffer),
    				     mxGetPi(out), n_Y, n_Z, mxGetPr(out), mxGetPi(
    					     out));*/
    			
    			

    
    
    
    
    
    

                }

           results.out_real = as2v_boxDoubleArray(out_real, n_X, n_Y, n_Z);
           results.out_complex = as2v_boxDoubleArray(out_complex, n_X, n_Y, n_Z);
           results.buffer_real = as2v_boxDoubleArray(buffer_real, n_AScan, n_AScan_block, 1);
           results.buffer_complex = as2v_boxDoubleArray(buffer_complex, n_AScan, n_AScan_block, 1);

    	return results;
    }
    
    return results;
}

void printIntro(void){
    PRINT("\naddSig2Vol_2 SSE1 Assembler Optimized 64bit LINUX&Windows v3.1 (Multiple Rec-AScan Vers.)\n\n Calculate the ellip. backprojection.\nUses SSE. Features: Win&Linux, 32&64bit version, SSE1-Assembler&C-Implementation, Multithreaded by PosixThreads (under windows pthreadsVC2.dll might be needed)\n\t %s. M.Zapf KIT-IPE\n\n",__DATE__);
    
    #ifdef addsig2vol_debug
        PRINT("addsig2vol_debug: 1");
        #ifdef C_CODE
            PRINT("C_CODE: 1, use c implementation");
        #elif
            PRINT("C_CODE: 0, use asm implementation");
        #endif

    #endif
}


int main(){

    unsigned int seed = time(NULL);
    srand(seed);

    unsigned int count = 100;
    unsigned int n_AScan = 1000;
    unsigned int x=100;
    float speedz[1] ={1500};


    upstreamDoubleArray(NULL);
     setUpstreamDoubleArray_callback(NULL);
     setUpstreamFloatArray_callback(NULL);

    as2v_doubleArray AScan    = as2v_mallocDoubleArray(n_AScan, count, 1);
    as2v_floatArray pix_vect  = as2v_mallocFloatArray(3, 1, 1);
    as2v_floatArray rec_posz  = as2v_mallocFloatArray(3, count, 1);
    as2v_floatArray send_posz = as2v_mallocFloatArray(3, count, 1);
    as2v_floatArray speed     = as2v_boxFloatArray(speedz, 1, 1, 1);

    as2v_doubleArray IMAGE_SUM = as2v_mallocDoubleArray(x, x, x);

    float timeintz[1] ={1e-7};
    float resz[1] ={0.001};

    for(int i =0; i<n_AScan;i++ ){
        AScan.data[i] = (int)(rand()%2);
    }
    for(int j=1; j < count; j++){
        for(int i =0; i<n_AScan;i++){
            AScan.data[i+n_AScan*j] = AScan.data[i];
        }
    }
    for(int i =0; i<count*3;i++)
    {   rec_posz.data[i] = 0.01*rand()/RAND_MAX;
        send_posz.data[i] = 0.01*rand()/RAND_MAX;
    }




    as2v_results results = as2v_addsig2vol_3(&AScan, NULL,
        &pix_vect, &rec_posz, &send_posz, &speed, resz, timeintz,
        &IMAGE_SUM, NULL);

    as2v_doubleArray outz = results.out_real;

    briefDoubleArray(&AScan);
    briefFloatArray(&send_posz);
    saveDoubleArray(&AScan, "AScan");

    PRINT("data size: %i\n”", outz.len);
    PRINT("x  y  z: %i  %i  %i\n", outz.x, outz.y, outz.z);

    saveDoubleArray(&outz, "outz");

    free(AScan.data);
    free(pix_vect.data);
    free(rec_posz.data);
    free(send_posz.data);
    free(IMAGE_SUM.data);

    return 0;
}
