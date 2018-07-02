#include "addsig2vol_3.h"
#include "as2v_array.h"
#include "mex.h"
#include "addsig2vol_3_unittests.h"
#include "threadstats.h"

//define matlab in and out
#define out       (plhs[0])
#define out2      (plhs[1])

#define tt1       (plhs[2])
#define tt2       (plhs[3])
#define tt3       (plhs[4])
#define tt4       (plhs[5])

#define AScan     (prhs[0])
#define pix_vect  (prhs[1])
#define rec_pos   (prhs[2])
#define send_pos  (prhs[3])
#define speed     (prhs[4])
#define res       (prhs[5])
#define timeint   (prhs[6])
#define IMAGE_XYZ (prhs[7])
#define IMAGE_SUM (prhs[8])

static char doThreadMeasuring = 0;

// Octave NOTE
// - printf druckt nach Octave, und zwar zur Ausführzeit. mexprint druckt am Ende

void upstreamCopiedDoubleArray(cArray* array)
{
    if(array){
    //data gets copied, so freeing mxArray done by the MATLAB memory manager dosn't destroy C data
    mxArray* mxa = caNewMxarrayFromCArray(array, NULL);
    // NOTE: this does a segfault work in octave
    //mexCallMATLAB(0, NULL, 1, mxa, "disp");
    }
}

void upstreamCopiedFloatArray(cArray* array)
{
    if(array){
    //data gets copied, so freeing mxArray done by the MATLAB memory manager dosn't destroy C data
    mxArray* mxa = caNewMxarrayFromCArray(array, NULL);
    // NOTE: this does a segfault work in octave
    //mexCallMATLAB(0, NULL, 1, mxa, "disp");
    }
}

//Jumping point for mex
void mexFunction (int nlhs, mxArray*plhs[], int nrhs, const mxArray*prhs[]) {
    //use mexCallMATLAB for debugging data
    //setUpstreamDoubleArray_callback(upstreamCopiedDoubleArray);
    //setUpstreamFloatArray_callback(upstreamCopiedFloatArray);

    // TODO unterscheiden zwischen Octave und Matlab bei Print Funktionen
    //as2v_setPrintCallback(&printf);
    //caSetPrintCallback(&printf);


    //if (nlhs > 2) mexErrMsgTxt("Too many output arguments.");

    switch (nrhs) {
        default:
            mexErrMsgTxt("Incorrect number of arguments.\n");
            break;

        case 0:;
            // Say hello
            printIntro();
            // Give plhs and prhs definitions
            mexPrintf("\n\n#define out       plhs[0] (Double(1:end))\n#define out2      plhs[1] DEBUG (Double(1:end))\n#define AScan     prhs[0] (Double(NxM))\n#define pix_vect  prhs[1] (Single(1:3))\n#define rec_pos   prhs[2] (Single(1:3xM) or Single(1:3x1))\n#define send_pos  prhs[3] (Single(1:3xM) or Single(1:3x1))\n#define speed     prhs[4] (Single (1x1 or 1xM))\n#define res       prhs[5] (Single)\n#define timeint   prhs[6] (Single)\n#define IMAGE_XYZ prhs[7] (UINT32(1:3))\n#define IMAGE_SUM prhs[8] (Double(1:end))\n");
            // Run benchmarks
            as2v_benchLocal();
            break;

        case 1:;
            // Force code to run on prhs[0] threads
            as2v_overwriteBenchresultToThreadcount_n((uint32_t) ceil(*((double*)mxGetPr(prhs[0]))));

            //runTests();
            break;

        case 2:;
            // Force code to run on prhs[0] threads
            // Parameter mode
            unsigned int param = (unsigned int) *(mxGetPr(prhs[0]));
            double value = (double) *(mxGetPr(prhs[1]));
            switch(param){
                case 0:;
                    as2v_overwriteBenchresultToThreadcount_n((uint32_t) value);
                break;
                case 1:
                    setSegmentation((unsigned int) value);
                break;
                case 2:;
                    if (((unsigned int) value) > 0) doThreadMeasuring = 0xff;
                    else doThreadMeasuring = 0;
                break;
                default:;break;
            }

            break;

        case 9:;
            // Pre C99: No declarations next after label
            ;
            // printf("Ascan %p\n", AScan);
            // printf("IMAGE_SUM %p\n", IMAGE_SUM);
            // printf("pix_vectz %p\n", pix_vect);
            // printf("rec_posz %p\n", rec_pos);
            // printf("send_posz %p\n", send_pos);
            // printf("speedz %p\n", speed);


            tsclearAll();
            tsclock(90);

            // mwSize* dim = (mwSize*) mxGetDimensions(IMAGE_SUM);
            // mwSize numberOfDimensions = mxGetNumberOfDimensions(IMAGE_SUM);
            // out = mxCreateNumericArray(numberOfDimensions, dim, mxDOUBLE_CLASS, mxREAL);

            // dim = (mwSize*) mxGetDimensions(AScan);
            // dim[0] = dim[0] *5; // TODO get INTERP_RATIO
            // dim[1] = dim[1]; // Buffer für gesamten Ascan
            // numberOfDimensions = 2;
            // out2 = mxCreateNumericArray(numberOfDimensions, dim, mxDOUBLE_CLASS, mxREAL);

            // ANNahme: out und out 2 haben die richtige Größe für den Output.
            // cArrayDouble out_image = caNewDoubleArrayFromMxarray(out, mxREAL);
            // cArrayDouble out_buffer = caNewDoubleArrayFromMxarray(out2, mxREAL);



            //   as2v_addsig2vol_3_mex(AScan, NULL, pix_vect, rec_pos, send_pos, speed, res, timeint, IMAGE_SUM, NULL, out, out2);

            // //printf("seg %i\n",(unsigned int) speedz.data[0] );
            // setSegmentation((unsigned int) speedz.data[0] );
            //
            //
            // //Run algorithm
            // speedz.data[0] = 1500; // fixed hack

            // version B
            // // 0x0 resizing is much faster than full allocation

            if (doThreadMeasuring > 0)
            {
                mwSize* dimS = (mwSize*) mxGetDimensions(AScan);
                mwSize* dimI = (mwSize*) mxGetDimensions(IMAGE_SUM);

                unsigned int tasks = as2v_determineSegmentation(dimI[0], dimI[1], dimI[2]);
                unsigned int runs = dimS[1]; //Number of Ascans

                double* _dataTimeStartTask = mxMalloc(tasks*runs*sizeof(double));
                double* _dataTimeEndTask = mxMalloc(tasks*runs*sizeof(double));
                double* _dataMoveToTask = mxMalloc(tasks*runs*sizeof(double));
                double* _dataThreadnumber = mxMalloc(tasks*runs*sizeof(double));

                threadstats_init_mexed(as2v_getThreadsNumber(), tasks, runs, _dataTimeStartTask, _dataTimeEndTask, _dataMoveToTask, _dataThreadnumber);

                mwSize sizing[2] = {tasks, runs};
                mwSize nod = 2;
                tt1 = mxCreateDoubleMatrix( 0, 0, mxREAL );
                tt2 = mxCreateDoubleMatrix( 0, 0, mxREAL );
                tt3 = mxCreateDoubleMatrix( 0, 0, mxREAL );
                tt4 = mxCreateDoubleMatrix( 0, 0, mxREAL );

                mxSetDimensions(tt1, sizing, nod);
                mxSetDimensions(tt2, sizing, nod);
                mxSetDimensions(tt3, sizing, nod);
                mxSetDimensions(tt4, sizing, nod);

                mxSetPr(tt1, _dataTimeStartTask);
                mxSetPr(tt2, _dataTimeEndTask);
                mxSetPr(tt3, _dataMoveToTask);
                mxSetPr(tt4, _dataThreadnumber);
            }

            #ifdef BUILDMEX

                mwSize* dim = (mwSize*) mxGetDimensions(IMAGE_SUM);
                mwSize numberOfDimensions = mxGetNumberOfDimensions(IMAGE_SUM);
                double* outdata = mxMalloc(dim[0]*dim[1]*dim[2]*sizeof(double));
                //cArrayDouble out_image = caNewArrayDoubleFromData(outdata, dim[0],dim[1],dim[2]);
                out = mxCreateDoubleMatrix( 0, 0, mxREAL );
                mxSetDimensions(out, dim, numberOfDimensions);
                mxSetPr(out, outdata);

                tsclock(90);
                tsprint(90, TS_MILI);


                tsclock(91);
                as2v_addsig2vol_3_mex(AScan, AScan, pix_vect, rec_pos, send_pos, speed, mxGetPr(res), mxGetPr(timeint), IMAGE_SUM, IMAGE_SUM, out, NULL);
                tsclock(91);
                tsprint(91, TS_MILI);
            #endif

            #ifndef BUILDMEX
                //Box mxarrays into input arrays (data linked, not copied) -> Freeing will be done by MATLAB memory management
                cArrayDouble AScan_realz = caNewDoubleArrayFromMxarray(AScan, mxREAL);
                cArrayDouble AScan_complexz = caNewDoubleArrayFromMxarray(AScan, mxCOMPLEX);
                cArrayDouble IMAGE_SUM_realz = caNewDoubleArrayFromMxarray(IMAGE_SUM, mxREAL);
                cArrayDouble IMAGE_SUM_complexz = caNewDoubleArrayFromMxarray(IMAGE_SUM, mxCOMPLEX);
                cArrayFloat pix_vectz = caNewFloatArrayFromMxarray(pix_vect,  mxREAL);
                cArrayFloat rec_posz = caNewFloatArrayFromMxarray(rec_pos, mxREAL);
                cArrayFloat send_posz = caNewFloatArrayFromMxarray(send_pos,  mxREAL);
                cArrayFloat speedz = caNewFloatArrayFromMxarray(speed, mxREAL);

                mwSize* dim = (mwSize*) mxGetDimensions(IMAGE_SUM);
                mwSize numberOfDimensions = mxGetNumberOfDimensions(IMAGE_SUM);
                double* outdata = mxMalloc(dim[0]*dim[1]*dim[2]*sizeof(double));
                cArrayDouble out_image = caNewArrayDoubleFromData(outdata, dim[0],dim[1],dim[2]);
                out = mxCreateDoubleMatrix( 0, 0, mxREAL );
                mxSetDimensions(out, dim, numberOfDimensions);
                mxSetPr(out, outdata);

                tsclock(90);
                tsprint(90, TS_MILI);

                tsclock(91);
                as2v_addsig2vol_3(&AScan_realz, &AScan_complexz, &pix_vectz, &rec_posz, &send_posz, &speedz, mxGetPr(res), mxGetPr(timeint), &IMAGE_SUM_realz, &IMAGE_SUM_complexz, &out_image, NULL);
                tsclock(91);
                tsprint(91, TS_MILI);

            #endif


            //
            // if (doThreadMeasuring)
            // {
            //     threadstats_free_mexed();
            // }
        return;
    }
}
