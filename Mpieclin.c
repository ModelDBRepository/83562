/* general purpose sigmoid:
 * use as template for single input, single output function system (no states)
 */

#define RAND_FACT 2147483647.0	/* largest random number generated */

#define S_FUNCTION_NAME Msigmoid2 /* same as filename */

#include "simstruc.h"
extern double exp();
extern double log();
/* extern long random(); */

/* input arguments (parameters) */

#define SLOPE    ssGetArg(S, 0)
#define THRESH  ssGetArg(S, 1)

/*
 * mdlInitializeSizes - initialize the sizes array
 *
 * for single input, single output function system (no states)
 * notice use of ssSetNumInputArgs() to set number of input pars 
 */

static void mdlInitializeSizes(S)
    SimStruct *S;
{
    ssSetNumContStates(    S, 0);      /* number of continuous states */
    ssSetNumDiscStates(    S, 0);      /* number of discrete states */
    ssSetNumInputs(        S, 1);      /* number of inputs */
    ssSetNumOutputs(       S, 1);      /* number of outputs */
    ssSetDirectFeedThrough(S, 1);      /* direct feedthrough flag */
    ssSetNumSampleTimes(   S, 1);      /* number of sample times */
    ssSetNumInputArgs(     S, 2);      /* number of input arguments */
    ssSetNumRWork(         S, 0);      /* number of real work vector elements */
    ssSetNumIWork(         S, 0);      /* number of integer work vector elements */
    ssSetNumPWork(         S, 0);      /* number of pointer work vector elements */
}

/*
 * mdlInitializeSampleTimes - initialize the sample times array
 *
 * This function is used to specify the sample time(s) for your S-function.
 * If your S-function is continuous, you must specify a sample time of 0.0.
 * Sample times must be registered in ascending order.
 */

static void mdlInitializeSampleTimes(S)
    SimStruct *S;
{
    ssSetSampleTimeEvent(S, 0, 0.0);
    ssSetOffsetTimeEvent(S, 0, 0.0);
}

/* have to have th efiollowing functions even though they are empty 
 * otherwise get segmentation violations */

static void mdlInitializeConditions(x0, S)
    double *x0;
    SimStruct *S;
{
}

static void mdlUpdate(x, u, S, tid)
    double *x, *u;
    SimStruct *S;
    int tid;
{
}

static void mdlDerivatives(dx, x, u, S, tid)
    double *dx, *x, *u;
    SimStruct *S;
    int tid;
{
}

static void mdlTerminate(S)
    SimStruct *S;
{
}

/***************** main function calculation **********/
/*
 * mdlOutputs - compute the outputs
 *
 */

static void mdlOutputs(y, x, u, S, tid)
    double *y, *x, *u;
    SimStruct *S;
    int tid;
{
     double slope, thresh;
     double act;
     double y_out;

     slope = mxGetPr(SLOPE)[0];
     thresh = mxGetPr(THRESH)[0];

     act = u[0];
     if (act <= thresh) {
	  y_out = 0;
     }
     else if ((act > thresh) && (act < 1/slope + thresh)) {
	  y_out = slope * (act - thresh);
     }
     else {
	  y_out = 1;
     }
     y[0] = y_out;

}


#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
