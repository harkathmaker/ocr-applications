#define __OCR__
#include "ocr.h"

/* STREAM Headers */
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <sys/time.h> 

#ifndef STREAM_ARRAY_SIZE
#	define STREAM_ARRAY_SIZE	100
#endif

#ifdef NTIMES
#if NTIMES <= 1
#	define NTIMES	10
#endif
#endif

#ifndef NTIMES
#	define NTIMES	10
#endif

#ifndef OFFSET
#	define OFFSET	0
#endif

#define HLINE "-------------------------------------------------------------\n"

#ifndef MIN
#	define MIN(x,y) ((x)<(y)?(x):(y))
#endif

#ifndef MAX
#	define MAX(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

static STREAM_TYPE	a[STREAM_ARRAY_SIZE + OFFSET],
					b[STREAM_ARRAY_SIZE + OFFSET],
					c[STREAM_ARRAY_SIZE + OFFSET];

static double	avgtime[4] = {0}, maxtime[4] = {0},
				mintime[4] = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};


static char *label[4] = {"Copy:      ", "Scale:     ", "Add:       ", "Triad:     "};

static double	bytes[4] = {
	2 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
	2 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
	3 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
	3 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE
};

extern double mysecond();
extern void checkSTREAMresults();
#ifdef TUNED
extern void tuned_STREAM_Copy();
extern void tuned_STREAM_Scale(STREAM_TYPE scalar);
extern void tuned_STREAM_Add();
extern void tuned_STREAM_Triad(STREAM_TYPE scalar);
#endif
#ifdef _OPENMP
extern int omp_get_num_threads();
#endif





ocrGuid_t copy(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	ssize_t j;
	STREAM_TYPE ** data = depv[0].ptr;
	// data[1] = c, data[0] = a
	for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		data[1][j] = data[0][j];
	PRINTF("Finished COPY\n");
	//c[STREAM_ARRAY_SIZE - 1] = %f\na[STREAM_ARRAY_SIZE - 1] = %f\n", data[1][STREAM_ARRAY_SIZE - 1], data[0][STREAM_ARRAY_SIZE - 1]);  
	return NULL_GUID;
}

ocrGuid_t scale(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	ssize_t j;
	STREAM_TYPE scalar = 3.0;
	STREAM_TYPE ** data = depv[0].ptr;
	// data[0] = b, data[1] = c, b = 3 * c
	for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		data[0][j] = scalar*data[1][j];
	PRINTF("Finished SCALE\n");
	//b[STREAM_ARRAY_SIZE - 1] = %f\nc[STREAM_ARRAY_SIZE - 1] = %f\n", data[0][STREAM_ARRAY_SIZE - 1], data[1][STREAM_ARRAY_SIZE - 1]); 
	return NULL_GUID;
}

ocrGuid_t add(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	ssize_t j;
	STREAM_TYPE ** data = depv[0].ptr;
	// data[0] = a, data[1] = b, data[2] = c,    c = a + b
	for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		data[2][j] = data[0][j] + data[1][j];
	PRINTF("Finished ADD\n");
	//c = %f + %f = %f\n", data[0][1], data[1][1], data[2][1]); 
	return NULL_GUID;
}

ocrGuid_t triad(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	ssize_t j;
	STREAM_TYPE scalar = 3.0;
	STREAM_TYPE ** data = depv[0].ptr;
	// data[0] = a, data[1] = b, data[2] = c,    a = b + scalar * c
	for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		data[0][j] = data[1][j] + scalar * data[2][j];
	PRINTF("Finished TRIAD\n");
	//a = %f + 3.0 * %f = %f\n", data[1][1], data[2][1], data[0][1]); 
	return NULL_GUID;
}

ocrGuid_t pipelineEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]){
	ocrGuid_t * guids = (ocrGuid_t *) depv[0].ptr;
	ocrGuid_t dataGuid = guids[0];
	ocrGuid_t doneEventGuid = guids[1];


	// COPY
	// EDT Tempalte for copy
	ocrGuid_t copyEdtTemplateGuid;
	ocrEdtTemplateCreate(&copyEdtTemplateGuid, copy, 0 /*paramc*/, EDT_PARAM_UNK /*depc*/);

	// Datablock for copy
	STREAM_TYPE ** copyDataArray;
	ocrGuid_t copyDataGuid;
	DBCREATE(&copyDataGuid, (void **) &copyDataArray, sizeof(STREAM_TYPE)*STREAM_ARRAY_SIZE, 0, NULL_GUID, NO_ALLOC);  

	// Formatting datablock for copy
	copyDataArray[0] = a;
	copyDataArray[1] = c; 
	//PRINTF("copyDataArray[0][0] = %f, copyDataArray[1][0] = %f\n", copyDataArray[0][0], copyDataArray[1][0]);

	// EDT for copy
	ocrGuid_t copyEdtGuid;
	ocrGuid_t copyOutput;
	ocrEdtCreate(&copyEdtGuid, copyEdtTemplateGuid, EDT_PARAM_DEF, NULL, 2, NULL/*&copyDataGuid*/,
				/*prop=*/EDT_PROP_NONE, NULL_GUID, &copyOutput);



	// SCALE
	// EDT Template for scale
	ocrGuid_t scaleEdtTemplateGuid;
	ocrEdtTemplateCreate(&scaleEdtTemplateGuid, scale, 0 /*paramc*/, EDT_PARAM_UNK /*depc*/);

	// Datablock for scale
	STREAM_TYPE ** scaleDataArray;
	ocrGuid_t scaleDataGuid;
	DBCREATE(&scaleDataGuid, (void **) &scaleDataArray, sizeof(STREAM_TYPE)*STREAM_ARRAY_SIZE, 0, NULL_GUID, NO_ALLOC);  

	// Formatting datablock for scale
	scaleDataArray[0] = b;
	scaleDataArray[1] = c; 
	//PRINTF("scaleDataArray[0][0] = %f, scaleDataArray[1][0] = %f\n", scaleDataArray[0][0], scaleDataArray[1][0]);

	// EDT for scale
	ocrGuid_t scaleEdtGuid;
	ocrGuid_t scaleOutput;
	ocrEdtCreate(&scaleEdtGuid, scaleEdtTemplateGuid, EDT_PARAM_DEF, NULL, 2, NULL,
				/*prop=*/EDT_PROP_NONE, NULL_GUID, &scaleOutput);



	// ADD
	// EDT Template for add
	ocrGuid_t addEdtTemplateGuid;
	ocrEdtTemplateCreate(&addEdtTemplateGuid, add, 0 /*paramc*/, EDT_PARAM_UNK /*depc*/);

	// Datablock for add
	STREAM_TYPE ** addDataArray;
	ocrGuid_t addDataGuid;
	DBCREATE(&addDataGuid, (void **) &addDataArray, sizeof(STREAM_TYPE)*STREAM_ARRAY_SIZE, 0, NULL_GUID, NO_ALLOC);  

	// Formatting datablock for add
	addDataArray[0] = a;
	addDataArray[1] = b;
	addDataArray[2] = c;
	//PRINTF("addDataArray[0][0] = %f, addDataArray[1][0] = %f\n", addDataArray[0][0], addDataArray[1][0]);

	// EDT for add
	ocrGuid_t addEdtGuid;
	ocrGuid_t addOutput;
	ocrEdtCreate(&addEdtGuid, addEdtTemplateGuid, EDT_PARAM_DEF, NULL, 2, NULL,
				/*prop=*/EDT_PROP_NONE, NULL_GUID, &addOutput);



	// TRIAD
	// EDT Template for triad
	ocrGuid_t triadEdtTemplateGuid;
	ocrEdtTemplateCreate(&triadEdtTemplateGuid, triad, 0 /*paramc*/, EDT_PARAM_UNK /*depc*/);

	// Datablock for triad
	STREAM_TYPE ** triadDataArray;
	ocrGuid_t triadDataGuid;
	DBCREATE(&triadDataGuid, (void **) &triadDataArray, sizeof(STREAM_TYPE)*STREAM_ARRAY_SIZE, 0, NULL_GUID, NO_ALLOC);  

	// Formatting datablock for triad
	triadDataArray[0] = a;
	triadDataArray[1] = b;
	triadDataArray[2] = c;
	//PRINTF("triadDataArray[0][0] = %f, triadDataArray[1][0] = %f\n", triadDataArray[0][0], triadDataArray[1][0]);

	// EDT for triad
	ocrGuid_t triadEdtGuid;
	ocrGuid_t triadOutput;
	ocrEdtCreate(&triadEdtGuid, triadEdtTemplateGuid, EDT_PARAM_DEF, NULL, 2, NULL,
				/*prop=*/EDT_PROP_NONE, NULL_GUID, &triadOutput);





	ocrGuid_t copyDone;
	ocrGuid_t scaleDone;
	ocrGuid_t addDone;
	ocrGuid_t triadDone;

	ocrEventCreate(&copyDone, OCR_EVENT_IDEM_T, false);
	ocrEventCreate(&scaleDone, OCR_EVENT_IDEM_T, false);
	ocrEventCreate(&addDone, OCR_EVENT_IDEM_T, false);
	ocrEventCreate(&triadDone, OCR_EVENT_IDEM_T, false);

	ocrAddDependence(copyDataGuid, copyEdtGuid, 0, DB_MODE_EW);
	// dataGuid dependence on first test in pipeline
	ocrAddDependence(dataGuid, copyEdtGuid, 1, DB_MODE_EW);
	ocrAddDependence(copyOutput, copyDone, 0, DB_MODE_RO); 

	ocrAddDependence(scaleDataGuid, scaleEdtGuid, 0, DB_MODE_EW);
	ocrAddDependence(copyDone, scaleEdtGuid, 1, DB_MODE_RO);

	ocrAddDependence(scaleOutput, scaleDone, 0, DB_MODE_RO);
	ocrAddDependence(addDataGuid, addEdtGuid, 0, DB_MODE_EW);
	ocrAddDependence(scaleDone, addEdtGuid, 1, DB_MODE_RO);

	ocrAddDependence(addOutput, addDone, 0, DB_MODE_RO);
	ocrAddDependence(triadDataGuid, triadEdtGuid, 0, DB_MODE_EW);
	ocrAddDependence(addDone, triadEdtGuid, 1, DB_MODE_RO);
	ocrAddDependence(triadOutput, triadDone, 0, DB_MODE_RO); 

	// Output dependence signalling end of pipeline
	ocrAddDependence(triadDone, doneEventGuid, 0, DB_MODE_RO);
	return NULL_GUID;
}



ocrGuid_t iterEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	// Read the data-block guid, passed through the event that scheduled this 'iterEdt' instance
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;

	// In this implementation we pass down the event the pipeline 
	// need to satisfy on completion. We can make the next iteration
	// to depend on this event being satisfied.

	// Create the event the pipeline EDT satisfies on pipeline completion
	ocrGuid_t iterDoneEventGuid;
	ocrEventCreate(&iterDoneEventGuid, OCR_EVENT_ONCE_T, false);

	// Create a new datablock to hold data and event guids to synchronize on
	ocrGuid_t * guidArray;
	ocrGuid_t guidArrayGuid;
	DBCREATE(&guidArrayGuid,(void **) &guidArray, sizeof(ocrGuid_t)*2, /*flags=*/0, /*loc=*/NULL_GUID, NO_ALLOC);
	guidArray[0] = dataGuid;
	guidArray[1] = iterDoneEventGuid;

	// Setup the pipeline EDT
	ocrGuid_t pipeTemplateGuid;
	ocrEdtTemplateCreate(&pipeTemplateGuid, pipelineEdt, 3 /*paramc*/, 1 /*depc*/);

	// Create this iteration of the pipeline EDT
	ocrGuid_t pipeGuid;
	ocrEdtCreate(&pipeGuid, pipeTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, &guidArrayGuid,
				/*prop=*/EDT_PROP_NONE, NULL_GUID, NULL);

	u64 i = paramv[0];
	u64 nbIt = paramv[1];
	i++;
	if (i <= nbIt) {
		PRINTF("Pushing iteration %lu\n", i);
		// Setup next iteration of 'iterEdt' that depends on 'iterDoneEventGuid' being satisfied with some data-block
		ocrGuid_t ndepv[1];
		ndepv[0] = iterDoneEventGuid;
		// Setup next iteration arguments (can reuse paramv here since it's copied)
		paramv[0] = i;

		ocrGuid_t nextIterTemplateGuid;
		ocrEdtTemplateCreate(&nextIterTemplateGuid, iterEdt, 3 /*paramc*/, 1/*depc*/);

		ocrGuid_t nextIterGuid;
		ocrEdtCreate(&nextIterGuid, nextIterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, ndepv,
					EDT_PROP_NONE, NULL_GUID, NULL);        
	} else {
		ocrShutdown();   
	}
	return NULL_GUID;
}



ocrGuid_t terminateEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	ocrShutdown(); // This is the last EDT to execute, terminate
	return NULL_GUID;
}




ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	int			quantum, checktick();
	int			BytesPerWord;
	int			k;
	ssize_t		j;
	STREAM_TYPE		scalar;
	double		t, times[4][NTIMES];

	/* --- SETUP --- determine precision and check timing --- */
	// printf(HLINE);
	// printf("STREAM version $Revision: 5.10 $\n");
	// printf(HLINE);
	// BytesPerWord = sizeof(STREAM_TYPE);
	// printf("This system uses %d bytes per array element.\n", BytesPerWord);
	// printf(HLINE);

// #ifdef N
	// printf("*****  WARNING: ******\n");
	// printf("      It appears that you set the preprocessor variable N when compiling this code.\n");
	// printf("      This version of the code uses the preprocesor variable STREAM_ARRAY_SIZE to control the array size\n");
	// printf("      Reverting to default value of STREAM_ARRAY_SIZE=%llu\n",(unsigned long long) STREAM_ARRAY_SIZE);
	// printf("*****  WARNING: ******\n");
// #endif

	// printf("Array size = %llu (elements), Offset = %d (elements)\n" , (unsigned long long) STREAM_ARRAY_SIZE, OFFSET);
	// printf("Memory per array = %.1f MiB (= %.1f GiB).\n", 
	// BytesPerWord * ( (double) STREAM_ARRAY_SIZE / 1024.0/1024.0),
	// BytesPerWord * ( (double) STREAM_ARRAY_SIZE / 1024.0/1024.0/1024.0));
	// printf("Total memory required = %.1f MiB (= %.1f GiB).\n",
			// (3.0 * BytesPerWord) * ( (double) STREAM_ARRAY_SIZE / 1024.0/1024.),
			// (3.0 * BytesPerWord) * ( (double) STREAM_ARRAY_SIZE / 1024.0/1024./1024.));
	// printf("Each kernel will be executed %d times.\n", NTIMES);
	// printf(" The *best* time for each kernel (excluding the first iteration)\n"); 
	// printf(" will be used to compute the reported bandwidth.\n");

#ifdef _OPENMP
	printf(HLINE);
#pragma omp parallel 
	{
#pragma omp master
		{
		k = omp_get_num_threads();
		printf ("Number of Threads requested = %i\n",k);
		}
	}
#endif

#ifdef _OPENMP
	k = 0;
#pragma omp parallel
#pragma omp atomic 
		k++;
		printf ("Number of Threads counted = %i\n",k);
#endif

	/* Get initial value for system clock. */
#pragma omp parallel for
	for (j=0; j<STREAM_ARRAY_SIZE; j++) {
		a[j] = 1.0;
		b[j] = 2.0;
		c[j] = 0.0;
	}
	// printf(HLINE);

	// if ( (quantum = checktick()) >= 1) 
		// printf("Your clock granularity/precision appears to be %d microseconds.\n", quantum);
	// else {
		// printf("Your clock granularity appears to be less than one microsecond.\n");
		// quantum = 1;
	// }

	t = mysecond();
#pragma omp parallel for
	for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		a[j] = 2.0E0 * a[j];
	t = 1.0E6 * (mysecond() - t);

	// printf("Each test below will take on the order of %d microseconds.\n", (int) t  );
	// printf("   (= %d clock ticks)\n", (int) (t/quantum) );
	// printf("Increase the size of the arrays if this shows that\n");
	// printf("you are not getting at least 20 clock ticks per test.\n");
	// printf(HLINE);
	// printf("WARNING -- The above is only a rough guideline.\n");
	// printf("For best results, please be sure you know the\n");
	// printf("precision of your system timer.\n");
	// printf(HLINE);






	/* --- MAIN LOOP --- repeat test cases NTIMES times --- */
	u64 i = 0;
	u64 nbIt = NTIMES;

	u64 nparamv[3];
	nparamv[0] = i;
	nparamv[1] = nbIt;
	nparamv[2] = STREAM_ARRAY_SIZE;

	// Setup datablock
	u64 * dataArray;
	ocrGuid_t dataGuid;
	DBCREATE(&dataGuid,(void **) &dataArray, sizeof(u64)*STREAM_ARRAY_SIZE, /*flags=*/0, /*loc=*/NULL_GUID, NO_ALLOC);

	ocrGuid_t iterTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, iterEdt, 3 /*paramc*/, 1 /*depc*/);

	ocrGuid_t startEventGuid;
	// Create a ONCE event
	ocrEventCreate(&startEventGuid, OCR_EVENT_ONCE_T, true);

	// Create pipeline EDT that depends on the ONCE event
	ocrGuid_t iterGuid;
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, &startEventGuid,
				/*prop=*/EDT_PROP_NONE, NULL_GUID, NULL);

	// Satisfy the ONCE event with the datablock guid to schedule the iter EDT
	ocrEventSatisfy(startEventGuid, dataGuid);






	// Pipeline Version
	// ocrGuid_t pipeTemplateGuid;
	// ocrEdtTemplateCreate(&pipeTemplateGuid, pipelineEdt, 0 /*paramc*/, 0 /*depc*/);

	// ocrGuid_t outputEventGuid;
	// ocrGuid_t pipeGuid;
	// ocrEdtCreate(&pipeGuid, pipeTemplateGuid, EDT_PARAM_DEF, NULL, EDT_PARAM_DEF, NULL,
				// /*prop=*/EDT_PROP_FINISH, NULL_GUID, &outputEventGuid);

	// ocrGuid_t terminateTemplateGuid;
	// ocrEdtTemplateCreate(&terminateTemplateGuid, terminateEdt, 0 /*paramc*/, 1 /*depc*/);

	// ocrGuid_t terminateGuid;
	// ocrEdtCreate(&terminateGuid, terminateTemplateGuid, EDT_PARAM_DEF, 0, EDT_PARAM_DEF, &outputEventGuid,
				// /*prop=*/0, NULL_GUID, NULL);
	return NULL_GUID;
}












# define	M	20

int
checktick()
    {
    int		i, minDelta, Delta;
    double	t1, t2, timesfound[M];

/*  Collect a sequence of M unique time values from the system. */

    for (i = 0; i < M; i++) {
	t1 = mysecond();
	while( ((t2=mysecond()) - t1) < 1.0E-6 )
	    ;
	timesfound[i] = t1 = t2;
	}

/*
 * Determine the minimum difference between these M values.
 * This result will be our estimate (in microseconds) for the
 * clock granularity.
 */

    minDelta = 1000000;
    for (i = 1; i < M; i++) {
	Delta = (int)( 1.0E6 * (timesfound[i]-timesfound[i-1]));
	minDelta = MIN(minDelta, MAX(Delta,0));
	}

   return(minDelta);
    }



/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */

#include <sys/time.h>

double mysecond()
{
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

#ifndef abs
#define abs(a) ((a) >= 0 ? (a) : -(a))
#endif
void checkSTREAMresults ()
{
	STREAM_TYPE aj,bj,cj,scalar;
	STREAM_TYPE aSumErr,bSumErr,cSumErr;
	STREAM_TYPE aAvgErr,bAvgErr,cAvgErr;
	double epsilon;
	ssize_t	j;
	int	k,ierr,err;

    /* reproduce initialization */
	aj = 1.0;
	bj = 2.0;
	cj = 0.0;
    /* a[] is modified during timing check */
	aj = 2.0E0 * aj;
    /* now execute timing loop */
	scalar = 3.0;
	for (k=0; k<NTIMES; k++)
        {
            cj = aj;
            bj = scalar*cj;
            cj = aj+bj;
            aj = bj+scalar*cj;
        }

    /* accumulate deltas between observed and expected results */
	aSumErr = 0.0;
	bSumErr = 0.0;
	cSumErr = 0.0;
	for (j=0; j<STREAM_ARRAY_SIZE; j++) {
		aSumErr += abs(a[j] - aj);
		bSumErr += abs(b[j] - bj);
		cSumErr += abs(c[j] - cj);
		// if (j == 417) printf("Index 417: c[j]: %f, cj: %f\n",c[j],cj);	// MCCALPIN
	}
	aAvgErr = aSumErr / (STREAM_TYPE) STREAM_ARRAY_SIZE;
	bAvgErr = bSumErr / (STREAM_TYPE) STREAM_ARRAY_SIZE;
	cAvgErr = cSumErr / (STREAM_TYPE) STREAM_ARRAY_SIZE;

	if (sizeof(STREAM_TYPE) == 4) {
		epsilon = 1.e-6;
	}
	else if (sizeof(STREAM_TYPE) == 8) {
		epsilon = 1.e-13;
	}
	else {
		printf("WEIRD: sizeof(STREAM_TYPE) = %lu\n",sizeof(STREAM_TYPE));
		epsilon = 1.e-6;
	}

	err = 0;
	if (abs(aAvgErr/aj) > epsilon) {
		err++;
		printf ("Failed Validation on array a[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
		printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",aj,aAvgErr,abs(aAvgErr)/aj);
		ierr = 0;
		for (j=0; j<STREAM_ARRAY_SIZE; j++) {
			if (abs(a[j]/aj-1.0) > epsilon) {
				ierr++;
#ifdef VERBOSE
				if (ierr < 10) {
					printf("         array a: index: %ld, expected: %e, observed: %e, relative error: %e\n",
						j,aj,a[j],abs((aj-a[j])/aAvgErr));
				}
#endif
			}
		}
		printf("     For array a[], %d errors were found.\n",ierr);
	}
	if (abs(bAvgErr/bj) > epsilon) {
		err++;
		printf ("Failed Validation on array b[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
		printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",bj,bAvgErr,abs(bAvgErr)/bj);
		printf ("     AvgRelAbsErr > Epsilon (%e)\n",epsilon);
		ierr = 0;
		for (j=0; j<STREAM_ARRAY_SIZE; j++) {
			if (abs(b[j]/bj-1.0) > epsilon) {
				ierr++;
#ifdef VERBOSE
				if (ierr < 10) {
					printf("         array b: index: %ld, expected: %e, observed: %e, relative error: %e\n",
						j,bj,b[j],abs((bj-b[j])/bAvgErr));
				}
#endif
			}
		}
		printf("     For array b[], %d errors were found.\n",ierr);
	}
	if (abs(cAvgErr/cj) > epsilon) {
		err++;
		printf ("Failed Validation on array c[], AvgRelAbsErr > epsilon (%e)\n",epsilon);
		printf ("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n",cj,cAvgErr,abs(cAvgErr)/cj);
		printf ("     AvgRelAbsErr > Epsilon (%e)\n",epsilon);
		ierr = 0;
		for (j=0; j<STREAM_ARRAY_SIZE; j++) {
			if (abs(c[j]/cj-1.0) > epsilon) {
				ierr++;
#ifdef VERBOSE
				if (ierr < 10) {
					printf("         array c: index: %ld, expected: %e, observed: %e, relative error: %e\n",
						j,cj,c[j],abs((cj-c[j])/cAvgErr));
				}
#endif
			}
		}
		printf("     For array c[], %d errors were found.\n",ierr);
	}
	if (err == 0) {
		printf ("Solution Validates: avg error less than %e on all three arrays\n",epsilon);
	}
#ifdef VERBOSE
	printf ("Results Validation Verbose Results: \n");
	printf ("    Expected a(1), b(1), c(1): %f %f %f \n",aj,bj,cj);
	printf ("    Observed a(1), b(1), c(1): %f %f %f \n",a[1],b[1],c[1]);
	printf ("    Rel Errors on a, b, c:     %e %e %e \n",abs(aAvgErr/aj),abs(bAvgErr/bj),abs(cAvgErr/cj));
#endif
}
