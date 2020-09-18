#include "typedef.h"

void testStruct(radarDataType   *radarData,
                stormStructType *stormStruct,
                retParamType    *retParam,  
                int *nmu, radarRetType    *radarRet,
                int *iflag, float *rms1,
                float *rms2, float *dn, int *iit, 
                float *xscalev,
                float *randemiss, float *localZAngle, 
                float *wfractPix, long *ichunk, 
                int *i0, int *j0, float *dZms, int *msFlag)
{
  int i;
  printf("%3i %3i %3i\n", stormStruct->nodes[0],stormStruct->nodes[1],
         stormStruct->nodes[2]);
  printf("%3i %3i %6.2f\n", stormStruct->rainType,
         stormStruct->iSurf,stormStruct->freezH);
  printf("%6.2f \n",radarData->xlong);
  for(i=0;i<88;i++)
    printf("%6.2f ",radarData->z13obs[i]);
  printf("%3i %3i\n",radarData->ngates,*nmu);
  printf("%i %li %i\n",*msFlag,*ichunk,*iit);
  printf("%g %g %g\n",*dn,*wfractPix,*localZAngle);
}


