/*
 * rt_look.h
 *
 * Code generation for model "actor".
 *
 * Model version              : 1.1483
 * Simulink Coder version : 8.11 (R2016b) 25-Aug-2016
 * C++ source code generated on : Fri Jan 19 11:20:06 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_rt_look_h_
#define RTW_HEADER_rt_look_h_
#include "rtwtypes.h"
#ifdef DOINTERPSEARCH
#include <float.h>
#endif

#ifdef __cplusplus

extern "C" {

#endif

#ifndef INTERP
# define INTERP(x,x1,x2,y1,y2)         ( (y1)+(((y2) - (y1))/((x2) - (x1)))*((x)-(x1)) )
#endif

#ifndef ZEROTECHNIQUE
#define ZEROTECHNIQUE

  typedef enum {
    NORMAL_INTERP,
    AVERAGE_VALUE,
    MIDDLE_VALUE
  } ZeroTechnique;

#endif

  extern int_T rt_GetLookupIndex(const real_T *x, int_T xlen, real_T u) ;

#ifdef __cplusplus

}                                      /* extern "C" */
#endif
#endif                                 /* RTW_HEADER_rt_look_h_ */
