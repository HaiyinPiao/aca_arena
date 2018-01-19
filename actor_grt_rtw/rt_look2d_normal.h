/*
 * rt_look2d_normal.h
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

#ifndef RTW_HEADER_rt_look2d_normal_h_
#define RTW_HEADER_rt_look2d_normal_h_
#include "rtwtypes.h"
#include "rt_look.h"
#ifdef __cplusplus

extern "C" {

#endif

  extern real_T rt_Lookup2D_Normal (const real_T *xVals, const int_T numX,
    const real_T *yVals, const int_T numY,
    const real_T *zVals,
    const real_T x, const real_T y);

#ifdef __cplusplus

}                                      /* extern "C" */
#endif
#endif                                 /* RTW_HEADER_rt_look2d_normal_h_ */
