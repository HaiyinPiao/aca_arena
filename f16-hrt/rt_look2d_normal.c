/*
 * rt_look2d_normal.c
 *
 * Code generation for model "f16".
 *
 * Model version              : 1.1162
 * Simulink Coder version : 8.11 (R2016b) 25-Aug-2016
 * C source code generated on : Tue Jan  9 21:21:58 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#include "rt_look2d_normal.h"

/* 2D normal lookup routine for data type of real_T. */
real_T rt_Lookup2D_Normal(const real_T *xVals, const int_T numX,
  const real_T *yVals, const int_T numY,
  const real_T *zVals,
  const real_T x, const real_T y)
{
  int_T xIdx, yIdx;
  real_T ylo, yhi;
  real_T Zx0yhi, Zx0ylo, xlo, xhi;
  real_T corner1, corner2;
  xIdx = rt_GetLookupIndex(xVals,numX,x);
  xlo = xVals[xIdx];
  xhi = xVals[xIdx+1];
  yIdx = rt_GetLookupIndex(yVals,numY,y);
  ylo = yVals[yIdx];
  yhi = yVals[yIdx+1];
  corner1 = *(zVals + xIdx + (numX * yIdx));
  corner2 = *(zVals + (xIdx+1) + (numX * yIdx));
  Zx0ylo = INTERP(x, xlo, xhi, corner1, corner2);
  corner1 = *(zVals + xIdx + (numX * (yIdx+1)));
  corner2 = *(zVals + (xIdx+1) + (numX*(yIdx+1)));
  Zx0yhi = INTERP(x, xlo, xhi, corner1, corner2);
  return (INTERP(y,ylo,yhi,Zx0ylo,Zx0yhi));
}
