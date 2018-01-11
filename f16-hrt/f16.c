/*
 * f16.c
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

#include "f16.h"
#include "f16_private.h"

/* Block signals (auto storage) */
B_f16_T f16_B;

/* Continuous states */
X_f16_T f16_X;

/* Block states (auto storage) */
DW_f16_T f16_DW;

/* External inputs (root inport signals with auto storage) */
ExtU_f16_T f16_U;

/* External outputs (root outports fed by signals with auto storage) */
ExtY_f16_T f16_Y;

/* Real-time model */
RT_MODEL_f16_T f16_M_;
RT_MODEL_f16_T *const f16_M = &f16_M_;
real_T look2_binlxpw(real_T u0, real_T u1, const real_T bp0[], const real_T bp1[],
                     const real_T table[], const uint32_T maxIndex[], uint32_T
                     stride)
{
  real_T frac;
  uint32_T bpIndices[2];
  real_T fractions[2];
  real_T yL_1d;
  uint32_T iRght;
  uint32_T bpIdx;
  uint32_T iLeft;

  /* Lookup 2-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex[0U]]) {
    /* Binary Search */
    bpIdx = maxIndex[0U] >> 1U;
    iLeft = 0U;
    iRght = maxIndex[0U];
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex[0U] - 1U;
    frac = (u0 - bp0[maxIndex[0U] - 1U]) / (bp0[maxIndex[0U]] - bp0[maxIndex[0U]
      - 1U]);
  }

  fractions[0U] = frac;
  bpIndices[0U] = iLeft;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u1 <= bp1[0U]) {
    iLeft = 0U;
    frac = (u1 - bp1[0U]) / (bp1[1U] - bp1[0U]);
  } else if (u1 < bp1[maxIndex[1U]]) {
    /* Binary Search */
    bpIdx = maxIndex[1U] >> 1U;
    iLeft = 0U;
    iRght = maxIndex[1U];
    while (iRght - iLeft > 1U) {
      if (u1 < bp1[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u1 - bp1[iLeft]) / (bp1[iLeft + 1U] - bp1[iLeft]);
  } else {
    iLeft = maxIndex[1U] - 1U;
    frac = (u1 - bp1[maxIndex[1U] - 1U]) / (bp1[maxIndex[1U]] - bp1[maxIndex[1U]
      - 1U]);
  }

  /* Interpolation 2-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  bpIdx = iLeft * stride + bpIndices[0U];
  yL_1d = (table[bpIdx + 1U] - table[bpIdx]) * fractions[0U] + table[bpIdx];
  bpIdx += stride;
  return (((table[bpIdx + 1U] - table[bpIdx]) * fractions[0U] + table[bpIdx]) -
          yL_1d) * frac + yL_1d;
}

real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T bpIdx;

  /* Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  /* Interpolation 1-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  return (table[iLeft + 1U] - table[iLeft]) * frac + table[iLeft];
}

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = (ODE4_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 13;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  f16_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  f16_step();
  f16_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  f16_step();
  f16_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  f16_step();
  f16_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(u0_0, u1_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = (rtNaN);
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Model step function */
void f16_step(void)
{
  real_T alf2CL[21];
  real_T y[189];
  real_T x[9];
  int32_T low_i;
  int32_T low_ip1;
  int32_T high_i;
  int32_T mid_i;
  real_T b_y[21];
  real_T rtb_Switch2_i;
  real_T rtb_clp_inv;
  real_T rtb_Product3_m;
  real_T rtb_Product1_lf;
  real_T rtb_Product1_pa;
  real_T rtb_Integrator1_idx_0;
  real_T rtb_sincos_o1_idx_0;
  real_T rtb_sincos_o1_idx_1;
  real_T rtb_sincos_o1_idx_2;
  real_T rtb_sincos_o2_idx_0;
  real_T rtb_sincos_o2_idx_1;
  int32_T exitg1;
  int32_T exitg2;
  if (rtmIsMajorTimeStep(f16_M)) {
    /* set solver stop time */
    if (!(f16_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&f16_M->solverInfo, ((f16_M->Timing.clockTickH0 + 1)
        * f16_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&f16_M->solverInfo, ((f16_M->Timing.clockTick0 + 1) *
        f16_M->Timing.stepSize0 + f16_M->Timing.clockTickH0 *
        f16_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(f16_M)) {
    f16_M->Timing.t[0] = rtsiGetT(&f16_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(f16_M)) {
    /* Gain: '<S31>/1//2' incorporates:
     *  Constant: '<S15>/Initial Euler Angles'
     */
    rtb_sincos_o1_idx_0 = f16_P.u2_Gain * f16_P.InitialEulerAngles_Value[2];
    rtb_sincos_o1_idx_1 = f16_P.u2_Gain * f16_P.InitialEulerAngles_Value[1];
    rtb_sincos_o1_idx_2 = f16_P.u2_Gain * f16_P.InitialEulerAngles_Value[0];

    /* Trigonometry: '<S31>/sincos' */
    rtb_sincos_o2_idx_0 = cos(rtb_sincos_o1_idx_0);
    rtb_sincos_o1_idx_0 = sin(rtb_sincos_o1_idx_0);
    rtb_sincos_o2_idx_1 = cos(rtb_sincos_o1_idx_1);
    rtb_sincos_o1_idx_1 = sin(rtb_sincos_o1_idx_1);
    rtb_Switch2_i = cos(rtb_sincos_o1_idx_2);
    rtb_sincos_o1_idx_2 = sin(rtb_sincos_o1_idx_2);

    /* Fcn: '<S31>/q0' incorporates:
     *  Trigonometry: '<S31>/sincos'
     */
    f16_B.q0 = rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_1 * rtb_Switch2_i +
      rtb_sincos_o1_idx_0 * rtb_sincos_o1_idx_1 * rtb_sincos_o1_idx_2;

    /* Fcn: '<S31>/q1' incorporates:
     *  Trigonometry: '<S31>/sincos'
     */
    f16_B.q1 = rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_1 * rtb_sincos_o1_idx_2 -
      rtb_sincos_o1_idx_0 * rtb_sincos_o1_idx_1 * rtb_Switch2_i;

    /* Fcn: '<S31>/q2' incorporates:
     *  Trigonometry: '<S31>/sincos'
     */
    f16_B.q2 = rtb_sincos_o2_idx_0 * rtb_sincos_o1_idx_1 * rtb_Switch2_i +
      rtb_sincos_o1_idx_0 * rtb_sincos_o2_idx_1 * rtb_sincos_o1_idx_2;

    /* Fcn: '<S31>/q3' incorporates:
     *  Trigonometry: '<S31>/sincos'
     */
    f16_B.q3 = rtb_sincos_o1_idx_0 * rtb_sincos_o2_idx_1 * rtb_Switch2_i -
      rtb_sincos_o2_idx_0 * rtb_sincos_o1_idx_1 * rtb_sincos_o1_idx_2;
  }

  /* Integrator: '<S15>/q0 q1 q2 q3' */
  if (f16_DW.q0q1q2q3_IWORK != 0) {
    f16_X.q0q1q2q3_CSTATE[0] = f16_B.q0;
    f16_X.q0q1q2q3_CSTATE[1] = f16_B.q1;
    f16_X.q0q1q2q3_CSTATE[2] = f16_B.q2;
    f16_X.q0q1q2q3_CSTATE[3] = f16_B.q3;
  }

  /* Sqrt: '<S34>/sqrt' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  Product: '<S35>/Product'
   *  Product: '<S35>/Product1'
   *  Product: '<S35>/Product2'
   *  Product: '<S35>/Product3'
   *  Sum: '<S35>/Sum'
   */
  rtb_sincos_o2_idx_0 = sqrt(((f16_X.q0q1q2q3_CSTATE[0] * f16_X.q0q1q2q3_CSTATE
    [0] + f16_X.q0q1q2q3_CSTATE[1] * f16_X.q0q1q2q3_CSTATE[1]) +
    f16_X.q0q1q2q3_CSTATE[2] * f16_X.q0q1q2q3_CSTATE[2]) +
    f16_X.q0q1q2q3_CSTATE[3] * f16_X.q0q1q2q3_CSTATE[3]);

  /* Product: '<S33>/Product' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  rtb_Switch2_i = f16_X.q0q1q2q3_CSTATE[0] / rtb_sincos_o2_idx_0;

  /* Product: '<S33>/Product1' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  rtb_clp_inv = f16_X.q0q1q2q3_CSTATE[1] / rtb_sincos_o2_idx_0;

  /* Product: '<S33>/Product2' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  rtb_Product3_m = f16_X.q0q1q2q3_CSTATE[2] / rtb_sincos_o2_idx_0;

  /* Product: '<S33>/Product3' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  rtb_sincos_o2_idx_0 = f16_X.q0q1q2q3_CSTATE[3] / rtb_sincos_o2_idx_0;

  /* Trigonometry: '<S30>/Trigonometric Function1' incorporates:
   *  Fcn: '<S30>/fcn1'
   *  Fcn: '<S30>/fcn2'
   */
  rtb_Integrator1_idx_0 = rt_atan2d_snf((rtb_clp_inv * rtb_Product3_m +
    rtb_Switch2_i * rtb_sincos_o2_idx_0) * 2.0, ((rtb_Switch2_i * rtb_Switch2_i
    + rtb_clp_inv * rtb_clp_inv) - rtb_Product3_m * rtb_Product3_m) -
    rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_0);

  /* Fcn: '<S30>/fcn3' */
  rtb_sincos_o1_idx_0 = (rtb_clp_inv * rtb_sincos_o2_idx_0 - rtb_Switch2_i *
    rtb_Product3_m) * -2.0;

  /* Trigonometry: '<S30>/trigFcn' */
  if (rtb_sincos_o1_idx_0 > 1.0) {
    rtb_sincos_o1_idx_0 = 1.0;
  } else {
    if (rtb_sincos_o1_idx_0 < -1.0) {
      rtb_sincos_o1_idx_0 = -1.0;
    }
  }

  rtb_sincos_o1_idx_0 = asin(rtb_sincos_o1_idx_0);

  /* End of Trigonometry: '<S30>/trigFcn' */

  /* Fcn: '<S30>/fcn4' */
  rtb_sincos_o2_idx_1 = (rtb_Product3_m * rtb_sincos_o2_idx_0 + rtb_Switch2_i *
    rtb_clp_inv) * 2.0;

  /* Fcn: '<S30>/fcn5' */
  rtb_Switch2_i = ((rtb_Switch2_i * rtb_Switch2_i - rtb_clp_inv * rtb_clp_inv) -
                   rtb_Product3_m * rtb_Product3_m) + rtb_sincos_o2_idx_0 *
    rtb_sincos_o2_idx_0;

  /* Trigonometry: '<S30>/Trigonometric Function3' */
  rtb_Switch2_i = rt_atan2d_snf(rtb_sincos_o2_idx_1, rtb_Switch2_i);

  /* Outport: '<Root>/attitude_g' */
  f16_Y.attitude_g[0] = rtb_Switch2_i;
  f16_Y.attitude_g[1] = rtb_sincos_o1_idx_0;
  f16_Y.attitude_g[2] = rtb_Integrator1_idx_0;
  if (rtmIsMajorTimeStep(f16_M)) {
    /* Constant: '<S1>/xg_0' */
    f16_B.xg_0[0] = f16_P.xg_0_Value[0];
    f16_B.xg_0[1] = f16_P.xg_0_Value[1];
    f16_B.xg_0[2] = f16_P.xg_0_Value[2];
  }

  /* Integrator: '<S1>/Integrator2' */
  if (f16_DW.Integrator2_IWORK != 0) {
    f16_X.Integrator2_CSTATE[0] = f16_B.xg_0[0];
    f16_X.Integrator2_CSTATE[1] = f16_B.xg_0[1];
    f16_X.Integrator2_CSTATE[2] = f16_B.xg_0[2];
  }

  f16_B.Integrator2[0] = f16_X.Integrator2_CSTATE[0];

  /* Outport: '<Root>/Xg' */
  f16_Y.Xg[0] = f16_B.Integrator2[0];

  /* Integrator: '<S1>/Integrator2' */
  f16_B.Integrator2[1] = f16_X.Integrator2_CSTATE[1];

  /* Outport: '<Root>/Xg' */
  f16_Y.Xg[1] = f16_B.Integrator2[1];

  /* Integrator: '<S1>/Integrator2' */
  f16_B.Integrator2[2] = f16_X.Integrator2_CSTATE[2];

  /* Outport: '<Root>/Xg' */
  f16_Y.Xg[2] = f16_B.Integrator2[2];

  /* UnitConversion: '<S2>/Unit Conversion' */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  f16_B.UnitConversion[0] = 57.295779513082323 * rtb_Switch2_i;
  f16_B.UnitConversion[1] = 57.295779513082323 * rtb_sincos_o1_idx_0;
  f16_B.UnitConversion[2] = 57.295779513082323 * rtb_Integrator1_idx_0;
  if (rtmIsMajorTimeStep(f16_M)) {
    /* Constant: '<S1>/Vk0' */
    f16_B.Vk0[0] = f16_P.Vk0_Value[0];
    f16_B.Vk0[1] = f16_P.Vk0_Value[1];
    f16_B.Vk0[2] = f16_P.Vk0_Value[2];
  }

  /* Integrator: '<S1>/Integrator1' */
  if (f16_DW.Integrator1_IWORK != 0) {
    f16_X.Integrator1_CSTATE[0] = f16_B.Vk0[0];
    f16_X.Integrator1_CSTATE[1] = f16_B.Vk0[1];
    f16_X.Integrator1_CSTATE[2] = f16_B.Vk0[2];
  }

  /* Sqrt: '<S28>/sqrt' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  Product: '<S29>/Product'
   *  Product: '<S29>/Product1'
   *  Product: '<S29>/Product2'
   *  Product: '<S29>/Product3'
   *  Sum: '<S29>/Sum'
   *  UnaryMinus: '<S12>/Unary Minus'
   *  UnaryMinus: '<S12>/Unary Minus1'
   *  UnaryMinus: '<S12>/Unary Minus2'
   */
  rtb_clp_inv = sqrt(((f16_X.q0q1q2q3_CSTATE[0] * f16_X.q0q1q2q3_CSTATE[0] +
                       -f16_X.q0q1q2q3_CSTATE[1] * -f16_X.q0q1q2q3_CSTATE[1]) +
                      -f16_X.q0q1q2q3_CSTATE[2] * -f16_X.q0q1q2q3_CSTATE[2]) +
                     -f16_X.q0q1q2q3_CSTATE[3] * -f16_X.q0q1q2q3_CSTATE[3]);

  /* Product: '<S24>/Product2' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus1'
   */
  rtb_Switch2_i = -f16_X.q0q1q2q3_CSTATE[2] / rtb_clp_inv;

  /* Product: '<S24>/Product3' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus2'
   */
  rtb_sincos_o2_idx_0 = -f16_X.q0q1q2q3_CSTATE[3] / rtb_clp_inv;

  /* Product: '<S24>/Product1' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus'
   */
  rtb_sincos_o2_idx_1 = -f16_X.q0q1q2q3_CSTATE[1] / rtb_clp_inv;

  /* Product: '<S24>/Product' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  rtb_clp_inv = f16_X.q0q1q2q3_CSTATE[0] / rtb_clp_inv;

  /* Sum: '<S25>/Sum' incorporates:
   *  Constant: '<S25>/Constant'
   *  Gain: '<S25>/Gain'
   *  Gain: '<S25>/Gain1'
   *  Gain: '<S25>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S25>/Product'
   *  Product: '<S25>/Product1'
   *  Product: '<S25>/Product2'
   *  Product: '<S25>/Product3'
   *  Product: '<S25>/Product4'
   *  Product: '<S25>/Product5'
   *  Product: '<S25>/Product6'
   *  Product: '<S25>/Product7'
   *  Product: '<S25>/Product8'
   *  Sum: '<S25>/Sum1'
   *  Sum: '<S25>/Sum2'
   *  Sum: '<S25>/Sum3'
   */
  f16_B.Sum = (((f16_P.Constant_Value - rtb_Switch2_i * rtb_Switch2_i) -
                rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_0) * f16_P.Gain2_Gain *
               f16_X.Integrator1_CSTATE[0] + (rtb_sincos_o2_idx_1 *
    rtb_Switch2_i + rtb_clp_inv * rtb_sincos_o2_idx_0) * f16_P.Gain_Gain *
               f16_X.Integrator1_CSTATE[1]) + (rtb_sincos_o2_idx_1 *
    rtb_sincos_o2_idx_0 - rtb_clp_inv * rtb_Switch2_i) * f16_P.Gain1_Gain *
    f16_X.Integrator1_CSTATE[2];

  /* Sum: '<S26>/Sum' incorporates:
   *  Constant: '<S26>/Constant'
   *  Gain: '<S26>/Gain'
   *  Gain: '<S26>/Gain1'
   *  Gain: '<S26>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S26>/Product'
   *  Product: '<S26>/Product1'
   *  Product: '<S26>/Product2'
   *  Product: '<S26>/Product3'
   *  Product: '<S26>/Product4'
   *  Product: '<S26>/Product5'
   *  Product: '<S26>/Product6'
   *  Product: '<S26>/Product7'
   *  Product: '<S26>/Product8'
   *  Sum: '<S26>/Sum1'
   *  Sum: '<S26>/Sum2'
   *  Sum: '<S26>/Sum3'
   */
  f16_B.Sum_n = (((f16_P.Constant_Value_f - rtb_sincos_o2_idx_1 *
                   rtb_sincos_o2_idx_1) - rtb_sincos_o2_idx_0 *
                  rtb_sincos_o2_idx_0) * f16_P.Gain2_Gain_e *
                 f16_X.Integrator1_CSTATE[1] + (rtb_sincos_o2_idx_1 *
    rtb_Switch2_i - rtb_clp_inv * rtb_sincos_o2_idx_0) * f16_P.Gain_Gain_h *
                 f16_X.Integrator1_CSTATE[0]) + (rtb_clp_inv *
    rtb_sincos_o2_idx_1 + rtb_Switch2_i * rtb_sincos_o2_idx_0) *
    f16_P.Gain1_Gain_b * f16_X.Integrator1_CSTATE[2];

  /* Sum: '<S27>/Sum' incorporates:
   *  Constant: '<S27>/Constant'
   *  Gain: '<S27>/Gain'
   *  Gain: '<S27>/Gain1'
   *  Gain: '<S27>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S27>/Product'
   *  Product: '<S27>/Product1'
   *  Product: '<S27>/Product2'
   *  Product: '<S27>/Product3'
   *  Product: '<S27>/Product4'
   *  Product: '<S27>/Product5'
   *  Product: '<S27>/Product6'
   *  Product: '<S27>/Product7'
   *  Product: '<S27>/Product8'
   *  Sum: '<S27>/Sum1'
   *  Sum: '<S27>/Sum2'
   *  Sum: '<S27>/Sum3'
   */
  f16_B.Sum_nr = ((rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_0 + rtb_clp_inv *
                   rtb_Switch2_i) * f16_P.Gain_Gain_b *
                  f16_X.Integrator1_CSTATE[0] + (rtb_Switch2_i *
    rtb_sincos_o2_idx_0 - rtb_clp_inv * rtb_sincos_o2_idx_1) *
                  f16_P.Gain1_Gain_p * f16_X.Integrator1_CSTATE[1]) +
    ((f16_P.Constant_Value_o - rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_1) -
     rtb_Switch2_i * rtb_Switch2_i) * f16_P.Gain2_Gain_eb *
    f16_X.Integrator1_CSTATE[2];

  /* Gain: '<Root>/Gain4' incorporates:
   *  DotProduct: '<Root>/Dot Product'
   *  SignalConversion: '<Root>/TmpSignal ConversionAtDot ProductInport1'
   *  Sqrt: '<Root>/Sqrt'
   */
  rtb_sincos_o2_idx_1 = sqrt((f16_B.Sum * f16_B.Sum + f16_B.Sum_n * f16_B.Sum_n)
    + f16_B.Sum_nr * f16_B.Sum_nr) * f16_P.Gain4_Gain;
  if (rtmIsMajorTimeStep(f16_M)) {
    /* Product: '<S95>/Product1' incorporates:
     *  Constant: '<S7>/Constant'
     */
    f16_B.Product1 = f16_P.Constant_Value_h * f16_P.Constant_Value_h;

    /* Product: '<S95>/Product2' incorporates:
     *  Constant: '<S7>/Constant'
     */
    f16_B.Product2 = f16_P.Constant_Value_h * f16_P.Constant_Value_h;
  }

  /* Gain: '<Root>/Gain' */
  rtb_Switch2_i = f16_P.Gain_Gain_l * f16_B.Integrator2[2];

  /* Saturate: '<S52>/Limit  altitude  to troposhere' */
  if (rtb_Switch2_i > f16_P.Limitaltitudetotroposhere_Upper) {
    rtb_sincos_o2_idx_0 = f16_P.Limitaltitudetotroposhere_Upper;
  } else if (rtb_Switch2_i < f16_P.Limitaltitudetotroposhere_Lower) {
    rtb_sincos_o2_idx_0 = f16_P.Limitaltitudetotroposhere_Lower;
  } else {
    rtb_sincos_o2_idx_0 = rtb_Switch2_i;
  }

  /* End of Saturate: '<S52>/Limit  altitude  to troposhere' */

  /* Sum: '<S52>/Sum1' incorporates:
   *  Constant: '<S52>/Sea Level  Temperature'
   *  Gain: '<S52>/Lapse Rate'
   */
  rtb_sincos_o2_idx_0 = f16_P.SeaLevelTemperature_Value - f16_P.LapseRate_Gain *
    rtb_sincos_o2_idx_0;

  /* Product: '<S54>/Product1' incorporates:
   *  Gain: '<S52>/gamma*R'
   *  Product: '<S95>/Product'
   *  Sqrt: '<S52>/a'
   *  Sqrt: '<S54>/vt'
   *  Sum: '<S95>/Sum'
   */
  rtb_Product1_pa = sqrt((rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_1 +
    f16_B.Product1) + f16_B.Product2) / sqrt(f16_P.gammaR_Gain *
    rtb_sincos_o2_idx_0);
  if (rtmIsMajorTimeStep(f16_M)) {
    /* Product: '<S55>/Product1' incorporates:
     *  Constant: '<S7>/Constant'
     */
    f16_B.Product1_a = f16_P.Constant_Value_h * f16_P.Constant_Value_h;

    /* Product: '<S55>/Product2' incorporates:
     *  Constant: '<S7>/Constant'
     */
    f16_B.Product2_f = f16_P.Constant_Value_h * f16_P.Constant_Value_h;
  }

  /* Gain: '<S52>/1//T0' */
  rtb_Integrator1_idx_0 = f16_P.uT0_Gain * rtb_sincos_o2_idx_0;

  /* Sum: '<S52>/Sum' incorporates:
   *  Constant: '<S52>/Altitude of Troposphere'
   */
  rtb_sincos_o1_idx_0 = f16_P.AltitudeofTroposphere_Value - rtb_Switch2_i;

  /* Math: '<S52>/(T//T0)^(g//LR) ' incorporates:
   *  Constant: '<S52>/Constant'
   */
  if ((rtb_Integrator1_idx_0 < 0.0) && (f16_P.Constant_Value_p > floor
       (f16_P.Constant_Value_p))) {
    rtb_clp_inv = -rt_powd_snf(-rtb_Integrator1_idx_0, f16_P.Constant_Value_p);
  } else {
    rtb_clp_inv = rt_powd_snf(rtb_Integrator1_idx_0, f16_P.Constant_Value_p);
  }

  /* End of Math: '<S52>/(T//T0)^(g//LR) ' */

  /* Saturate: '<S52>/Limit  altitude  to Stratosphere' */
  if (rtb_sincos_o1_idx_0 > f16_P.LimitaltitudetoStratosphere_Upp) {
    rtb_sincos_o1_idx_0 = f16_P.LimitaltitudetoStratosphere_Upp;
  } else {
    if (rtb_sincos_o1_idx_0 < f16_P.LimitaltitudetoStratosphere_Low) {
      rtb_sincos_o1_idx_0 = f16_P.LimitaltitudetoStratosphere_Low;
    }
  }

  /* End of Saturate: '<S52>/Limit  altitude  to Stratosphere' */

  /* Product: '<S7>/Product' incorporates:
   *  Constant: '<S7>/Constant1'
   *  Gain: '<S51>/1//2rhoV^2'
   *  Gain: '<S52>/g//R'
   *  Gain: '<S52>/rho0'
   *  Math: '<S52>/Stratosphere Model'
   *  Product: '<S51>/Product2'
   *  Product: '<S52>/Product'
   *  Product: '<S52>/Product1'
   *  Product: '<S52>/Product3'
   *  Product: '<S55>/Product'
   *  Sum: '<S55>/Sum'
   *
   * About '<S52>/Stratosphere Model':
   *  Operator: exp
   */
  rtb_Product3_m = ((rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_1 +
                     f16_B.Product1_a) + f16_B.Product2_f) * (rtb_clp_inv /
    rtb_Integrator1_idx_0 * f16_P.rho0_Gain * exp(1.0 / rtb_sincos_o2_idx_0 *
    (f16_P.gR_Gain * rtb_sincos_o1_idx_0))) * f16_P.u2rhoV2_Gain * f16_P.S;

  /* Product: '<S98>/Divide' incorporates:
   *  Constant: '<S10>/alf_up_lim'
   *  Constant: '<S98>/W'
   *  Lookup2D: '<S98>/LIFT'
   *  Product: '<S98>/Lmax'
   */
  rtb_Integrator1_idx_0 = rt_Lookup2D_Normal(f16_P.MACH, 9, f16_P.ALPHA, 21,
    f16_P.LIFT, rtb_Product1_pa, f16_P.alf_up_lim_Value) * rtb_Product3_m /
    (f16_P.m0 * f16_P.g);

  /* Switch: '<S96>/Switch2' incorporates:
   *  Integrator: '<S10>/Integrator'
   *  RelationalOperator: '<S96>/LowerRelop1'
   */
  if (!(f16_X.Integrator_CSTATE > rtb_Integrator1_idx_0)) {
    /* Product: '<S97>/Divide' incorporates:
     *  Constant: '<S10>/alf_lo_lim'
     *  Constant: '<S97>/W'
     *  Lookup2D: '<S97>/LIFT'
     *  Product: '<S97>/Lmax'
     */
    rtb_Integrator1_idx_0 = rt_Lookup2D_Normal(f16_P.MACH, 9, f16_P.ALPHA, 21,
      f16_P.LIFT, rtb_Product1_pa, f16_P.alf_lo_lim_Value) * rtb_Product3_m /
      (f16_P.m0 * f16_P.g);

    /* Switch: '<S96>/Switch' incorporates:
     *  RelationalOperator: '<S96>/UpperRelop'
     */
    if (!(f16_X.Integrator_CSTATE < rtb_Integrator1_idx_0)) {
      rtb_Integrator1_idx_0 = f16_X.Integrator_CSTATE;
    }

    /* End of Switch: '<S96>/Switch' */
  }

  /* End of Switch: '<S96>/Switch2' */

  /* Product: '<S4>/Product' incorporates:
   *  Constant: '<S4>/W'
   */
  rtb_sincos_o2_idx_1 = f16_P.m0 * f16_P.g * rtb_Integrator1_idx_0;

  /* Product: '<S4>/Divide1' */
  rtb_clp_inv = rtb_sincos_o2_idx_1 / rtb_Product3_m;

  /* MATLAB Function: '<S4>/MATLAB Function' incorporates:
   *  Constant: '<S4>/n_nc1'
   *  Constant: '<S4>/n_nc2'
   *  Constant: '<S4>/n_nc3'
   */
  /* MATLAB Function 'L_and_alpha_dynamics/MATLAB Function': '<S42>:1' */
  /* '<S42>:1:2' */
  memcpy(&y[0], &f16_P.LIFT[0], 189U * sizeof(real_T));
  memcpy(&x[0], &f16_P.MACH[0], 9U * sizeof(real_T));
  memset(&alf2CL[0], 0, 21U * sizeof(real_T));
  low_i = 1;
  do {
    exitg2 = 0;
    if (low_i < 10) {
      if (rtIsNaN(f16_P.MACH[low_i - 1])) {
        exitg2 = 1;
      } else {
        low_i++;
      }
    } else {
      if (f16_P.MACH[1] < f16_P.MACH[0]) {
        rtb_sincos_o2_idx_0 = x[0];
        x[0] = x[8];
        x[8] = rtb_sincos_o2_idx_0;
        rtb_sincos_o2_idx_0 = x[1];
        x[1] = x[7];
        x[7] = rtb_sincos_o2_idx_0;
        rtb_sincos_o2_idx_0 = x[2];
        x[2] = x[6];
        x[6] = rtb_sincos_o2_idx_0;
        rtb_sincos_o2_idx_0 = x[3];
        x[3] = x[5];
        x[5] = rtb_sincos_o2_idx_0;
        for (low_i = 0; low_i < 21; low_i++) {
          low_ip1 = low_i * 9;
          rtb_sincos_o2_idx_0 = y[low_ip1];
          y[low_ip1] = y[low_ip1 + 8];
          y[low_ip1 + 8] = rtb_sincos_o2_idx_0;
          rtb_sincos_o2_idx_0 = y[low_ip1 + 1];
          y[low_ip1 + 1] = y[low_ip1 + 7];
          y[low_ip1 + 7] = rtb_sincos_o2_idx_0;
          rtb_sincos_o2_idx_0 = y[low_ip1 + 2];
          y[low_ip1 + 2] = y[low_ip1 + 6];
          y[low_ip1 + 6] = rtb_sincos_o2_idx_0;
          rtb_sincos_o2_idx_0 = y[low_ip1 + 3];
          y[low_ip1 + 3] = y[low_ip1 + 5];
          y[low_ip1 + 5] = rtb_sincos_o2_idx_0;
        }
      }

      if (rtIsNaN(rtb_Product1_pa)) {
        for (low_i = 0; low_i < 21; low_i++) {
          alf2CL[low_i] = (rtNaN);
        }
      } else if (rtb_Product1_pa > x[8]) {
        rtb_sincos_o2_idx_0 = (rtb_Product1_pa - x[8]) / (x[8] - x[7]);
        for (low_i = 0; low_i < 21; low_i++) {
          alf2CL[low_i] = (y[low_i * 9 + 8] - y[low_i * 9 + 7]) *
            rtb_sincos_o2_idx_0 + y[low_i * 9 + 8];
        }
      } else if (rtb_Product1_pa < x[0]) {
        rtb_sincos_o2_idx_0 = (rtb_Product1_pa - x[0]) / (x[1] - x[0]);
        for (low_i = 0; low_i < 21; low_i++) {
          alf2CL[low_i] = (y[low_i * 9 + 1] - y[low_i * 9]) *
            rtb_sincos_o2_idx_0 + y[low_i * 9];
        }
      } else {
        low_i = 1;
        low_ip1 = 2;
        high_i = 9;
        while (high_i > low_ip1) {
          mid_i = (low_i + high_i) >> 1;
          if (rtb_Product1_pa >= x[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }

        rtb_sincos_o2_idx_0 = (rtb_Product1_pa - x[low_i - 1]) / (x[low_i] -
          x[low_i - 1]);
        if (rtb_sincos_o2_idx_0 == 0.0) {
          for (low_ip1 = 0; low_ip1 < 21; low_ip1++) {
            alf2CL[low_ip1] = y[(low_ip1 * 9 + low_i) - 1];
          }
        } else if (rtb_sincos_o2_idx_0 == 1.0) {
          for (low_ip1 = 0; low_ip1 < 21; low_ip1++) {
            alf2CL[low_ip1] = y[low_ip1 * 9 + low_i];
          }
        } else {
          for (low_ip1 = 0; low_ip1 < 21; low_ip1++) {
            if (y[(low_ip1 * 9 + low_i) - 1] == y[low_ip1 * 9 + low_i]) {
              alf2CL[low_ip1] = y[(low_ip1 * 9 + low_i) - 1];
            } else {
              alf2CL[low_ip1] = y[(low_ip1 * 9 + low_i) - 1] * (1.0 -
                rtb_sincos_o2_idx_0) + y[low_ip1 * 9 + low_i] *
                rtb_sincos_o2_idx_0;
            }
          }
        }
      }

      exitg2 = 1;
    }
  } while (exitg2 == 0);

  /* '<S42>:1:3' */
  memcpy(&b_y[0], &f16_P.ALPHA[0], 21U * sizeof(real_T));
  rtb_sincos_o2_idx_0 = 0.0;
  low_i = 1;
  do {
    exitg1 = 0;
    if (low_i < 22) {
      if (rtIsNaN(alf2CL[low_i - 1])) {
        exitg1 = 1;
      } else {
        low_i++;
      }
    } else {
      if (alf2CL[1] < alf2CL[0]) {
        for (low_i = 0; low_i < 10; low_i++) {
          rtb_sincos_o2_idx_0 = alf2CL[low_i];
          alf2CL[low_i] = alf2CL[20 - low_i];
          alf2CL[20 - low_i] = rtb_sincos_o2_idx_0;
        }

        for (low_i = 0; low_i < 10; low_i++) {
          rtb_sincos_o2_idx_0 = b_y[low_i];
          b_y[low_i] = b_y[20 - low_i];
          b_y[20 - low_i] = rtb_sincos_o2_idx_0;
        }
      }

      if (rtIsNaN(rtb_clp_inv)) {
        rtb_sincos_o2_idx_0 = (rtNaN);
      } else if (rtb_clp_inv > alf2CL[20]) {
        rtb_sincos_o2_idx_0 = (rtb_clp_inv - alf2CL[20]) / (alf2CL[20] - alf2CL
          [19]) * (b_y[20] - b_y[19]) + b_y[20];
      } else if (rtb_clp_inv < alf2CL[0]) {
        rtb_sincos_o2_idx_0 = (rtb_clp_inv - alf2CL[0]) / (alf2CL[1] - alf2CL[0])
          * (b_y[1] - b_y[0]) + b_y[0];
      } else {
        low_i = 1;
        low_ip1 = 2;
        high_i = 21;
        while (high_i > low_ip1) {
          mid_i = (low_i + high_i) >> 1;
          if (rtb_clp_inv >= alf2CL[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }

        rtb_clp_inv = (rtb_clp_inv - alf2CL[low_i - 1]) / (alf2CL[low_i] -
          alf2CL[low_i - 1]);
        if (rtb_clp_inv == 0.0) {
          rtb_sincos_o2_idx_0 = b_y[low_i - 1];
        } else if (rtb_clp_inv == 1.0) {
          rtb_sincos_o2_idx_0 = b_y[low_i];
        } else if (b_y[low_i - 1] == b_y[low_i]) {
          rtb_sincos_o2_idx_0 = b_y[low_i - 1];
        } else {
          rtb_sincos_o2_idx_0 = (1.0 - rtb_clp_inv) * b_y[low_i - 1] +
            rtb_clp_inv * b_y[low_i];
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /* UnitConversion: '<S41>/Unit Conversion' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_clp_inv = 0.017453292519943295 * rtb_sincos_o2_idx_0;

  /* Product: '<S3>/Product' incorporates:
   *  Lookup2D: '<S3>/DRAG'
   *  UnitConversion: '<S40>/Unit Conversion'
   */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  f16_B.Product = rt_Lookup2D_Normal(f16_P.MACH, 9, f16_P.ALPHA, 21, f16_P.DRAG,
    rtb_Product1_pa, 57.295779513082323 * rtb_clp_inv) * rtb_Product3_m;
  if (rtmIsMajorTimeStep(f16_M)) {
  }

  /* Lookup_n-D: '<S100>/THRUST_AB_MAX_RATIO' */
  rtb_Switch2_i = look2_binlxpw(rtb_Product1_pa, rtb_Switch2_i,
    f16_P.THRUST_MA_AB_IDX, f16_P.THRUST_H_IDX, f16_P.THRUST_AB_TAB,
    f16_P.THRUST_AB_MAX_RATIO_maxIndex, 14U);

  /* Product: '<S100>/Divide' incorporates:
   *  Constant: '<S100>/Constant4'
   *  Constant: '<S100>/Constant5'
   *  Constant: '<S100>/W'
   *  Product: '<S100>/Product'
   *  Sum: '<S100>/Add'
   */
  rtb_Switch2_i = (rtb_Switch2_i * f16_P.THRUST_AB * f16_P.g - f16_B.Product) /
    (f16_P.m0 * f16_P.g);

  /* Switch: '<S99>/Switch2' incorporates:
   *  Integrator: '<S11>/Integrator'
   *  RelationalOperator: '<S99>/LowerRelop1'
   */
  if (!(f16_X.Integrator_CSTATE_b > rtb_Switch2_i)) {
    /* Switch: '<S99>/Switch' incorporates:
     *  Constant: '<S11>/Constant6'
     *  RelationalOperator: '<S99>/UpperRelop'
     */
    if (f16_X.Integrator_CSTATE_b < f16_P.Constant6_Value) {
      rtb_Switch2_i = f16_P.Constant6_Value;
    } else {
      rtb_Switch2_i = f16_X.Integrator_CSTATE_b;
    }

    /* End of Switch: '<S99>/Switch' */
  }

  /* End of Switch: '<S99>/Switch2' */

  /* Gain: '<S6>/Gain1' */
  rtb_sincos_o1_idx_2 = f16_P.Gain1_Gain_h * rtb_sincos_o2_idx_1;

  /* Sqrt: '<S49>/sqrt' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  Product: '<S50>/Product'
   *  Product: '<S50>/Product1'
   *  Product: '<S50>/Product2'
   *  Product: '<S50>/Product3'
   *  Sum: '<S50>/Sum'
   *  UnaryMinus: '<S12>/Unary Minus'
   *  UnaryMinus: '<S12>/Unary Minus1'
   *  UnaryMinus: '<S12>/Unary Minus2'
   *  UnaryMinus: '<S43>/Unary Minus'
   *  UnaryMinus: '<S43>/Unary Minus1'
   *  UnaryMinus: '<S43>/Unary Minus2'
   */
  rtb_Product3_m = sqrt(((f16_X.q0q1q2q3_CSTATE[0] * f16_X.q0q1q2q3_CSTATE[0] +
    -(-f16_X.q0q1q2q3_CSTATE[1]) * -(-f16_X.q0q1q2q3_CSTATE[1])) +
    -(-f16_X.q0q1q2q3_CSTATE[2]) * -(-f16_X.q0q1q2q3_CSTATE[2])) +
                        -(-f16_X.q0q1q2q3_CSTATE[3]) * -(-f16_X.q0q1q2q3_CSTATE
    [3]));

  /* Product: '<S45>/Product2' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus1'
   *  UnaryMinus: '<S43>/Unary Minus1'
   */
  rtb_sincos_o2_idx_1 = -(-f16_X.q0q1q2q3_CSTATE[2]) / rtb_Product3_m;

  /* Product: '<S45>/Product3' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus2'
   *  UnaryMinus: '<S43>/Unary Minus2'
   */
  rtb_Product1_pa = -(-f16_X.q0q1q2q3_CSTATE[3]) / rtb_Product3_m;

  /* Product: '<S45>/Product1' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus'
   *  UnaryMinus: '<S43>/Unary Minus'
   */
  rtb_sincos_o2_idx_0 = -(-f16_X.q0q1q2q3_CSTATE[1]) / rtb_Product3_m;

  /* Product: '<S45>/Product' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  rtb_Product3_m = f16_X.q0q1q2q3_CSTATE[0] / rtb_Product3_m;

  /* Sum: '<S46>/Sum' incorporates:
   *  Constant: '<S46>/Constant'
   *  Constant: '<S6>/W in global axes'
   *  Gain: '<S46>/Gain'
   *  Gain: '<S46>/Gain1'
   *  Gain: '<S46>/Gain2'
   *  Product: '<S46>/Product'
   *  Product: '<S46>/Product1'
   *  Product: '<S46>/Product2'
   *  Product: '<S46>/Product3'
   *  Product: '<S46>/Product4'
   *  Product: '<S46>/Product5'
   *  Product: '<S46>/Product6'
   *  Product: '<S46>/Product7'
   *  Product: '<S46>/Product8'
   *  Sum: '<S46>/Sum1'
   *  Sum: '<S46>/Sum2'
   *  Sum: '<S46>/Sum3'
   */
  f16_B.Sum_e = (((f16_P.Constant_Value_c - rtb_sincos_o2_idx_1 *
                   rtb_sincos_o2_idx_1) - rtb_Product1_pa * rtb_Product1_pa) *
                 f16_P.Gain2_Gain_m * f16_P.Winglobalaxes_Value[0] +
                 (rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_1 + rtb_Product3_m *
                  rtb_Product1_pa) * f16_P.Gain_Gain_c *
                 f16_P.Winglobalaxes_Value[1]) + (rtb_sincos_o2_idx_0 *
    rtb_Product1_pa - rtb_Product3_m * rtb_sincos_o2_idx_1) * f16_P.Gain1_Gain_k
    * f16_P.Winglobalaxes_Value[2];

  /* Sum: '<S47>/Sum' incorporates:
   *  Constant: '<S47>/Constant'
   *  Constant: '<S6>/W in global axes'
   *  Gain: '<S47>/Gain'
   *  Gain: '<S47>/Gain1'
   *  Gain: '<S47>/Gain2'
   *  Product: '<S47>/Product'
   *  Product: '<S47>/Product1'
   *  Product: '<S47>/Product2'
   *  Product: '<S47>/Product3'
   *  Product: '<S47>/Product4'
   *  Product: '<S47>/Product5'
   *  Product: '<S47>/Product6'
   *  Product: '<S47>/Product7'
   *  Product: '<S47>/Product8'
   *  Sum: '<S47>/Sum1'
   *  Sum: '<S47>/Sum2'
   *  Sum: '<S47>/Sum3'
   */
  f16_B.Sum_p = (((f16_P.Constant_Value_g - rtb_sincos_o2_idx_0 *
                   rtb_sincos_o2_idx_0) - rtb_Product1_pa * rtb_Product1_pa) *
                 f16_P.Gain2_Gain_mr * f16_P.Winglobalaxes_Value[1] +
                 (rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_1 - rtb_Product3_m *
                  rtb_Product1_pa) * f16_P.Gain_Gain_j *
                 f16_P.Winglobalaxes_Value[0]) + (rtb_Product3_m *
    rtb_sincos_o2_idx_0 + rtb_sincos_o2_idx_1 * rtb_Product1_pa) *
    f16_P.Gain1_Gain_g * f16_P.Winglobalaxes_Value[2];

  /* Sum: '<S48>/Sum' incorporates:
   *  Constant: '<S48>/Constant'
   *  Constant: '<S6>/W in global axes'
   *  Gain: '<S48>/Gain'
   *  Gain: '<S48>/Gain1'
   *  Gain: '<S48>/Gain2'
   *  Product: '<S48>/Product'
   *  Product: '<S48>/Product1'
   *  Product: '<S48>/Product2'
   *  Product: '<S48>/Product3'
   *  Product: '<S48>/Product4'
   *  Product: '<S48>/Product5'
   *  Product: '<S48>/Product6'
   *  Product: '<S48>/Product7'
   *  Product: '<S48>/Product8'
   *  Sum: '<S48>/Sum1'
   *  Sum: '<S48>/Sum2'
   *  Sum: '<S48>/Sum3'
   */
  f16_B.Sum_o = ((rtb_sincos_o2_idx_0 * rtb_Product1_pa + rtb_Product3_m *
                  rtb_sincos_o2_idx_1) * f16_P.Gain_Gain_d *
                 f16_P.Winglobalaxes_Value[0] + (rtb_sincos_o2_idx_1 *
    rtb_Product1_pa - rtb_Product3_m * rtb_sincos_o2_idx_0) * f16_P.Gain1_Gain_f
                 * f16_P.Winglobalaxes_Value[1]) + ((f16_P.Constant_Value_pf -
    rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_0) - rtb_sincos_o2_idx_1 *
    rtb_sincos_o2_idx_1) * f16_P.Gain2_Gain_c * f16_P.Winglobalaxes_Value[2];

  /* Gain: '<S6>/Gain2' incorporates:
   *  Constant: '<S5>/W'
   *  Constant: '<S6>/Constant'
   *  Constant: '<S6>/Constant1'
   *  Constant: '<S6>/Constant2'
   *  Constant: '<S6>/Constant3'
   *  Constant: '<S6>/Constant5'
   *  Constant: '<S6>/Constant6'
   *  Gain: '<S6>/Gain'
   *  Gain: '<S6>/Gain1'
   *  Product: '<S5>/Product'
   *  Sum: '<S5>/Add'
   *  Sum: '<S6>/Add'
   */
  rtb_sincos_o2_idx_1 = 1.0 / f16_P.m0;
  f16_B.Gain2[0] = ((((f16_P.m0 * f16_P.g * rtb_Switch2_i + f16_B.Product) +
                      f16_P.Gain_Gain_lr * f16_B.Product) + f16_P.Gain1_Gain_h *
                     f16_P.Constant6_Value_g) + f16_B.Sum_e) *
    rtb_sincos_o2_idx_1;
  f16_B.Gain2[1] = (((f16_P.Gain_Gain_lr * f16_P.Constant3_Value +
                      f16_P.Constant1_Value) + f16_P.Gain1_Gain_h *
                     f16_P.Constant5_Value) + f16_B.Sum_p) * rtb_sincos_o2_idx_1;
  f16_B.Gain2[2] = (((f16_P.Gain_Gain_lr * f16_P.Constant2_Value +
                      f16_P.Constant_Value_k) + rtb_sincos_o1_idx_2) +
                    f16_B.Sum_o) * rtb_sincos_o2_idx_1;

  /* DotProduct: '<S1>/Dot Product' incorporates:
   *  Integrator: '<S1>/Integrator1'
   */
  rtb_sincos_o2_idx_0 = (f16_X.Integrator1_CSTATE[0] * f16_X.Integrator1_CSTATE
    [0] + f16_X.Integrator1_CSTATE[1] * f16_X.Integrator1_CSTATE[1]) +
    f16_X.Integrator1_CSTATE[2] * f16_X.Integrator1_CSTATE[2];

  /* Product: '<S1>/Divide' incorporates:
   *  DotProduct: '<S1>/Dot Product'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S36>/i x j'
   *  Product: '<S36>/j x k'
   *  Product: '<S36>/k x i'
   *  Product: '<S37>/i x k'
   *  Product: '<S37>/j x i'
   *  Product: '<S37>/k x j'
   *  Sum: '<S16>/Sum'
   */
  rtb_sincos_o1_idx_0 = (f16_X.Integrator1_CSTATE[1] * f16_B.Gain2[2] -
    f16_X.Integrator1_CSTATE[2] * f16_B.Gain2[1]) / rtb_sincos_o2_idx_0;
  rtb_sincos_o1_idx_1 = (f16_X.Integrator1_CSTATE[2] * f16_B.Gain2[0] -
    f16_X.Integrator1_CSTATE[0] * f16_B.Gain2[2]) / rtb_sincos_o2_idx_0;
  rtb_sincos_o1_idx_2 = (f16_X.Integrator1_CSTATE[0] * f16_B.Gain2[1] -
    f16_X.Integrator1_CSTATE[1] * f16_B.Gain2[0]) / rtb_sincos_o2_idx_0;

  /* Sum: '<S1>/Add' incorporates:
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S38>/i x j'
   *  Product: '<S38>/j x k'
   *  Product: '<S38>/k x i'
   *  Product: '<S39>/i x k'
   *  Product: '<S39>/j x i'
   *  Product: '<S39>/k x j'
   *  Sum: '<S17>/Sum'
   */
  f16_B.Add[0] = f16_B.Gain2[0] - (rtb_sincos_o1_idx_1 *
    f16_X.Integrator1_CSTATE[2] - rtb_sincos_o1_idx_2 *
    f16_X.Integrator1_CSTATE[1]);
  f16_B.Add[1] = f16_B.Gain2[1] - (rtb_sincos_o1_idx_2 *
    f16_X.Integrator1_CSTATE[0] - rtb_sincos_o1_idx_0 *
    f16_X.Integrator1_CSTATE[2]);
  f16_B.Add[2] = f16_B.Gain2[2] - (rtb_sincos_o1_idx_0 *
    f16_X.Integrator1_CSTATE[1] - rtb_sincos_o1_idx_1 *
    f16_X.Integrator1_CSTATE[0]);

  /* Saturate: '<S8>/p_limit' incorporates:
   *  Integrator: '<S8>/Integrator'
   */
  rtb_sincos_o2_idx_0 = -100.0 / f16_P.RAD2DEG;
  rtb_sincos_o2_idx_1 = 100.0 / f16_P.RAD2DEG;
  if (f16_X.Integrator_CSTATE_d > rtb_sincos_o2_idx_1) {
    rtb_sincos_o2_idx_0 = rtb_sincos_o2_idx_1;
  } else {
    if (!(f16_X.Integrator_CSTATE_d < rtb_sincos_o2_idx_0)) {
      rtb_sincos_o2_idx_0 = f16_X.Integrator_CSTATE_d;
    }
  }

  /* End of Saturate: '<S8>/p_limit' */

  /* Sum: '<S1>/Add1' incorporates:
   *  Constant: '<S1>/Constant'
   *  Constant: '<S1>/Constant1'
   */
  rtb_sincos_o1_idx_0 += rtb_sincos_o2_idx_0;
  rtb_sincos_o1_idx_1 += f16_P.Constant1_Value_k;
  rtb_sincos_o1_idx_2 += f16_P.Constant_Value_n;

  /* Sqrt: '<S22>/sqrt' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  Product: '<S23>/Product'
   *  Product: '<S23>/Product1'
   *  Product: '<S23>/Product2'
   *  Product: '<S23>/Product3'
   *  Sum: '<S23>/Sum'
   *  UnaryMinus: '<S12>/Unary Minus'
   *  UnaryMinus: '<S12>/Unary Minus1'
   *  UnaryMinus: '<S12>/Unary Minus2'
   */
  rtb_Product3_m = sqrt(((f16_X.q0q1q2q3_CSTATE[0] * f16_X.q0q1q2q3_CSTATE[0] +
    -f16_X.q0q1q2q3_CSTATE[1] * -f16_X.q0q1q2q3_CSTATE[1]) +
    -f16_X.q0q1q2q3_CSTATE[2] * -f16_X.q0q1q2q3_CSTATE[2]) +
                        -f16_X.q0q1q2q3_CSTATE[3] * -f16_X.q0q1q2q3_CSTATE[3]);

  /* Product: '<S18>/Product' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  rtb_Product1_pa = f16_X.q0q1q2q3_CSTATE[0] / rtb_Product3_m;

  /* Product: '<S18>/Product1' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus'
   */
  rtb_Product1_lf = -f16_X.q0q1q2q3_CSTATE[1] / rtb_Product3_m;

  /* Product: '<S18>/Product2' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus1'
   */
  rtb_sincos_o2_idx_1 = -f16_X.q0q1q2q3_CSTATE[2] / rtb_Product3_m;

  /* Product: '<S18>/Product3' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  UnaryMinus: '<S12>/Unary Minus2'
   */
  rtb_Product3_m = -f16_X.q0q1q2q3_CSTATE[3] / rtb_Product3_m;

  /* Sum: '<S19>/Sum' incorporates:
   *  Constant: '<S19>/Constant'
   *  Gain: '<S19>/Gain'
   *  Gain: '<S19>/Gain1'
   *  Gain: '<S19>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S19>/Product'
   *  Product: '<S19>/Product1'
   *  Product: '<S19>/Product2'
   *  Product: '<S19>/Product3'
   *  Product: '<S19>/Product4'
   *  Product: '<S19>/Product5'
   *  Product: '<S19>/Product6'
   *  Product: '<S19>/Product7'
   *  Product: '<S19>/Product8'
   *  Sum: '<S19>/Sum1'
   *  Sum: '<S19>/Sum2'
   *  Sum: '<S19>/Sum3'
   */
  f16_B.Sum_pg = (((f16_P.Constant_Value_m - rtb_sincos_o2_idx_1 *
                    rtb_sincos_o2_idx_1) - rtb_Product3_m * rtb_Product3_m) *
                  f16_P.Gain2_Gain_k * f16_X.Integrator1_CSTATE[0] +
                  (rtb_Product1_lf * rtb_sincos_o2_idx_1 + rtb_Product1_pa *
                   rtb_Product3_m) * f16_P.Gain_Gain_dr *
                  f16_X.Integrator1_CSTATE[1]) + (rtb_Product1_lf *
    rtb_Product3_m - rtb_Product1_pa * rtb_sincos_o2_idx_1) *
    f16_P.Gain1_Gain_py * f16_X.Integrator1_CSTATE[2];

  /* Sum: '<S20>/Sum' incorporates:
   *  Constant: '<S20>/Constant'
   *  Gain: '<S20>/Gain'
   *  Gain: '<S20>/Gain1'
   *  Gain: '<S20>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S20>/Product'
   *  Product: '<S20>/Product1'
   *  Product: '<S20>/Product2'
   *  Product: '<S20>/Product3'
   *  Product: '<S20>/Product4'
   *  Product: '<S20>/Product5'
   *  Product: '<S20>/Product6'
   *  Product: '<S20>/Product7'
   *  Product: '<S20>/Product8'
   *  Sum: '<S20>/Sum1'
   *  Sum: '<S20>/Sum2'
   *  Sum: '<S20>/Sum3'
   */
  f16_B.Sum_l = (((f16_P.Constant_Value_a - rtb_Product1_lf * rtb_Product1_lf) -
                  rtb_Product3_m * rtb_Product3_m) * f16_P.Gain2_Gain_cm *
                 f16_X.Integrator1_CSTATE[1] + (rtb_Product1_lf *
    rtb_sincos_o2_idx_1 - rtb_Product1_pa * rtb_Product3_m) * f16_P.Gain_Gain_a *
                 f16_X.Integrator1_CSTATE[0]) + (rtb_Product1_pa *
    rtb_Product1_lf + rtb_sincos_o2_idx_1 * rtb_Product3_m) *
    f16_P.Gain1_Gain_p5 * f16_X.Integrator1_CSTATE[2];

  /* Sum: '<S21>/Sum' incorporates:
   *  Constant: '<S21>/Constant'
   *  Gain: '<S21>/Gain'
   *  Gain: '<S21>/Gain1'
   *  Gain: '<S21>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S21>/Product'
   *  Product: '<S21>/Product1'
   *  Product: '<S21>/Product2'
   *  Product: '<S21>/Product3'
   *  Product: '<S21>/Product4'
   *  Product: '<S21>/Product5'
   *  Product: '<S21>/Product6'
   *  Product: '<S21>/Product7'
   *  Product: '<S21>/Product8'
   *  Sum: '<S21>/Sum1'
   *  Sum: '<S21>/Sum2'
   *  Sum: '<S21>/Sum3'
   */
  f16_B.Sum_c = ((rtb_Product1_lf * rtb_Product3_m + rtb_Product1_pa *
                  rtb_sincos_o2_idx_1) * f16_P.Gain_Gain_o *
                 f16_X.Integrator1_CSTATE[0] + (rtb_sincos_o2_idx_1 *
    rtb_Product3_m - rtb_Product1_pa * rtb_Product1_lf) * f16_P.Gain1_Gain_o *
                 f16_X.Integrator1_CSTATE[1]) + ((f16_P.Constant_Value_b -
    rtb_Product1_lf * rtb_Product1_lf) - rtb_sincos_o2_idx_1 *
    rtb_sincos_o2_idx_1) * f16_P.Gain2_Gain_h * f16_X.Integrator1_CSTATE[2];

  /* Gain: '<S32>/High Gain Quaternion Normalization' incorporates:
   *  Constant: '<S32>/Constant'
   *  DotProduct: '<S32>/Dot Product'
   *  Integrator: '<S15>/q0 q1 q2 q3'
   *  Sum: '<S32>/Sum'
   */
  rtb_Product3_m = (f16_P.Constant_Value_hd - (((f16_X.q0q1q2q3_CSTATE[0] *
    f16_X.q0q1q2q3_CSTATE[0] + f16_X.q0q1q2q3_CSTATE[1] * f16_X.q0q1q2q3_CSTATE
    [1]) + f16_X.q0q1q2q3_CSTATE[2] * f16_X.q0q1q2q3_CSTATE[2]) +
    f16_X.q0q1q2q3_CSTATE[3] * f16_X.q0q1q2q3_CSTATE[3])) *
    f16_P.HighGainQuaternionNormalization;

  /* Fcn: '<S32>/q0dot' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  f16_B.q0dot = ((f16_X.q0q1q2q3_CSTATE[1] * rtb_sincos_o1_idx_0 +
                  f16_X.q0q1q2q3_CSTATE[2] * rtb_sincos_o1_idx_1) +
                 f16_X.q0q1q2q3_CSTATE[3] * rtb_sincos_o1_idx_2) * -0.5 +
    rtb_Product3_m * f16_X.q0q1q2q3_CSTATE[0];

  /* Fcn: '<S32>/q1dot' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  f16_B.q1dot = ((f16_X.q0q1q2q3_CSTATE[0] * rtb_sincos_o1_idx_0 +
                  f16_X.q0q1q2q3_CSTATE[2] * rtb_sincos_o1_idx_2) -
                 f16_X.q0q1q2q3_CSTATE[3] * rtb_sincos_o1_idx_1) * 0.5 +
    rtb_Product3_m * f16_X.q0q1q2q3_CSTATE[1];

  /* Fcn: '<S32>/q2dot' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  f16_B.q2dot = ((f16_X.q0q1q2q3_CSTATE[0] * rtb_sincos_o1_idx_1 +
                  f16_X.q0q1q2q3_CSTATE[3] * rtb_sincos_o1_idx_0) -
                 f16_X.q0q1q2q3_CSTATE[1] * rtb_sincos_o1_idx_2) * 0.5 +
    rtb_Product3_m * f16_X.q0q1q2q3_CSTATE[2];

  /* Fcn: '<S32>/q3dot' incorporates:
   *  Integrator: '<S15>/q0 q1 q2 q3'
   */
  f16_B.q3dot = ((f16_X.q0q1q2q3_CSTATE[0] * rtb_sincos_o1_idx_2 +
                  f16_X.q0q1q2q3_CSTATE[1] * rtb_sincos_o1_idx_1) -
                 f16_X.q0q1q2q3_CSTATE[2] * rtb_sincos_o1_idx_0) * 0.5 +
    rtb_Product3_m * f16_X.q0q1q2q3_CSTATE[3];
  if (rtmIsMajorTimeStep(f16_M)) {
  }

  /* Lookup_n-D: '<S8>/clp_inv' */
  rtb_clp_inv = look1_binlxpw(rtb_clp_inv, f16_P.clp_inv_alf, f16_P.clp_inv, 11U);

  /* Product: '<S8>/Product' incorporates:
   *  Inport: '<Root>/mu_dotc'
   *  Sum: '<S8>/Sum'
   */
  f16_B.Product_n = (f16_U.mu_dotc - rtb_sincos_o2_idx_0) * rtb_clp_inv;

  /* Saturate: '<S10>/n_n_c_lim' incorporates:
   *  Inport: '<Root>/n_nc'
   */
  if (f16_U.n_nc > f16_P.n_n_c_lim_UpperSat) {
    rtb_sincos_o2_idx_0 = f16_P.n_n_c_lim_UpperSat;
  } else if (f16_U.n_nc < f16_P.n_n_c_lim_LowerSat) {
    rtb_sincos_o2_idx_0 = f16_P.n_n_c_lim_LowerSat;
  } else {
    rtb_sincos_o2_idx_0 = f16_U.n_nc;
  }

  /* End of Saturate: '<S10>/n_n_c_lim' */

  /* Gain: '<S10>/Gain' incorporates:
   *  Sum: '<S10>/Sum'
   */
  f16_B.Gain = (rtb_sincos_o2_idx_0 - rtb_Integrator1_idx_0) *
    f16_P.Gain_Gain_a5;

  /* Saturate: '<S11>/n_xc_lim' incorporates:
   *  Inport: '<Root>/n_xc'
   */
  if (f16_U.n_xc > f16_P.n_xc_lim_UpperSat) {
    rtb_sincos_o2_idx_0 = f16_P.n_xc_lim_UpperSat;
  } else if (f16_U.n_xc < f16_P.n_xc_lim_LowerSat) {
    rtb_sincos_o2_idx_0 = f16_P.n_xc_lim_LowerSat;
  } else {
    rtb_sincos_o2_idx_0 = f16_U.n_xc;
  }

  /* End of Saturate: '<S11>/n_xc_lim' */

  /* Gain: '<S11>/Gain' incorporates:
   *  Sum: '<S11>/Sum'
   */
  f16_B.Gain_e = (rtb_sincos_o2_idx_0 - rtb_Switch2_i) * f16_P.Gain_Gain_dw;
  if (rtmIsMajorTimeStep(f16_M)) {
    /* Matfile logging */
    /*rt_UpdateTXYLogVars(f16_M->rtwLogInfo, (f16_M->Timing.t));*/
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(f16_M)) {
    /* Update for Integrator: '<S15>/q0 q1 q2 q3' */
    f16_DW.q0q1q2q3_IWORK = 0;

    /* Update for Integrator: '<S1>/Integrator2' */
    f16_DW.Integrator2_IWORK = 0;

    /* Update for Integrator: '<S1>/Integrator1' */
    f16_DW.Integrator1_IWORK = 0;
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(f16_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(f16_M)!=-1) &&
          !((rtmGetTFinal(f16_M)-(((f16_M->Timing.clockTick1+
               f16_M->Timing.clockTickH1* 4294967296.0)) * 0.033333333333333333))
            > (((f16_M->Timing.clockTick1+f16_M->Timing.clockTickH1*
                 4294967296.0)) * 0.033333333333333333) * (DBL_EPSILON))) {
        rtmSetErrorStatus(f16_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&f16_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++f16_M->Timing.clockTick0)) {
      ++f16_M->Timing.clockTickH0;
    }

    f16_M->Timing.t[0] = rtsiGetSolverStopTime(&f16_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.033333333333333333s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.033333333333333333, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      f16_M->Timing.clockTick1++;
      if (!f16_M->Timing.clockTick1) {
        f16_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void f16_derivatives(void)
{
  XDot_f16_T *_rtXdot;
  _rtXdot = ((XDot_f16_T *) f16_M->derivs);

  /* Derivatives for Integrator: '<S15>/q0 q1 q2 q3' */
  _rtXdot->q0q1q2q3_CSTATE[0] = f16_B.q0dot;
  _rtXdot->q0q1q2q3_CSTATE[1] = f16_B.q1dot;
  _rtXdot->q0q1q2q3_CSTATE[2] = f16_B.q2dot;
  _rtXdot->q0q1q2q3_CSTATE[3] = f16_B.q3dot;

  /* Derivatives for Integrator: '<S1>/Integrator2' */
  _rtXdot->Integrator2_CSTATE[0] = f16_B.Sum_pg;
  _rtXdot->Integrator2_CSTATE[1] = f16_B.Sum_l;
  _rtXdot->Integrator2_CSTATE[2] = f16_B.Sum_c;

  /* Derivatives for Integrator: '<S1>/Integrator1' */
  _rtXdot->Integrator1_CSTATE[0] = f16_B.Add[0];
  _rtXdot->Integrator1_CSTATE[1] = f16_B.Add[1];
  _rtXdot->Integrator1_CSTATE[2] = f16_B.Add[2];

  /* Derivatives for Integrator: '<S10>/Integrator' */
  _rtXdot->Integrator_CSTATE = f16_B.Gain;

  /* Derivatives for Integrator: '<S11>/Integrator' */
  _rtXdot->Integrator_CSTATE_b = f16_B.Gain_e;

  /* Derivatives for Integrator: '<S8>/Integrator' */
  _rtXdot->Integrator_CSTATE_d = f16_B.Product_n;
}

/* Model initialize function */
void f16_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)f16_M, 0,
                sizeof(RT_MODEL_f16_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&f16_M->solverInfo, &f16_M->Timing.simTimeStep);
    rtsiSetTPtr(&f16_M->solverInfo, &rtmGetTPtr(f16_M));
    rtsiSetStepSizePtr(&f16_M->solverInfo, &f16_M->Timing.stepSize0);
    rtsiSetdXPtr(&f16_M->solverInfo, &f16_M->derivs);
    rtsiSetContStatesPtr(&f16_M->solverInfo, (real_T **) &f16_M->contStates);
    rtsiSetNumContStatesPtr(&f16_M->solverInfo, &f16_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&f16_M->solverInfo,
      &f16_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&f16_M->solverInfo,
      &f16_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&f16_M->solverInfo,
      &f16_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&f16_M->solverInfo, (&rtmGetErrorStatus(f16_M)));
    rtsiSetRTModelPtr(&f16_M->solverInfo, f16_M);
  }

  rtsiSetSimTimeStep(&f16_M->solverInfo, MAJOR_TIME_STEP);
  f16_M->intgData.y = f16_M->odeY;
  f16_M->intgData.f[0] = f16_M->odeF[0];
  f16_M->intgData.f[1] = f16_M->odeF[1];
  f16_M->intgData.f[2] = f16_M->odeF[2];
  f16_M->intgData.f[3] = f16_M->odeF[3];
  f16_M->contStates = ((X_f16_T *) &f16_X);
  rtsiSetSolverData(&f16_M->solverInfo, (void *)&f16_M->intgData);
  rtsiSetSolverName(&f16_M->solverInfo,"ode4");
  rtmSetTPtr(f16_M, &f16_M->Timing.tArray[0]);
  rtmSetTFinal(f16_M, -1);
  f16_M->Timing.stepSize0 = 0.033333333333333333;
  rtmSetFirstInitCond(f16_M, 1);

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    f16_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(f16_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(f16_M->rtwLogInfo, (NULL));
    rtliSetLogT(f16_M->rtwLogInfo, "tout");
    rtliSetLogX(f16_M->rtwLogInfo, "");
    rtliSetLogXFinal(f16_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(f16_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(f16_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(f16_M->rtwLogInfo, 0);
    rtliSetLogDecimation(f16_M->rtwLogInfo, 1);
    rtliSetLogY(f16_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(f16_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(f16_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &f16_B), 0,
                sizeof(B_f16_T));

  /* states (continuous) */
  {
    (void) memset((void *)&f16_X, 0,
                  sizeof(X_f16_T));
  }

  /* states (dwork) */
  (void) memset((void *)&f16_DW, 0,
                sizeof(DW_f16_T));

  /* external inputs */
  (void)memset((void *)&f16_U, 0, sizeof(ExtU_f16_T));

  /* external outputs */
  (void) memset((void *)&f16_Y, 0,
                sizeof(ExtY_f16_T));

  /* Matfile logging */
  /*rt_StartDataLoggingWithStartTime(f16_M->rtwLogInfo, 0.0, rtmGetTFinal(f16_M),
    f16_M->Timing.stepSize0, (&rtmGetErrorStatus(f16_M)));*/

  /* Start for Constant: '<S1>/xg_0' */
  f16_B.xg_0[0] = f16_P.xg_0_Value[0];
  f16_B.xg_0[1] = f16_P.xg_0_Value[1];
  f16_B.xg_0[2] = f16_P.xg_0_Value[2];

  /* Start for Constant: '<S1>/Vk0' */
  f16_B.Vk0[0] = f16_P.Vk0_Value[0];
  f16_B.Vk0[1] = f16_P.Vk0_Value[1];
  f16_B.Vk0[2] = f16_P.Vk0_Value[2];

  /* InitializeConditions for Integrator: '<S15>/q0 q1 q2 q3' incorporates:
   *  InitializeConditions for Integrator: '<S1>/Integrator2'
   */
  if (rtmIsFirstInitCond(f16_M)) {
    f16_X.q0q1q2q3_CSTATE[0] = 0.0;
    f16_X.q0q1q2q3_CSTATE[1] = 0.0;
    f16_X.q0q1q2q3_CSTATE[2] = 0.0;
    f16_X.q0q1q2q3_CSTATE[3] = 0.0;
    f16_X.Integrator2_CSTATE[0] = 0.0;
    f16_X.Integrator2_CSTATE[1] = 0.0;
    f16_X.Integrator2_CSTATE[2] = -1000.0;
  }

  f16_DW.q0q1q2q3_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S15>/q0 q1 q2 q3' */

  /* InitializeConditions for Integrator: '<S1>/Integrator2' */
  f16_DW.Integrator2_IWORK = 1;

  /* InitializeConditions for Integrator: '<S1>/Integrator1' */
  if (rtmIsFirstInitCond(f16_M)) {
    f16_X.Integrator1_CSTATE[0] = 200.0;
    f16_X.Integrator1_CSTATE[1] = 0.0;
    f16_X.Integrator1_CSTATE[2] = 0.0;
  }

  f16_DW.Integrator1_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S1>/Integrator1' */

  /* InitializeConditions for Integrator: '<S10>/Integrator' */
  f16_X.Integrator_CSTATE = f16_P.Integrator_IC;

  /* InitializeConditions for Integrator: '<S11>/Integrator' */
  f16_X.Integrator_CSTATE_b = f16_P.Integrator_IC_e;

  /* InitializeConditions for Integrator: '<S8>/Integrator' */
  f16_X.Integrator_CSTATE_d = f16_P.Integrator_IC_e2;

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(f16_M)) {
    rtmSetFirstInitCond(f16_M, 0);
  }
}

/* Model terminate function */
void f16_terminate(void)
{
  /* (no terminate code required) */
}
