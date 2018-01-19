/*
 * actor.cpp
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

#include "actor.h"
#include "actor_private.h"

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
void actorModelClass::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
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
  actor_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  this->step();
  actor_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  this->step();
  actor_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  this->step();
  actor_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = std::abs(u0);
    tmp_0 = std::abs(u1);
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
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
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

    y = atan2((real_T)u0_0, (real_T)u1_0);
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

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (std::abs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = std::floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = std::ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

real_T rt_modd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  if (u1 == 0.0) {
    y = u0;
  } else if (!((!rtIsNaN(u0)) && (!rtIsInf(u0)) && ((!rtIsNaN(u1)) && (!rtIsInf
                (u1))))) {
    y = (rtNaN);
  } else {
    tmp = u0 / u1;
    if (u1 <= std::floor(u1)) {
      y = u0 - std::floor(tmp) * u1;
    } else if (std::abs(tmp - rt_roundd_snf(tmp)) <= DBL_EPSILON * std::abs(tmp))
    {
      y = 0.0;
    } else {
      y = (tmp - std::floor(tmp)) * u1;
    }
  }

  return y;
}

/* Model step function */
void actorModelClass::step()
{
  real_T alf2CL[21];
  real_T y[189];
  real_T x[9];
  int32_T low_i;
  int32_T low_ip1;
  int32_T high_i;
  int32_T mid_i;
  real_T b_y[21];
  real_T rtb_Product1_jz;
  real_T rtb_Switch2;
  real_T rtb_jxi;
  real_T rtb_ixk;
  real_T rtb_QS;
  real_T rtb_UnitConversion_b;
  real_T rtb_sqrt;
  boolean_T rtb_Compare_g;
  real_T rtb_Sum_n_idx_0;
  real_T rtb_sincos_o1_idx_0;
  real_T rtb_sincos_o1_idx_1;
  real_T rtb_sincos_o1_idx_2;
  real_T rtb_sincos_o2_idx_0;
  real_T rtb_sincos_o2_idx_1;
  real_T u0;
  int32_T exitg1;
  int32_T exitg2;
  if (rtmIsMajorTimeStep((&actor_M))) {
    /* set solver stop time */
    if (!((&actor_M)->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&(&actor_M)->solverInfo, (((&actor_M)
        ->Timing.clockTickH0 + 1) * (&actor_M)->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&(&actor_M)->solverInfo, (((&actor_M)
        ->Timing.clockTick0 + 1) * (&actor_M)->Timing.stepSize0 + (&actor_M)
        ->Timing.clockTickH0 * (&actor_M)->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep((&actor_M))) {
    (&actor_M)->Timing.t[0] = rtsiGetT(&(&actor_M)->solverInfo);
  }

  if (rtmIsMajorTimeStep((&actor_M))) {
    /* Constant: '<S1>/Vk0' */
    actor_B.Vk0[0] = actor_P.Vk0_Value[0];
    actor_B.Vk0[1] = actor_P.Vk0_Value[1];
    actor_B.Vk0[2] = actor_P.Vk0_Value[2];

    /* Gain: '<S27>/1//2' incorporates:
     *  Constant: '<S11>/att_g0'
     */
    rtb_sincos_o1_idx_0 = actor_P.u2_Gain * actor_P.att_g0_Value[2];
    rtb_sincos_o1_idx_1 = actor_P.u2_Gain * actor_P.att_g0_Value[1];
    rtb_sincos_o1_idx_2 = actor_P.u2_Gain * actor_P.att_g0_Value[0];

    /* Trigonometry: '<S27>/sincos' */
    rtb_sincos_o2_idx_0 = std::cos(rtb_sincos_o1_idx_0);
    rtb_sincos_o1_idx_0 = std::sin(rtb_sincos_o1_idx_0);
    rtb_sincos_o2_idx_1 = std::cos(rtb_sincos_o1_idx_1);
    rtb_sincos_o1_idx_1 = std::sin(rtb_sincos_o1_idx_1);
    rtb_jxi = std::cos(rtb_sincos_o1_idx_2);
    rtb_QS = std::sin(rtb_sincos_o1_idx_2);

    /* Fcn: '<S27>/q0' incorporates:
     *  Trigonometry: '<S27>/sincos'
     */
    actor_B.q0 = rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_1 * rtb_jxi +
      rtb_sincos_o1_idx_0 * rtb_sincos_o1_idx_1 * rtb_QS;

    /* Fcn: '<S27>/q1' incorporates:
     *  Trigonometry: '<S27>/sincos'
     */
    actor_B.q1 = rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_1 * rtb_QS -
      rtb_sincos_o1_idx_0 * rtb_sincos_o1_idx_1 * rtb_jxi;

    /* Fcn: '<S27>/q2' incorporates:
     *  Trigonometry: '<S27>/sincos'
     */
    actor_B.q2 = rtb_sincos_o2_idx_0 * rtb_sincos_o1_idx_1 * rtb_jxi +
      rtb_sincos_o1_idx_0 * rtb_sincos_o2_idx_1 * rtb_QS;

    /* Fcn: '<S27>/q3' incorporates:
     *  Trigonometry: '<S27>/sincos'
     */
    actor_B.q3 = rtb_sincos_o1_idx_0 * rtb_sincos_o2_idx_1 * rtb_jxi -
      rtb_sincos_o2_idx_0 * rtb_sincos_o1_idx_1 * rtb_QS;
  }

  /* Integrator: '<S1>/Integrator1' */
  if (actor_DW.Integrator1_IWORK != 0) {
    actor_X.Vk[0] = actor_B.Vk0[0];
    actor_X.Vk[1] = actor_B.Vk0[1];
    actor_X.Vk[2] = actor_B.Vk0[2];
  }

  /* Integrator: '<S11>/q0 q1 q2 q3' */
  if (actor_DW.q0q1q2q3_IWORK != 0) {
    actor_X.q0[0] = actor_B.q0;
    actor_X.q0[1] = actor_B.q1;
    actor_X.q0[2] = actor_B.q2;
    actor_X.q0[3] = actor_B.q3;
  }

  /* Sqrt: '<S24>/sqrt' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Product: '<S25>/Product'
   *  Product: '<S25>/Product1'
   *  Product: '<S25>/Product2'
   *  Product: '<S25>/Product3'
   *  Sum: '<S25>/Sum'
   *  UnaryMinus: '<S8>/Unary Minus'
   *  UnaryMinus: '<S8>/Unary Minus1'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_Product1_jz = std::sqrt(((actor_X.q0[0] * actor_X.q0[0] + -actor_X.q0[1] *
    -actor_X.q0[1]) + -actor_X.q0[2] * -actor_X.q0[2]) + -actor_X.q0[3] *
    -actor_X.q0[3]);

  /* Product: '<S20>/Product2' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S8>/Unary Minus1'
   */
  rtb_Switch2 = -actor_X.q0[2] / rtb_Product1_jz;

  /* Product: '<S20>/Product3' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_jxi = -actor_X.q0[3] / rtb_Product1_jz;

  /* Product: '<S20>/Product1' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S8>/Unary Minus'
   */
  rtb_ixk = -actor_X.q0[1] / rtb_Product1_jz;

  /* Product: '<S20>/Product' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  rtb_Product1_jz = actor_X.q0[0] / rtb_Product1_jz;

  /* SignalConversion: '<Root>/TmpSignal ConversionAtDot ProductInport1' incorporates:
   *  Constant: '<S21>/Constant'
   *  Constant: '<S22>/Constant'
   *  Constant: '<S23>/Constant'
   *  Gain: '<S21>/Gain'
   *  Gain: '<S21>/Gain1'
   *  Gain: '<S21>/Gain2'
   *  Gain: '<S22>/Gain'
   *  Gain: '<S22>/Gain1'
   *  Gain: '<S22>/Gain2'
   *  Gain: '<S23>/Gain'
   *  Gain: '<S23>/Gain1'
   *  Gain: '<S23>/Gain2'
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
   *  Product: '<S22>/Product'
   *  Product: '<S22>/Product1'
   *  Product: '<S22>/Product2'
   *  Product: '<S22>/Product3'
   *  Product: '<S22>/Product4'
   *  Product: '<S22>/Product5'
   *  Product: '<S22>/Product6'
   *  Product: '<S22>/Product7'
   *  Product: '<S22>/Product8'
   *  Product: '<S23>/Product'
   *  Product: '<S23>/Product1'
   *  Product: '<S23>/Product2'
   *  Product: '<S23>/Product3'
   *  Product: '<S23>/Product4'
   *  Product: '<S23>/Product5'
   *  Product: '<S23>/Product6'
   *  Product: '<S23>/Product7'
   *  Product: '<S23>/Product8'
   *  Sum: '<S21>/Sum'
   *  Sum: '<S21>/Sum1'
   *  Sum: '<S21>/Sum2'
   *  Sum: '<S21>/Sum3'
   *  Sum: '<S22>/Sum'
   *  Sum: '<S22>/Sum1'
   *  Sum: '<S22>/Sum2'
   *  Sum: '<S22>/Sum3'
   *  Sum: '<S23>/Sum'
   *  Sum: '<S23>/Sum1'
   *  Sum: '<S23>/Sum2'
   *  Sum: '<S23>/Sum3'
   */
  rtb_sincos_o1_idx_0 = (((actor_P.Constant_Value_k - rtb_Switch2 * rtb_Switch2)
    - rtb_jxi * rtb_jxi) * actor_P.Gain2_Gain * actor_X.Vk[0] + (rtb_ixk *
    rtb_Switch2 + rtb_Product1_jz * rtb_jxi) * actor_P.Gain_Gain_jt *
    actor_X.Vk[1]) + (rtb_ixk * rtb_jxi - rtb_Product1_jz * rtb_Switch2) *
    actor_P.Gain1_Gain * actor_X.Vk[2];
  rtb_sincos_o1_idx_1 = (((actor_P.Constant_Value_f - rtb_ixk * rtb_ixk) -
    rtb_jxi * rtb_jxi) * actor_P.Gain2_Gain_e * actor_X.Vk[1] + (rtb_ixk *
    rtb_Switch2 - rtb_Product1_jz * rtb_jxi) * actor_P.Gain_Gain_h * actor_X.Vk
    [0]) + (rtb_Product1_jz * rtb_ixk + rtb_Switch2 * rtb_jxi) *
    actor_P.Gain1_Gain_b * actor_X.Vk[2];
  rtb_sincos_o1_idx_2 = ((rtb_ixk * rtb_jxi + rtb_Product1_jz * rtb_Switch2) *
    actor_P.Gain_Gain_b * actor_X.Vk[0] + (rtb_Switch2 * rtb_jxi -
    rtb_Product1_jz * rtb_ixk) * actor_P.Gain1_Gain_p * actor_X.Vk[1]) +
    ((actor_P.Constant_Value_o - rtb_ixk * rtb_ixk) - rtb_Switch2 * rtb_Switch2)
    * actor_P.Gain2_Gain_eb * actor_X.Vk[2];

  /* Gain: '<Root>/Gain4' incorporates:
   *  DotProduct: '<Root>/Dot Product'
   *  Sqrt: '<Root>/Sqrt'
   */
  rtb_sincos_o2_idx_1 = std::sqrt((rtb_sincos_o1_idx_0 * rtb_sincos_o1_idx_0 +
    rtb_sincos_o1_idx_1 * rtb_sincos_o1_idx_1) + rtb_sincos_o1_idx_2 *
    rtb_sincos_o1_idx_2) * actor_P.Gain4_Gain;
  if (rtmIsMajorTimeStep((&actor_M))) {
    /* Product: '<S88>/Product1' incorporates:
     *  Constant: '<S4>/Constant'
     */
    actor_B.Product1 = actor_P.Constant_Value_h * actor_P.Constant_Value_h;

    /* Product: '<S88>/Product2' incorporates:
     *  Constant: '<S4>/Constant'
     */
    actor_B.Product2 = actor_P.Constant_Value_h * actor_P.Constant_Value_h;

    /* Constant: '<S1>/xg0' */
    actor_B.xg0[0] = actor_P.xg0_Value[0];
    actor_B.xg0[1] = actor_P.xg0_Value[1];
    actor_B.xg0[2] = actor_P.xg0_Value[2];
  }

  /* Integrator: '<S1>/Integrator2' */
  if (actor_DW.Integrator2_IWORK != 0) {
    actor_X.xg[0] = actor_B.xg0[0];
    actor_X.xg[1] = actor_B.xg0[1];
    actor_X.xg[2] = actor_B.xg0[2];
  }

  /* Gain: '<Root>/Gain' incorporates:
   *  Integrator: '<S1>/Integrator2'
   */
  rtb_ixk = actor_P.Gain_Gain_l * actor_X.xg[2];

  /* Saturate: '<S45>/Limit  altitude  to troposhere' */
  if (rtb_ixk > actor_P.Limitaltitudetotroposhere_Upper) {
    rtb_jxi = actor_P.Limitaltitudetotroposhere_Upper;
  } else if (rtb_ixk < actor_P.Limitaltitudetotroposhere_Lower) {
    rtb_jxi = actor_P.Limitaltitudetotroposhere_Lower;
  } else {
    rtb_jxi = rtb_ixk;
  }

  /* End of Saturate: '<S45>/Limit  altitude  to troposhere' */

  /* Sum: '<S45>/Sum1' incorporates:
   *  Constant: '<S45>/Sea Level  Temperature'
   *  Gain: '<S45>/Lapse Rate'
   */
  rtb_jxi = actor_P.SeaLevelTemperature_Value - actor_P.LapseRate_Gain * rtb_jxi;

  /* Product: '<S47>/Product1' incorporates:
   *  Gain: '<S45>/gamma*R'
   *  Product: '<S88>/Product'
   *  Sqrt: '<S45>/a'
   *  Sqrt: '<S47>/vt'
   *  Sum: '<S88>/Sum'
   */
  rtb_Switch2 = std::sqrt((rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_1 +
    actor_B.Product1) + actor_B.Product2) / std::sqrt(actor_P.gammaR_Gain *
    rtb_jxi);
  if (rtmIsMajorTimeStep((&actor_M))) {
    /* Product: '<S48>/Product1' incorporates:
     *  Constant: '<S4>/Constant'
     */
    actor_B.Product1_a = actor_P.Constant_Value_h * actor_P.Constant_Value_h;

    /* Product: '<S48>/Product2' incorporates:
     *  Constant: '<S4>/Constant'
     */
    actor_B.Product2_f = actor_P.Constant_Value_h * actor_P.Constant_Value_h;
  }

  /* Gain: '<S45>/1//T0' */
  rtb_sincos_o2_idx_0 = actor_P.uT0_Gain * rtb_jxi;

  /* Sum: '<S45>/Sum' incorporates:
   *  Constant: '<S45>/Altitude of Troposphere'
   */
  u0 = actor_P.AltitudeofTroposphere_Value - rtb_ixk;

  /* Math: '<S45>/(T//T0)^(g//LR) ' incorporates:
   *  Constant: '<S45>/Constant'
   */
  if ((rtb_sincos_o2_idx_0 < 0.0) && (actor_P.Constant_Value_p > std::floor
       (actor_P.Constant_Value_p))) {
    rtb_QS = -rt_powd_snf(-rtb_sincos_o2_idx_0, actor_P.Constant_Value_p);
  } else {
    rtb_QS = rt_powd_snf(rtb_sincos_o2_idx_0, actor_P.Constant_Value_p);
  }

  /* End of Math: '<S45>/(T//T0)^(g//LR) ' */

  /* Saturate: '<S45>/Limit  altitude  to Stratosphere' */
  if (u0 > actor_P.LimitaltitudetoStratosphere_Upp) {
    u0 = actor_P.LimitaltitudetoStratosphere_Upp;
  } else {
    if (u0 < actor_P.LimitaltitudetoStratosphere_Low) {
      u0 = actor_P.LimitaltitudetoStratosphere_Low;
    }
  }

  /* End of Saturate: '<S45>/Limit  altitude  to Stratosphere' */

  /* Product: '<S4>/Product' incorporates:
   *  Constant: '<S4>/Constant1'
   *  Gain: '<S44>/1//2rhoV^2'
   *  Gain: '<S45>/g//R'
   *  Gain: '<S45>/rho0'
   *  Math: '<S45>/Stratosphere Model'
   *  Product: '<S44>/Product2'
   *  Product: '<S45>/Product'
   *  Product: '<S45>/Product1'
   *  Product: '<S45>/Product3'
   *  Product: '<S48>/Product'
   *  Sum: '<S48>/Sum'
   *
   * About '<S45>/Stratosphere Model':
   *  Operator: exp
   */
  rtb_QS = ((rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_1 + actor_B.Product1_a) +
            actor_B.Product2_f) * (rtb_QS / rtb_sincos_o2_idx_0 *
    actor_P.rho0_Gain * std::exp(1.0 / rtb_jxi * (actor_P.gR_Gain * u0))) *
    actor_P.u2rhoV2_Gain * actor_P.S;
  if (rtmIsMajorTimeStep((&actor_M))) {
    /* Memory: '<S7>/Memory1' */
    rtb_UnitConversion_b = actor_DW.Memory1_PreviousInput;

    /* Memory: '<S7>/Memory' */
    rtb_sqrt = actor_DW.Memory_PreviousInput;

    /* Product: '<S138>/Product1' incorporates:
     *  Trigonometry: '<S138>/Trigonometric Function'
     */
    actor_B.Product1_e = rtb_UnitConversion_b * std::sin(rtb_sqrt);

    /* Product: '<S137>/Product1' incorporates:
     *  Trigonometry: '<S137>/Trigonometric Function'
     */
    actor_B.Product1_ey = rtb_UnitConversion_b * std::sin(rtb_sqrt);
  }

  /* Product: '<S138>/Divide' incorporates:
   *  Constant: '<S131>/alf_up_lim'
   *  Constant: '<S138>/W'
   *  Lookup2D: '<S138>/LIFT'
   *  Product: '<S138>/Lmax'
   *  Sum: '<S138>/Add'
   */
  rtb_sincos_o2_idx_0 = (rt_Lookup2D_Normal(*(real_T (*)[9])&actor_P.MACH[0], 9,
    *(real_T (*)[21])&actor_P.ALPHA[0], 21, *(real_T (*)[189])&actor_P.LIFT[0],
    rtb_Switch2, actor_P.alf_up_lim_Value) * rtb_QS + actor_B.Product1_e) /
    (actor_P.m0 * actor_P.g);

  /* Switch: '<S136>/Switch2' incorporates:
   *  Integrator: '<S131>/n_n0'
   *  RelationalOperator: '<S136>/LowerRelop1'
   */
  if (!(actor_X.n_n0_CSTATE > rtb_sincos_o2_idx_0)) {
    /* Product: '<S137>/Divide' incorporates:
     *  Constant: '<S131>/alf_lo_lim'
     *  Constant: '<S137>/W'
     *  Lookup2D: '<S137>/LIFT'
     *  Product: '<S137>/Lmax'
     *  Sum: '<S137>/Add'
     */
    rtb_sincos_o2_idx_0 = (rt_Lookup2D_Normal(*(real_T (*)[9])&actor_P.MACH[0],
      9, *(real_T (*)[21])&actor_P.ALPHA[0], 21, *(real_T (*)[189])&
      actor_P.LIFT[0], rtb_Switch2, actor_P.alf_lo_lim_Value) * rtb_QS +
      actor_B.Product1_ey) / (actor_P.m0 * actor_P.g);

    /* Switch: '<S136>/Switch' incorporates:
     *  RelationalOperator: '<S136>/UpperRelop'
     */
    if (!(actor_X.n_n0_CSTATE < rtb_sincos_o2_idx_0)) {
      rtb_sincos_o2_idx_0 = actor_X.n_n0_CSTATE;
    }

    /* End of Switch: '<S136>/Switch' */
  }

  /* End of Switch: '<S136>/Switch2' */
  if (rtmIsMajorTimeStep((&actor_M))) {
    /* Product: '<S128>/Product1' incorporates:
     *  Memory: '<S128>/Memory'
     *  Memory: '<S7>/Memory2'
     *  Trigonometry: '<S128>/Trigonometric Function'
     */
    actor_B.Product1_j = actor_DW.Memory2_PreviousInput * std::sin
      (actor_DW.Memory_PreviousInput_p);
  }

  /* Sum: '<S128>/Add' incorporates:
   *  Constant: '<S128>/W'
   *  Product: '<S128>/Product'
   */
  rtb_jxi = actor_P.m0 * actor_P.g * rtb_sincos_o2_idx_0 - actor_B.Product1_j;

  /* Product: '<S128>/Divide1' */
  rtb_sqrt = rtb_jxi / rtb_QS;

  /* MATLAB Function: '<S128>/MATLAB Function' incorporates:
   *  Constant: '<S128>/n_nc1'
   *  Constant: '<S128>/n_nc2'
   *  Constant: '<S128>/n_nc3'
   */
  /* MATLAB Function 'load factor dynamics/L_and_alpha_dynamics/MATLAB Function': '<S135>:1' */
  /* '<S135>:1:2' */
  memcpy(&y[0], &actor_P.LIFT[0], 189U * sizeof(real_T));
  memcpy(&x[0], &actor_P.MACH[0], 9U * sizeof(real_T));
  memset(&alf2CL[0], 0, 21U * sizeof(real_T));
  low_i = 1;
  do {
    exitg2 = 0;
    if (low_i < 10) {
      if (rtIsNaN(actor_P.MACH[low_i - 1])) {
        exitg2 = 1;
      } else {
        low_i++;
      }
    } else {
      if (actor_P.MACH[1] < actor_P.MACH[0]) {
        rtb_Product1_jz = x[0];
        x[0] = x[8];
        x[8] = rtb_Product1_jz;
        rtb_Product1_jz = x[1];
        x[1] = x[7];
        x[7] = rtb_Product1_jz;
        rtb_Product1_jz = x[2];
        x[2] = x[6];
        x[6] = rtb_Product1_jz;
        rtb_Product1_jz = x[3];
        x[3] = x[5];
        x[5] = rtb_Product1_jz;
        for (low_i = 0; low_i < 21; low_i++) {
          low_ip1 = low_i * 9;
          rtb_Product1_jz = y[low_ip1];
          y[low_ip1] = y[low_ip1 + 8];
          y[low_ip1 + 8] = rtb_Product1_jz;
          rtb_Product1_jz = y[low_ip1 + 1];
          y[low_ip1 + 1] = y[low_ip1 + 7];
          y[low_ip1 + 7] = rtb_Product1_jz;
          rtb_Product1_jz = y[low_ip1 + 2];
          y[low_ip1 + 2] = y[low_ip1 + 6];
          y[low_ip1 + 6] = rtb_Product1_jz;
          rtb_Product1_jz = y[low_ip1 + 3];
          y[low_ip1 + 3] = y[low_ip1 + 5];
          y[low_ip1 + 5] = rtb_Product1_jz;
        }
      }

      if (rtIsNaN(rtb_Switch2)) {
        for (low_i = 0; low_i < 21; low_i++) {
          alf2CL[low_i] = (rtNaN);
        }
      } else if (rtb_Switch2 > x[8]) {
        rtb_Product1_jz = (rtb_Switch2 - x[8]) / (x[8] - x[7]);
        for (low_i = 0; low_i < 21; low_i++) {
          alf2CL[low_i] = (y[low_i * 9 + 8] - y[low_i * 9 + 7]) *
            rtb_Product1_jz + y[low_i * 9 + 8];
        }
      } else if (rtb_Switch2 < x[0]) {
        rtb_Product1_jz = (rtb_Switch2 - x[0]) / (x[1] - x[0]);
        for (low_i = 0; low_i < 21; low_i++) {
          alf2CL[low_i] = (y[low_i * 9 + 1] - y[low_i * 9]) * rtb_Product1_jz +
            y[low_i * 9];
        }
      } else {
        low_i = 1;
        low_ip1 = 2;
        high_i = 9;
        while (high_i > low_ip1) {
          mid_i = (low_i + high_i) >> 1;
          if (rtb_Switch2 >= x[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }

        rtb_Product1_jz = (rtb_Switch2 - x[low_i - 1]) / (x[low_i] - x[low_i - 1]);
        if (rtb_Product1_jz == 0.0) {
          for (low_ip1 = 0; low_ip1 < 21; low_ip1++) {
            alf2CL[low_ip1] = y[(low_ip1 * 9 + low_i) - 1];
          }
        } else if (rtb_Product1_jz == 1.0) {
          for (low_ip1 = 0; low_ip1 < 21; low_ip1++) {
            alf2CL[low_ip1] = y[low_ip1 * 9 + low_i];
          }
        } else {
          for (low_ip1 = 0; low_ip1 < 21; low_ip1++) {
            if (y[(low_ip1 * 9 + low_i) - 1] == y[low_ip1 * 9 + low_i]) {
              alf2CL[low_ip1] = y[(low_ip1 * 9 + low_i) - 1];
            } else {
              alf2CL[low_ip1] = y[(low_ip1 * 9 + low_i) - 1] * (1.0 -
                rtb_Product1_jz) + y[low_ip1 * 9 + low_i] * rtb_Product1_jz;
            }
          }
        }
      }

      exitg2 = 1;
    }
  } while (exitg2 == 0);

  /* '<S135>:1:3' */
  memcpy(&b_y[0], &actor_P.ALPHA[0], 21U * sizeof(real_T));
  rtb_Product1_jz = 0.0;
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
          rtb_Product1_jz = alf2CL[low_i];
          alf2CL[low_i] = alf2CL[20 - low_i];
          alf2CL[20 - low_i] = rtb_Product1_jz;
        }

        for (low_i = 0; low_i < 10; low_i++) {
          rtb_Product1_jz = b_y[low_i];
          b_y[low_i] = b_y[20 - low_i];
          b_y[20 - low_i] = rtb_Product1_jz;
        }
      }

      if (rtIsNaN(rtb_sqrt)) {
        rtb_Product1_jz = (rtNaN);
      } else if (rtb_sqrt > alf2CL[20]) {
        rtb_Product1_jz = (rtb_sqrt - alf2CL[20]) / (alf2CL[20] - alf2CL[19]) *
          (b_y[20] - b_y[19]) + b_y[20];
      } else if (rtb_sqrt < alf2CL[0]) {
        rtb_Product1_jz = (rtb_sqrt - alf2CL[0]) / (alf2CL[1] - alf2CL[0]) *
          (b_y[1] - b_y[0]) + b_y[0];
      } else {
        low_i = 1;
        low_ip1 = 2;
        high_i = 21;
        while (high_i > low_ip1) {
          mid_i = (low_i + high_i) >> 1;
          if (rtb_sqrt >= alf2CL[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }

        rtb_sqrt = (rtb_sqrt - alf2CL[low_i - 1]) / (alf2CL[low_i] -
          alf2CL[low_i - 1]);
        if (rtb_sqrt == 0.0) {
          rtb_Product1_jz = b_y[low_i - 1];
        } else if (rtb_sqrt == 1.0) {
          rtb_Product1_jz = b_y[low_i];
        } else if (b_y[low_i - 1] == b_y[low_i]) {
          rtb_Product1_jz = b_y[low_i - 1];
        } else {
          rtb_Product1_jz = (1.0 - rtb_sqrt) * b_y[low_i - 1] + rtb_sqrt *
            b_y[low_i];
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /* UnitConversion: '<S134>/Unit Conversion' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  actor_B.UnitConversion = 0.017453292519943295 * rtb_Product1_jz;

  /* Outport: '<Root>/alpha' */
  actor_Y.alpha = actor_B.UnitConversion;

  /* Lookup_n-D: '<S140>/THRUST_AB_MAX_RATIO' */
  rtb_ixk = look2_binlxpw(rtb_Switch2, rtb_ixk, actor_P.THRUST_MA_AB_IDX,
    actor_P.THRUST_H_IDX, actor_P.THRUST_AB_TAB,
    actor_P.THRUST_AB_MAX_RATIO_maxIndex, 14U);

  /* Product: '<S140>/Product1' incorporates:
   *  Constant: '<S140>/Constant4'
   *  Constant: '<S140>/Constant5'
   *  Product: '<S140>/Product'
   *  Trigonometry: '<S140>/Trigonometric Function'
   */
  rtb_sqrt = rtb_ixk * actor_P.THRUST_AB * actor_P.g * std::cos
    (actor_B.UnitConversion);

  /* Product: '<S127>/Product' incorporates:
   *  Lookup2D: '<S127>/DRAG'
   *  UnitConversion: '<S133>/Unit Conversion'
   */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  rtb_ixk = rt_Lookup2D_Normal(*(real_T (*)[9])&actor_P.MACH[0], 9, *(real_T (*)
    [21])&actor_P.ALPHA[0], 21, *(real_T (*)[189])&actor_P.DRAG[0], rtb_Switch2,
    57.295779513082323 * actor_B.UnitConversion) * rtb_QS;

  /* Product: '<S140>/Divide' incorporates:
   *  Constant: '<S140>/W'
   *  Sum: '<S140>/Add'
   */
  rtb_Switch2 = (rtb_sqrt - rtb_ixk) / (actor_P.m0 * actor_P.g);

  /* Switch: '<S139>/Switch2' incorporates:
   *  Integrator: '<S132>/n_x0'
   *  RelationalOperator: '<S139>/LowerRelop1'
   */
  if (!(actor_X.n_x0_CSTATE > rtb_Switch2)) {
    /* Switch: '<S139>/Switch' incorporates:
     *  Constant: '<S132>/Constant6'
     *  RelationalOperator: '<S139>/UpperRelop'
     */
    if (actor_X.n_x0_CSTATE < actor_P.Constant6_Value) {
      rtb_Switch2 = actor_P.Constant6_Value;
    } else {
      rtb_Switch2 = actor_X.n_x0_CSTATE;
    }

    /* End of Switch: '<S139>/Switch' */
  }

  /* End of Switch: '<S139>/Switch2' */

  /* Product: '<S129>/Divide' incorporates:
   *  Constant: '<S129>/W'
   *  Product: '<S129>/Product'
   *  Sum: '<S129>/Add'
   *  Trigonometry: '<S129>/Trigonometric Function'
   */
  actor_B.Divide = (actor_P.m0 * actor_P.g * rtb_Switch2 + rtb_ixk) / std::cos
    (actor_B.UnitConversion);

  /* Outport: '<Root>/T' */
  actor_Y.T = actor_B.Divide;

  /* Outport: '<Root>/L' */
  actor_Y.L = rtb_jxi;

  /* Outport: '<Root>/D' */
  actor_Y.D = rtb_ixk;

  /* Sqrt: '<S30>/sqrt' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Product: '<S31>/Product'
   *  Product: '<S31>/Product1'
   *  Product: '<S31>/Product2'
   *  Product: '<S31>/Product3'
   *  Sum: '<S31>/Sum'
   */
  rtb_ixk = std::sqrt(((actor_X.q0[0] * actor_X.q0[0] + actor_X.q0[1] *
                        actor_X.q0[1]) + actor_X.q0[2] * actor_X.q0[2]) +
                      actor_X.q0[3] * actor_X.q0[3]);

  /* Product: '<S29>/Product' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  rtb_jxi = actor_X.q0[0] / rtb_ixk;

  /* Product: '<S29>/Product1' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  rtb_Product1_jz = actor_X.q0[1] / rtb_ixk;

  /* Product: '<S29>/Product2' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  rtb_UnitConversion_b = actor_X.q0[2] / rtb_ixk;

  /* Product: '<S29>/Product3' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  rtb_ixk = actor_X.q0[3] / rtb_ixk;

  /* Trigonometry: '<S26>/Trigonometric Function1' incorporates:
   *  Fcn: '<S26>/fcn1'
   *  Fcn: '<S26>/fcn2'
   */
  rtb_Sum_n_idx_0 = rt_atan2d_snf((rtb_Product1_jz * rtb_UnitConversion_b +
    rtb_jxi * rtb_ixk) * 2.0, ((rtb_jxi * rtb_jxi + rtb_Product1_jz *
    rtb_Product1_jz) - rtb_UnitConversion_b * rtb_UnitConversion_b) - rtb_ixk *
    rtb_ixk);

  /* Fcn: '<S26>/fcn3' */
  u0 = (rtb_Product1_jz * rtb_ixk - rtb_jxi * rtb_UnitConversion_b) * -2.0;

  /* Fcn: '<S26>/fcn4' */
  rtb_QS = (rtb_UnitConversion_b * rtb_ixk + rtb_jxi * rtb_Product1_jz) * 2.0;

  /* Fcn: '<S26>/fcn5' */
  rtb_jxi = ((rtb_jxi * rtb_jxi - rtb_Product1_jz * rtb_Product1_jz) -
             rtb_UnitConversion_b * rtb_UnitConversion_b) + rtb_ixk * rtb_ixk;

  /* Outport: '<Root>/attk_g' incorporates:
   *  Trigonometry: '<S26>/Trigonometric Function3'
   */
  actor_Y.attk_g[0] = rt_atan2d_snf(rtb_QS, rtb_jxi);

  /* Trigonometry: '<S26>/trigFcn' */
  if (u0 > 1.0) {
    u0 = 1.0;
  } else {
    if (u0 < -1.0) {
      u0 = -1.0;
    }
  }

  /* Outport: '<Root>/attk_g' incorporates:
   *  Trigonometry: '<S26>/trigFcn'
   */
  actor_Y.attk_g[1] = std::asin(u0);
  actor_Y.attk_g[2] = rtb_Sum_n_idx_0;

  /* Gain: '<S93>/1//2' incorporates:
   *  Constant: '<S5>/Constant'
   *  Constant: '<S5>/Constant1'
   */
  rtb_Sum_n_idx_0 = actor_P.u2_Gain_j * actor_P.Constant1_Value_h;
  rtb_ixk = actor_P.u2_Gain_j * actor_B.UnitConversion;
  rtb_sqrt = actor_P.u2_Gain_j * actor_P.Constant_Value_c;

  /* Trigonometry: '<S93>/sincos' */
  rtb_Product1_jz = std::cos(rtb_Sum_n_idx_0);

  /* Outport: '<Root>/Xg' incorporates:
   *  Integrator: '<S1>/Integrator2'
   */
  actor_Y.Xg[0] = actor_X.xg[0];

  /* Outport: '<Root>/Vg' */
  actor_Y.Vg[0] = rtb_sincos_o1_idx_0;

  /* Trigonometry: '<S93>/sincos' */
  rtb_Sum_n_idx_0 = std::sin(rtb_Sum_n_idx_0);
  rtb_sincos_o1_idx_0 = std::cos(rtb_ixk);

  /* Outport: '<Root>/Xg' incorporates:
   *  Integrator: '<S1>/Integrator2'
   */
  actor_Y.Xg[1] = actor_X.xg[1];

  /* Outport: '<Root>/Vg' */
  actor_Y.Vg[1] = rtb_sincos_o1_idx_1;

  /* Trigonometry: '<S93>/sincos' */
  rtb_ixk = std::sin(rtb_ixk);
  rtb_jxi = std::cos(rtb_sqrt);
  rtb_sincos_o1_idx_1 = std::sin(rtb_sqrt);

  /* Outport: '<Root>/Xg' incorporates:
   *  Integrator: '<S1>/Integrator2'
   */
  actor_Y.Xg[2] = actor_X.xg[2];

  /* Outport: '<Root>/Vg' */
  actor_Y.Vg[2] = rtb_sincos_o1_idx_2;

  /* Fcn: '<S93>/q0' incorporates:
   *  Trigonometry: '<S93>/sincos'
   */
  rtb_QS = rtb_Product1_jz * rtb_sincos_o1_idx_0 * rtb_jxi + rtb_Sum_n_idx_0 *
    rtb_ixk * rtb_sincos_o1_idx_1;

  /* UnaryMinus: '<S89>/Unary Minus' incorporates:
   *  Fcn: '<S93>/q1'
   *  Trigonometry: '<S93>/sincos'
   */
  rtb_sqrt = -(rtb_Product1_jz * rtb_sincos_o1_idx_0 * rtb_sincos_o1_idx_1 -
               rtb_Sum_n_idx_0 * rtb_ixk * rtb_jxi);

  /* UnaryMinus: '<S89>/Unary Minus1' incorporates:
   *  Fcn: '<S93>/q2'
   *  Trigonometry: '<S93>/sincos'
   */
  rtb_UnitConversion_b = -(rtb_Product1_jz * rtb_ixk * rtb_jxi + rtb_Sum_n_idx_0
    * rtb_sincos_o1_idx_0 * rtb_sincos_o1_idx_1);

  /* UnaryMinus: '<S89>/Unary Minus2' incorporates:
   *  Fcn: '<S93>/q3'
   *  Trigonometry: '<S93>/sincos'
   */
  rtb_ixk = -(rtb_Sum_n_idx_0 * rtb_sincos_o1_idx_0 * rtb_jxi - rtb_Product1_jz *
              rtb_ixk * rtb_sincos_o1_idx_1);

  /* Sum: '<S94>/Sum' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Product: '<S94>/Product'
   *  Product: '<S94>/Product1'
   *  Product: '<S94>/Product2'
   *  Product: '<S94>/Product3'
   *  UnaryMinus: '<S8>/Unary Minus'
   *  UnaryMinus: '<S8>/Unary Minus1'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_jxi = ((actor_X.q0[0] * rtb_QS - -actor_X.q0[1] * rtb_sqrt) - -actor_X.q0
             [2] * rtb_UnitConversion_b) - -actor_X.q0[3] * rtb_ixk;

  /* Sum: '<S95>/Sum' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Product: '<S95>/Product'
   *  Product: '<S95>/Product1'
   *  Product: '<S95>/Product2'
   *  Product: '<S95>/Product3'
   *  UnaryMinus: '<S8>/Unary Minus'
   *  UnaryMinus: '<S8>/Unary Minus1'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_Product1_jz = ((actor_X.q0[0] * rtb_sqrt + -actor_X.q0[1] * rtb_QS) +
                     -actor_X.q0[2] * rtb_ixk) - -actor_X.q0[3] *
    rtb_UnitConversion_b;

  /* Sum: '<S96>/Sum' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Product: '<S96>/Product'
   *  Product: '<S96>/Product1'
   *  Product: '<S96>/Product2'
   *  Product: '<S96>/Product3'
   *  UnaryMinus: '<S8>/Unary Minus'
   *  UnaryMinus: '<S8>/Unary Minus1'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_sincos_o1_idx_1 = ((actor_X.q0[0] * rtb_UnitConversion_b - -actor_X.q0[1] *
    rtb_ixk) + -actor_X.q0[2] * rtb_QS) + -actor_X.q0[3] * rtb_sqrt;

  /* Sum: '<S97>/Sum' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Product: '<S97>/Product'
   *  Product: '<S97>/Product1'
   *  Product: '<S97>/Product2'
   *  Product: '<S97>/Product3'
   *  UnaryMinus: '<S8>/Unary Minus'
   *  UnaryMinus: '<S8>/Unary Minus1'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_ixk = ((actor_X.q0[0] * rtb_ixk + -actor_X.q0[1] * rtb_UnitConversion_b) -
             -actor_X.q0[2] * rtb_sqrt) + -actor_X.q0[3] * rtb_QS;

  /* Sqrt: '<S99>/sqrt' incorporates:
   *  Product: '<S100>/Product'
   *  Product: '<S100>/Product1'
   *  Product: '<S100>/Product2'
   *  Product: '<S100>/Product3'
   *  Sum: '<S100>/Sum'
   *  UnaryMinus: '<S90>/Unary Minus'
   *  UnaryMinus: '<S90>/Unary Minus1'
   *  UnaryMinus: '<S90>/Unary Minus2'
   */
  rtb_sqrt = std::sqrt(((rtb_jxi * rtb_jxi + -rtb_Product1_jz * -rtb_Product1_jz)
                        + -rtb_sincos_o1_idx_1 * -rtb_sincos_o1_idx_1) +
                       -rtb_ixk * -rtb_ixk);

  /* Product: '<S98>/Product' */
  rtb_jxi /= rtb_sqrt;

  /* Product: '<S98>/Product1' incorporates:
   *  UnaryMinus: '<S90>/Unary Minus'
   */
  rtb_Product1_jz = -rtb_Product1_jz / rtb_sqrt;

  /* Product: '<S98>/Product2' incorporates:
   *  UnaryMinus: '<S90>/Unary Minus1'
   */
  rtb_UnitConversion_b = -rtb_sincos_o1_idx_1 / rtb_sqrt;

  /* Product: '<S98>/Product3' incorporates:
   *  UnaryMinus: '<S90>/Unary Minus2'
   */
  rtb_QS = -rtb_ixk / rtb_sqrt;

  /* Trigonometry: '<S92>/Trigonometric Function1' incorporates:
   *  Fcn: '<S92>/fcn1'
   *  Fcn: '<S92>/fcn2'
   */
  rtb_sincos_o1_idx_1 = rt_atan2d_snf((rtb_Product1_jz * rtb_UnitConversion_b +
    rtb_jxi * rtb_QS) * 2.0, ((rtb_jxi * rtb_jxi + rtb_Product1_jz *
    rtb_Product1_jz) - rtb_UnitConversion_b * rtb_UnitConversion_b) - rtb_QS *
    rtb_QS);

  /* Fcn: '<S92>/fcn3' */
  u0 = (rtb_Product1_jz * rtb_QS - rtb_jxi * rtb_UnitConversion_b) * -2.0;

  /* Fcn: '<S92>/fcn4' */
  rtb_sqrt = (rtb_UnitConversion_b * rtb_QS + rtb_jxi * rtb_Product1_jz) * 2.0;

  /* Fcn: '<S92>/fcn5' */
  rtb_jxi = ((rtb_jxi * rtb_jxi - rtb_Product1_jz * rtb_Product1_jz) -
             rtb_UnitConversion_b * rtb_UnitConversion_b) + rtb_QS * rtb_QS;

  /* Outport: '<Root>/att_g' incorporates:
   *  Trigonometry: '<S92>/Trigonometric Function3'
   */
  actor_Y.att_g[0] = rt_atan2d_snf(rtb_sqrt, rtb_jxi);

  /* Trigonometry: '<S92>/trigFcn' */
  if (u0 > 1.0) {
    u0 = 1.0;
  } else {
    if (u0 < -1.0) {
      u0 = -1.0;
    }
  }

  /* Outport: '<Root>/att_g' incorporates:
   *  Trigonometry: '<S92>/trigFcn'
   */
  actor_Y.att_g[1] = std::asin(u0);
  actor_Y.att_g[2] = rtb_sincos_o1_idx_1;

  /* Outport: '<Root>/TAS' */
  actor_Y.TAS = rtb_sincos_o2_idx_1;

  /* Outport: '<Root>/nxk' */
  actor_Y.nxk = rtb_Switch2;

  /* Outport: '<Root>/nnk' */
  actor_Y.nnk = rtb_sincos_o2_idx_0;

  /* Saturate: '<S130>/mudot_limit' incorporates:
   *  Integrator: '<S130>/mudot0'
   */
  rtb_QS = -100.0 / actor_P.RAD2DEG;
  rtb_sqrt = 100.0 / actor_P.RAD2DEG;
  if (actor_X.mudot0_CSTATE > rtb_sqrt) {
    rtb_QS = rtb_sqrt;
  } else {
    if (!(actor_X.mudot0_CSTATE < rtb_QS)) {
      rtb_QS = actor_X.mudot0_CSTATE;
    }
  }

  /* End of Saturate: '<S130>/mudot_limit' */

  /* Outport: '<Root>/mudot' */
  actor_Y.mudot = rtb_QS;
  if (rtmIsMajorTimeStep((&actor_M))) {
    /* UnitConversion: '<S119>/Unit Conversion' incorporates:
     *  Constant: '<S105>/ref_pos'
     */
    /* Unit Conversion - from: deg to: rad
       Expression: output = (0.0174533*input) + (0) */
    rtb_sqrt = 0.017453292519943295 * actor_P.FlatEarthtoLLA_psi;

    /* Trigonometry: '<S105>/SinCos' */
    actor_B.SinCos_o1 = std::sin(rtb_sqrt);
    actor_B.SinCos_o2 = std::cos(rtb_sqrt);

    /* Sum: '<S123>/Sum' incorporates:
     *  Constant: '<S123>/Constant'
     *  Constant: '<S123>/f'
     */
    rtb_sqrt = actor_P.f_Value - actor_P.Constant_Value_pz;

    /* Sqrt: '<S124>/sqrt' incorporates:
     *  Constant: '<S124>/Constant'
     *  Product: '<S124>/Product1'
     *  Sum: '<S124>/Sum1'
     */
    rtb_sqrt = std::sqrt(actor_P.Constant_Value_i - rtb_sqrt * rtb_sqrt);

    /* Switch: '<S116>/Switch' incorporates:
     *  Abs: '<S116>/Abs'
     *  Bias: '<S116>/Bias'
     *  Bias: '<S116>/Bias1'
     *  Constant: '<S101>/initial_pos'
     *  Constant: '<S116>/Constant2'
     *  Constant: '<S117>/Constant'
     *  Math: '<S116>/Math Function1'
     *  RelationalOperator: '<S117>/Compare'
     */
    if (std::abs(actor_P.FlatEarthtoLLA_LL0[0]) >
        actor_P.CompareToConstant_const) {
      rtb_sincos_o2_idx_1 = rt_modd_snf(actor_P.FlatEarthtoLLA_LL0[0] +
        actor_P.Bias_Bias_h, actor_P.Constant2_Value) + actor_P.Bias1_Bias_h;
    } else {
      rtb_sincos_o2_idx_1 = actor_P.FlatEarthtoLLA_LL0[0];
    }

    /* End of Switch: '<S116>/Switch' */

    /* Abs: '<S113>/Abs1' */
    rtb_jxi = std::abs(rtb_sincos_o2_idx_1);

    /* RelationalOperator: '<S115>/Compare' incorporates:
     *  Constant: '<S115>/Constant'
     */
    rtb_Compare_g = (rtb_jxi > actor_P.CompareToConstant_const_p);

    /* Switch: '<S113>/Switch' incorporates:
     *  Bias: '<S113>/Bias'
     *  Bias: '<S113>/Bias1'
     *  Gain: '<S113>/Gain'
     *  Product: '<S113>/Divide1'
     */
    if (rtb_Compare_g) {
      /* Signum: '<S113>/Sign1' */
      if (rtb_sincos_o2_idx_1 < 0.0) {
        rtb_sincos_o2_idx_1 = -1.0;
      } else if (rtb_sincos_o2_idx_1 > 0.0) {
        rtb_sincos_o2_idx_1 = 1.0;
      } else {
        if (rtb_sincos_o2_idx_1 == 0.0) {
          rtb_sincos_o2_idx_1 = 0.0;
        }
      }

      /* End of Signum: '<S113>/Sign1' */
      actor_B.Switch = ((rtb_jxi + actor_P.Bias_Bias_e) * actor_P.Gain_Gain_j +
                        actor_P.Bias1_Bias_b) * rtb_sincos_o2_idx_1;
    } else {
      actor_B.Switch = rtb_sincos_o2_idx_1;
    }

    /* End of Switch: '<S113>/Switch' */

    /* UnitConversion: '<S121>/Unit Conversion' */
    /* Unit Conversion - from: deg to: rad
       Expression: output = (0.0174533*input) + (0) */
    rtb_UnitConversion_b = 0.017453292519943295 * actor_B.Switch;

    /* Trigonometry: '<S122>/Trigonometric Function1' */
    rtb_sincos_o2_idx_1 = std::sin(rtb_UnitConversion_b);

    /* Sum: '<S122>/Sum1' incorporates:
     *  Constant: '<S122>/Constant'
     *  Product: '<S122>/Product1'
     */
    rtb_sincos_o2_idx_1 = actor_P.Constant_Value_d - rtb_sqrt * rtb_sqrt *
      rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_1;

    /* Product: '<S120>/Product1' incorporates:
     *  Constant: '<S120>/Constant1'
     *  Sqrt: '<S120>/sqrt'
     */
    rtb_jxi = actor_P.Constant1_Value_a / std::sqrt(rtb_sincos_o2_idx_1);

    /* Trigonometry: '<S120>/Trigonometric Function1' incorporates:
     *  Constant: '<S120>/Constant'
     *  Constant: '<S120>/Constant2'
     *  Product: '<S120>/Product2'
     *  Product: '<S120>/Product3'
     *  Sum: '<S120>/Sum1'
     */
    actor_B.TrigonometricFunction1 = rt_atan2d_snf(actor_P.Constant2_Value_h,
      (actor_P.Constant_Value_kq - rtb_sqrt * rtb_sqrt) * rtb_jxi /
      rtb_sincos_o2_idx_1);

    /* Trigonometry: '<S120>/Trigonometric Function2' incorporates:
     *  Constant: '<S120>/Constant3'
     *  Product: '<S120>/Product4'
     *  Trigonometry: '<S120>/Trigonometric Function'
     */
    actor_B.TrigonometricFunction2 = rt_atan2d_snf(actor_P.Constant3_Value,
      rtb_jxi * std::cos(rtb_UnitConversion_b));

    /* Switch: '<S104>/Switch1' incorporates:
     *  Constant: '<S104>/Constant'
     *  Constant: '<S104>/Constant1'
     */
    if (rtb_Compare_g) {
      rtb_sqrt = actor_P.Constant_Value;
    } else {
      rtb_sqrt = actor_P.Constant1_Value;
    }

    /* End of Switch: '<S104>/Switch1' */

    /* Sum: '<S104>/Sum' incorporates:
     *  Constant: '<S101>/initial_pos'
     */
    rtb_sincos_o2_idx_1 = rtb_sqrt + actor_P.FlatEarthtoLLA_LL0[1];

    /* Switch: '<S114>/Switch' incorporates:
     *  Abs: '<S114>/Abs'
     *  Bias: '<S114>/Bias'
     *  Bias: '<S114>/Bias1'
     *  Constant: '<S114>/Constant2'
     *  Constant: '<S118>/Constant'
     *  Math: '<S114>/Math Function1'
     *  RelationalOperator: '<S118>/Compare'
     */
    if (std::abs(rtb_sincos_o2_idx_1) > actor_P.CompareToConstant_const_i) {
      actor_B.Switch_j = rt_modd_snf(rtb_sincos_o2_idx_1 + actor_P.Bias_Bias_b,
        actor_P.Constant2_Value_j) + actor_P.Bias1_Bias_k;
    } else {
      actor_B.Switch_j = rtb_sincos_o2_idx_1;
    }

    /* End of Switch: '<S114>/Switch' */
  }

  /* Gain: '<S3>/Gain4' */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  rtb_jxi = actor_P.m0 * actor_P.g;

  /* Gain: '<S3>/Gain3' incorporates:
   *  Constant: '<S3>/Constant5'
   *  Constant: '<S3>/Constant6'
   *  Gain: '<S3>/Gain1'
   */
  rtb_ixk = actor_P.m0 * actor_P.g;
  rtb_Product1_jz = actor_P.Gain1_Gain_h * actor_P.Constant6_Value_g * rtb_ixk;
  rtb_sincos_o1_idx_0 = actor_P.Gain1_Gain_h * actor_P.Constant5_Value * rtb_ixk;
  rtb_sincos_o1_idx_2 = actor_P.Gain1_Gain_h * rtb_sincos_o2_idx_0 * rtb_ixk;

  /* Sqrt: '<S42>/sqrt' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Product: '<S43>/Product'
   *  Product: '<S43>/Product1'
   *  Product: '<S43>/Product2'
   *  Product: '<S43>/Product3'
   *  Sum: '<S43>/Sum'
   *  UnaryMinus: '<S36>/Unary Minus'
   *  UnaryMinus: '<S36>/Unary Minus1'
   *  UnaryMinus: '<S36>/Unary Minus2'
   *  UnaryMinus: '<S8>/Unary Minus'
   *  UnaryMinus: '<S8>/Unary Minus1'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_sincos_o2_idx_0 = std::sqrt(((actor_X.q0[0] * actor_X.q0[0] +
    -(-actor_X.q0[1]) * -(-actor_X.q0[1])) + -(-actor_X.q0[2]) * -(-actor_X.q0[2]))
    + -(-actor_X.q0[3]) * -(-actor_X.q0[3]));

  /* Product: '<S38>/Product2' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S36>/Unary Minus1'
   *  UnaryMinus: '<S8>/Unary Minus1'
   */
  rtb_sqrt = -(-actor_X.q0[2]) / rtb_sincos_o2_idx_0;

  /* Product: '<S38>/Product3' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S36>/Unary Minus2'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_sincos_o2_idx_1 = -(-actor_X.q0[3]) / rtb_sincos_o2_idx_0;

  /* Product: '<S38>/Product1' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S36>/Unary Minus'
   *  UnaryMinus: '<S8>/Unary Minus'
   */
  rtb_UnitConversion_b = -(-actor_X.q0[1]) / rtb_sincos_o2_idx_0;

  /* Product: '<S38>/Product' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  rtb_sincos_o2_idx_0 = actor_X.q0[0] / rtb_sincos_o2_idx_0;

  /* Gain: '<S3>/Gain2' incorporates:
   *  Constant: '<S39>/Constant'
   *  Constant: '<S3>/Constant'
   *  Constant: '<S3>/Constant1'
   *  Constant: '<S3>/W in global axes'
   *  Constant: '<S40>/Constant'
   *  Constant: '<S41>/Constant'
   *  Gain: '<S39>/Gain'
   *  Gain: '<S39>/Gain1'
   *  Gain: '<S39>/Gain2'
   *  Gain: '<S3>/Gain4'
   *  Gain: '<S40>/Gain'
   *  Gain: '<S40>/Gain1'
   *  Gain: '<S40>/Gain2'
   *  Gain: '<S41>/Gain'
   *  Gain: '<S41>/Gain1'
   *  Gain: '<S41>/Gain2'
   *  Product: '<S39>/Product'
   *  Product: '<S39>/Product1'
   *  Product: '<S39>/Product2'
   *  Product: '<S39>/Product3'
   *  Product: '<S39>/Product4'
   *  Product: '<S39>/Product5'
   *  Product: '<S39>/Product6'
   *  Product: '<S39>/Product7'
   *  Product: '<S39>/Product8'
   *  Product: '<S40>/Product'
   *  Product: '<S40>/Product1'
   *  Product: '<S40>/Product2'
   *  Product: '<S40>/Product3'
   *  Product: '<S40>/Product4'
   *  Product: '<S40>/Product5'
   *  Product: '<S40>/Product6'
   *  Product: '<S40>/Product7'
   *  Product: '<S40>/Product8'
   *  Product: '<S41>/Product'
   *  Product: '<S41>/Product1'
   *  Product: '<S41>/Product2'
   *  Product: '<S41>/Product3'
   *  Product: '<S41>/Product4'
   *  Product: '<S41>/Product5'
   *  Product: '<S41>/Product6'
   *  Product: '<S41>/Product7'
   *  Product: '<S41>/Product8'
   *  Sum: '<S39>/Sum'
   *  Sum: '<S39>/Sum1'
   *  Sum: '<S39>/Sum2'
   *  Sum: '<S39>/Sum3'
   *  Sum: '<S3>/Add'
   *  Sum: '<S40>/Sum'
   *  Sum: '<S40>/Sum1'
   *  Sum: '<S40>/Sum2'
   *  Sum: '<S40>/Sum3'
   *  Sum: '<S41>/Sum'
   *  Sum: '<S41>/Sum1'
   *  Sum: '<S41>/Sum2'
   *  Sum: '<S41>/Sum3'
   */
  rtb_ixk = 1.0 / actor_P.m0;
  rtb_sincos_o1_idx_1 = (((((actor_P.Constant_Value_cc - rtb_sqrt * rtb_sqrt) -
    rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_1) * actor_P.Gain2_Gain_m *
    actor_P.Winglobalaxes_Value[0] + (rtb_UnitConversion_b * rtb_sqrt +
    rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_1) * actor_P.Gain_Gain_c *
    actor_P.Winglobalaxes_Value[1]) + (rtb_UnitConversion_b *
    rtb_sincos_o2_idx_1 - rtb_sincos_o2_idx_0 * rtb_sqrt) * actor_P.Gain1_Gain_k
    * actor_P.Winglobalaxes_Value[2]) + (rtb_jxi * rtb_Switch2 + rtb_Product1_jz))
    * rtb_ixk;
  rtb_Switch2 = (((((actor_P.Constant_Value_g - rtb_UnitConversion_b *
                     rtb_UnitConversion_b) - rtb_sincos_o2_idx_1 *
                    rtb_sincos_o2_idx_1) * actor_P.Gain2_Gain_mr *
                   actor_P.Winglobalaxes_Value[1] + (rtb_UnitConversion_b *
    rtb_sqrt - rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_1) * actor_P.Gain_Gain_j3
                   * actor_P.Winglobalaxes_Value[0]) + (rtb_sincos_o2_idx_0 *
    rtb_UnitConversion_b + rtb_sqrt * rtb_sincos_o2_idx_1) *
                  actor_P.Gain1_Gain_g * actor_P.Winglobalaxes_Value[2]) +
                 (rtb_jxi * actor_P.Constant1_Value_i + rtb_sincos_o1_idx_0)) *
    rtb_ixk;
  rtb_jxi = ((((rtb_UnitConversion_b * rtb_sincos_o2_idx_1 + rtb_sincos_o2_idx_0
                * rtb_sqrt) * actor_P.Gain_Gain_d * actor_P.Winglobalaxes_Value
               [0] + (rtb_sqrt * rtb_sincos_o2_idx_1 - rtb_sincos_o2_idx_0 *
                      rtb_UnitConversion_b) * actor_P.Gain1_Gain_f *
               actor_P.Winglobalaxes_Value[1]) + ((actor_P.Constant_Value_pf -
    rtb_UnitConversion_b * rtb_UnitConversion_b) - rtb_sqrt * rtb_sqrt) *
              actor_P.Gain2_Gain_c * actor_P.Winglobalaxes_Value[2]) + (rtb_jxi *
              actor_P.Constant_Value_k2 + rtb_sincos_o1_idx_2)) * rtb_ixk;

  /* DotProduct: '<S1>/Dot Product' incorporates:
   *  Integrator: '<S1>/Integrator1'
   */
  rtb_sqrt = (actor_X.Vk[0] * actor_X.Vk[0] + actor_X.Vk[1] * actor_X.Vk[1]) +
    actor_X.Vk[2] * actor_X.Vk[2];

  /* Product: '<S1>/Divide' incorporates:
   *  DotProduct: '<S1>/Dot Product'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S32>/i x j'
   *  Product: '<S32>/j x k'
   *  Product: '<S32>/k x i'
   *  Product: '<S33>/i x k'
   *  Product: '<S33>/j x i'
   *  Product: '<S33>/k x j'
   *  Sum: '<S12>/Sum'
   */
  rtb_Product1_jz = (actor_X.Vk[1] * rtb_jxi - actor_X.Vk[2] * rtb_Switch2) /
    rtb_sqrt;
  rtb_sincos_o1_idx_0 = (actor_X.Vk[2] * rtb_sincos_o1_idx_1 - actor_X.Vk[0] *
    rtb_jxi) / rtb_sqrt;
  rtb_sincos_o1_idx_2 = (actor_X.Vk[0] * rtb_Switch2 - actor_X.Vk[1] *
    rtb_sincos_o1_idx_1) / rtb_sqrt;

  /* Sum: '<S1>/Add' incorporates:
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S34>/i x j'
   *  Product: '<S34>/j x k'
   *  Product: '<S34>/k x i'
   *  Product: '<S35>/i x k'
   *  Product: '<S35>/j x i'
   *  Product: '<S35>/k x j'
   *  Sum: '<S13>/Sum'
   */
  actor_B.Add[0] = rtb_sincos_o1_idx_1 - (rtb_sincos_o1_idx_0 * actor_X.Vk[2] -
    rtb_sincos_o1_idx_2 * actor_X.Vk[1]);
  actor_B.Add[1] = rtb_Switch2 - (rtb_sincos_o1_idx_2 * actor_X.Vk[0] -
    rtb_Product1_jz * actor_X.Vk[2]);
  actor_B.Add[2] = rtb_jxi - (rtb_Product1_jz * actor_X.Vk[1] -
    rtb_sincos_o1_idx_0 * actor_X.Vk[0]);

  /* Sum: '<S1>/Add1' incorporates:
   *  Constant: '<S1>/Constant'
   *  Constant: '<S1>/Constant1'
   */
  rtb_Product1_jz += rtb_QS;
  rtb_sincos_o1_idx_0 += actor_P.Constant1_Value_k;
  rtb_UnitConversion_b = actor_P.Constant_Value_n + rtb_sincos_o1_idx_2;

  /* Sqrt: '<S18>/sqrt' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Product: '<S19>/Product'
   *  Product: '<S19>/Product1'
   *  Product: '<S19>/Product2'
   *  Product: '<S19>/Product3'
   *  Sum: '<S19>/Sum'
   *  UnaryMinus: '<S8>/Unary Minus'
   *  UnaryMinus: '<S8>/Unary Minus1'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_QS = std::sqrt(((actor_X.q0[0] * actor_X.q0[0] + -actor_X.q0[1] *
                       -actor_X.q0[1]) + -actor_X.q0[2] * -actor_X.q0[2]) +
                     -actor_X.q0[3] * -actor_X.q0[3]);

  /* Product: '<S14>/Product' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  rtb_sqrt = actor_X.q0[0] / rtb_QS;

  /* Product: '<S14>/Product1' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S8>/Unary Minus'
   */
  rtb_sincos_o2_idx_0 = -actor_X.q0[1] / rtb_QS;

  /* Product: '<S14>/Product2' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S8>/Unary Minus1'
   */
  rtb_sincos_o2_idx_1 = -actor_X.q0[2] / rtb_QS;

  /* Product: '<S14>/Product3' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  UnaryMinus: '<S8>/Unary Minus2'
   */
  rtb_jxi = -actor_X.q0[3] / rtb_QS;

  /* Sum: '<S15>/Sum' incorporates:
   *  Constant: '<S15>/Constant'
   *  Gain: '<S15>/Gain'
   *  Gain: '<S15>/Gain1'
   *  Gain: '<S15>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S15>/Product'
   *  Product: '<S15>/Product1'
   *  Product: '<S15>/Product2'
   *  Product: '<S15>/Product3'
   *  Product: '<S15>/Product4'
   *  Product: '<S15>/Product5'
   *  Product: '<S15>/Product6'
   *  Product: '<S15>/Product7'
   *  Product: '<S15>/Product8'
   *  Sum: '<S15>/Sum1'
   *  Sum: '<S15>/Sum2'
   *  Sum: '<S15>/Sum3'
   */
  actor_B.Sum = (((actor_P.Constant_Value_m - rtb_sincos_o2_idx_1 *
                   rtb_sincos_o2_idx_1) - rtb_jxi * rtb_jxi) *
                 actor_P.Gain2_Gain_k * actor_X.Vk[0] + (rtb_sincos_o2_idx_0 *
    rtb_sincos_o2_idx_1 + rtb_sqrt * rtb_jxi) * actor_P.Gain_Gain_dr *
                 actor_X.Vk[1]) + (rtb_sincos_o2_idx_0 * rtb_jxi - rtb_sqrt *
    rtb_sincos_o2_idx_1) * actor_P.Gain1_Gain_py * actor_X.Vk[2];

  /* Sum: '<S16>/Sum' incorporates:
   *  Constant: '<S16>/Constant'
   *  Gain: '<S16>/Gain'
   *  Gain: '<S16>/Gain1'
   *  Gain: '<S16>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S16>/Product'
   *  Product: '<S16>/Product1'
   *  Product: '<S16>/Product2'
   *  Product: '<S16>/Product3'
   *  Product: '<S16>/Product4'
   *  Product: '<S16>/Product5'
   *  Product: '<S16>/Product6'
   *  Product: '<S16>/Product7'
   *  Product: '<S16>/Product8'
   *  Sum: '<S16>/Sum1'
   *  Sum: '<S16>/Sum2'
   *  Sum: '<S16>/Sum3'
   */
  actor_B.Sum_l = (((actor_P.Constant_Value_ao - rtb_sincos_o2_idx_0 *
                     rtb_sincos_o2_idx_0) - rtb_jxi * rtb_jxi) *
                   actor_P.Gain2_Gain_cm * actor_X.Vk[1] + (rtb_sincos_o2_idx_0 *
    rtb_sincos_o2_idx_1 - rtb_sqrt * rtb_jxi) * actor_P.Gain_Gain_a *
                   actor_X.Vk[0]) + (rtb_sqrt * rtb_sincos_o2_idx_0 +
    rtb_sincos_o2_idx_1 * rtb_jxi) * actor_P.Gain1_Gain_p5 * actor_X.Vk[2];

  /* Sum: '<S17>/Sum' incorporates:
   *  Constant: '<S17>/Constant'
   *  Gain: '<S17>/Gain'
   *  Gain: '<S17>/Gain1'
   *  Gain: '<S17>/Gain2'
   *  Integrator: '<S1>/Integrator1'
   *  Product: '<S17>/Product'
   *  Product: '<S17>/Product1'
   *  Product: '<S17>/Product2'
   *  Product: '<S17>/Product3'
   *  Product: '<S17>/Product4'
   *  Product: '<S17>/Product5'
   *  Product: '<S17>/Product6'
   *  Product: '<S17>/Product7'
   *  Product: '<S17>/Product8'
   *  Sum: '<S17>/Sum1'
   *  Sum: '<S17>/Sum2'
   *  Sum: '<S17>/Sum3'
   */
  actor_B.Sum_c = ((rtb_sincos_o2_idx_0 * rtb_jxi + rtb_sqrt *
                    rtb_sincos_o2_idx_1) * actor_P.Gain_Gain_o * actor_X.Vk[0] +
                   (rtb_sincos_o2_idx_1 * rtb_jxi - rtb_sqrt *
                    rtb_sincos_o2_idx_0) * actor_P.Gain1_Gain_o * actor_X.Vk[1])
    + ((actor_P.Constant_Value_b - rtb_sincos_o2_idx_0 * rtb_sincos_o2_idx_0) -
       rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_1) * actor_P.Gain2_Gain_h *
    actor_X.Vk[2];

  /* Gain: '<S28>/High Gain Quaternion Normalization' incorporates:
   *  Constant: '<S28>/Constant'
   *  DotProduct: '<S28>/Dot Product'
   *  Integrator: '<S11>/q0 q1 q2 q3'
   *  Sum: '<S28>/Sum'
   */
  rtb_QS = (actor_P.Constant_Value_hd - (((actor_X.q0[0] * actor_X.q0[0] +
    actor_X.q0[1] * actor_X.q0[1]) + actor_X.q0[2] * actor_X.q0[2]) +
             actor_X.q0[3] * actor_X.q0[3])) *
    actor_P.HighGainQuaternionNormalization;

  /* Fcn: '<S28>/q0dot' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  actor_B.q0dot = ((actor_X.q0[1] * rtb_Product1_jz + actor_X.q0[2] *
                    rtb_sincos_o1_idx_0) + actor_X.q0[3] * rtb_UnitConversion_b)
    * -0.5 + rtb_QS * actor_X.q0[0];

  /* Fcn: '<S28>/q1dot' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  actor_B.q1dot = ((actor_X.q0[0] * rtb_Product1_jz + actor_X.q0[2] *
                    rtb_UnitConversion_b) - actor_X.q0[3] * rtb_sincos_o1_idx_0)
    * 0.5 + rtb_QS * actor_X.q0[1];

  /* Fcn: '<S28>/q2dot' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  actor_B.q2dot = ((actor_X.q0[0] * rtb_sincos_o1_idx_0 + actor_X.q0[3] *
                    rtb_Product1_jz) - actor_X.q0[1] * rtb_UnitConversion_b) *
    0.5 + rtb_QS * actor_X.q0[2];

  /* Fcn: '<S28>/q3dot' incorporates:
   *  Integrator: '<S11>/q0 q1 q2 q3'
   */
  actor_B.q3dot = ((actor_X.q0[0] * rtb_UnitConversion_b + actor_X.q0[1] *
                    rtb_sincos_o1_idx_0) - actor_X.q0[2] * rtb_Product1_jz) *
    0.5 + rtb_QS * actor_X.q0[3];

  /* Product: '<S130>/Product' incorporates:
   *  Inport: '<Root>/mudot_c'
   *  Integrator: '<S130>/mudot0'
   *  Lookup_n-D: '<S130>/clp_inv'
   *  Sum: '<S130>/Sum'
   */
  actor_B.Product = (actor_U.mudot_c - actor_X.mudot0_CSTATE) * look1_binlxpw
    (actor_B.UnitConversion, actor_P.clp_inv_alf, actor_P.clp_inv, 11U);

  /* Saturate: '<S131>/n_n_c_lim' incorporates:
   *  Inport: '<Root>/nnk_c'
   */
  if (actor_U.nnk_c > actor_P.n_n_c_lim_UpperSat) {
    rtb_sqrt = actor_P.n_n_c_lim_UpperSat;
  } else if (actor_U.nnk_c < actor_P.n_n_c_lim_LowerSat) {
    rtb_sqrt = actor_P.n_n_c_lim_LowerSat;
  } else {
    rtb_sqrt = actor_U.nnk_c;
  }

  /* End of Saturate: '<S131>/n_n_c_lim' */

  /* Gain: '<S131>/Gain' incorporates:
   *  Integrator: '<S131>/n_n0'
   *  Sum: '<S131>/Sum'
   */
  actor_B.Gain = (rtb_sqrt - actor_X.n_n0_CSTATE) * actor_P.Gain_Gain_a5;

  /* Saturate: '<S132>/n_xc_lim' incorporates:
   *  Inport: '<Root>/nxk_c'
   */
  if (actor_U.nxk_c > actor_P.n_xc_lim_UpperSat) {
    rtb_sqrt = actor_P.n_xc_lim_UpperSat;
  } else if (actor_U.nxk_c < actor_P.n_xc_lim_LowerSat) {
    rtb_sqrt = actor_P.n_xc_lim_LowerSat;
  } else {
    rtb_sqrt = actor_U.nxk_c;
  }

  /* End of Saturate: '<S132>/n_xc_lim' */

  /* Gain: '<S132>/Gain' incorporates:
   *  Integrator: '<S132>/n_x0'
   *  Sum: '<S132>/Sum'
   */
  actor_B.Gain_e = (rtb_sqrt - actor_X.n_x0_CSTATE) * actor_P.Gain_Gain_dw;
  if (rtmIsMajorTimeStep((&actor_M))) {
    /* Update for Integrator: '<S1>/Integrator1' */
    actor_DW.Integrator1_IWORK = 0;

    /* Update for Integrator: '<S11>/q0 q1 q2 q3' */
    actor_DW.q0q1q2q3_IWORK = 0;

    /* Update for Integrator: '<S1>/Integrator2' */
    actor_DW.Integrator2_IWORK = 0;
    if (rtmIsMajorTimeStep((&actor_M))) {
      /* Update for Memory: '<S7>/Memory1' */
      actor_DW.Memory1_PreviousInput = actor_B.Divide;

      /* Update for Memory: '<S7>/Memory' */
      actor_DW.Memory_PreviousInput = actor_B.UnitConversion;

      /* Update for Memory: '<S7>/Memory2' */
      actor_DW.Memory2_PreviousInput = actor_B.Divide;

      /* Update for Memory: '<S128>/Memory' */
      actor_DW.Memory_PreviousInput_p = actor_B.UnitConversion;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep((&actor_M))) {
    rt_ertODEUpdateContinuousStates(&(&actor_M)->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++(&actor_M)->Timing.clockTick0)) {
      ++(&actor_M)->Timing.clockTickH0;
    }

    (&actor_M)->Timing.t[0] = rtsiGetSolverStopTime(&(&actor_M)->solverInfo);

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
      (&actor_M)->Timing.clockTick1++;
      if (!(&actor_M)->Timing.clockTick1) {
        (&actor_M)->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void actorModelClass::actor_derivatives()
{
  XDot_actor_T *_rtXdot;
  _rtXdot = ((XDot_actor_T *) (&actor_M)->derivs);

  /* Derivatives for Integrator: '<S1>/Integrator1' */
  _rtXdot->Vk[0] = actor_B.Add[0];
  _rtXdot->Vk[1] = actor_B.Add[1];
  _rtXdot->Vk[2] = actor_B.Add[2];

  /* Derivatives for Integrator: '<S11>/q0 q1 q2 q3' */
  _rtXdot->q0[0] = actor_B.q0dot;
  _rtXdot->q0[1] = actor_B.q1dot;
  _rtXdot->q0[2] = actor_B.q2dot;
  _rtXdot->q0[3] = actor_B.q3dot;

  /* Derivatives for Integrator: '<S1>/Integrator2' */
  _rtXdot->xg[0] = actor_B.Sum;
  _rtXdot->xg[1] = actor_B.Sum_l;
  _rtXdot->xg[2] = actor_B.Sum_c;

  /* Derivatives for Integrator: '<S131>/n_n0' */
  _rtXdot->n_n0_CSTATE = actor_B.Gain;

  /* Derivatives for Integrator: '<S132>/n_x0' */
  _rtXdot->n_x0_CSTATE = actor_B.Gain_e;

  /* Derivatives for Integrator: '<S130>/mudot0' */
  _rtXdot->mudot0_CSTATE = actor_B.Product;
}

/* Model initialize function */
void actorModelClass::initialize()
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)(&actor_M), 0,
                sizeof(RT_MODEL_actor_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&(&actor_M)->solverInfo, &(&actor_M)
                          ->Timing.simTimeStep);
    rtsiSetTPtr(&(&actor_M)->solverInfo, &rtmGetTPtr((&actor_M)));
    rtsiSetStepSizePtr(&(&actor_M)->solverInfo, &(&actor_M)->Timing.stepSize0);
    rtsiSetdXPtr(&(&actor_M)->solverInfo, &(&actor_M)->derivs);
    rtsiSetContStatesPtr(&(&actor_M)->solverInfo, (real_T **) &(&actor_M)
                         ->contStates);
    rtsiSetNumContStatesPtr(&(&actor_M)->solverInfo, &(&actor_M)
      ->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&actor_M)->solverInfo, &(&actor_M)
      ->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&actor_M)->solverInfo, &(&actor_M)
      ->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&actor_M)->solverInfo, &(&actor_M)
      ->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&(&actor_M)->solverInfo, (&rtmGetErrorStatus((&actor_M))));
    rtsiSetRTModelPtr(&(&actor_M)->solverInfo, (&actor_M));
  }

  rtsiSetSimTimeStep(&(&actor_M)->solverInfo, MAJOR_TIME_STEP);
  (&actor_M)->intgData.y = (&actor_M)->odeY;
  (&actor_M)->intgData.f[0] = (&actor_M)->odeF[0];
  (&actor_M)->intgData.f[1] = (&actor_M)->odeF[1];
  (&actor_M)->intgData.f[2] = (&actor_M)->odeF[2];
  (&actor_M)->intgData.f[3] = (&actor_M)->odeF[3];
  (&actor_M)->contStates = ((X_actor_T *) &actor_X);
  rtsiSetSolverData(&(&actor_M)->solverInfo, (void *)&(&actor_M)->intgData);
  rtsiSetSolverName(&(&actor_M)->solverInfo,"ode4");
  rtmSetTPtr((&actor_M), &(&actor_M)->Timing.tArray[0]);
  (&actor_M)->Timing.stepSize0 = 0.033333333333333333;
  rtmSetFirstInitCond((&actor_M), 1);

  /* block I/O */
  (void) memset(((void *) &actor_B), 0,
                sizeof(B_actor_T));

  /* states (continuous) */
  {
    (void) memset((void *)&actor_X, 0,
                  sizeof(X_actor_T));
  }

  /* states (dwork) */
  (void) memset((void *)&actor_DW, 0,
                sizeof(DW_actor_T));

  /* external inputs */
  (void)memset((void *)&actor_U, 0, sizeof(ExtU_actor_T));

  /* external outputs */
  (void) memset((void *)&actor_Y, 0,
                sizeof(ExtY_actor_T));

  /* Start for Constant: '<S1>/Vk0' */
  actor_B.Vk0[0] = actor_P.Vk0_Value[0];

  /* Start for Constant: '<S1>/xg0' */
  actor_B.xg0[0] = actor_P.xg0_Value[0];

  /* Start for Constant: '<S1>/Vk0' */
  actor_B.Vk0[1] = actor_P.Vk0_Value[1];

  /* Start for Constant: '<S1>/xg0' */
  actor_B.xg0[1] = actor_P.xg0_Value[1];

  /* Start for Constant: '<S1>/Vk0' */
  actor_B.Vk0[2] = actor_P.Vk0_Value[2];

  /* Start for Constant: '<S1>/xg0' */
  actor_B.xg0[2] = actor_P.xg0_Value[2];

  /* InitializeConditions for Integrator: '<S1>/Integrator1' incorporates:
   *  InitializeConditions for Integrator: '<S11>/q0 q1 q2 q3'
   */
  if (rtmIsFirstInitCond((&actor_M))) {
    actor_X.Vk[0] = 150.0;
    actor_X.Vk[1] = 0.0;
    actor_X.Vk[2] = 0.0;
    actor_X.q0[0] = 0.0;
    actor_X.q0[1] = 0.0;
    actor_X.q0[2] = 0.0;
    actor_X.q0[3] = 0.0;
  }

  actor_DW.Integrator1_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S1>/Integrator1' */

  /* InitializeConditions for Integrator: '<S11>/q0 q1 q2 q3' */
  actor_DW.q0q1q2q3_IWORK = 1;

  /* InitializeConditions for Integrator: '<S1>/Integrator2' */
  if (rtmIsFirstInitCond((&actor_M))) {
    actor_X.xg[0] = 0.0;
    actor_X.xg[1] = 0.0;
    actor_X.xg[2] = -3000.0;
  }

  actor_DW.Integrator2_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S1>/Integrator2' */

  /* InitializeConditions for Memory: '<S7>/Memory1' */
  actor_DW.Memory1_PreviousInput = actor_P.Memory1_X0;

  /* InitializeConditions for Memory: '<S7>/Memory' */
  actor_DW.Memory_PreviousInput = actor_P.Memory_X0;

  /* InitializeConditions for Integrator: '<S131>/n_n0' */
  actor_X.n_n0_CSTATE = actor_P.n_n0_IC;

  /* InitializeConditions for Memory: '<S7>/Memory2' */
  actor_DW.Memory2_PreviousInput = actor_P.Memory2_X0;

  /* InitializeConditions for Memory: '<S128>/Memory' */
  actor_DW.Memory_PreviousInput_p = actor_P.Memory_X0_e;

  /* InitializeConditions for Integrator: '<S132>/n_x0' */
  actor_X.n_x0_CSTATE = actor_P.n_x0_IC;

  /* InitializeConditions for Integrator: '<S130>/mudot0' */
  actor_X.mudot0_CSTATE = actor_P.mudot0_IC;

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond((&actor_M))) {
    rtmSetFirstInitCond((&actor_M), 0);
  }
}

/* Model terminate function */
void actorModelClass::terminate()
{
  /* (no terminate code required) */
}

/* Constructor */
actorModelClass::actorModelClass()
{
  static const P_actor_T actor_P_temp = {
    { -4.0, -3.0, -1.0, 0.0, 2.0, 3.0, 4.0, 6.0, 7.0, 9.0, 10.0, 11.0, 13.0,
      15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 32.0 },

    { 0.0319, 0.023785609, 0.026183285, 0.028264766, 0.0423745, 0.051810209,
      0.042255636, 0.04307142, 0.042129381, 0.0266, 0.022685014, 0.024634346,
      0.026741452, 0.039078291, 0.049220323, 0.042670171, 0.042978838,
      0.042199574, 0.0208, 0.020686567, 0.021766258, 0.023941285, 0.033098236,
      0.045001383, 0.04349924, 0.043657328, 0.042689599, 0.021, 0.021014735,
      0.022032502, 0.024489408, 0.034918732, 0.045287169, 0.043913775,
      0.04498432, 0.045223593, 0.0251, 0.026761642, 0.028134601, 0.031552311,
      0.047567606, 0.055171989, 0.053744408, 0.054114619, 0.053348427, 0.0287,
      0.031554579, 0.033976233, 0.039258225, 0.056586035, 0.066872897,
      0.062663674, 0.06082523, 0.059851526, 0.0341, 0.037766032, 0.042082066,
      0.050600871, 0.069033834, 0.080743601, 0.072991253, 0.070959774,
      0.067465323, 0.0483, 0.056767839, 0.065785216, 0.081894831, 0.104523528,
      0.115528274, 0.100405633, 0.095616852, 0.089530681, 0.058, 0.070336874,
      0.082183635, 0.100659502, 0.125373825, 0.136551161, 0.118150894,
      0.110090442, 0.10056336, 0.0819, 0.105959113, 0.124139874, 0.139716419,
      0.172665466, 0.184856919, 0.158989735, 0.139037621, 0.122628718, 0.0979,
      0.127475591, 0.147328986, 0.161577537, 0.200095224, 0.210860381,
      0.180719114, 0.153511211, 0.133661397, 0.119, 0.151806379, 0.171393607,
      0.183667805, 0.228726178, 0.238866074, 0.202448493, 0.1679848, 0.144694076,
      0.1666, 0.203307167, 0.225034444, 0.234046927, 0.288441273, 0.293669986,
      0.24590725, 0.19693198, 0.166759434, 0.2165, 0.26047946, 0.286310527,
      0.29298709, 0.34833429, 0.343906732, 0.289366008, 0.225879159, 0.188824791,
      0.28, 0.327532867, 0.348580697, 0.357265132, 0.408227307, 0.390412611,
      0.332824765, 0.254826339, 0.210890149, 0.3563, 0.410726042, 0.41880884,
      0.425575612, 0.468120324, 0.436918489, 0.376283522, 0.283773518,
      0.232955507, 0.44, 0.506757529, 0.493806398, 0.493886092, 0.528013341,
      0.483424368, 0.41974228, 0.312720697, 0.255020865, 0.528, 0.595238583,
      0.57052363, 0.562196572, 0.587906358, 0.529930247, 0.463201037,
      0.341667877, 0.277086223, 0.6201, 0.680343189, 0.647240862, 0.630507052,
      0.647799375, 0.576436126, 0.506659795, 0.370615056, 0.299151581, 0.7084,
      0.765447795, 0.723958093, 0.698817531, 0.707692392, 0.622942005,
      0.550118552, 0.399562236, 0.321216939, 0.912809347, 0.97820931,
      0.915751172, 0.869593731, 0.857424934, 0.739206702, 0.658765446,
      0.471930184, 0.376380333 },

    { -0.2027, -0.233008047, -0.255361196, -0.289057361, -0.321313984,
      -0.174786201, -0.124258473, -0.123011285, -0.111455189, -0.1297,
      -0.151482541, -0.166341746, -0.188175655, -0.20881539, -0.09883645,
      -0.069714407, -0.071576776, -0.0724588, 0.015, 0.011357815, 0.011935364,
      0.013317842, 0.015190782, 0.055105019, 0.039373725, 0.028936504,
      0.005983346, 0.0856, 0.091398796, 0.102598654, 0.111739846, 0.119067858,
      0.134528943, 0.093917791, 0.07649893, 0.048146305, 0.2257, 0.252674007,
      0.275046579, 0.306238523, 0.311470781, 0.300015468, 0.209805707,
      0.166467602, 0.13050199, 0.2965, 0.331830821, 0.363051245, 0.409726333,
      0.406149678, 0.391002311, 0.26465974, 0.210780122, 0.170106735, 0.3666,
      0.410816826, 0.450683737, 0.509724241, 0.49557352, 0.478966463,
      0.317654153, 0.255643088, 0.209515734, 0.5057, 0.571140841, 0.622493291,
      0.686523077, 0.656953717, 0.646940047, 0.417049162, 0.34272286,
      0.287128671, 0.5754, 0.655422419, 0.706259548, 0.765805736, 0.73868157,
      0.724614208, 0.467092199, 0.384969181, 0.325935139, 0.7084, 0.821713179,
      0.858577802, 0.894334416, 0.896983004, 0.854172851, 0.571687725,
      0.469461824, 0.403548075, 0.7803, 0.897077376, 0.925210069, 0.945216936,
      0.964849851, 0.909473794, 0.625089983, 0.511708146, 0.442354544, 0.8649,
      0.9617838, 0.985401719, 0.99511589, 1.026873924, 0.956764732, 0.678492241,
      0.553954467, 0.481161012, 1.0049, 1.07904034, 1.103881787, 1.098141517,
      1.133751049, 1.051444654, 0.785296756, 0.63844711, 0.558773948, 1.1142,
      1.182871693, 1.214991785, 1.199788633, 1.230743789, 1.175337979,
      0.892101272, 0.722939753, 0.636386885, 1.2213, 1.285787066, 1.300515898,
      1.301251935, 1.327736529, 1.32984256, 0.998905788, 0.807432395,
      0.713999821, 1.3157, 1.401969851, 1.379158835, 1.388953154, 1.424729269,
      1.484347141, 1.105710303, 0.891925038, 0.791612758, 1.4044, 1.509254639,
      1.440975509, 1.476654374, 1.521722009, 1.638851722, 1.212514819,
      0.976417681, 0.869225695, 1.4854, 1.572870819, 1.49672526, 1.564355593,
      1.618714749, 1.793356303, 1.319319335, 1.060910324, 0.946838631, 1.5575,
      1.625069424, 1.55247501, 1.652056813, 1.715707489, 1.947860884,
      1.426123851, 1.145402967, 1.024451568, 1.5983, 1.677268029, 1.608224761,
      1.739758032, 1.812700229, 2.102365465, 1.532928366, 1.229895609,
      1.102064504, 1.648215844, 1.807764543, 1.747599137, 1.959011081,
      2.05518208, 2.488626918, 1.799939655, 1.441127216, 1.296096845 },

    { 0.2, 0.6, 0.8, 0.9, 1.0, 1.1, 1.515, 1.816, 2.201 },
    57.3,
    27.87,
    13226.75,

    { 1.1816, 1.1308, 1.115, 1.1284, 1.1707, 1.2411, 1.3287, 1.4365, 1.5711,
      1.7301, 1.8314, 1.97, 2.07, 2.2, 1.0, 0.9599, 0.9474, 0.9589, 0.9942,
      1.0529, 1.1254, 1.2149, 1.326, 1.4579, 1.57, 1.69, 1.8, 1.92, 0.8184,
      0.789, 0.7798, 0.7894, 0.8177, 0.8648, 0.9221, 0.9933, 1.0809, 1.1857,
      1.3086, 1.41, 1.53, 1.64, 0.6627, 0.6406, 0.634, 0.642, 0.6647, 0.7017,
      0.7462, 0.8021, 0.87, 0.9512, 1.0474, 1.24, 1.34, 1.44, 0.528, 0.5116,
      0.507, 0.5134, 0.5309, 0.5596, 0.5936, 0.636, 0.6874, 0.7495, 0.8216, 0.91,
      1.0, 1.1, 0.3756, 0.3645, 0.3615, 0.3661, 0.3784, 0.3983, 0.4219, 0.4509,
      0.486, 0.5289, 0.5786, 0.6359, 0.72, 0.8, 0.2327, 0.2258, 0.224, 0.2268,
      0.2345, 0.2467, 0.2614, 0.2794, 0.3011, 0.3277, 0.3585, 0.394, 0.46, 0.52
    },

    { -3047.99999536704, 0.0, 3047.99999536704, 6095.99999073408,
      9143.99998610112, 12191.99998146816, 15239.999976835199 },

    { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6 },

    { 2.7777777777777777, 2.785515320334262, 2.2573363431151243,
      2.3809523809523809, 2.6109660574412534, 2.6666666666666665,
      3.03951367781155, 3.4013605442176873, 4.3478260869565215,
      4.7619047619047619, 8.3333333333333339, 10.0 },

    { -0.175, -0.087, 0.0, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611,
      0.698, 0.785 },
    9.8,
    11773.296,

    { 44.94, -73.097 },
    180.0,
    90.0,
    180.0,
    180.0,
    90.0,
    180.0,
    0.0,
    -90.0,
    -1.0,
    90.0,
    180.0,
    -180.0,
    180.0,
    -180.0,
    0.0,
    180.0,
    -90.0,
    -1.0,
    90.0,
    360.0,
    180.0,
    -180.0,
    360.0,
    180.0,
    -180.0,

    { 150.0, 0.0, 0.0 },

    { 0.0, 0.0, 0.0 },
    0.5,
    0.5,
    2.0,
    2.0,
    2.0,
    2.0,
    0.5,
    2.0,
    2.0,
    2.0,
    2.0,
    0.5,
    2.0,
    1.0,
    0.0,
    288.15,

    { 0.0, 0.0, -3000.0 },
    -1.0,
    11000.0,
    0.0,
    0.0065,
    401.87433999999996,
    27.0,
    0.00347041471455839,
    5.2558756014667134,
    1.225,
    11000.0,
    0.0,
    -9000.0,
    0.034163191409533639,
    0.5,
    0.0,
    0.0,
    0.0,
    -10.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.5,
    0.0,
    180.0,
    1.0,
    6.378137E+6,
    1.0,
    1.0,
    0.0033528106647474805,
    1.0,
    1.0,
    1.0,
    360.0,
    0.0,
    360.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    -1.0,

    { 0.0, 0.0, 115378.30080000001 },
    0.5,
    2.0,
    2.0,
    2.0,
    2.0,
    0.5,
    2.0,
    2.0,
    2.0,
    2.0,
    0.5,
    2.0,
    0.0,
    0.0,
    0.5,
    2.0,
    2.0,
    2.0,
    0.5,
    2.0,
    2.0,
    2.0,
    0.5,
    2.0,
    2.0,
    2.0,
    1.0,
    1.0,
    9.0,
    -3.0,
    2.0,
    2.0,
    0.0,
    2.0,

    { 13U, 6U }
  };                                   /* Modifiable parameters */

  /* Initialize tunable parameters */
  actor_P = actor_P_temp;
}

/* Destructor */
actorModelClass::~actorModelClass()
{
  /* Currently there is no destructor body generated.*/
}

/* Real-Time Model get method */
RT_MODEL_actor_T * actorModelClass::getRTM()
{
  return (&actor_M);
}
