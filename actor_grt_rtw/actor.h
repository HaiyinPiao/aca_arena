/*
 * actor.h
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

#ifndef RTW_HEADER_actor_h_
#define RTW_HEADER_actor_h_
#include <cmath>
#include <float.h>
#include <math.h>
#include <string.h>
#ifndef actor_COMMON_INCLUDES_
# define actor_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* actor_COMMON_INCLUDES_ */

#include "actor_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rt_look.h"
#include "rt_look2d_normal.h"
#include "rtGetInf.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlkStateChangeFlag
# define rtmGetBlkStateChangeFlag(rtm) ((rtm)->blkStateChange)
#endif

#ifndef rtmSetBlkStateChangeFlag
# define rtmSetBlkStateChangeFlag(rtm, val) ((rtm)->blkStateChange = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

/* Block signals (auto storage) */
typedef struct {
  real_T Vk0[3];                       /* '<S1>/Vk0' */
  real_T q0;                           /* '<S27>/q0' */
  real_T q1;                           /* '<S27>/q1' */
  real_T q2;                           /* '<S27>/q2' */
  real_T q3;                           /* '<S27>/q3' */
  real_T Product1;                     /* '<S88>/Product1' */
  real_T Product2;                     /* '<S88>/Product2' */
  real_T xg0[3];                       /* '<S1>/xg0' */
  real_T Product1_a;                   /* '<S48>/Product1' */
  real_T Product2_f;                   /* '<S48>/Product2' */
  real_T Product1_e;                   /* '<S138>/Product1' */
  real_T Product1_ey;                  /* '<S137>/Product1' */
  real_T Product1_j;                   /* '<S128>/Product1' */
  real_T UnitConversion;               /* '<S134>/Unit Conversion' */
  real_T Divide;                       /* '<S129>/Divide' */
  real_T SinCos_o1;                    /* '<S105>/SinCos' */
  real_T SinCos_o2;                    /* '<S105>/SinCos' */
  real_T Switch;                       /* '<S113>/Switch' */
  real_T TrigonometricFunction1;       /* '<S120>/Trigonometric Function1' */
  real_T TrigonometricFunction2;       /* '<S120>/Trigonometric Function2' */
  real_T Switch_j;                     /* '<S114>/Switch' */
  real_T Add[3];                       /* '<S1>/Add' */
  real_T Sum;                          /* '<S15>/Sum' */
  real_T Sum_l;                        /* '<S16>/Sum' */
  real_T Sum_c;                        /* '<S17>/Sum' */
  real_T q0dot;                        /* '<S28>/q0dot' */
  real_T q1dot;                        /* '<S28>/q1dot' */
  real_T q2dot;                        /* '<S28>/q2dot' */
  real_T q3dot;                        /* '<S28>/q3dot' */
  real_T Product;                      /* '<S130>/Product' */
  real_T Gain;                         /* '<S131>/Gain' */
  real_T Gain_e;                       /* '<S132>/Gain' */
} B_actor_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T Memory1_PreviousInput;        /* '<S7>/Memory1' */
  real_T Memory_PreviousInput;         /* '<S7>/Memory' */
  real_T Memory2_PreviousInput;        /* '<S7>/Memory2' */
  real_T Memory_PreviousInput_p;       /* '<S128>/Memory' */
  int_T Integrator1_IWORK;             /* '<S1>/Integrator1' */
  int_T q0q1q2q3_IWORK;                /* '<S11>/q0 q1 q2 q3' */
  int_T Integrator2_IWORK;             /* '<S1>/Integrator2' */
} DW_actor_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T Vk[3];                        /* '<S1>/Integrator1' */
  real_T q0[4];                        /* '<S11>/q0 q1 q2 q3' */
  real_T xg[3];                        /* '<S1>/Integrator2' */
  real_T n_n0_CSTATE;                  /* '<S131>/n_n0' */
  real_T n_x0_CSTATE;                  /* '<S132>/n_x0' */
  real_T mudot0_CSTATE;                /* '<S130>/mudot0' */
} X_actor_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T Vk[3];                        /* '<S1>/Integrator1' */
  real_T q0[4];                        /* '<S11>/q0 q1 q2 q3' */
  real_T xg[3];                        /* '<S1>/Integrator2' */
  real_T n_n0_CSTATE;                  /* '<S131>/n_n0' */
  real_T n_x0_CSTATE;                  /* '<S132>/n_x0' */
  real_T mudot0_CSTATE;                /* '<S130>/mudot0' */
} XDot_actor_T;

/* State disabled  */
typedef struct {
  boolean_T Vk[3];                     /* '<S1>/Integrator1' */
  boolean_T q0[4];                     /* '<S11>/q0 q1 q2 q3' */
  boolean_T xg[3];                     /* '<S1>/Integrator2' */
  boolean_T n_n0_CSTATE;               /* '<S131>/n_n0' */
  boolean_T n_x0_CSTATE;               /* '<S132>/n_x0' */
  boolean_T mudot0_CSTATE;             /* '<S130>/mudot0' */
} XDis_actor_T;

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
} ODE4_IntgData;

#endif

/* External inputs (root inport signals with auto storage) */
typedef struct {
  real_T nxk_c;                        /* '<Root>/nxk_c' */
  real_T nnk_c;                        /* '<Root>/nnk_c' */
  real_T mudot_c;                      /* '<Root>/mudot_c' */
} ExtU_actor_T;

/* External outputs (root outports fed by signals with auto storage) */
typedef struct {
  real_T alpha;                        /* '<Root>/alpha' */
  real_T T;                            /* '<Root>/T' */
  real_T L;                            /* '<Root>/L' */
  real_T D;                            /* '<Root>/D' */
  real_T attk_g[3];                    /* '<Root>/attk_g' */
  real_T att_g[3];                     /* '<Root>/att_g' */
  real_T Xg[3];                        /* '<Root>/Xg' */
  real_T Vg[3];                        /* '<Root>/Vg' */
  real_T TAS;                          /* '<Root>/TAS' */
  real_T nxk;                          /* '<Root>/nxk' */
  real_T nnk;                          /* '<Root>/nnk' */
  real_T mudot;                        /* '<Root>/mudot' */
} ExtY_actor_T;

/* Parameters (auto storage) */
struct P_actor_T_ {
  real_T ALPHA[21];                    /* Variable: ALPHA
                                        * Referenced by:
                                        *   '<S127>/DRAG'
                                        *   '<S128>/n_nc2'
                                        *   '<S137>/LIFT'
                                        *   '<S138>/LIFT'
                                        */
  real_T DRAG[189];                    /* Variable: DRAG
                                        * Referenced by: '<S127>/DRAG'
                                        */
  real_T LIFT[189];                    /* Variable: LIFT
                                        * Referenced by:
                                        *   '<S128>/n_nc3'
                                        *   '<S137>/LIFT'
                                        *   '<S138>/LIFT'
                                        */
  real_T MACH[9];                      /* Variable: MACH
                                        * Referenced by:
                                        *   '<S127>/DRAG'
                                        *   '<S128>/n_nc1'
                                        *   '<S137>/LIFT'
                                        *   '<S138>/LIFT'
                                        */
  real_T RAD2DEG;                      /* Variable: RAD2DEG
                                        * Referenced by: '<S130>/mudot_limit'
                                        */
  real_T S;                            /* Variable: S
                                        * Referenced by: '<S4>/Constant1'
                                        */
  real_T THRUST_AB;                    /* Variable: THRUST_AB
                                        * Referenced by: '<S140>/Constant4'
                                        */
  real_T THRUST_AB_TAB[98];            /* Variable: THRUST_AB_TAB
                                        * Referenced by: '<S140>/THRUST_AB_MAX_RATIO'
                                        */
  real_T THRUST_H_IDX[7];              /* Variable: THRUST_H_IDX
                                        * Referenced by: '<S140>/THRUST_AB_MAX_RATIO'
                                        */
  real_T THRUST_MA_AB_IDX[14];         /* Variable: THRUST_MA_AB_IDX
                                        * Referenced by: '<S140>/THRUST_AB_MAX_RATIO'
                                        */
  real_T clp_inv[12];                  /* Variable: clp_inv
                                        * Referenced by: '<S130>/clp_inv'
                                        */
  real_T clp_inv_alf[12];              /* Variable: clp_inv_alf
                                        * Referenced by: '<S130>/clp_inv'
                                        */
  real_T g;                            /* Variable: g
                                        * Referenced by:
                                        *   '<S3>/Gain3'
                                        *   '<S3>/Gain4'
                                        *   '<S128>/W'
                                        *   '<S129>/W'
                                        *   '<S137>/W'
                                        *   '<S138>/W'
                                        *   '<S140>/Constant5'
                                        *   '<S140>/W'
                                        */
  real_T m0;                           /* Variable: m0
                                        * Referenced by:
                                        *   '<S3>/Gain2'
                                        *   '<S3>/Gain3'
                                        *   '<S3>/Gain4'
                                        *   '<S128>/W'
                                        *   '<S129>/W'
                                        *   '<S137>/W'
                                        *   '<S138>/W'
                                        *   '<S140>/W'
                                        */
  real_T FlatEarthtoLLA_LL0[2];        /* Mask Parameter: FlatEarthtoLLA_LL0
                                        * Referenced by: '<S101>/initial_pos'
                                        */
  real_T CompareToConstant_const;      /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S117>/Constant'
                                        */
  real_T CompareToConstant_const_p;    /* Mask Parameter: CompareToConstant_const_p
                                        * Referenced by: '<S115>/Constant'
                                        */
  real_T CompareToConstant_const_i;    /* Mask Parameter: CompareToConstant_const_i
                                        * Referenced by: '<S118>/Constant'
                                        */
  real_T CompareToConstant_const_j;    /* Mask Parameter: CompareToConstant_const_j
                                        * Referenced by: '<S111>/Constant'
                                        */
  real_T CompareToConstant_const_k;    /* Mask Parameter: CompareToConstant_const_k
                                        * Referenced by: '<S109>/Constant'
                                        */
  real_T CompareToConstant_const_c;    /* Mask Parameter: CompareToConstant_const_c
                                        * Referenced by: '<S112>/Constant'
                                        */
  real_T FlatEarthtoLLA_psi;           /* Mask Parameter: FlatEarthtoLLA_psi
                                        * Referenced by: '<S105>/ref_pos'
                                        */
  real_T Bias_Bias;                    /* Expression: -90
                                        * Referenced by: '<S107>/Bias'
                                        */
  real_T Gain_Gain;                    /* Expression: -1
                                        * Referenced by: '<S107>/Gain'
                                        */
  real_T Bias1_Bias;                   /* Expression: +90
                                        * Referenced by: '<S107>/Bias1'
                                        */
  real_T Bias_Bias_f;                  /* Expression: 180
                                        * Referenced by: '<S110>/Bias'
                                        */
  real_T Bias1_Bias_l;                 /* Expression: -180
                                        * Referenced by: '<S110>/Bias1'
                                        */
  real_T Bias_Bias_l;                  /* Expression: 180
                                        * Referenced by: '<S108>/Bias'
                                        */
  real_T Bias1_Bias_j;                 /* Expression: -180
                                        * Referenced by: '<S108>/Bias1'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S104>/Constant1'
                                        */
  real_T Constant_Value;               /* Expression: 180
                                        * Referenced by: '<S104>/Constant'
                                        */
  real_T Bias_Bias_e;                  /* Expression: -90
                                        * Referenced by: '<S113>/Bias'
                                        */
  real_T Gain_Gain_j;                  /* Expression: -1
                                        * Referenced by: '<S113>/Gain'
                                        */
  real_T Bias1_Bias_b;                 /* Expression: +90
                                        * Referenced by: '<S113>/Bias1'
                                        */
  real_T Constant2_Value;              /* Expression: 360
                                        * Referenced by: '<S116>/Constant2'
                                        */
  real_T Bias_Bias_h;                  /* Expression: 180
                                        * Referenced by: '<S116>/Bias'
                                        */
  real_T Bias1_Bias_h;                 /* Expression: -180
                                        * Referenced by: '<S116>/Bias1'
                                        */
  real_T Constant2_Value_j;            /* Expression: 360
                                        * Referenced by: '<S114>/Constant2'
                                        */
  real_T Bias_Bias_b;                  /* Expression: 180
                                        * Referenced by: '<S114>/Bias'
                                        */
  real_T Bias1_Bias_k;                 /* Expression: -180
                                        * Referenced by: '<S114>/Bias1'
                                        */
  real_T Vk0_Value[3];                 /* Expression: [150 0 0]
                                        * Referenced by: '<S1>/Vk0'
                                        */
  real_T att_g0_Value[3];              /* Expression: [0 0 0]
                                        * Referenced by: '<S11>/att_g0'
                                        */
  real_T u2_Gain;                      /* Expression: 0.5
                                        * Referenced by: '<S27>/1//2'
                                        */
  real_T Constant_Value_k;             /* Expression: 0.5
                                        * Referenced by: '<S21>/Constant'
                                        */
  real_T Gain2_Gain;                   /* Expression: 2
                                        * Referenced by: '<S21>/Gain2'
                                        */
  real_T Gain_Gain_jt;                 /* Expression: 2
                                        * Referenced by: '<S21>/Gain'
                                        */
  real_T Gain1_Gain;                   /* Expression: 2
                                        * Referenced by: '<S21>/Gain1'
                                        */
  real_T Gain_Gain_h;                  /* Expression: 2
                                        * Referenced by: '<S22>/Gain'
                                        */
  real_T Constant_Value_f;             /* Expression: 0.5
                                        * Referenced by: '<S22>/Constant'
                                        */
  real_T Gain2_Gain_e;                 /* Expression: 2
                                        * Referenced by: '<S22>/Gain2'
                                        */
  real_T Gain1_Gain_b;                 /* Expression: 2
                                        * Referenced by: '<S22>/Gain1'
                                        */
  real_T Gain_Gain_b;                  /* Expression: 2
                                        * Referenced by: '<S23>/Gain'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: 2
                                        * Referenced by: '<S23>/Gain1'
                                        */
  real_T Constant_Value_o;             /* Expression: 0.5
                                        * Referenced by: '<S23>/Constant'
                                        */
  real_T Gain2_Gain_eb;                /* Expression: 2
                                        * Referenced by: '<S23>/Gain2'
                                        */
  real_T Gain4_Gain;                   /* Expression: 1.0
                                        * Referenced by: '<Root>/Gain4'
                                        */
  real_T Constant_Value_h;             /* Expression: 0.0
                                        * Referenced by: '<S4>/Constant'
                                        */
  real_T SeaLevelTemperature_Value;    /* Expression: T0
                                        * Referenced by: '<S45>/Sea Level  Temperature'
                                        */
  real_T xg0_Value[3];                 /* Expression: [0 0 -3000]
                                        * Referenced by: '<S1>/xg0'
                                        */
  real_T Gain_Gain_l;                  /* Expression: -1
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Limitaltitudetotroposhere_Upper;/* Expression: h_trop
                                          * Referenced by: '<S45>/Limit  altitude  to troposhere'
                                          */
  real_T Limitaltitudetotroposhere_Lower;/* Expression: h0
                                          * Referenced by: '<S45>/Limit  altitude  to troposhere'
                                          */
  real_T LapseRate_Gain;               /* Expression: L
                                        * Referenced by: '<S45>/Lapse Rate'
                                        */
  real_T gammaR_Gain;                  /* Expression: gamma*R
                                        * Referenced by: '<S45>/gamma*R'
                                        */
  real_T alf_up_lim_Value;             /* Expression: 27.0
                                        * Referenced by: '<S131>/alf_up_lim'
                                        */
  real_T uT0_Gain;                     /* Expression: 1/T0
                                        * Referenced by: '<S45>/1//T0'
                                        */
  real_T Constant_Value_p;             /* Expression: g/(L*R)
                                        * Referenced by: '<S45>/Constant'
                                        */
  real_T rho0_Gain;                    /* Expression: rho0
                                        * Referenced by: '<S45>/rho0'
                                        */
  real_T AltitudeofTroposphere_Value;  /* Expression: h_trop
                                        * Referenced by: '<S45>/Altitude of Troposphere'
                                        */
  real_T LimitaltitudetoStratosphere_Upp;/* Expression: 0
                                          * Referenced by: '<S45>/Limit  altitude  to Stratosphere'
                                          */
  real_T LimitaltitudetoStratosphere_Low;/* Expression: h_trop-h_strat
                                          * Referenced by: '<S45>/Limit  altitude  to Stratosphere'
                                          */
  real_T gR_Gain;                      /* Expression: g/R
                                        * Referenced by: '<S45>/g//R'
                                        */
  real_T u2rhoV2_Gain;                 /* Expression: 1/2
                                        * Referenced by: '<S44>/1//2rhoV^2'
                                        */
  real_T Memory1_X0;                   /* Expression: 0
                                        * Referenced by: '<S7>/Memory1'
                                        */
  real_T Memory_X0;                    /* Expression: 0
                                        * Referenced by: '<S7>/Memory'
                                        */
  real_T n_n0_IC;                      /* Expression: 0
                                        * Referenced by: '<S131>/n_n0'
                                        */
  real_T alf_lo_lim_Value;             /* Expression: -10.0
                                        * Referenced by: '<S131>/alf_lo_lim'
                                        */
  real_T Memory2_X0;                   /* Expression: 0
                                        * Referenced by: '<S7>/Memory2'
                                        */
  real_T Memory_X0_e;                  /* Expression: 0
                                        * Referenced by: '<S128>/Memory'
                                        */
  real_T n_x0_IC;                      /* Expression: 0
                                        * Referenced by: '<S132>/n_x0'
                                        */
  real_T Constant6_Value;              /* Expression: 0.0
                                        * Referenced by: '<S132>/Constant6'
                                        */
  real_T Constant1_Value_h;            /* Expression: 0
                                        * Referenced by: '<S5>/Constant1'
                                        */
  real_T Constant_Value_c;             /* Expression: 0
                                        * Referenced by: '<S5>/Constant'
                                        */
  real_T u2_Gain_j;                    /* Expression: 0.5
                                        * Referenced by: '<S93>/1//2'
                                        */
  real_T mudot0_IC;                    /* Expression: 0
                                        * Referenced by: '<S130>/mudot0'
                                        */
  real_T Constant_Value_a;             /* Expression: 180
                                        * Referenced by: '<S103>/Constant'
                                        */
  real_T Constant2_Value_h;            /* Expression: 1
                                        * Referenced by: '<S120>/Constant2'
                                        */
  real_T Constant1_Value_a;            /* Expression: R
                                        * Referenced by: '<S120>/Constant1'
                                        */
  real_T Constant_Value_d;             /* Expression: 1
                                        * Referenced by: '<S122>/Constant'
                                        */
  real_T Constant_Value_i;             /* Expression: 1
                                        * Referenced by: '<S124>/Constant'
                                        */
  real_T Constant_Value_pz;            /* Expression: F
                                        * Referenced by: '<S123>/Constant'
                                        */
  real_T f_Value;                      /* Expression: 1
                                        * Referenced by: '<S123>/f'
                                        */
  real_T Constant_Value_kq;            /* Expression: 1
                                        * Referenced by: '<S120>/Constant'
                                        */
  real_T Constant3_Value;              /* Expression: 1
                                        * Referenced by: '<S120>/Constant3'
                                        */
  real_T Constant2_Value_i;            /* Expression: 360
                                        * Referenced by: '<S110>/Constant2'
                                        */
  real_T Constant1_Value_b;            /* Expression: 0
                                        * Referenced by: '<S103>/Constant1'
                                        */
  real_T Constant2_Value_l;            /* Expression: 360
                                        * Referenced by: '<S108>/Constant2'
                                        */
  real_T Constant8_Value;              /* Expression: 0
                                        * Referenced by: '<S6>/Constant8'
                                        */
  real_T Constant1_Value_i;            /* Expression: 0
                                        * Referenced by: '<S3>/Constant1'
                                        */
  real_T Constant_Value_k2;            /* Expression: 0
                                        * Referenced by: '<S3>/Constant'
                                        */
  real_T Constant6_Value_g;            /* Expression: 0
                                        * Referenced by: '<S3>/Constant6'
                                        */
  real_T Constant5_Value;              /* Expression: 0
                                        * Referenced by: '<S3>/Constant5'
                                        */
  real_T Gain1_Gain_h;                 /* Expression: -1
                                        * Referenced by: '<S3>/Gain1'
                                        */
  real_T Winglobalaxes_Value[3];       /* Expression: [0 0 m0*g]
                                        * Referenced by: '<S3>/W in global axes'
                                        */
  real_T Constant_Value_cc;            /* Expression: 0.5
                                        * Referenced by: '<S39>/Constant'
                                        */
  real_T Gain2_Gain_m;                 /* Expression: 2
                                        * Referenced by: '<S39>/Gain2'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 2
                                        * Referenced by: '<S39>/Gain'
                                        */
  real_T Gain1_Gain_k;                 /* Expression: 2
                                        * Referenced by: '<S39>/Gain1'
                                        */
  real_T Gain_Gain_j3;                 /* Expression: 2
                                        * Referenced by: '<S40>/Gain'
                                        */
  real_T Constant_Value_g;             /* Expression: 0.5
                                        * Referenced by: '<S40>/Constant'
                                        */
  real_T Gain2_Gain_mr;                /* Expression: 2
                                        * Referenced by: '<S40>/Gain2'
                                        */
  real_T Gain1_Gain_g;                 /* Expression: 2
                                        * Referenced by: '<S40>/Gain1'
                                        */
  real_T Gain_Gain_d;                  /* Expression: 2
                                        * Referenced by: '<S41>/Gain'
                                        */
  real_T Gain1_Gain_f;                 /* Expression: 2
                                        * Referenced by: '<S41>/Gain1'
                                        */
  real_T Constant_Value_pf;            /* Expression: 0.5
                                        * Referenced by: '<S41>/Constant'
                                        */
  real_T Gain2_Gain_c;                 /* Expression: 2
                                        * Referenced by: '<S41>/Gain2'
                                        */
  real_T Constant1_Value_k;            /* Expression: 0
                                        * Referenced by: '<S1>/Constant1'
                                        */
  real_T Constant_Value_n;             /* Expression: 0
                                        * Referenced by: '<S1>/Constant'
                                        */
  real_T Constant_Value_m;             /* Expression: 0.5
                                        * Referenced by: '<S15>/Constant'
                                        */
  real_T Gain_Gain_dr;                 /* Expression: 2
                                        * Referenced by: '<S15>/Gain'
                                        */
  real_T Gain1_Gain_py;                /* Expression: 2
                                        * Referenced by: '<S15>/Gain1'
                                        */
  real_T Gain2_Gain_k;                 /* Expression: 2
                                        * Referenced by: '<S15>/Gain2'
                                        */
  real_T Constant_Value_ao;            /* Expression: 0.5
                                        * Referenced by: '<S16>/Constant'
                                        */
  real_T Gain_Gain_a;                  /* Expression: 2
                                        * Referenced by: '<S16>/Gain'
                                        */
  real_T Gain1_Gain_p5;                /* Expression: 2
                                        * Referenced by: '<S16>/Gain1'
                                        */
  real_T Gain2_Gain_cm;                /* Expression: 2
                                        * Referenced by: '<S16>/Gain2'
                                        */
  real_T Constant_Value_b;             /* Expression: 0.5
                                        * Referenced by: '<S17>/Constant'
                                        */
  real_T Gain_Gain_o;                  /* Expression: 2
                                        * Referenced by: '<S17>/Gain'
                                        */
  real_T Gain1_Gain_o;                 /* Expression: 2
                                        * Referenced by: '<S17>/Gain1'
                                        */
  real_T Gain2_Gain_h;                 /* Expression: 2
                                        * Referenced by: '<S17>/Gain2'
                                        */
  real_T Constant_Value_hd;            /* Expression: 1
                                        * Referenced by: '<S28>/Constant'
                                        */
  real_T HighGainQuaternionNormalization;/* Expression: 1.0
                                          * Referenced by: '<S28>/High Gain Quaternion Normalization'
                                          */
  real_T n_n_c_lim_UpperSat;           /* Expression: 9
                                        * Referenced by: '<S131>/n_n_c_lim'
                                        */
  real_T n_n_c_lim_LowerSat;           /* Expression: -3
                                        * Referenced by: '<S131>/n_n_c_lim'
                                        */
  real_T Gain_Gain_a5;                 /* Expression: 2
                                        * Referenced by: '<S131>/Gain'
                                        */
  real_T n_xc_lim_UpperSat;            /* Expression: 2
                                        * Referenced by: '<S132>/n_xc_lim'
                                        */
  real_T n_xc_lim_LowerSat;            /* Expression: 0
                                        * Referenced by: '<S132>/n_xc_lim'
                                        */
  real_T Gain_Gain_dw;                 /* Expression: 2
                                        * Referenced by: '<S132>/Gain'
                                        */
  uint32_T THRUST_AB_MAX_RATIO_maxIndex[2];/* Computed Parameter: THRUST_AB_MAX_RATIO_maxIndex
                                            * Referenced by: '<S140>/THRUST_AB_MAX_RATIO'
                                            */
};

/* Real-time Model Data Structure */
struct tag_RTM_actor_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_actor_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T blkStateChange;
  real_T odeY[13];
  real_T odeF[4][13];
  ODE4_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    boolean_T firstInitCondFlag;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

#ifdef __cplusplus

extern "C" {

#endif

#ifdef __cplusplus

}
#endif

/* Class declaration for model actor */
class actorModelClass {
  /* public data and function members */
 public:
  /* External inputs */
  ExtU_actor_T actor_U;

  /* External outputs */
  ExtY_actor_T actor_Y;

  /* model initialize function */
  void initialize();

  /* model step function */
  void step();

  /* model terminate function */
  void terminate();

  /* Constructor */
  actorModelClass();

  /* Destructor */
  ~actorModelClass();

  /* Real-Time Model get method */
  RT_MODEL_actor_T * getRTM();

  /* Block signals */
  B_actor_T actor_B;

  /* Tunable parameters */
  P_actor_T actor_P;

  /* private data and function members */
 private:


  /* Block states */
  DW_actor_T actor_DW;
  X_actor_T actor_X;                   /* Block continuous states */

  /* Real-Time Model */
  RT_MODEL_actor_T actor_M;

  /* Continuous states update member function*/
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  /* Derivatives member function */
  void actor_derivatives();
};

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'actor'
 * '<S1>'   : 'actor/3DoF plant'
 * '<S2>'   : 'actor/Angle Conversion'
 * '<S3>'   : 'actor/acck Synthesise'
 * '<S4>'   : 'actor/atmosphere'
 * '<S5>'   : 'actor/attitude calibration by adding alpha'
 * '<S6>'   : 'actor/flight gear viewport'
 * '<S7>'   : 'actor/load factor dynamics'
 * '<S8>'   : 'actor/3DoF plant/Quaternion Conjugate1'
 * '<S9>'   : 'actor/3DoF plant/Quaternion Rotation1'
 * '<S10>'  : 'actor/3DoF plant/Quaternion Rotation2'
 * '<S11>'  : 'actor/3DoF plant/anguler_integration'
 * '<S12>'  : 'actor/3DoF plant/wx(Iw)'
 * '<S13>'  : 'actor/3DoF plant/wx(Iw)1'
 * '<S14>'  : 'actor/3DoF plant/Quaternion Rotation1/Quaternion Normalize'
 * '<S15>'  : 'actor/3DoF plant/Quaternion Rotation1/V1'
 * '<S16>'  : 'actor/3DoF plant/Quaternion Rotation1/V2'
 * '<S17>'  : 'actor/3DoF plant/Quaternion Rotation1/V3'
 * '<S18>'  : 'actor/3DoF plant/Quaternion Rotation1/Quaternion Normalize/Quaternion Modulus'
 * '<S19>'  : 'actor/3DoF plant/Quaternion Rotation1/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S20>'  : 'actor/3DoF plant/Quaternion Rotation2/Quaternion Normalize'
 * '<S21>'  : 'actor/3DoF plant/Quaternion Rotation2/V1'
 * '<S22>'  : 'actor/3DoF plant/Quaternion Rotation2/V2'
 * '<S23>'  : 'actor/3DoF plant/Quaternion Rotation2/V3'
 * '<S24>'  : 'actor/3DoF plant/Quaternion Rotation2/Quaternion Normalize/Quaternion Modulus'
 * '<S25>'  : 'actor/3DoF plant/Quaternion Rotation2/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S26>'  : 'actor/3DoF plant/anguler_integration/Quaternions to Rotation Angles'
 * '<S27>'  : 'actor/3DoF plant/anguler_integration/Rotation Angles to Quaternions1'
 * '<S28>'  : 'actor/3DoF plant/anguler_integration/qdot'
 * '<S29>'  : 'actor/3DoF plant/anguler_integration/Quaternions to Rotation Angles/Quaternion Normalize'
 * '<S30>'  : 'actor/3DoF plant/anguler_integration/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus'
 * '<S31>'  : 'actor/3DoF plant/anguler_integration/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S32>'  : 'actor/3DoF plant/wx(Iw)/Subsystem'
 * '<S33>'  : 'actor/3DoF plant/wx(Iw)/Subsystem1'
 * '<S34>'  : 'actor/3DoF plant/wx(Iw)1/Subsystem'
 * '<S35>'  : 'actor/3DoF plant/wx(Iw)1/Subsystem1'
 * '<S36>'  : 'actor/acck Synthesise/Quaternion Conjugate1'
 * '<S37>'  : 'actor/acck Synthesise/Quaternion Rotation1'
 * '<S38>'  : 'actor/acck Synthesise/Quaternion Rotation1/Quaternion Normalize'
 * '<S39>'  : 'actor/acck Synthesise/Quaternion Rotation1/V1'
 * '<S40>'  : 'actor/acck Synthesise/Quaternion Rotation1/V2'
 * '<S41>'  : 'actor/acck Synthesise/Quaternion Rotation1/V3'
 * '<S42>'  : 'actor/acck Synthesise/Quaternion Rotation1/Quaternion Normalize/Quaternion Modulus'
 * '<S43>'  : 'actor/acck Synthesise/Quaternion Rotation1/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S44>'  : 'actor/atmosphere/Dynamic Pressure'
 * '<S45>'  : 'actor/atmosphere/ISA Atmosphere Model'
 * '<S46>'  : 'actor/atmosphere/Ideal Airspeed Correction'
 * '<S47>'  : 'actor/atmosphere/Mach Number'
 * '<S48>'  : 'actor/atmosphere/Dynamic Pressure/dot'
 * '<S49>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2EAS'
 * '<S50>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2TAS'
 * '<S51>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2CAS'
 * '<S52>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2TAS'
 * '<S53>'  : 'actor/atmosphere/Ideal Airspeed Correction/Mach <= 1.0'
 * '<S54>'  : 'actor/atmosphere/Ideal Airspeed Correction/Pressure Conversion'
 * '<S55>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2CAS'
 * '<S56>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2EAS'
 * '<S57>'  : 'actor/atmosphere/Ideal Airspeed Correction/Velocity Conversion'
 * '<S58>'  : 'actor/atmosphere/Ideal Airspeed Correction/Velocity Conversion1'
 * '<S59>'  : 'actor/atmosphere/Ideal Airspeed Correction/Velocity Conversion2'
 * '<S60>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach'
 * '<S61>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach/Relative Ratio'
 * '<S62>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach/Relative Ratio/Density Conversion'
 * '<S63>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach/Relative Ratio/Pressure Conversion'
 * '<S64>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach/Relative Ratio/Temperature Conversion'
 * '<S65>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach'
 * '<S66>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach/Relative Ratio'
 * '<S67>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach/Relative Ratio/Density Conversion'
 * '<S68>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach/Relative Ratio/Pressure Conversion'
 * '<S69>'  : 'actor/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach/Relative Ratio/Temperature Conversion'
 * '<S70>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach'
 * '<S71>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach/Relative Ratio'
 * '<S72>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach/Relative Ratio/Density Conversion'
 * '<S73>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach/Relative Ratio/Pressure Conversion'
 * '<S74>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach/Relative Ratio/Temperature Conversion'
 * '<S75>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach'
 * '<S76>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach/Relative Ratio'
 * '<S77>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach/Relative Ratio/Density Conversion'
 * '<S78>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach/Relative Ratio/Pressure Conversion'
 * '<S79>'  : 'actor/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach/Relative Ratio/Temperature Conversion'
 * '<S80>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2CAS/Relative Ratio'
 * '<S81>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2CAS/Relative Ratio/Density Conversion'
 * '<S82>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2CAS/Relative Ratio/Pressure Conversion'
 * '<S83>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2CAS/Relative Ratio/Temperature Conversion'
 * '<S84>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2EAS/Relative Ratio'
 * '<S85>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2EAS/Relative Ratio/Density Conversion'
 * '<S86>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2EAS/Relative Ratio/Pressure Conversion'
 * '<S87>'  : 'actor/atmosphere/Ideal Airspeed Correction/TAS2EAS/Relative Ratio/Temperature Conversion'
 * '<S88>'  : 'actor/atmosphere/Mach Number/dot'
 * '<S89>'  : 'actor/attitude calibration by adding alpha/Quaternion Conjugate1'
 * '<S90>'  : 'actor/attitude calibration by adding alpha/Quaternion Conjugate2'
 * '<S91>'  : 'actor/attitude calibration by adding alpha/Quaternion Multiplication'
 * '<S92>'  : 'actor/attitude calibration by adding alpha/Quaternions to Rotation Angles'
 * '<S93>'  : 'actor/attitude calibration by adding alpha/Rotation Angles to Quaternions1'
 * '<S94>'  : 'actor/attitude calibration by adding alpha/Quaternion Multiplication/q0'
 * '<S95>'  : 'actor/attitude calibration by adding alpha/Quaternion Multiplication/q1'
 * '<S96>'  : 'actor/attitude calibration by adding alpha/Quaternion Multiplication/q2'
 * '<S97>'  : 'actor/attitude calibration by adding alpha/Quaternion Multiplication/q3'
 * '<S98>'  : 'actor/attitude calibration by adding alpha/Quaternions to Rotation Angles/Quaternion Normalize'
 * '<S99>'  : 'actor/attitude calibration by adding alpha/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus'
 * '<S100>' : 'actor/attitude calibration by adding alpha/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S101>' : 'actor/flight gear viewport/Flat Earth to LLA'
 * '<S102>' : 'actor/flight gear viewport/FlightGear Preconfigured 6DoF Animation'
 * '<S103>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap'
 * '<S104>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap1'
 * '<S105>' : 'actor/flight gear viewport/Flat Earth to LLA/LongLat_offset'
 * '<S106>' : 'actor/flight gear viewport/Flat Earth to LLA/pos_deg'
 * '<S107>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90'
 * '<S108>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap/Wrap Longitude'
 * '<S109>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Compare To Constant'
 * '<S110>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Wrap Angle 180'
 * '<S111>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S112>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap/Wrap Longitude/Compare To Constant'
 * '<S113>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90'
 * '<S114>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap1/Wrap Longitude'
 * '<S115>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Compare To Constant'
 * '<S116>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Wrap Angle 180'
 * '<S117>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S118>' : 'actor/flight gear viewport/Flat Earth to LLA/LatLong wrap1/Wrap Longitude/Compare To Constant'
 * '<S119>' : 'actor/flight gear viewport/Flat Earth to LLA/LongLat_offset/Angle Conversion2'
 * '<S120>' : 'actor/flight gear viewport/Flat Earth to LLA/LongLat_offset/Find Radian//Distance'
 * '<S121>' : 'actor/flight gear viewport/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/Angle Conversion2'
 * '<S122>' : 'actor/flight gear viewport/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/denom'
 * '<S123>' : 'actor/flight gear viewport/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/e'
 * '<S124>' : 'actor/flight gear viewport/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/e^4'
 * '<S125>' : 'actor/flight gear viewport/FlightGear Preconfigured 6DoF Animation/Angle Conversion'
 * '<S126>' : 'actor/flight gear viewport/FlightGear Preconfigured 6DoF Animation/Angle Conversion1'
 * '<S127>' : 'actor/load factor dynamics/D_dynamics'
 * '<S128>' : 'actor/load factor dynamics/L_and_alpha_dynamics'
 * '<S129>' : 'actor/load factor dynamics/T_dynamics'
 * '<S130>' : 'actor/load factor dynamics/mudot_dynamics'
 * '<S131>' : 'actor/load factor dynamics/n_n_dynamics'
 * '<S132>' : 'actor/load factor dynamics/n_x_dynamics'
 * '<S133>' : 'actor/load factor dynamics/D_dynamics/Angle Conversion'
 * '<S134>' : 'actor/load factor dynamics/L_and_alpha_dynamics/Angle Conversion'
 * '<S135>' : 'actor/load factor dynamics/L_and_alpha_dynamics/MATLAB Function'
 * '<S136>' : 'actor/load factor dynamics/n_n_dynamics/n_n_lim'
 * '<S137>' : 'actor/load factor dynamics/n_n_dynamics/n_n_lo_lim'
 * '<S138>' : 'actor/load factor dynamics/n_n_dynamics/n_n_up_lim'
 * '<S139>' : 'actor/load factor dynamics/n_x_dynamics/n_x_lim'
 * '<S140>' : 'actor/load factor dynamics/n_x_dynamics/n_x_up_lim'
 */
#endif                                 /* RTW_HEADER_actor_h_ */
