/*
 * f16.h
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

#ifndef RTW_HEADER_f16_h_
#define RTW_HEADER_f16_h_
#include <math.h>
#include <string.h>
#include <float.h>
#include <stddef.h>
#ifndef f16_COMMON_INCLUDES_
# define f16_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#endif                                 /* f16_COMMON_INCLUDES_ */

#include "f16_types.h"

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

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
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

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
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

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

/* Block signals (auto storage) */
typedef struct {
  real_T q0;                           /* '<S31>/q0' */
  real_T q1;                           /* '<S31>/q1' */
  real_T q2;                           /* '<S31>/q2' */
  real_T q3;                           /* '<S31>/q3' */
  real_T xg_0[3];                      /* '<S1>/xg_0' */
  real_T Integrator2[3];               /* '<S1>/Integrator2' */
  real_T UnitConversion[3];            /* '<S2>/Unit Conversion' */
  real_T Vk0[3];                       /* '<S1>/Vk0' */
  real_T Sum;                          /* '<S25>/Sum' */
  real_T Sum_n;                        /* '<S26>/Sum' */
  real_T Sum_nr;                       /* '<S27>/Sum' */
  real_T Product1;                     /* '<S95>/Product1' */
  real_T Product2;                     /* '<S95>/Product2' */
  real_T Product1_a;                   /* '<S55>/Product1' */
  real_T Product2_f;                   /* '<S55>/Product2' */
  real_T Product;                      /* '<S3>/Product' */
  real_T Sum_e;                        /* '<S46>/Sum' */
  real_T Sum_p;                        /* '<S47>/Sum' */
  real_T Sum_o;                        /* '<S48>/Sum' */
  real_T Gain2[3];                     /* '<S6>/Gain2' */
  real_T Add[3];                       /* '<S1>/Add' */
  real_T Sum_pg;                       /* '<S19>/Sum' */
  real_T Sum_l;                        /* '<S20>/Sum' */
  real_T Sum_c;                        /* '<S21>/Sum' */
  real_T q0dot;                        /* '<S32>/q0dot' */
  real_T q1dot;                        /* '<S32>/q1dot' */
  real_T q2dot;                        /* '<S32>/q2dot' */
  real_T q3dot;                        /* '<S32>/q3dot' */
  real_T Product_n;                    /* '<S8>/Product' */
  real_T Gain;                         /* '<S10>/Gain' */
  real_T Gain_e;                       /* '<S11>/Gain' */
} B_f16_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  struct {
    void *LoggedData;
  } Euler_PWORK;                       /* '<Root>/Euler' */

  struct {
    void *LoggedData;
  } Euler1_PWORK;                      /* '<Root>/Euler1' */

  struct {
    void *LoggedData;
  } Euler2_PWORK;                      /* '<Root>/Euler2' */

  struct {
    void *LoggedData;
  } Euler3_PWORK;                      /* '<Root>/Euler3' */

  struct {
    void *LoggedData;
  } Scope_PWORK;                       /* '<S6>/Scope' */

  struct {
    void *LoggedData;
  } Scope1_PWORK;                      /* '<S6>/Scope1' */

  int_T q0q1q2q3_IWORK;                /* '<S15>/q0 q1 q2 q3' */
  int_T Integrator2_IWORK;             /* '<S1>/Integrator2' */
  int_T Integrator1_IWORK;             /* '<S1>/Integrator1' */
} DW_f16_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T q0q1q2q3_CSTATE[4];           /* '<S15>/q0 q1 q2 q3' */
  real_T Integrator2_CSTATE[3];        /* '<S1>/Integrator2' */
  real_T Integrator1_CSTATE[3];        /* '<S1>/Integrator1' */
  real_T Integrator_CSTATE;            /* '<S10>/Integrator' */
  real_T Integrator_CSTATE_b;          /* '<S11>/Integrator' */
  real_T Integrator_CSTATE_d;          /* '<S8>/Integrator' */
} X_f16_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T q0q1q2q3_CSTATE[4];           /* '<S15>/q0 q1 q2 q3' */
  real_T Integrator2_CSTATE[3];        /* '<S1>/Integrator2' */
  real_T Integrator1_CSTATE[3];        /* '<S1>/Integrator1' */
  real_T Integrator_CSTATE;            /* '<S10>/Integrator' */
  real_T Integrator_CSTATE_b;          /* '<S11>/Integrator' */
  real_T Integrator_CSTATE_d;          /* '<S8>/Integrator' */
} XDot_f16_T;

/* State disabled  */
typedef struct {
  boolean_T q0q1q2q3_CSTATE[4];        /* '<S15>/q0 q1 q2 q3' */
  boolean_T Integrator2_CSTATE[3];     /* '<S1>/Integrator2' */
  boolean_T Integrator1_CSTATE[3];     /* '<S1>/Integrator1' */
  boolean_T Integrator_CSTATE;         /* '<S10>/Integrator' */
  boolean_T Integrator_CSTATE_b;       /* '<S11>/Integrator' */
  boolean_T Integrator_CSTATE_d;       /* '<S8>/Integrator' */
} XDis_f16_T;

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
  real_T n_nc;                         /* '<Root>/n_nc' */
  real_T n_xc;                         /* '<Root>/n_xc' */
  real_T mu_dotc;                      /* '<Root>/mu_dotc' */
} ExtU_f16_T;

/* External outputs (root outports fed by signals with auto storage) */
typedef struct {
  real_T attitude_g[3];                /* '<Root>/attitude_g' */
  real_T Xg[3];                        /* '<Root>/Xg' */
} ExtY_f16_T;

/* Parameters (auto storage) */
struct P_f16_T_ {
  real_T ALPHA[21];                    /* Variable: ALPHA
                                        * Referenced by:
                                        *   '<S3>/DRAG'
                                        *   '<S4>/n_nc2'
                                        *   '<S97>/LIFT'
                                        *   '<S98>/LIFT'
                                        */
  real_T DRAG[189];                    /* Variable: DRAG
                                        * Referenced by: '<S3>/DRAG'
                                        */
  real_T LIFT[189];                    /* Variable: LIFT
                                        * Referenced by:
                                        *   '<S4>/n_nc3'
                                        *   '<S97>/LIFT'
                                        *   '<S98>/LIFT'
                                        */
  real_T MACH[9];                      /* Variable: MACH
                                        * Referenced by:
                                        *   '<S3>/DRAG'
                                        *   '<S4>/n_nc1'
                                        *   '<S97>/LIFT'
                                        *   '<S98>/LIFT'
                                        */
  real_T RAD2DEG;                      /* Variable: RAD2DEG
                                        * Referenced by: '<S8>/p_limit'
                                        */
  real_T S;                            /* Variable: S
                                        * Referenced by: '<S7>/Constant1'
                                        */
  real_T THRUST_AB;                    /* Variable: THRUST_AB
                                        * Referenced by: '<S100>/Constant4'
                                        */
  real_T THRUST_AB_TAB[98];            /* Variable: THRUST_AB_TAB
                                        * Referenced by: '<S100>/THRUST_AB_MAX_RATIO'
                                        */
  real_T THRUST_H_IDX[7];              /* Variable: THRUST_H_IDX
                                        * Referenced by: '<S100>/THRUST_AB_MAX_RATIO'
                                        */
  real_T THRUST_MA_AB_IDX[14];         /* Variable: THRUST_MA_AB_IDX
                                        * Referenced by: '<S100>/THRUST_AB_MAX_RATIO'
                                        */
  real_T clp_inv[12];                  /* Variable: clp_inv
                                        * Referenced by: '<S8>/clp_inv'
                                        */
  real_T clp_inv_alf[12];              /* Variable: clp_inv_alf
                                        * Referenced by: '<S8>/clp_inv'
                                        */
  real_T g;                            /* Variable: g
                                        * Referenced by:
                                        *   '<S4>/W'
                                        *   '<S5>/W'
                                        *   '<S97>/W'
                                        *   '<S98>/W'
                                        *   '<S100>/Constant5'
                                        *   '<S100>/W'
                                        */
  real_T m0;                           /* Variable: m0
                                        * Referenced by:
                                        *   '<S4>/W'
                                        *   '<S5>/W'
                                        *   '<S6>/Gain2'
                                        *   '<S97>/W'
                                        *   '<S98>/W'
                                        *   '<S100>/W'
                                        */
  real_T InitialEulerAngles_Value[3];  /* Expression: [0 0 0]
                                        * Referenced by: '<S15>/Initial Euler Angles'
                                        */
  real_T u2_Gain;                      /* Expression: 0.5
                                        * Referenced by: '<S31>/1//2'
                                        */
  real_T xg_0_Value[3];                /* Expression: [0 0 -1000]
                                        * Referenced by: '<S1>/xg_0'
                                        */
  real_T Vk0_Value[3];                 /* Expression: [200 0 0]
                                        * Referenced by: '<S1>/Vk0'
                                        */
  real_T Constant_Value;               /* Expression: 0.5
                                        * Referenced by: '<S25>/Constant'
                                        */
  real_T Gain2_Gain;                   /* Expression: 2
                                        * Referenced by: '<S25>/Gain2'
                                        */
  real_T Gain_Gain;                    /* Expression: 2
                                        * Referenced by: '<S25>/Gain'
                                        */
  real_T Gain1_Gain;                   /* Expression: 2
                                        * Referenced by: '<S25>/Gain1'
                                        */
  real_T Gain_Gain_h;                  /* Expression: 2
                                        * Referenced by: '<S26>/Gain'
                                        */
  real_T Constant_Value_f;             /* Expression: 0.5
                                        * Referenced by: '<S26>/Constant'
                                        */
  real_T Gain2_Gain_e;                 /* Expression: 2
                                        * Referenced by: '<S26>/Gain2'
                                        */
  real_T Gain1_Gain_b;                 /* Expression: 2
                                        * Referenced by: '<S26>/Gain1'
                                        */
  real_T Gain_Gain_b;                  /* Expression: 2
                                        * Referenced by: '<S27>/Gain'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: 2
                                        * Referenced by: '<S27>/Gain1'
                                        */
  real_T Constant_Value_o;             /* Expression: 0.5
                                        * Referenced by: '<S27>/Constant'
                                        */
  real_T Gain2_Gain_eb;                /* Expression: 2
                                        * Referenced by: '<S27>/Gain2'
                                        */
  real_T Gain4_Gain;                   /* Expression: 1.0
                                        * Referenced by: '<Root>/Gain4'
                                        */
  real_T Constant_Value_h;             /* Expression: 0.0
                                        * Referenced by: '<S7>/Constant'
                                        */
  real_T SeaLevelTemperature_Value;    /* Expression: T0
                                        * Referenced by: '<S52>/Sea Level  Temperature'
                                        */
  real_T Gain_Gain_l;                  /* Expression: -1
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Limitaltitudetotroposhere_Upper;/* Expression: h_trop
                                          * Referenced by: '<S52>/Limit  altitude  to troposhere'
                                          */
  real_T Limitaltitudetotroposhere_Lower;/* Expression: h0
                                          * Referenced by: '<S52>/Limit  altitude  to troposhere'
                                          */
  real_T LapseRate_Gain;               /* Expression: L
                                        * Referenced by: '<S52>/Lapse Rate'
                                        */
  real_T gammaR_Gain;                  /* Expression: gamma*R
                                        * Referenced by: '<S52>/gamma*R'
                                        */
  real_T alf_up_lim_Value;             /* Expression: 27.0
                                        * Referenced by: '<S10>/alf_up_lim'
                                        */
  real_T uT0_Gain;                     /* Expression: 1/T0
                                        * Referenced by: '<S52>/1//T0'
                                        */
  real_T Constant_Value_p;             /* Expression: g/(L*R)
                                        * Referenced by: '<S52>/Constant'
                                        */
  real_T rho0_Gain;                    /* Expression: rho0
                                        * Referenced by: '<S52>/rho0'
                                        */
  real_T AltitudeofTroposphere_Value;  /* Expression: h_trop
                                        * Referenced by: '<S52>/Altitude of Troposphere'
                                        */
  real_T LimitaltitudetoStratosphere_Upp;/* Expression: 0
                                          * Referenced by: '<S52>/Limit  altitude  to Stratosphere'
                                          */
  real_T LimitaltitudetoStratosphere_Low;/* Expression: h_trop-h_strat
                                          * Referenced by: '<S52>/Limit  altitude  to Stratosphere'
                                          */
  real_T gR_Gain;                      /* Expression: g/R
                                        * Referenced by: '<S52>/g//R'
                                        */
  real_T u2rhoV2_Gain;                 /* Expression: 1/2
                                        * Referenced by: '<S51>/1//2rhoV^2'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S10>/Integrator'
                                        */
  real_T alf_lo_lim_Value;             /* Expression: -10.0
                                        * Referenced by: '<S10>/alf_lo_lim'
                                        */
  real_T Integrator_IC_e;              /* Expression: 0
                                        * Referenced by: '<S11>/Integrator'
                                        */
  real_T Constant6_Value;              /* Expression: 0.0
                                        * Referenced by: '<S11>/Constant6'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S6>/Constant1'
                                        */
  real_T Constant_Value_k;             /* Expression: 0
                                        * Referenced by: '<S6>/Constant'
                                        */
  real_T Constant3_Value;              /* Expression: 0
                                        * Referenced by: '<S6>/Constant3'
                                        */
  real_T Constant2_Value;              /* Expression: 0
                                        * Referenced by: '<S6>/Constant2'
                                        */
  real_T Gain_Gain_lr;                 /* Expression: -1
                                        * Referenced by: '<S6>/Gain'
                                        */
  real_T Constant6_Value_g;            /* Expression: 0
                                        * Referenced by: '<S6>/Constant6'
                                        */
  real_T Constant5_Value;              /* Expression: 0
                                        * Referenced by: '<S6>/Constant5'
                                        */
  real_T Gain1_Gain_h;                 /* Expression: -1
                                        * Referenced by: '<S6>/Gain1'
                                        */
  real_T Winglobalaxes_Value[3];       /* Expression: [0 0 m0*g]
                                        * Referenced by: '<S6>/W in global axes'
                                        */
  real_T Constant_Value_c;             /* Expression: 0.5
                                        * Referenced by: '<S46>/Constant'
                                        */
  real_T Gain2_Gain_m;                 /* Expression: 2
                                        * Referenced by: '<S46>/Gain2'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 2
                                        * Referenced by: '<S46>/Gain'
                                        */
  real_T Gain1_Gain_k;                 /* Expression: 2
                                        * Referenced by: '<S46>/Gain1'
                                        */
  real_T Gain_Gain_j;                  /* Expression: 2
                                        * Referenced by: '<S47>/Gain'
                                        */
  real_T Constant_Value_g;             /* Expression: 0.5
                                        * Referenced by: '<S47>/Constant'
                                        */
  real_T Gain2_Gain_mr;                /* Expression: 2
                                        * Referenced by: '<S47>/Gain2'
                                        */
  real_T Gain1_Gain_g;                 /* Expression: 2
                                        * Referenced by: '<S47>/Gain1'
                                        */
  real_T Gain_Gain_d;                  /* Expression: 2
                                        * Referenced by: '<S48>/Gain'
                                        */
  real_T Gain1_Gain_f;                 /* Expression: 2
                                        * Referenced by: '<S48>/Gain1'
                                        */
  real_T Constant_Value_pf;            /* Expression: 0.5
                                        * Referenced by: '<S48>/Constant'
                                        */
  real_T Gain2_Gain_c;                 /* Expression: 2
                                        * Referenced by: '<S48>/Gain2'
                                        */
  real_T Integrator_IC_e2;             /* Expression: 0
                                        * Referenced by: '<S8>/Integrator'
                                        */
  real_T Constant1_Value_k;            /* Expression: 0
                                        * Referenced by: '<S1>/Constant1'
                                        */
  real_T Constant_Value_n;             /* Expression: 0
                                        * Referenced by: '<S1>/Constant'
                                        */
  real_T Constant_Value_m;             /* Expression: 0.5
                                        * Referenced by: '<S19>/Constant'
                                        */
  real_T Gain_Gain_dr;                 /* Expression: 2
                                        * Referenced by: '<S19>/Gain'
                                        */
  real_T Gain1_Gain_py;                /* Expression: 2
                                        * Referenced by: '<S19>/Gain1'
                                        */
  real_T Gain2_Gain_k;                 /* Expression: 2
                                        * Referenced by: '<S19>/Gain2'
                                        */
  real_T Constant_Value_a;             /* Expression: 0.5
                                        * Referenced by: '<S20>/Constant'
                                        */
  real_T Gain_Gain_a;                  /* Expression: 2
                                        * Referenced by: '<S20>/Gain'
                                        */
  real_T Gain1_Gain_p5;                /* Expression: 2
                                        * Referenced by: '<S20>/Gain1'
                                        */
  real_T Gain2_Gain_cm;                /* Expression: 2
                                        * Referenced by: '<S20>/Gain2'
                                        */
  real_T Constant_Value_b;             /* Expression: 0.5
                                        * Referenced by: '<S21>/Constant'
                                        */
  real_T Gain_Gain_o;                  /* Expression: 2
                                        * Referenced by: '<S21>/Gain'
                                        */
  real_T Gain1_Gain_o;                 /* Expression: 2
                                        * Referenced by: '<S21>/Gain1'
                                        */
  real_T Gain2_Gain_h;                 /* Expression: 2
                                        * Referenced by: '<S21>/Gain2'
                                        */
  real_T Constant_Value_hd;            /* Expression: 1
                                        * Referenced by: '<S32>/Constant'
                                        */
  real_T HighGainQuaternionNormalization;/* Expression: 1.0
                                          * Referenced by: '<S32>/High Gain Quaternion Normalization'
                                          */
  real_T n_n_c_lim_UpperSat;           /* Expression: 9
                                        * Referenced by: '<S10>/n_n_c_lim'
                                        */
  real_T n_n_c_lim_LowerSat;           /* Expression: -3
                                        * Referenced by: '<S10>/n_n_c_lim'
                                        */
  real_T Gain_Gain_a5;                 /* Expression: 2
                                        * Referenced by: '<S10>/Gain'
                                        */
  real_T n_xc_lim_UpperSat;            /* Expression: 2
                                        * Referenced by: '<S11>/n_xc_lim'
                                        */
  real_T n_xc_lim_LowerSat;            /* Expression: 0
                                        * Referenced by: '<S11>/n_xc_lim'
                                        */
  real_T Gain_Gain_dw;                 /* Expression: 2
                                        * Referenced by: '<S11>/Gain'
                                        */
  uint32_T THRUST_AB_MAX_RATIO_maxIndex[2];/* Computed Parameter: THRUST_AB_MAX_RATIO_maxIndex
                                            * Referenced by: '<S100>/THRUST_AB_MAX_RATIO'
                                            */
};

/* Real-time Model Data Structure */
struct tag_RTM_f16_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_f16_T *contStates;
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
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_f16_T f16_P;

/* Block signals (auto storage) */
extern B_f16_T f16_B;

/* Continuous states (auto storage) */
extern X_f16_T f16_X;

/* Block states (auto storage) */
extern DW_f16_T f16_DW;

/* External inputs (root inport signals with auto storage) */
extern ExtU_f16_T f16_U;

/* External outputs (root outports fed by signals with auto storage) */
extern ExtY_f16_T f16_Y;

/* Model entry point functions */
extern void f16_initialize(void);
extern void f16_step(void);
extern void f16_terminate(void);

/* Real-time Model object */
extern RT_MODEL_f16_T *const f16_M;

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
 * '<Root>' : 'f16'
 * '<S1>'   : 'f16/3DoF plant'
 * '<S2>'   : 'f16/Angle Conversion'
 * '<S3>'   : 'f16/D_dynamics'
 * '<S4>'   : 'f16/L_and_alpha_dynamics'
 * '<S5>'   : 'f16/T_dynamics'
 * '<S6>'   : 'f16/acck Synthesise'
 * '<S7>'   : 'f16/atmosphere'
 * '<S8>'   : 'f16/mu_dot_dynamics'
 * '<S9>'   : 'f16/mu_dynamics'
 * '<S10>'  : 'f16/n_n_dynamics'
 * '<S11>'  : 'f16/n_x_dynamics'
 * '<S12>'  : 'f16/3DoF plant/Quaternion Conjugate1'
 * '<S13>'  : 'f16/3DoF plant/Quaternion Rotation1'
 * '<S14>'  : 'f16/3DoF plant/Quaternion Rotation2'
 * '<S15>'  : 'f16/3DoF plant/anguler_integration'
 * '<S16>'  : 'f16/3DoF plant/wx(Iw)'
 * '<S17>'  : 'f16/3DoF plant/wx(Iw)1'
 * '<S18>'  : 'f16/3DoF plant/Quaternion Rotation1/Quaternion Normalize'
 * '<S19>'  : 'f16/3DoF plant/Quaternion Rotation1/V1'
 * '<S20>'  : 'f16/3DoF plant/Quaternion Rotation1/V2'
 * '<S21>'  : 'f16/3DoF plant/Quaternion Rotation1/V3'
 * '<S22>'  : 'f16/3DoF plant/Quaternion Rotation1/Quaternion Normalize/Quaternion Modulus'
 * '<S23>'  : 'f16/3DoF plant/Quaternion Rotation1/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S24>'  : 'f16/3DoF plant/Quaternion Rotation2/Quaternion Normalize'
 * '<S25>'  : 'f16/3DoF plant/Quaternion Rotation2/V1'
 * '<S26>'  : 'f16/3DoF plant/Quaternion Rotation2/V2'
 * '<S27>'  : 'f16/3DoF plant/Quaternion Rotation2/V3'
 * '<S28>'  : 'f16/3DoF plant/Quaternion Rotation2/Quaternion Normalize/Quaternion Modulus'
 * '<S29>'  : 'f16/3DoF plant/Quaternion Rotation2/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S30>'  : 'f16/3DoF plant/anguler_integration/Quaternions to Rotation Angles'
 * '<S31>'  : 'f16/3DoF plant/anguler_integration/Rotation Angles to Quaternions1'
 * '<S32>'  : 'f16/3DoF plant/anguler_integration/qdot'
 * '<S33>'  : 'f16/3DoF plant/anguler_integration/Quaternions to Rotation Angles/Quaternion Normalize'
 * '<S34>'  : 'f16/3DoF plant/anguler_integration/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus'
 * '<S35>'  : 'f16/3DoF plant/anguler_integration/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S36>'  : 'f16/3DoF plant/wx(Iw)/Subsystem'
 * '<S37>'  : 'f16/3DoF plant/wx(Iw)/Subsystem1'
 * '<S38>'  : 'f16/3DoF plant/wx(Iw)1/Subsystem'
 * '<S39>'  : 'f16/3DoF plant/wx(Iw)1/Subsystem1'
 * '<S40>'  : 'f16/D_dynamics/Angle Conversion'
 * '<S41>'  : 'f16/L_and_alpha_dynamics/Angle Conversion'
 * '<S42>'  : 'f16/L_and_alpha_dynamics/MATLAB Function'
 * '<S43>'  : 'f16/acck Synthesise/Quaternion Conjugate1'
 * '<S44>'  : 'f16/acck Synthesise/Quaternion Rotation1'
 * '<S45>'  : 'f16/acck Synthesise/Quaternion Rotation1/Quaternion Normalize'
 * '<S46>'  : 'f16/acck Synthesise/Quaternion Rotation1/V1'
 * '<S47>'  : 'f16/acck Synthesise/Quaternion Rotation1/V2'
 * '<S48>'  : 'f16/acck Synthesise/Quaternion Rotation1/V3'
 * '<S49>'  : 'f16/acck Synthesise/Quaternion Rotation1/Quaternion Normalize/Quaternion Modulus'
 * '<S50>'  : 'f16/acck Synthesise/Quaternion Rotation1/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S51>'  : 'f16/atmosphere/Dynamic Pressure'
 * '<S52>'  : 'f16/atmosphere/ISA Atmosphere Model'
 * '<S53>'  : 'f16/atmosphere/Ideal Airspeed Correction'
 * '<S54>'  : 'f16/atmosphere/Mach Number'
 * '<S55>'  : 'f16/atmosphere/Dynamic Pressure/dot'
 * '<S56>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2EAS'
 * '<S57>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2TAS'
 * '<S58>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2CAS'
 * '<S59>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2TAS'
 * '<S60>'  : 'f16/atmosphere/Ideal Airspeed Correction/Mach <= 1.0'
 * '<S61>'  : 'f16/atmosphere/Ideal Airspeed Correction/Pressure Conversion'
 * '<S62>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2CAS'
 * '<S63>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2EAS'
 * '<S64>'  : 'f16/atmosphere/Ideal Airspeed Correction/Velocity Conversion'
 * '<S65>'  : 'f16/atmosphere/Ideal Airspeed Correction/Velocity Conversion1'
 * '<S66>'  : 'f16/atmosphere/Ideal Airspeed Correction/Velocity Conversion2'
 * '<S67>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach'
 * '<S68>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach/Relative Ratio'
 * '<S69>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach/Relative Ratio/Density Conversion'
 * '<S70>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach/Relative Ratio/Pressure Conversion'
 * '<S71>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2EAS/Calculate Mach/Relative Ratio/Temperature Conversion'
 * '<S72>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach'
 * '<S73>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach/Relative Ratio'
 * '<S74>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach/Relative Ratio/Density Conversion'
 * '<S75>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach/Relative Ratio/Pressure Conversion'
 * '<S76>'  : 'f16/atmosphere/Ideal Airspeed Correction/CAS2TAS/Calculate Mach/Relative Ratio/Temperature Conversion'
 * '<S77>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach'
 * '<S78>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach/Relative Ratio'
 * '<S79>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach/Relative Ratio/Density Conversion'
 * '<S80>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach/Relative Ratio/Pressure Conversion'
 * '<S81>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2CAS/Calculate Mach/Relative Ratio/Temperature Conversion'
 * '<S82>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach'
 * '<S83>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach/Relative Ratio'
 * '<S84>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach/Relative Ratio/Density Conversion'
 * '<S85>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach/Relative Ratio/Pressure Conversion'
 * '<S86>'  : 'f16/atmosphere/Ideal Airspeed Correction/EAS2TAS/Calculate Mach/Relative Ratio/Temperature Conversion'
 * '<S87>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2CAS/Relative Ratio'
 * '<S88>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2CAS/Relative Ratio/Density Conversion'
 * '<S89>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2CAS/Relative Ratio/Pressure Conversion'
 * '<S90>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2CAS/Relative Ratio/Temperature Conversion'
 * '<S91>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2EAS/Relative Ratio'
 * '<S92>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2EAS/Relative Ratio/Density Conversion'
 * '<S93>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2EAS/Relative Ratio/Pressure Conversion'
 * '<S94>'  : 'f16/atmosphere/Ideal Airspeed Correction/TAS2EAS/Relative Ratio/Temperature Conversion'
 * '<S95>'  : 'f16/atmosphere/Mach Number/dot'
 * '<S96>'  : 'f16/n_n_dynamics/n_n_lim'
 * '<S97>'  : 'f16/n_n_dynamics/n_n_lo_lim'
 * '<S98>'  : 'f16/n_n_dynamics/n_n_up_lim'
 * '<S99>'  : 'f16/n_x_dynamics/n_x_lim'
 * '<S100>' : 'f16/n_x_dynamics/n_x_up_lim'
 */
#endif                                 /* RTW_HEADER_f16_h_ */
