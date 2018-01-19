TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


INCLUDEPATH += /home/haiyin/software/matlab/extern/include\
INCLUDEPATH += /home/haiyin/software/matlab/simulink/include\
INCLUDEPATH += /home/haiyin/software/matlab/rtw/c/src\
INCLUDEPATH += /home/haiyin/code_repo/aca_arena/actor_grt_rtw

#INCLUDEPATH += /usr/local/MATLAB/R2016b/extern/include\
#INCLUDEPATH += /usr/local/MATLAB/R2016b/simulink/include\
#INCLUDEPATH += /usr/local/MATLAB/R2016b/rtw/c/src\
#INCLUDEPATH += /home/haiyinpiao/code_repo/aca_arena/actor_grt_rtw

INCLUDEPATH += ../irrlicht-1.8.4/include
LIBS += -L../irrlicht-1.8.4/lib/Linux \
        -lIrrlicht -lGL -lXxf86vm -lXext -lX11

HEADERS += \
    actor_grt_rtw/actor.h \
    actor_grt_rtw/actor_private.h \
    actor_grt_rtw/actor_types.h \
    actor_grt_rtw/multiword_types.h \
    actor_grt_rtw/rt_defines.h \
    actor_grt_rtw/rtGetInf.h \
    actor_grt_rtw/rtGetNaN.h \
    actor_grt_rtw/rt_look.h \
    actor_grt_rtw/rt_look2d_normal.h \
    actor_grt_rtw/rtmodel.h \
    actor_grt_rtw/rt_nonfinite.h \
    actor_grt_rtw/rtwtypes.h

SOURCES += \
    main.cpp \
    actor_grt_rtw/actor.cpp \
    actor_grt_rtw/rtGetInf.cpp \
    actor_grt_rtw/rtGetNaN.cpp \
    actor_grt_rtw/rt_look.cpp \
    actor_grt_rtw/rt_look2d_normal.cpp \
    actor_grt_rtw/rt_nonfinite.cpp

