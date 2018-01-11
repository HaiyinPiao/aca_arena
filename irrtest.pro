TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


#INCLUDEPATH += /home/haiyin/software/matlab/extern/include\
#INCLUDEPATH += /home/haiyin/software/matlab/simulink/include\
#INCLUDEPATH += /home/haiyin/software/matlab/rtw/c/src\
#INCLUDEPATH += /home/haiyin/code_repo/aca_fighter/f16-hrt

INCLUDEPATH += /usr/local/MATLAB/R2016b/extern/include\
INCLUDEPATH += /usr/local/MATLAB/R2016b/simulink/include\
INCLUDEPATH += /usr/local/MATLAB/R2016b/rtw/c/src\
INCLUDEPATH += /home/haiyinpiao/code_repo/aca_arena/f16-hrt

INCLUDEPATH += ../irrlicht-1.8.4/include
LIBS += -L../irrlicht-1.8.4/lib/Linux \
        -lIrrlicht -lGL -lXxf86vm -lXext -lX11

HEADERS += \
    f16-hrt/rtwtypes.h \
    f16-hrt/rtmodel.h \
    f16-hrt/rtGetNaN.h \
    f16-hrt/rtGetInf.h \
    f16-hrt/rt_nonfinite.h \
    f16-hrt/rt_look2d_normal.h \
    f16-hrt/rt_look.h \
    f16-hrt/rt_defines.h \
    f16-hrt/multiword_types.h \
    f16-hrt/f16.h \
    f16-hrt/f16_types.h \
    f16-hrt/f16_private.h \
    f16-hrt/builtin_typeid_types.h

SOURCES += \
    f16-hrt/rtGetNaN.c \
    f16-hrt/rtGetInf.c \
    f16-hrt/rt_nonfinite.c \
    f16-hrt/rt_look2d_normal.c \
    f16-hrt/rt_look.c \
    f16-hrt/f16.c \
    f16-hrt/f16_data.c \
    main.cpp

