TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += $$PWD \
    $$PWD/src/ \
    $$PWD/src/Limiter/ \
    $$PWD/src/Flux/ \
    $$PWD/src/Problem/ \
    $$PWD/src/Timestep/ \
    $$PWD/src/Indicator/ \
    $$PWD/src/Boundary/ \
    $$PWD/src/Integrator/ \
    $$PWD/src/Parameters/ \
    $$PWD/src/Mesh/ \
    $$PWD/src/Definitions/

SOURCES += main.cpp \
    src/Boundary/Boundary.cpp \
    src/Boundary/BoundaryPeriodic.cpp \
    src/Boundary/BoundarySoft.cpp \
    src/Boundary/BoundaryWall.cpp \
    src/Definitions/defs.cpp \
    src/Flux/Flux.cpp \
    src/Flux/FluxGodunovType.cpp \
    src/Flux/FluxHLL.cpp \
    src/Flux/FluxHLLC.cpp \
    src/Flux/FluxLaxFriedrichs.cpp \
    src/Indicator/Indicator.cpp \
    src/Indicator/IndicatorEverywhere.cpp \
    src/Indicator/IndicatorHarten.cpp \
    src/Indicator/IndicatorKrivodonova.cpp \
    src/Indicator/IndicatorNowhere.cpp \
    src/Limiter/Limiter.cpp \
    src/Limiter/LimiterFinDiff.cpp \
    src/Limiter/LimiterHWENO_SC.cpp \
    src/Limiter/LimiterHWENO_SC_Char.cpp \
    src/Limiter/LimiterWENO.cpp \
    src/Limiter/LimiterWENO_S.cpp \
    src/Mesh/Mesh1D.cpp \
    src/Parameters/Params.cpp \
    src/Problem/Problem.cpp \
    src/Problem/ProblemGas1D.cpp \
    src/Problem/ProblemMHD1D.cpp \
    src/Problem/ProblemTransfer1D.cpp \
    src/Timestep/Timestep.cpp \
    src/Timestep/TimestepEuler.cpp \
    src/Timestep/TimestepRK2.cpp \
    src/Timestep/TimestepRK2TVD.cpp \
    src/Timestep/TimestepRK3TVD.cpp \
    src/Limiter/LimiterHWENO.cpp

HEADERS += \
    src/Boundary/Boundary.h \
    src/Boundary/BoundaryPeriodic.h \
    src/Boundary/BoundarySoft.h \
    src/Boundary/BoundaryWall.h \
    src/Definitions/defs.h \
    src/Flux/Flux.h \
    src/Flux/FluxGodunovType.h \
    src/Flux/FluxHLL.h \
    src/Flux/FluxHLLC.h \
    src/Flux/FluxLaxFriedrichs.h \
    src/Indicator/Indicator.h \
    src/Indicator/IndicatorEverywhere.h \
    src/Indicator/IndicatorHarten.h \
    src/Indicator/IndicatorKrivodonova.h \
    src/Indicator/IndicatorNowhere.h \
    src/Integrator/Integrator.h \
    src/Limiter/Limiter.h \
    src/Limiter/LimiterFinDiff.h \
    src/Limiter/LimiterHWENO_SC.h \
    src/Limiter/LimiterHWENO_SC_Char.h \
    src/Limiter/LimiterWENO_S.h \
    src/Mesh/Mesh1D.h \
    src/Parameters/Params.h \
    src/Problem/Problem.h \
    src/Problem/ProblemGas1D.h \
    src/Problem/ProblemMHD1D.h \
    src/Problem/ProblemTransfer1D.h \
    src/Timestep/Timestep.h \
    src/Timestep/TimestepEuler.h \
    src/Timestep/TimestepRK2.h \
    src/Timestep/TimestepRK2TVD.h \
    src/Timestep/TimestepRK3TVD.h \
    src/Limiter/LimiterHWENO.h \
    src/Limiter/LimiterWENO.h \
    src/numvector/numvector.h
