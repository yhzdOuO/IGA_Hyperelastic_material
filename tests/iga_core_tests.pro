QT -= gui
CONFIG += console c++17
CONFIG -= app_bundle
TEMPLATE = app
TARGET = iga_core_tests

QMAKE_CXXFLAGS += -O2
INCLUDEPATH += .. "E:/code/lib/eigen/eigen-3.4.1"

SOURCES += \
    iga_core_tests.cpp \
    ../nurbs/NurbsAlgorithm.cpp \
    ../nurbs/TensorProductNurbsSpace.cpp \
    ../mechanics/SaintVenantKirchhoffMaterial.cpp \
    ../mechanics/ElementCommon.cpp \
    ../mechanics/HyperelasticMaterial.cpp \
    ../mechanics/UpdatedLagrangianElement.cpp \
    ../mechanics/TotalLagrangianElement.cpp \
    ../mechanics/GlobalNonlinearSystem.cpp

HEADERS += \
    ../nurbs/ControlPoints.h \
    ../nurbs/KnotVector.h \
    ../nurbs/NurbsAlgorithm.h \
    ../nurbs/TensorProductNurbsSpace.h \
    ../mechanics/HyperelasticMaterial.h \
    ../mechanics/SaintVenantKirchhoffMaterial.h \
    ../mechanics/ElementCommon.h \
    ../mechanics/UpdatedLagrangianElement.h \
    ../mechanics/TotalLagrangianElement.h \
    ../mechanics/GlobalNonlinearSystem.h
