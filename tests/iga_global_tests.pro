QT -= gui
CONFIG += console c++17
CONFIG -= app_bundle
TEMPLATE = app
TARGET = iga_global_tests

QMAKE_CXXFLAGS += -O2
INCLUDEPATH += .. "E:/code/lib/eigen/eigen-3.4.1"

SOURCES += \
    iga_global_tests.cpp \
    ../nurbs/NurbsAlgorithm.cpp \
    ../nurbs/TensorProductNurbsSpace.cpp \
    ../mechanics/NeumannBoundaryLoad.cpp \
    ../mechanics/SaintVenantKirchhoffMaterial.cpp \
    ../mechanics/TotalLagrangianElement.cpp \
    ../mechanics/GlobalNonlinearSystem.cpp

HEADERS += \
    ../nurbs/ControlPoints.h \
    ../nurbs/KnotVector.h \
    ../nurbs/NurbsAlgorithm.h \
    ../nurbs/TensorProductNurbsSpace.h \
    ../mechanics/ExternalLoad.h \
    ../mechanics/HyperelasticMaterial.h \
    ../mechanics/NeumannBoundaryLoad.h \
    ../mechanics/SaintVenantKirchhoffMaterial.h \
    ../mechanics/TotalLagrangianElement.h \
    ../mechanics/GlobalNonlinearSystem.h
