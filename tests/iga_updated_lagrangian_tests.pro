QT -= gui
CONFIG += console c++17
CONFIG -= app_bundle
TEMPLATE = app
TARGET = iga_updated_lagrangian_tests

QMAKE_CXXFLAGS += -O2
INCLUDEPATH += .. "E:/code/lib/eigen/eigen-3.4.1"

SOURCES += \
    iga_updated_lagrangian_tests.cpp \
    ../nurbs/NurbsAlgorithm.cpp \
    ../nurbs/TensorProductNurbsSpace.cpp \
    ../mechanics/ElementCommon.cpp \
    ../mechanics/HyperelasticMaterial.cpp \
    ../mechanics/SaintVenantKirchhoffMaterial.cpp \
    ../mechanics/NeoHookeanMaterial.cpp \
    ../mechanics/PrincipalStretchMaterial.cpp \
    ../mechanics/TotalLagrangianElement.cpp \
    ../mechanics/UpdatedLagrangianElement.cpp \
    ../mechanics/GlobalNonlinearSystem.cpp \
    ../mechanics/GeneralizedAlphaDynamics.cpp

HEADERS += \
    ../mechanics/ElementCommon.h \
    ../mechanics/HyperelasticMaterial.h \
    ../mechanics/SaintVenantKirchhoffMaterial.h \
    ../mechanics/NeoHookeanMaterial.h \
    ../mechanics/PrincipalStretchMaterial.h \
    ../mechanics/TotalLagrangianElement.h \
    ../mechanics/UpdatedLagrangianElement.h \
    ../mechanics/GlobalNonlinearSystem.h \
    ../mechanics/GeneralizedAlphaDynamics.h
