QT       += core gui \
            openglwidgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# Eigen-heavy validation builds still need optimization (see the project spec).
QMAKE_CXXFLAGS_DEBUG += -O2

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    camera/camera.cpp \
    camera/cameraAnimator.cpp \
    camera/cameraMouseCtrl.cpp \
    camera/cameraPose.cpp \
    camera/cameraProjection.cpp \
    nurbs/NurbsAlgorithm.cpp \
    nurbs/TensorProductNurbsSpace.cpp \
    mechanics/NeumannBoundaryLoad.cpp \
    mechanics/SaintVenantKirchhoffMaterial.cpp \
    mechanics/NeoHookeanMaterial.cpp \
    mechanics/PrincipalStretchMaterial.cpp \
    mechanics/ElementCommon.cpp \
    mechanics/HyperelasticMaterial.cpp \
    mechanics/UpdatedLagrangianElement.cpp \
    mechanics/TotalLagrangianElement.cpp \
    mechanics/GlobalNonlinearSystem.cpp \
    mechanics/GeneralizedAlphaDynamics.cpp \
    visualization/BoundarySampler.cpp \
    visualization/BoundarySamplingCache.cpp \
    visualization/CantileverExample.cpp \
    visualization/DynamicCantileverRunner.cpp \
    visualization/FieldEvaluator.cpp \
    visualization/IgaVisualizationWidget.cpp

HEADERS += \
    mainwindow.h \
    camera/camera.h \
    camera/cameraAnimator.h \
    camera/cameraMouseCtrl.h \
    camera/cameraPose.h \
    camera/cameraProjection.h \
    nurbs/ControlPoints.h \
    nurbs/KnotVector.h \
    nurbs/NurbsAlgorithm.h \
    nurbs/TensorProductNurbsSpace.h \
    mechanics/ExternalLoad.h \
    mechanics/HyperelasticMaterial.h \
    mechanics/NeumannBoundaryLoad.h \
    mechanics/SaintVenantKirchhoffMaterial.h \
    mechanics/NeoHookeanMaterial.h \
    mechanics/PrincipalStretchMaterial.h \
    mechanics/ElementCommon.h \
    mechanics/UpdatedLagrangianElement.h \
    mechanics/TotalLagrangianElement.h \
    mechanics/GlobalNonlinearSystem.h \
    mechanics/GeneralizedAlphaDynamics.h \
    visualization/BoundarySampler.h \
    visualization/BoundarySamplingCache.h \
    visualization/CantileverExample.h \
    visualization/DynamicCantileverRunner.h \
    visualization/FieldEvaluator.h \
    visualization/IgaVisualizationWidget.h \
    visualization/SurfaceMesh.h

FORMS += \
    mainwindow.ui

RESOURCES += visualization/visualization.qrc

INCLUDEPATH += "E:/code/lib/eigen/eigen-3.4.1"

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
