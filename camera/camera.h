#ifndef CAMERA_H
#define CAMERA_H

#include "cameraAnimator.h"
#include "cameraMouseCtrl.h"
#include "cameraProjection.h"

class Camera
{
public:
    Camera();
    void setProjectionMode(enum CameraProjection::ProjectionMode mode);
    void setVerticalAngle(float angle);
    void setViewportSize(int width, int height);
    void setNearAndFar(float nearPlane, float farPlane);
    void setPose(float yaw, float pitch, float distance, const QVector3D& center);
    QMatrix4x4 getViewMatrix() const;
    QMatrix4x4 getProjectionMatrix() const;
    void handleEvent(QEvent *event);
    void setAnimation(enum CameraAnimator::AnimationType animationType);
    void timeElapse(int t);

private:
    CameraProjection projection;
    CameraPose pose;
    CameraAnimator animator;
    CameraMouseController mouse;
};

#endif // CAMERA_H
