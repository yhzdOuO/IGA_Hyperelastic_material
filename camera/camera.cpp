#include "camera.h"

Camera::Camera()
    :animator(&pose),
    mouse(&pose)
{
    animator.setAnimation(CameraAnimator::AnimationType::Reset, 0);
}

void Camera::setProjectionMode(CameraProjection::ProjectionMode mode)
{
    projection.setProjectionMode(mode);
}

void Camera::setVerticalAngle(float angle)
{
    projection.setVerticalAngle(angle);
}

void Camera::handleEvent(QEvent *event)
{
    if (!animator.isRunning()) {
        mouse.handleEvent(event);
    }
}

void Camera::setAnimation(CameraAnimator::AnimationType animationType)
{
    animator.setAnimation(animationType, 750);
}

void Camera::timeElapse(int t)
{
    animator.updateTime(t);
}

void Camera::setViewportSize(int width, int height)
{
    mouse.setViewportSize(width, height);
    if (height > 0) {
        projection.setAspectRatio(1.0f * width / height);
    }
}

void Camera::setNearAndFar(float nearPlane, float farPlane)
{
    projection.setNearAndFar(nearPlane, farPlane);
}

void Camera::setPose(float yaw, float pitch, float distance, const QVector3D& center)
{
    pose.setYaw(yaw);
    pose.setPitch(pitch);
    pose.setDistance(distance);
    pose.setCenter(center);
}

QMatrix4x4 Camera::getViewMatrix() const
{
    return pose.getViewMatrix();
}

QMatrix4x4 Camera::getProjectionMatrix() const
{
    return projection.getProjectionMatrix();
}
