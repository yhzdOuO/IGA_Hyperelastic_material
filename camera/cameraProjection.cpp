#include "cameraProjection.h"

CameraProjection::CameraProjection() {
    verticalAngle = 45;
    aspectRatio = 1.0;
    nearPlane = 0.1;
    farPlane = 30;
    mode = ProjectionMode::Perspective;
}

void CameraProjection::setProjectionMode(ProjectionMode mode)
{
    this->mode = mode;
}

void CameraProjection::setVerticalAngle(float verticalAngle)
{
    this->verticalAngle = verticalAngle;
}

void CameraProjection::setAspectRatio(float aspectRatio)
{
    this->aspectRatio = aspectRatio;
}

void CameraProjection::setNearAndFar(float near, float far)
{
    nearPlane = near;
    farPlane = far;
}

QMatrix4x4 CameraProjection::getProjectionMatrix() const
{
    QMatrix4x4 projection;
    if (mode == ProjectionMode::Perspective) {
        projection.perspective(verticalAngle, aspectRatio, nearPlane, farPlane);
    } else {
        const float halfHeight = defaultDistance;
        projection.ortho(-halfHeight * aspectRatio, halfHeight * aspectRatio,
                         -halfHeight, halfHeight, nearPlane, farPlane);
    }
    return projection;
}
