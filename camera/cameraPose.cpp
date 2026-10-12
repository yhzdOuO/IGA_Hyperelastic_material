#include "CameraPose.h"

// ===== CameraPose =====
CameraPose::CameraPose()
    : yaw(0),
    pitch(0),
    distance(0),
    center(0, 0, 0) {}

CameraPose::CameraPose(float yaw, float pitch, float distance, QVector3D center)
    : yaw(yaw),
    pitch(pitch),
    distance(distance),
    center(center) {}

void CameraPose::setYaw(float yaw)
{
    if (yaw < 0 || yaw >= 2 * M_PI) {
        yaw -= floor(yaw / (2 * M_PI)) * (2 * M_PI);
    }
    this->yaw = yaw;
}

void CameraPose::setPitch(float pitch)
{
    pitch = qBound(-M_PI / 2, pitch, M_PI / 2); //这个限制的逻辑应该写入到cameraPose里面
    this->pitch = pitch;
}

void CameraPose::setDistance(float distance)
{
    if (distance < 0) {
        distance = 0;
    }
    this->distance = distance;
}

void CameraPose::setCenter(QVector3D center)
{
    this->center = center;
}

float CameraPose::getYaw() const
{
    return yaw;
}

float CameraPose::getPitch() const
{
    return pitch;
}

float CameraPose::getDistance() const
{
    return distance;
}

QVector3D CameraPose::getCenter() const
{
    return center;
}

QVector3D CameraPose::forward() const
{
    float cameraPitch = -pitch;
    float cameraYaw = M_PI + yaw;

    return QVector3D(cos(cameraPitch)*cos(cameraYaw), cos(cameraPitch)*sin(cameraYaw), sin(cameraPitch));
}

QVector3D CameraPose::right() const
{
    float cameraRight = M_PI / 2 + yaw;
    return QVector3D(cos(cameraRight), sin(cameraRight), 0);
}

QVector3D CameraPose::up() const
{
    QVector3D f = forward();
    QVector3D r = right();
    return QVector3D::crossProduct(r, f);
}

QVector3D CameraPose::position() const
{
    return center - forward() * distance;
}

QMatrix4x4 CameraPose::getViewMatrix() const
{
    QMatrix4x4 view;
    view.lookAt(position(), center, up());
    return view;
}

// QMatrix4x4 CameraPose::getCameraCoordinate() const
// {
//     QVector3D f = forward();
//     QVector3D r = right();
//     QVector3D u = up();
//     QMatrix4x4 camera;
//     camera.setColumn(0, QVector4D(r.x(), r.y(), r.z(), 0));
//     camera.setColumn(1, QVector4D(u.x(), u.y(), u.z(), 0));
//     camera.setColumn(2, QVector4D(f.x(), f.y(), f.z(), 0));
//     camera.setColumn(3, QVector4D(0,     0,     0,     1));
//     return camera;
// }
