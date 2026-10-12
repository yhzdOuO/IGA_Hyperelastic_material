#ifndef CAMERAPOSE_H
#define CAMERAPOSE_H

#include <QVector3D>
#include <QMatrix4x4>

// 相机姿态
class CameraPose {
public:
    CameraPose();
    // CameraPose(const CameraPose& pose);
    CameraPose(float yaw, float pitch, float distance, QVector3D center);

    void setYaw(float yaw);
    void setPitch(float pitch);
    void setDistance(float distance);
    void setCenter(QVector3D center);

    float getYaw() const;
    float getPitch() const;
    float getDistance() const;
    QVector3D getCenter() const;

    QVector3D forward() const;
    QVector3D right() const;
    QVector3D up() const;
    QVector3D position() const;
    QMatrix4x4 getViewMatrix() const;
    // QMatrix4x4 getCameraCoordinate() const;

    float yaw;
    float pitch;
    float distance;
    QVector3D center;
};

#endif // CAMERAPOSE_H
