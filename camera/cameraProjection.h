#ifndef CAMERAPROJECTION_H
#define CAMERAPROJECTION_H

#include <QMatrix4x4>

#define defaultDistance 3


class CameraProjection
{
public:
    /* 投影模式 */
    enum class ProjectionMode {
        Perspective,
        Orthogonal
    };

public:
    CameraProjection();
    void setProjectionMode(ProjectionMode mode);
    /* 操作width和height */
    void setVerticalAngle(float verticalAngle);
    void setAspectRatio(float aspectRatio);
    void setNearAndFar(float near, float far);
    /* 相机矩阵 */
    QMatrix4x4 getProjectionMatrix() const;

private:
    /* 相机属性 */
    float verticalAngle, aspectRatio, nearPlane, farPlane;
    ProjectionMode mode;
};

#endif // CAMERAPROJECTION_H
