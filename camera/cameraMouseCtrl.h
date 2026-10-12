#ifndef CAMERA_CONTROLLER_H
#define CAMERA_CONTROLLER_H

#include <QMouseEvent>
#include "cameraPose.h"

class CameraMouseController
{
public:
    explicit CameraMouseController(CameraPose* camera);

    bool isDragging() const;
    void setViewportSize(int width, int height);
    void handleEvent(QEvent* event);
    void bind(CameraPose* camera);

private:
    // Qt event handlers
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);

private:
    // Drag state
    Qt::MouseButton m_pressedButton = Qt::NoButton;

    // Viewport size
    int m_viewportWidth = 0;
    int m_viewportHeight = 0;

    // Previous mouse position
    int m_lastX = 0;
    int m_lastY = 0;

    // Sensitivity
    float m_sensitivity = 1.0f;

    // Pose caching
    CameraPose m_lastPose;

    // Bound camera
    CameraPose* m_camera = nullptr;
};

#endif // CAMERA_CONTROLLER_H
