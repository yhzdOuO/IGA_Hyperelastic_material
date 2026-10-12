#include "cameraMouseCtrl.h"

CameraMouseController::CameraMouseController(CameraPose* camera)
    : m_camera(camera)
{
}

bool CameraMouseController::isDragging() const
{
    return m_pressedButton != Qt::NoButton;
}

void CameraMouseController::setViewportSize(int width, int height)
{
    m_viewportWidth = width;
    m_viewportHeight = height;
}

void CameraMouseController::bind(CameraPose* camera)
{
    m_camera = camera;
}

void CameraMouseController::handleEvent(QEvent* event)
{
    switch (event->type()) {
    case QEvent::MouseButtonPress:
        mousePressEvent(static_cast<QMouseEvent*>(event));
        break;
    case QEvent::MouseButtonRelease:
        mouseReleaseEvent(static_cast<QMouseEvent*>(event));
        break;
    case QEvent::MouseMove:
        mouseMoveEvent(static_cast<QMouseEvent*>(event));
        break;
    case QEvent::Wheel:
        wheelEvent(static_cast<QWheelEvent*>(event));
        break;
    default:
        break;
    }
}

void CameraMouseController::mousePressEvent(QMouseEvent* event)
{
    m_pressedButton = event->button();

    if (m_pressedButton == Qt::LeftButton || m_pressedButton == Qt::RightButton) {
        m_lastX = event->pos().x();
        m_lastY = event->pos().y();
        m_lastPose = *m_camera; // store initial pose
    }
}

void CameraMouseController::mouseReleaseEvent(QMouseEvent* /*event*/)
{
    m_pressedButton = Qt::NoButton;
}

void CameraMouseController::mouseMoveEvent(QMouseEvent* event)
{
    if (m_pressedButton != Qt::LeftButton &&
        m_pressedButton != Qt::RightButton) {
        return;
    }

    float currX = event->pos().x();
    float currY = event->pos().y();

    // Normalized delta
    float dx = (currX - m_lastX) / float(m_viewportWidth);
    float dy = (currY - m_lastY) / float(m_viewportHeight);

    // Invert Y (screen coords: top->down)
    dy = -dy;

    if (m_pressedButton == Qt::LeftButton) {
        constexpr float kBaseSpeed = 4.0f;

        float yaw = m_lastPose.yaw + kBaseSpeed * m_sensitivity * -dx;
        float pitch = m_lastPose.pitch + kBaseSpeed * m_sensitivity * -dy;

        m_camera->setYaw(yaw);
        m_camera->setPitch(pitch);

    } else if (m_pressedButton == Qt::RightButton) {
        constexpr float kBaseSpeed = 3.0f;

        QVector3D delta = (kBaseSpeed * m_sensitivity) *
            (-dx * m_lastPose.right() + -dy * m_lastPose.up());

        m_camera->setCenter(m_lastPose.center + delta);
    }
}

void CameraMouseController::wheelEvent(QWheelEvent* event)
{
    constexpr float kBaseZoomSpeed = -2e-3f;

    float d = m_camera->getDistance() +
              kBaseZoomSpeed * m_sensitivity * event->angleDelta().y();

    m_camera->setDistance(d);
}
