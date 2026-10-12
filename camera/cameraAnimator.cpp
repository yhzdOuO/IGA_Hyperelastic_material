#include "cameraAnimator.h"
#include <QtMath>

CameraAnimator::CameraAnimator(CameraPose* camera)
    : m_running(false),
    m_totalDuration(0.0f),
    m_elapsedTime(0.0f),
    m_curvePower(3.0f),
    m_cameraPose(camera)
{
}

bool CameraAnimator::isRunning() const
{
    return m_running;
}

void CameraAnimator::setAnimation(AnimationType type, int duration)
{
    if (!m_cameraPose || m_running) return;

    const CameraPose& src = *m_cameraPose;
    CameraPose dst = src;  // make a copy

    switch (type) {
    case AnimationType::Reset:
        dst.yaw = 0;
        dst.pitch = M_PI / 8;
        dst.distance = 2.5f;
        dst.center = QVector3D(0, 0, 0.5f);
        break;

    case AnimationType::PerspectiveUp:
        dst.yaw = 0;
        dst.pitch = M_PI / 2;
        break;

    case AnimationType::PerspectiveDown:
        dst.yaw = 0;
        dst.pitch = -M_PI / 2;
        break;

    case AnimationType::PerspectiveLeft:
        dst.yaw = M_PI * 3 / 2;
        dst.pitch = 0;
        break;

    case AnimationType::PerspectiveRight:
        dst.yaw = M_PI / 2;
        dst.pitch = 0;
        break;

    case AnimationType::PerspectiveForward:
        dst.yaw = 0;
        dst.pitch = 0;
        break;

    case AnimationType::PerspectiveBack:
        dst.yaw = M_PI;
        dst.pitch = 0;
        break;
    }

    createAnimation(src, dst, duration);
}

void CameraAnimator::createAnimation(const CameraPose& src, const CameraPose& dst, int duration)
{
    m_srcPose = src;

    if (duration < 0)
        duration = 0;

    m_totalDuration = duration;
    m_elapsedTime = 0;
    m_running = true;

    // calculate delta (shortest yaw rotation)
    float yaw = dst.yaw - src.yaw;
    if (yaw < -M_PI) yaw += 2 * M_PI;
    else if (yaw > M_PI) yaw -= 2 * M_PI;

    m_deltaPose.yaw = yaw;
    m_deltaPose.pitch = dst.pitch - src.pitch;
    m_deltaPose.distance = dst.distance - src.distance;
    m_deltaPose.center = dst.center - src.center;

    if (duration == 0) {
        m_running = false;
        applyProgress(1.0f);
    }
}

void CameraAnimator::stop()
{
    m_running = false;
}

void CameraAnimator::updateTime(int dt)
{
    if (!m_running) return;

    m_elapsedTime += dt;
    if (m_elapsedTime >= m_totalDuration) {
        m_elapsedTime = m_totalDuration;
        m_running = false;
    }

    float progress = evaluateProgress(m_elapsedTime, m_totalDuration);
    applyProgress(progress);
}

float CameraAnimator::evaluateProgress(float elapsed, float totalDuration) const
{
    float t = 1.0f * elapsed / totalDuration;

    if (t <= 0.0f) return 0.0f;
    if (t >= 1.0f) return 1.0f;

    // classic ease-in-out curve
    if (t < 0.5f) {
        return qPow(2 * t, m_curvePower) / 2.0f;
    } else {
        return 1 - qPow(2 * (1 - t), m_curvePower) / 2.0f;
    }
}

void CameraAnimator::applyProgress(float progress)
{
    m_cameraPose->setYaw    (m_srcPose.yaw     + m_deltaPose.yaw     * progress);
    m_cameraPose->setPitch  (m_srcPose.pitch   + m_deltaPose.pitch   * progress);
    m_cameraPose->setDistance(m_srcPose.distance + m_deltaPose.distance * progress);
    m_cameraPose->setCenter (m_srcPose.center  + m_deltaPose.center  * progress);
}
