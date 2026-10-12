#ifndef CAMERA_ANIMATOR_H
#define CAMERA_ANIMATOR_H

#include "cameraPose.h"

class CameraAnimator
{
public:
    enum class AnimationType {
        Reset,
        PerspectiveUp,
        PerspectiveDown,
        PerspectiveForward,
        PerspectiveBack,
        PerspectiveLeft,
        PerspectiveRight
    };

public:
    CameraAnimator(CameraPose* camera);

    bool isRunning() const;

    void setAnimation(AnimationType type, int duration);
    void createAnimation(const CameraPose& src, const CameraPose& dst, int duration);
    void stop();

    void updateTime(int dt);     // 原 timeElapse()

private:
    float evaluateProgress(float elapsed, float totalDuration) const;
    void applyProgress(float progress); // 原 update()

private:
    // Animation timing
    bool m_running = false;
    int m_totalDuration = 0.0f;
    int m_elapsedTime = 0.0f;

    // Animation curve
    float m_curvePower = 3.0f;

    // Animation state
    CameraPose m_srcPose;
    CameraPose m_deltaPose;

    // Bound camera
    CameraPose* m_cameraPose = nullptr;
};

#endif // CAMERA_ANIMATOR_H
