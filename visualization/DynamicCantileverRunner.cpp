#include "visualization/DynamicCantileverRunner.h"

#include "mechanics/NeumannBoundaryLoad.h"
#include "mechanics/SaintVenantKirchhoffMaterial.h"

#include <QThread>
#include <algorithm>

namespace visualization {
namespace {

/**
 * @brief 返回悬臂梁上表面随时间变化的竖直参考牵引力。
 *
 * 载荷在 0.5 s 之前从零线性增长到接近峰值，并在 t=0.5 s 时突然撤销。
 * 返回值是单位参考面积上的牵引力，而不是整个上表面的合力。
 */
Eigen::Vector3d topSurfaceTraction(double time)
{
    constexpr double releaseTime = 0.5;
    constexpr double peakMagnitude = 100.0;

    // ---------- 功能块 1：卸载时刻及其以后保持零牵引力 ----------
    // 分支包含等号，使载荷在恰好 0.5 s 时从增长曲线瞬间跳到零。
    if (time >= releaseTime) {
        return Eigen::Vector3d::Zero();
    }

    // ---------- 功能块 2：卸载前按时间线性增加向下牵引力 ----------
    // clamp 同时让求解开始前的负时间安全地保持零载荷。
    const double ramp = std::clamp(time / releaseTime, 0.0, 1.0);
    return -peakMagnitude * ramp * Eigen::Vector3d::UnitZ();
}

} // namespace

void PlaybackHistory::append(PlaybackFrame frame)
{
    auto immutable = std::make_shared<const PlaybackFrame>(std::move(frame));
    frames_.push_back(immutable); // Only the solver thread writes history.
    std::atomic_store_explicit(&latestFrame_, std::move(immutable),
                               std::memory_order_release);
}

bool PlaybackHistory::interpolationFrames(double time, PlaybackFrame& lower,
                                          PlaybackFrame& upper,
                                          double& fraction) const
{
    // The GUI calls interpolation only after finished_ acquires the solver's release.
    if (frames_.empty()) return false;
    if (time <= frames_.front()->time) {
        lower = upper = *frames_.front();
        fraction = 0.0;
        return true;
    }
    if (time >= frames_.back()->time) {
        lower = upper = *frames_.back();
        fraction = 0.0;
        return true;
    }
    auto iterator = std::lower_bound(
        frames_.begin(), frames_.end(), time,
        [](const std::shared_ptr<const PlaybackFrame>& frame, double value) {
            return frame->time < value;
        });
    upper = **iterator;
    lower = **(iterator - 1);
    fraction = (time - lower.time) / (upper.time - lower.time);
    return true;
}

bool PlaybackHistory::latest(PlaybackFrame& frame) const
{
    const std::shared_ptr<const PlaybackFrame> snapshot =
        std::atomic_load_explicit(&latestFrame_, std::memory_order_acquire);
    if (!snapshot) return false;
    frame = *snapshot;
    return true;
}

double PlaybackHistory::endTime() const
{
    return frames_.empty() ? 0.0 : frames_.back()->time;
}

void PlaybackHistory::setFinished(QString error)
{
    error_ = std::move(error);
    finished_.store(true, std::memory_order_release);
}

bool PlaybackHistory::isFinished() const
{
    return finished_.load(std::memory_order_acquire);
}

QString PlaybackHistory::error() const
{
    finished_.load(std::memory_order_acquire);
    return error_;
}

DynamicCantileverRunner::DynamicCantileverRunner(
    nurbs::TensorProductNurbsSpace space, int samplesPerSpan,
    std::shared_ptr<PlaybackHistory> history)
    : space_(std::move(space)),
      samplingCache_(space_, samplesPerSpan),
      history_(std::move(history))
{
}

void DynamicCantileverRunner::run()
{
    try {
        // ---------- 动力学应用层入口 ----------
        // 本层只给出算例参数和回调；质量组装、初始加速度、时间步及
        // Newton 循环全部封装在 mechanics::solveDynamics() 中。
        const iga::SaintVenantKirchhoffMaterial material(0.0, 1.0e6);
        const int dofs = iga::numDofs(space_);
        // ---------- 功能块 1：在 zeta 最大参数面预积分时变 Neumann 载荷 ----------
        // 该面就是梁的上表面；空间积分只执行一次，时间步中只求值上面的牵引函数。
        const iga::TimeDependentLoad surfaceLoad =
            iga::makeReferenceFaceTractionLoad(space_, 2, 1, topSurfaceTraction);
        iga::DynamicsSettings settings;
        settings.timeStep = 4.9e-3;
        // ---------- 功能块 2：卸载后继续推进足够长时间，观察悬臂自由振动 ----------
        settings.numberOfSteps = 500;
        settings.spectralRadius = 0.8;
        settings.maxIterations = 20;
        settings.printIterations = false;
        int step = 0;
        iga::solveDynamics(
            space_, material, 1.0, Eigen::VectorXd::Zero(dofs),
            Eigen::VectorXd::Zero(dofs), iga::collectXiMinFaceDofs(space_),
            surfaceLoad, settings,
            [this, &material, &step](const iga::DynamicState& state) {
                // 回调收到的只可能是 t_n 的已收敛物理状态，不包含
                // d_(n+alpha_f) 或任何 Newton 中间构形。
                PlaybackFrame frame;
                frame.time = state.time;
                frame.displacements = state.displacements;
                frame.scalars = samplingCache_.evaluateScalars(
                    state.displacements, material);
                frame.kineticEnergy = state.kineticEnergy;
                frame.strainEnergy = state.strainEnergy;
                history_->append(std::move(frame));
                emit stepFinished(step++, state.time);
                return !QThread::currentThread()->isInterruptionRequested();
            });
        history_->setFinished();
    } catch (const std::exception& exception) {
        history_->setFinished(QString::fromUtf8(exception.what()));
    }
    emit simulationFinished();
}

} // namespace visualization
