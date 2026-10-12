#ifndef VISUALIZATION_DYNAMICCANTILEVERRUNNER_H
#define VISUALIZATION_DYNAMICCANTILEVERRUNNER_H

#include "mechanics/GeneralizedAlphaDynamics.h"
#include "visualization/BoundarySamplingCache.h"

#include <QObject>
#include <atomic>
#include <memory>

namespace visualization {

struct PlaybackFrame {
    double time = 0.0;
    Eigen::VectorXd displacements;
    std::vector<SampleScalars> scalars;
    double kineticEnergy = 0.0;
    double strainEnergy = 0.0;
};

/** @brief Immutable history plus an atomically exchanged live-frame snapshot. */
class PlaybackHistory {
public:
    void append(PlaybackFrame frame);
    bool interpolationFrames(double time, PlaybackFrame& lower,
                             PlaybackFrame& upper, double& fraction) const;
    bool latest(PlaybackFrame& frame) const;
    double endTime() const;
    void setFinished(QString error = {});
    bool isFinished() const;
    QString error() const;

private:
    std::vector<std::shared_ptr<const PlaybackFrame>> frames_;
    std::shared_ptr<const PlaybackFrame> latestFrame_;
    std::atomic<bool> finished_{false};
    QString error_;
};

/** @brief Runs the nonlinear cantilever time loop entirely outside the GUI thread. */
class DynamicCantileverRunner : public QObject {
    Q_OBJECT
public:
    DynamicCantileverRunner(nurbs::TensorProductNurbsSpace space,
                            int samplesPerSpan,
                            std::shared_ptr<PlaybackHistory> history);

public slots:
    void run();

signals:
    void stepFinished(int step, double time);
    void simulationFinished();

private:
    nurbs::TensorProductNurbsSpace space_;
    BoundarySamplingCache samplingCache_;
    std::shared_ptr<PlaybackHistory> history_;
};

} // namespace visualization

#endif // VISUALIZATION_DYNAMICCANTILEVERRUNNER_H
