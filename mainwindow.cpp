#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "mechanics/SaintVenantKirchhoffMaterial.h"
#include "visualization/BoundarySamplingCache.h"
#include "visualization/CantileverExample.h"
#include "visualization/DynamicCantileverRunner.h"
#include "visualization/IgaVisualizationWidget.h"

#include <QLabel>
#include <QKeySequence>
#include <QShortcut>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>
#include <algorithm>

MainWindow::MainWindow(int samplesPerSpan, SimulationMode mode, QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle(mode == SimulationMode::Static
                       ? tr("IGA cantilever statics")
                       : tr("IGA generalized-alpha cantilever dynamics"));
    resize(1180, 720);
    try {
        // 两种模式共用同一个 OpenGL 查看器；区别只在求解入口和数据更新方式。
        if (mode == SimulationMode::Static) {
            // ---------- 默认路径：启动时同步完成一次静力 Newton 求解 ----------
            visualization::CantileverVisualization result =
                visualization::solveCantileverVisualization(samplesPerSpan);
            viewer_ = new visualization::IgaVisualizationWidget(
                std::move(result.mesh), ui->centralwidget);
            auto* layout = new QVBoxLayout(ui->centralwidget);
            layout->setContentsMargins(0, 0, 0, 0);
            layout->addWidget(viewer_);
            connect(viewer_, &visualization::IgaVisualizationWidget::statusChanged,
                    this, [this](const QString& message) {
                        ui->statusbar->showMessage(message);
                    });
            ui->statusbar->showMessage(
                tr("Static solution: tip uz=%1, Newton iterations=%2")
                    .arg(result.tipDisplacement, 0, 'e', 4)
                    .arg(result.newtonIterations));
            return;
        }

        // ---------- 可选路径：GUI 线程准备动力学几何、采样缓存和初始帧 ----------
        const nurbs::TensorProductNurbsSpace space =
            visualization::makeCantileverSpace();
        const iga::SaintVenantKirchhoffMaterial material(0.0, 5.0e6);
        samplingCache_ = std::make_unique<visualization::BoundarySamplingCache>(
            space, samplesPerSpan);
        const Eigen::VectorXd zero = Eigen::VectorXd::Zero(iga::numDofs(space));
        visualization::SurfaceMesh initialMesh = samplingCache_->makeFrame(
            zero, samplingCache_->evaluateScalars(zero, material));
        viewer_ = new visualization::IgaVisualizationWidget(
            std::move(initialMesh), ui->centralwidget);
        auto* layout = new QVBoxLayout(ui->centralwidget);
        layout->setContentsMargins(0, 0, 0, 0);
        layout->addWidget(viewer_);
        connect(viewer_, &visualization::IgaVisualizationWidget::statusChanged,
                this, [this](const QString& message) {
                    ui->statusbar->showMessage(message);
                });

        // ---------- 调用链 2：求解对象移入后台线程 ----------
        // QThread::started 最终进入 DynamicCantileverRunner::run()，再调用
        // iga::solveDynamics()；GUI 线程自身不做单元装配或 Newton 迭代。
        playbackHistory_ = std::make_shared<visualization::PlaybackHistory>();
        solverThread_ = new QThread(this);
        auto* runner = new visualization::DynamicCantileverRunner(
            space, samplesPerSpan, playbackHistory_);
        runner->moveToThread(solverThread_);
        connect(solverThread_, &QThread::started,
                runner, &visualization::DynamicCantileverRunner::run);
        connect(runner, &visualization::DynamicCantileverRunner::stepFinished,
                this, [this](int step, double time) {
                    showLatestDynamicFrame();
                    ui->statusbar->showMessage(
                        tr("Solving dynamics: step %1, t=%2 s")
                            .arg(step).arg(time, 0, 'f', 4));
                });
        connect(runner, &visualization::DynamicCantileverRunner::simulationFinished,
                this, &MainWindow::startDynamicPlayback);
        connect(runner, &visualization::DynamicCantileverRunner::simulationFinished,
                solverThread_, &QThread::quit);
        connect(runner, &visualization::DynamicCantileverRunner::simulationFinished,
                runner, &QObject::deleteLater);

        playbackTimer_ = new QTimer(this);
        playbackTimer_->setInterval(16);
        connect(playbackTimer_, &QTimer::timeout, this, [this]() {
            if (playbackReady_ && playbackHistory_->isFinished()
                && playbackHistory_->error().isEmpty()) {
                showPlaybackFrame();
            }
        });
        playbackTimer_->start();

        // P 可在播放结束后重放，也可在播放过程中从 t=0 重新开始。
        auto* replayShortcut = new QShortcut(QKeySequence(Qt::Key_P), this);
        replayShortcut->setContext(Qt::WindowShortcut);
        connect(replayShortcut, &QShortcut::activated,
                this, &MainWindow::startDynamicPlayback);

        ui->statusbar->showMessage(tr("Dynamics solver started in background"));
        solverThread_->start();
    } catch (const std::exception& error) {
        auto* layout = new QVBoxLayout(ui->centralwidget);
        auto* label = new QLabel(tr("Unable to build cantilever visualization:\n%1")
                                     .arg(QString::fromUtf8(error.what())),
                                 ui->centralwidget);
        label->setAlignment(Qt::AlignCenter);
        layout->addWidget(label);
    }
}

MainWindow::~MainWindow()
{
    if (solverThread_ && solverThread_->isRunning()) {
        solverThread_->requestInterruption();
        solverThread_->quit();
        solverThread_->wait();
    }
    delete ui;
}

void MainWindow::startDynamicPlayback()
{
    if (!playbackHistory_ || !playbackTimer_) return;
    if (!playbackHistory_->isFinished()) {
        ui->statusbar->showMessage(
            tr("Replay is available after the dynamics simulation finishes"));
        return;
    }
    const QString error = playbackHistory_->error();
    if (!error.isEmpty()) {
        playbackReady_ = false;
        playbackTimer_->stop();
        ui->statusbar->showMessage(tr("Dynamics failed: %1").arg(error));
        return;
    }
    if (playbackHistory_->endTime() <= 0.0) {
        playbackReady_ = false;
        playbackTimer_->stop();
        ui->statusbar->showMessage(tr("Dynamics produced no replayable time interval"));
        return;
    }
    playbackClock_.restart();
    playbackReady_ = true;
    playbackTimer_->start();
    ui->statusbar->showMessage(
        tr("Playing dynamics at physical speed; press P to restart"));
}

void MainWindow::showLatestDynamicFrame()
{
    visualization::PlaybackFrame frame;
    if (!playbackHistory_ || !playbackHistory_->latest(frame)) return;
    viewer_->setDynamicFrame(
        samplingCache_->makeFrame(frame.displacements, frame.scalars), frame.time);
}

void MainWindow::showPlaybackFrame()
{
    const double endTime = playbackHistory_->endTime();
    if (endTime <= 0.0) return;
    // Do not wrap the last converged state directly to t=0. The endpoints of a
    // damped/non-periodic response do not match, and that artificial jump looked
    // like a solver blow-up just after the recorded 123 steps.
    const double elapsedTime = 1.0e-3 * playbackClock_.elapsed();
    const double time = std::min(elapsedTime, endTime);
    visualization::PlaybackFrame lower;
    visualization::PlaybackFrame upper;
    double fraction = 0.0;
    if (!playbackHistory_->interpolationFrames(time, lower, upper, fraction)) return;
    const Eigen::VectorXd displacement =
        (1.0 - fraction) * lower.displacements + fraction * upper.displacements;
    std::vector<visualization::SampleScalars> scalars(lower.scalars.size());
    for (std::size_t i = 0; i < scalars.size(); ++i) {
        for (std::size_t field = 0; field < scalars[i].size(); ++field) {
            scalars[i][field] = (1.0 - fraction) * lower.scalars[i][field]
                + fraction * upper.scalars[i][field];
        }
    }
    viewer_->setDynamicFrame(samplingCache_->makeFrame(displacement, scalars), time);
    if (elapsedTime >= endTime) {
        playbackReady_ = false;
        playbackTimer_->stop();
        ui->statusbar->showMessage(
            tr("Playback reached the final configured time step at t=%1 s; "
               "press P to replay")
                .arg(endTime, 0, 'f', 4));
    }
}
