#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QElapsedTimer>
#include <memory>

class QThread;
class QTimer;
namespace visualization {
class BoundarySamplingCache;
class IgaVisualizationWidget;
class PlaybackHistory;
}

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

/** @brief 程序启动后执行的求解类型。 */
enum class SimulationMode {
    Static,
    Dynamic
};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    /**
     * @brief 创建主窗口，并按 mode 启动静力或动力学求解。
     * @param samplesPerSpan 每个边界 knot span 的可视化采样数。
     * @param mode 求解模式；默认静力学。
     * @param parent 可选父控件。
     */
    MainWindow(int samplesPerSpan = 4,
               SimulationMode mode = SimulationMode::Static,
               QWidget *parent = nullptr);

    /** @brief Releases the generated user-interface object. */
    ~MainWindow();

private:
    void startDynamicPlayback();
    void showLatestDynamicFrame();
    void showPlaybackFrame();

    Ui::MainWindow *ui;
    visualization::IgaVisualizationWidget* viewer_ = nullptr;
    std::unique_ptr<visualization::BoundarySamplingCache> samplingCache_;
    std::shared_ptr<visualization::PlaybackHistory> playbackHistory_;
    QThread* solverThread_ = nullptr;
    QTimer* playbackTimer_ = nullptr;
    QElapsedTimer playbackClock_;
    bool playbackReady_ = false;
};
#endif // MAINWINDOW_H
