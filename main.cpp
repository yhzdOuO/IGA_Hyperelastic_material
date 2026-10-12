#include "mainwindow.h"

#include <QApplication>
#include <QCommandLineOption>
#include <QCommandLineParser>
#include <QDebug>

/**
 * @brief Starts the Qt event loop and displays the main window.
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 * @return Qt application exit code.
 */
int main(int argc, char *argv[])
{
    // 程序入口负责选择静力/动力学模式和显示采样密度；默认运行静力学。
    QApplication a(argc, argv);
    QCommandLineParser parser;
    parser.setApplicationDescription(
        "Interactive static/dynamic IGA cantilever visualization");
    parser.addHelpOption();
    QCommandLineOption samplesOption(
        {"s", "samples-per-span"},
        "Number of boundary samples per knot span and direction (>=2).",
        "count", "4");
    parser.addOption(samplesOption);
    QCommandLineOption modeOption(
        {"m", "mode"},
        "Simulation mode: static or dynamic (default: static).",
        "mode", "static");
    parser.addOption(modeOption); //&添加要process时查找的参数
    parser.process(a); //&查找添加的所有option
    bool validSamples = false;
    const int samplesPerSpan = parser.value(samplesOption).toInt(&validSamples);
    if (!validSamples || samplesPerSpan < 2) {
        qCritical() << "--samples-per-span must be an integer >= 2";
        return 2;
    }
    const QString modeName = parser.value(modeOption).trimmed().toLower();
    SimulationMode mode = SimulationMode::Static;
    if (modeName == "dynamic") {
        mode = SimulationMode::Dynamic;
    } else if (modeName != "static") {
        qCritical() << "--mode must be either 'static' or 'dynamic'";
        return 2;
    }
    MainWindow w(samplesPerSpan, mode);
    w.show();
    return a.exec();
}
