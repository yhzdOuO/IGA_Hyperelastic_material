#include "visualization/IgaVisualizationWidget.h"

#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QSurfaceFormat>
#include <QWheelEvent>
#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>

namespace visualization {
namespace {

// CPU 端 colorbar 使用与 fragment shader 相同语义的蓝-白-红发散色图。
QColor coolWarm(double t)
{
    t = std::clamp(t, 0.0, 1.0);
    if (t < 0.5) {
        const double a = 2.0 * t;
        return QColor::fromRgbF(0.23 + 0.77 * a, 0.30 + 0.70 * a,
                                0.75 + 0.25 * a);
    }
    const double a = 2.0 * (t - 0.5);
    return QColor::fromRgbF(1.0, 1.0 - 0.80 * a, 1.0 - 0.80 * a);
}

// 非负标量（默认 von Mises）使用近似 viridis 的顺序色图。
QColor viridis(double t)
{
    t = std::clamp(t, 0.0, 1.0);
    const double r = std::clamp(0.267 + t * (-0.40 + 1.12 * t), 0.0, 1.0);
    const double g = std::clamp(0.005 + t * (1.45 - 0.50 * t), 0.0, 1.0);
    const double b = std::clamp(0.329 + t * (0.85 - 1.10 * t), 0.0, 1.0);
    return QColor::fromRgbF(r, g, b);
}

} // namespace

IgaVisualizationWidget::IgaVisualizationWidget(SurfaceMesh mesh, QWidget* parent)
    : QOpenGLWidget(parent), mesh_(std::move(mesh))
{
    // 在窗口显示前请求 Core Profile 3.3、深度缓冲和多重采样。
    QSurfaceFormat format;
    format.setVersion(3, 3);
    format.setProfile(QSurfaceFormat::CoreProfile);
    format.setDepthBufferSize(24);
    format.setSamples(4);
    setFormat(format);
    setFocusPolicy(Qt::StrongFocus);
    setMinimumSize(720, 480);
}

IgaVisualizationWidget::~IgaVisualizationWidget()
{
    // OpenGL 对象必须在所属 context 为 current 时销毁。
    if (!initialized_) return;
    makeCurrent();
    surfaceVao_.destroy();
    referenceVao_.destroy();
    vertexBuffer_.destroy();
    referenceBuffer_.destroy();
    triangleBuffer_.destroy();
    wireBuffer_.destroy();
    doneCurrent();
}

void IgaVisualizationWidget::setDynamicFrame(SurfaceMesh mesh, double simulationTime)
{
    if (mesh.vertices.size() != mesh_.vertices.size()
        || mesh.triangleIndices.size() != mesh_.triangleIndices.size()) {
        qWarning("Ignoring dynamic frame with incompatible topology");
        return;
    }
    const ScalarField field = mesh_.activeField;
    const double scale = mesh_.exaggeration;
    mesh_ = std::move(mesh);
    mesh_.selectField(field);
    mesh_.updateGeometry(scale);
    simulationTime_ = simulationTime;
    uploadDynamicVertices();
    emitStatus();
    update();
}

void IgaVisualizationWidget::initializeGL()
{
    // 功能块 1：初始化函数表和基础渲染状态。边界绕序已通过测试，可剔除背面。
    initializeOpenGLFunctions();
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);
    glEnable(GL_MULTISAMPLE);

    // 功能块 2：从 Qt resource 编译/链接表面 shader，失败时立即给出日志。
    if (!surfaceProgram_.addShaderFromSourceFile(QOpenGLShader::Vertex,
                                                  ":/shaders/surface.vert")
        || !surfaceProgram_.addShaderFromSourceFile(QOpenGLShader::Fragment,
                                                     ":/shaders/surface.frag")
        || !surfaceProgram_.link()) {
        qFatal("Surface shader error: %s", qPrintable(surfaceProgram_.log()));
    }
    // 线框和参考构型只需位置与统一颜色，使用更简单的 shader。
    if (!lineProgram_.addShaderFromSourceFile(QOpenGLShader::Vertex,
                                               ":/shaders/wireframe.vert")
        || !lineProgram_.addShaderFromSourceFile(QOpenGLShader::Fragment,
                                                  ":/shaders/wireframe.frag")
        || !lineProgram_.link()) {
        qFatal("Line shader error: %s", qPrintable(lineProgram_.log()));
    }

    // 功能块 3：建立形变表面 VAO、交错动态 VBO 和三角形 EBO。
    surfaceVao_.create();
    surfaceVao_.bind();
    vertexBuffer_.create();
    vertexBuffer_.setUsagePattern(QOpenGLBuffer::DynamicDraw);
    vertexBuffer_.bind();
    std::vector<Vertex> vertices;
    vertices.reserve(mesh_.vertices.size());
    for (const SurfaceVertex& vertex : mesh_.vertices) vertices.push_back(vertex.draw);
    vertexBuffer_.allocate(vertices.data(), static_cast<int>(vertices.size() * sizeof(Vertex)));
    // Vertex 的三个成员按固定 location=0/1/2 映射到 shader 输入。
    surfaceProgram_.enableAttributeArray(0);
    surfaceProgram_.setAttributeBuffer(0, GL_FLOAT, offsetof(Vertex, position), 3,
                                       sizeof(Vertex));
    surfaceProgram_.enableAttributeArray(1);
    surfaceProgram_.setAttributeBuffer(1, GL_FLOAT, offsetof(Vertex, normal), 3,
                                       sizeof(Vertex));
    surfaceProgram_.enableAttributeArray(2);
    surfaceProgram_.setAttributeBuffer(2, GL_FLOAT, offsetof(Vertex, scalar), 1,
                                       sizeof(Vertex));
    // 三角形索引为静态数据；切换字段或倍率时无需重建。
    triangleBuffer_.create();
    triangleBuffer_.bind();
    triangleBuffer_.allocate(mesh_.triangleIndices.data(),
                             static_cast<int>(mesh_.triangleIndices.size()
                                              * sizeof(std::uint32_t)));
    surfaceVao_.release();

    // 功能块 4：单元边界线使用独立 EBO，避免显示三角剖分对角线。
    wireBuffer_.create();
    wireBuffer_.bind();
    wireBuffer_.allocate(mesh_.wireIndices.data(),
                         static_cast<int>(mesh_.wireIndices.size() * sizeof(std::uint32_t)));
    wireBuffer_.release();

    // 功能块 5：参考构型只上传 X，线框索引与当前构型共享。
    referenceVao_.create();
    referenceVao_.bind();
    referenceBuffer_.create();
    referenceBuffer_.setUsagePattern(QOpenGLBuffer::StaticDraw);
    referenceBuffer_.bind();
    std::vector<std::array<float, 3>> referencePositions(mesh_.vertices.size());
    for (std::size_t i = 0; i < mesh_.vertices.size(); ++i) {
        for (int component = 0; component < 3; ++component) {
            referencePositions[i][component] =
                static_cast<float>(mesh_.vertices[i].referencePosition(component));
        }
    }
    referenceBuffer_.allocate(referencePositions.data(),
                              static_cast<int>(referencePositions.size()
                                               * sizeof(referencePositions.front())));
    lineProgram_.enableAttributeArray(0);
    lineProgram_.setAttributeBuffer(0, GL_FLOAT, 0, 3,
                                    sizeof(referencePositions.front()));
    referenceVao_.release();

    // GPU 资源就绪后，按形变后包围盒自动设置初始视角。
    initialized_ = true;
    frameCamera();
    emitStatus();
}

void IgaVisualizationWidget::resizeGL(int width, int height)
{
    // 同步 viewport 宽高比；Camera 内部负责生成新的投影矩阵。
    camera_.setViewportSize(width, height);
}

void IgaVisualizationWidget::paintGL()
{
    // QPainter 绘制叠加层时可能改变 OpenGL 状态，因此每帧都显式恢复实体渲染状态。
    // 深度测试/写入保证近表面遮挡远表面；关闭混合确保表面完全不透明。
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);

    // 功能块 1：清理颜色和深度缓冲，组合当前相机的 MVP 矩阵。
    glClearColor(0.075f, 0.09f, 0.12f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    const QMatrix4x4 view = camera_.getViewMatrix();
    const QMatrix4x4 mvp = camera_.getProjectionMatrix() * view;

    // 功能块 2：绘制按当前标量场着色的形变表面。
    surfaceProgram_.bind();
    surfaceProgram_.setUniformValue("uMvp", mvp);
    surfaceProgram_.setUniformValue("uView", view);
    surfaceProgram_.setUniformValue("uNormalMatrix", view.normalMatrix());
    surfaceProgram_.setUniformValue("uMinimum", static_cast<float>(mesh_.scalarMinimum));
    surfaceProgram_.setUniformValue("uMaximum", static_cast<float>(mesh_.scalarMaximum));
    // 默认按物理量类型选色图；C 键可临时切换另一类色图进行比较。
    const bool diverging = alternatePalette_ ? !mesh_.usesDivergingColorMap()
                                              : mesh_.usesDivergingColorMap();
    surfaceProgram_.setUniformValue("uDiverging", diverging);
    surfaceVao_.bind();
    triangleBuffer_.bind();
    glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(mesh_.triangleIndices.size()),
                   GL_UNSIGNED_INT, nullptr);
    surfaceVao_.release();
    surfaceProgram_.release();

    // 功能块 3：按显示开关叠加当前单元边界与未形变参考边界。
    lineProgram_.bind();
    lineProgram_.setUniformValue("uMvp", mvp);
    wireBuffer_.bind();
    if (showWireframe_) {
        lineProgram_.setUniformValue("uColor", QVector4D(0.02f, 0.02f, 0.02f, 0.75f));
        surfaceVao_.bind();
        wireBuffer_.bind();
        glDrawElements(GL_LINES, static_cast<GLsizei>(mesh_.wireIndices.size()),
                       GL_UNSIGNED_INT, nullptr);
        surfaceVao_.release();
    }
    if (showReference_) {
        lineProgram_.setUniformValue("uColor", QVector4D(0.92f, 0.92f, 0.92f, 0.65f));
        referenceVao_.bind();
        wireBuffer_.bind();
        glDrawElements(GL_LINES, static_cast<GLsizei>(mesh_.wireIndices.size()),
                       GL_UNSIGNED_INT, nullptr);
        referenceVao_.release();
    }
    lineProgram_.release();
    // 功能块 4：二维叠加层不参与三维深度判断或背面剔除。
    // QPainter 会复用当前 OpenGL context，因此必须先切换到明确的 UI 绘制状态；
    // 否则启用的 GL_CULL_FACE 可能剔除 colorbar 的填充图元。
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // 最后用 QPainter 覆盖绘制 colorbar 和文字说明。下帧开头会恢复实体状态。
    drawOverlay();
}

// 鼠标事件统一转交已有 Camera；每次相机状态变化后请求重绘。
void IgaVisualizationWidget::mousePressEvent(QMouseEvent* event)
{
    camera_.handleEvent(event);
    update();
}

void IgaVisualizationWidget::mouseReleaseEvent(QMouseEvent* event)
{
    camera_.handleEvent(event);
    update();
}

void IgaVisualizationWidget::mouseMoveEvent(QMouseEvent* event)
{
    camera_.handleEvent(event);
    update();
}

void IgaVisualizationWidget::wheelEvent(QWheelEvent* event)
{
    camera_.handleEvent(event);
    update();
}

void IgaVisualizationWidget::keyPressEvent(QKeyEvent* event)
{
    // 功能块 1：切换纯显示选项，不触发重新求解或重新采样。
    bool changed = true;
    if (event->key() == Qt::Key_W) {
        showWireframe_ = !showWireframe_;
    } else if (event->key() == Qt::Key_R) {
        showReference_ = !showReference_;
    } else if (event->key() == Qt::Key_C) {
        alternatePalette_ = !alternatePalette_;
    } else if (event->key() == Qt::Key_Plus || event->key() == Qt::Key_Equal) {
        // +/- 只重算 x_draw 和解析法向，SurfaceMesh 中缓存的真实应力保持不变。
        mesh_.updateGeometry(mesh_.exaggeration * 2.0);
        uploadDynamicVertices();
    } else if (event->key() == Qt::Key_Minus) {
        mesh_.updateGeometry(std::max(0.125, mesh_.exaggeration * 0.5));
        uploadDynamicVertices();
    } else if (event->key() >= Qt::Key_1 && event->key() <= Qt::Key_8) {
        // 八个字段均已在采样时缓存，这里只选择数组分量并更新全局值域。
        mesh_.selectField(static_cast<ScalarField>(event->key() - Qt::Key_1));
        uploadDynamicVertices();
    } else if (event->key() == Qt::Key_Escape) {
        window()->close();
    } else {
        changed = false;
        QOpenGLWidget::keyPressEvent(event);
    }
    // 所有已识别操作统一刷新状态栏和画面。
    if (changed) {
        emitStatus();
        update();
    }
}

void IgaVisualizationWidget::uploadDynamicVertices()
{
    if (!initialized_) return;
    // 从完整 CPU 顶点中提取紧凑交错布局，并原位覆盖动态 VBO。
    std::vector<Vertex> vertices;
    vertices.reserve(mesh_.vertices.size());
    for (const SurfaceVertex& vertex : mesh_.vertices) vertices.push_back(vertex.draw);
    vertexBuffer_.bind();
    vertexBuffer_.write(0, vertices.data(), static_cast<int>(vertices.size() * sizeof(Vertex)));
    vertexBuffer_.release();
}

void IgaVisualizationWidget::frameCamera()
{
    // 功能块 1：统计当前显示倍率下形变几何的轴对齐包围盒。
    Eigen::Vector3d minimum = Eigen::Vector3d::Constant(
        std::numeric_limits<double>::infinity());
    Eigen::Vector3d maximum = Eigen::Vector3d::Constant(
        -std::numeric_limits<double>::infinity());
    for (const SurfaceVertex& vertex : mesh_.vertices) {
        const Eigen::Vector3d position = vertex.referencePosition
                                         + mesh_.exaggeration * vertex.displacement;
        minimum = minimum.cwiseMin(position);
        maximum = maximum.cwiseMax(position);
    }
    // 功能块 2：由包围盒中心和对角线确定斜视初始姿态，避免硬编码模型尺度。
    const Eigen::Vector3d center = 0.5 * (minimum + maximum);
    const double diagonal = std::max((maximum - minimum).norm(), 1.0e-3);
    camera_.setPose(0.68f, 0.35f, static_cast<float>(1.65 * diagonal),
                    QVector3D(center(0), center(1), center(2)));
    // 裁剪面也随模型尺度调整，兼顾深度精度并保证整个梁位于视锥内。
    camera_.setNearAndFar(static_cast<float>(std::max(1.0e-3, diagonal * 0.01)),
                          static_cast<float>(6.0 * diagonal));
    camera_.setViewportSize(width(), height());
}

void IgaVisualizationWidget::drawOverlay()
{
    // 功能块 1：逐屏幕行绘制色条，确保方向为顶部最大、底部最小。
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    const int barWidth = 22;
    const int barHeight = std::min(260, height() - 120);
    const QRect bar(width() - 78, 48, barWidth, barHeight);
    const bool diverging = alternatePalette_ ? !mesh_.usesDivergingColorMap()
                                              : mesh_.usesDivergingColorMap();
    for (int y = 0; y < bar.height(); ++y) {
        const double t = 1.0 - static_cast<double>(y) / std::max(1, bar.height() - 1);
        painter.setPen(diverging ? coolWarm(t) : viridis(t));
        painter.drawLine(bar.left(), bar.top() + y, bar.right(), bar.top() + y);
    }
    // 功能块 2：标注字段名、上下限和交互快捷键。
    painter.setPen(QColor(235, 238, 242));
    painter.drawRect(bar.adjusted(-1, -1, 1, 1));
    painter.drawText(width() - 150, 24, scalarFieldName(mesh_.activeField));
    painter.drawText(bar.left() - 46, bar.top() + 5,
                     QString::number(mesh_.scalarMaximum, 'e', 2));
    painter.drawText(bar.left() - 46, bar.bottom(),
                     QString::number(mesh_.scalarMinimum, 'e', 2));
    painter.drawText(14, height() - 16,
                     "LMB orbit | RMB pan | wheel zoom | W wire | R reference | P replay | +/- scale | 1-8 field | C palette");
    painter.drawText(14, 24, QString("physical state  t = %1 s").arg(simulationTime_, 0, 'f', 4));
}

void IgaVisualizationWidget::emitStatus()
{
    // 状态栏信息与 colorbar 使用相同的对称化值域。
    emit statusChanged(QString("t=%1 s | %2 | range [%3, %4] | deformation x%5")
                           .arg(simulationTime_, 0, 'f', 4)
                           .arg(scalarFieldName(mesh_.activeField))
                           .arg(mesh_.scalarMinimum, 0, 'e', 3)
                           .arg(mesh_.scalarMaximum, 0, 'e', 3)
                           .arg(mesh_.exaggeration, 0, 'g', 4));
}

} // namespace visualization
