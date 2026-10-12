#ifndef VISUALIZATION_IGAVISUALIZATIONWIDGET_H
#define VISUALIZATION_IGAVISUALIZATIONWIDGET_H

#include "camera/camera.h"
#include "visualization/SurfaceMesh.h"

#include <QOpenGLBuffer>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>

namespace visualization {

/**
 * @brief 使用 OpenGL 3.3 渲染 IGA 形变表面并处理相机、字段和显示交互。
 *
 * 数值求值已经在 SurfaceMesh 中完成；本控件只管理 GPU 资源和显示状态。
 */
class IgaVisualizationWidget : public QOpenGLWidget,
                               protected QOpenGLFunctions_3_3_Core {
    Q_OBJECT

public:
    /** @brief 接管一份 CPU 表面网格并创建可交互 OpenGL 控件。 */
    explicit IgaVisualizationWidget(SurfaceMesh mesh, QWidget* parent = nullptr);
    ~IgaVisualizationWidget() override;

    /** @brief Replaces only per-frame CPU data; static topology remains on the GPU. */
    void setDynamicFrame(SurfaceMesh mesh, double simulationTime);

signals:
    /** @brief 标量字段、值域或形变倍率变化时更新主窗口状态栏。 */
    void statusChanged(const QString& message);

protected:
    // QOpenGLWidget 生命周期和 Qt 输入事件入口。
    void initializeGL() override;
    void resizeGL(int width, int height) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;
    void keyPressEvent(QKeyEvent* event) override;

private:
    /** @brief 将位置、法向或当前标量的变化写回动态 VBO。 */
    void uploadDynamicVertices();
    /** @brief 根据当前形变后包围盒设置相机中心、距离和裁剪面。 */
    void frameCamera();
    /** @brief 使用 QPainter 绘制 colorbar、数值和快捷键说明。 */
    void drawOverlay();
    /** @brief 把当前字段、值域和倍率发送到状态栏。 */
    void emitStatus();

    // CPU 数据与相机状态。
    SurfaceMesh mesh_;
    Camera camera_;
    // 表面着色与纯色线框各使用一套 shader program。
    QOpenGLShaderProgram surfaceProgram_;
    QOpenGLShaderProgram lineProgram_;
    // 形变表面与参考构型使用独立 VAO/VBO，但共享线框索引。
    QOpenGLVertexArrayObject surfaceVao_;
    QOpenGLVertexArrayObject referenceVao_;
    QOpenGLBuffer vertexBuffer_{QOpenGLBuffer::VertexBuffer};
    QOpenGLBuffer referenceBuffer_{QOpenGLBuffer::VertexBuffer};
    QOpenGLBuffer triangleBuffer_{QOpenGLBuffer::IndexBuffer};
    QOpenGLBuffer wireBuffer_{QOpenGLBuffer::IndexBuffer};
    // 运行时显示开关；不影响求解结果。
    bool initialized_ = false;
    bool showWireframe_ = true;
    bool showReference_ = false;
    bool alternatePalette_ = false;
    double simulationTime_ = 0.0;
};

} // namespace visualization

#endif // VISUALIZATION_IGAVISUALIZATIONWIDGET_H
