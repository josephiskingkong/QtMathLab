#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QKeyEvent>
#include <QWidget>
#include <QStackedWidget>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    QStackedWidget* getStackedWidget();

protected:
    void keyPressEvent(QKeyEvent* event);

private slots:
    void on_pushButton_clicked();

    void makePlot(double a, double b, double minY, double maxY, const std::string& expression);
    void on_newton_clicked();
    void on_dichotomy_clicked();
    void on_newtonresult_clicked();
    void makePlotNewton(double a, double b, double minY, double maxY, const std::string &expression);
private:
    Ui::MainWindow *ui;
    QStackedWidget* stackedWidget;
};

#endif // MAINWINDOW_H
