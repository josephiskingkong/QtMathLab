#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QKeyEvent>
#include <QWidget>
#include <QStackedWidget>
#include <string>
#include <type_traits>


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
    void on_sorting_clicked();

    void on_addNewArrayButton_clicked();

    void on_sortButton_clicked();

    void makePlotArray(std::vector<int> arr);
    void on_clearArray_clicked();

    void on_randomArray_clicked();

    void on_goldRatioResult_clicked();

    void setButtonState(bool dichotomy, bool newton, bool sorting, bool coordinate);
    void on_coordinate_clicked();
    void on_coordinateBtn_clicked();

    void makePlotCoordinate(double a, double b, double minY, double maxY, const std::string &expression);
private:
    Ui::MainWindow *ui;
    QStackedWidget* stackedWidget;
};

#endif // MAINWINDOW_H
