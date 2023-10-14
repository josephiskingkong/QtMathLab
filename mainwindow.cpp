#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QRegularExpression>
#include <QValidator>
#include <QKeyEvent>
#include <cmath>
#include <QWidget>
#include "muParser.h"
#include <limits>

double evaluateMathExpression(const std::string& expression, double x) {
    mu::Parser parser;
    parser.DefineVar("x", &x);

    try {
        parser.SetExpr(expression);
        if (parser.Eval()) {
            return parser.Eval();
        } else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    } catch (mu::ParserError& e) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}


double function(double x) {
    return (27 - 18 * x + 2 * x * x) * exp(-x / 3);
}

double dichotomyMethodMin(double a, double b, double e, const std::string& expression)
{
    while (std::abs(a - b) > e)
    {
        double mid = (a + b) / 2.0;
        double shift = e / 2.0;

        double left = mid - shift;
        double right = mid + shift;

        if (evaluateMathExpression(expression, left) > evaluateMathExpression(expression, right))
        {
            a = mid;
        }
        else
        {
            b = mid;
        }
    }
    return (a + b) / 2.0;
}

double dichotomyMethodMax(double a, double b, double e, const std::string& expression)
{
    while (std::abs(a - b) > e)
    {
        double mid = (a + b) / 2.0;
        double shift = e / 2.0;

        double left = mid - shift;
        double right = mid + shift;

        if (evaluateMathExpression(expression, left) < evaluateMathExpression(expression, right))
        {
            a = mid;
        }
        else
        {
            b = mid;
        }
    }
    return (a + b) / 2.0;
}

double dichotomyMethodZero(double a, double b, double e, const std::string& expression)
{
    double leftValue = evaluateMathExpression(expression, a);
    double rightValue = evaluateMathExpression(expression, b);
    double midValue;

    while (b - a > e)
    {
        double mid = (a + b) / 2.0;
        midValue = evaluateMathExpression(expression, mid);

        if (midValue == 0.0)
        {
            return mid;
        }

        if (leftValue * midValue < 0)
        {
            b = mid;
            rightValue = midValue;
        }
        else
        {
            a = mid;
            leftValue = midValue;
        }
    }

    if (std::abs(leftValue) < 0.0001)
    {
        return a;
    }
    return 500.0;
}

double derivative(const std::string& expression, double x) {
    mu::Parser parser;
    parser.DefineVar("x", &x);
    parser.SetExpr(expression);

    if (parser.Eval()) {
        double f_x = parser.Eval();
        double h = 1e-8;
        double x_plus_h = x + h;

        parser.DefineVar("x", &x_plus_h);
        parser.SetExpr(expression);

        if (parser.Eval()) {
            double f_x_plus_h = parser.Eval();
            return (f_x_plus_h - f_x) / h;
        }
    }

    return std::numeric_limits<double>::quiet_NaN();
}

double newtonMax(double a, double b, double epsilon, const std::string& expression) {
    double x0 = b;
    double x = x0;

    int maxIterations = 1000;
    int iterations = 0;

    while (iterations < maxIterations) {
        double f_x = evaluateMathExpression(expression, x);
        double f_prime_x = derivative(expression, x);

        x = x - f_prime_x > epsilon ? x - f_x / f_prime_x : x - epsilon;

        if (x < a || x > b) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (std::abs(f_prime_x) < epsilon) {
            return x;
        }
        iterations++;
    }

    return x;
}

double newtonMin(double a, double b, double epsilon, const std::string& expression) {
    double x0 = b;
    double x = x0;

    int maxIterations = 1000;
    int iterations = 0;

    while (iterations < maxIterations) {
        double f_x = evaluateMathExpression(expression, x);
        double f_prime_x = derivative(expression, x);

        x = x + f_prime_x > epsilon ? x - f_x / f_prime_x : x + epsilon;

        if (x < a || x > b) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (std::abs(f_prime_x) < epsilon) {
            return x;
        }
        iterations++;
    }

    return x;
}


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QRegularExpression rxPositive("^(1[0-5]|[1-9])$");

    QValidator *validatorPositive = new QRegularExpressionValidator(rxPositive, this);

    ui->dichotomy->setEnabled(false);

//    QRegularExpression rx("^-?([0-9]|[1-8][0-9]|90)$");

//    QValidator *validator = new QRegularExpressionValidator(rx, this);

    ui->lineEdit->setValidator(validatorPositive);
    ui->lineEdit_6->setValidator(validatorPositive);
//    ui->lineEdit_2->setValidator(validator);
//    ui->lineEdit_3->setValidator(validator);
}



MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::keyPressEvent(QKeyEvent* event)
{
    QMainWindow::keyPressEvent(event);
}

void MainWindow::on_pushButton_clicked()
{
    QString resultString = "<b>Результат:<br>";

    double a = ui->lineEdit_2->text().toDouble();
    double b = ui->lineEdit_3->text().toDouble();
    double e = ui->lineEdit->text().toDouble();
    const std::string expression = ui->lineEdit_4->text().toStdString();

    e = pow(10, -e);

    if (ui->lineEdit->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите точность</font>";
        ui->textBrowser_4->setText(resultString);
        return;
    }

    if (ui->lineEdit_2->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите начальную точку</font>";
        ui->textBrowser_4->setText(resultString);
        return;
    }

    if (ui->lineEdit_3->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите конечную точку</font>";
        ui->textBrowser_4->setText(resultString);
        return;
    }

    if (a >= b) {
        resultString += "<font color='red'>Ошибка: Начальная точка равна или правее конечной.</font>";
        ui->textBrowser_4->setText(resultString);
        return;
    }

    if (std::isnan(evaluateMathExpression(expression, (a+b) / 2.0))) {
        resultString += "<font color='red'>Ошибка: Функция не может быть вычислена.</font>";
        ui->textBrowser_4->setText(resultString);
        return;
    }

    double minimumX = dichotomyMethodMin(a, b, e, expression);
    double maximumX = dichotomyMethodMax(a, b, e, expression);
    double zeroX = dichotomyMethodZero(a, b, e, expression);
    if (zeroX == 500)
    {
        resultString += "<font color='red'>Точка нуля не найдена<\b></font>";
    }
    else
    {
        resultString += "Точка нуля фукнции: " + QString::number(zeroX) + "\n<\b>";
    }

    ui->textBrowser_4->setText(resultString);
    makePlot(a, b, evaluateMathExpression(expression, minimumX), evaluateMathExpression(expression, maximumX), expression);
}

void MainWindow::makePlot(double a, double b, double minY, double maxY, const std::string& expression) {
    QCustomPlot* plotContainer = ui->plotContainer;

    plotContainer->xAxis->setRange(a, b);

    plotContainer->yAxis->setRange(minY, maxY);

    QVector<double> xData, yData;
    double step = 0.1;
    for (double x = a; x <= b; x += step) {
        xData.append(x);
        yData.append(evaluateMathExpression(expression, x));
    }

    plotContainer->clearGraphs();

    plotContainer->addGraph();
    plotContainer->graph(0)->setData(xData, yData);
    plotContainer->graph(0)->setPen(QPen(QColor(0, 0, 240), 2));


    plotContainer->replot();
}

void MainWindow::makePlotNewton(double a, double b, double minY, double maxY, const std::string& expression) {
    QCustomPlot* plotContainer = ui->plotContainer_2;

    plotContainer->xAxis->setRange(a, b);

    plotContainer->yAxis->setRange(minY, maxY);

    QVector<double> xData, yData;
    double step = 0.1;
    for (double x = a; x <= b; x += step) {
        xData.append(x);
        yData.append(evaluateMathExpression(expression, x));
    }

    plotContainer->clearGraphs();

    plotContainer->addGraph();
    plotContainer->graph(0)->setData(xData, yData);
    plotContainer->graph(0)->setPen(QPen(QColor(0, 0, 240), 2));


    plotContainer->replot();
}

QStackedWidget* MainWindow::getStackedWidget() {
    return ui->stackedWidget;
}

void MainWindow::on_dichotomy_clicked()
{
    ui->stackedWidget->setCurrentIndex(0);

    ui->dichotomy->setEnabled(false);

    ui->newton->setEnabled(true);
}

void MainWindow::on_newton_clicked()
{
    ui->stackedWidget->setCurrentIndex(1);

    ui->newton->setEnabled(false);

    ui->dichotomy->setEnabled(true);
}

void MainWindow::on_newtonresult_clicked()
{
    QString resultString = "<b>Результат:<br>";

    double a = ui->lineEdit_7->text().toDouble();
    double b = ui->lineEdit_8->text().toDouble();
    double e = ui->lineEdit_6->text().toDouble();
    const std::string expression = ui->lineEdit_5->text().toStdString();

    e = pow(10, -e);

    if (ui->lineEdit_6->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите точность</font>";
        ui->textBrowser_6->setText(resultString);
        return;
    }

    if (ui->lineEdit_7->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите начальную точку</font>";
        ui->textBrowser_6->setText(resultString);
        return;
    }

    if (ui->lineEdit_8->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите конечную точку</font>";
        ui->textBrowser_6->setText(resultString);
        return;
    }

    if (a >= b) {
        resultString += "<font color='red'>Ошибка: Начальная точка равна или правее конечной.</font>";
        ui->textBrowser_6->setText(resultString);
        return;
    }

    if (std::isnan(evaluateMathExpression(expression, (a+b) / 2.0))) {
        resultString += "<font color='red'>Ошибка: Функция не может быть вычислена.</font>";
        ui->textBrowser_6->setText(resultString);
        return;
    }

    double minimumX = newtonMin(a, b, e, expression);
    double maximumX = newtonMax(a, b, e, expression);

    if (!std::isnan(minimumX)) {
        resultString += "Минимум функции: " + QString::number(minimumX) + "\n";
        resultString += "Значение функции в этой точке: " + QString::number(evaluateMathExpression(expression, minimumX)) + "<br>";
    } else {
        resultString += "<font color='red'>Ошибка: Точка минимума не найдена</font><br>";
    }

    if (!std::isnan(maximumX)) {
        resultString += "Максимум функции: " + QString::number(maximumX) + "\n";
        resultString += "Значение функции в этой точке: " + QString::number(evaluateMathExpression(expression, maximumX)) + "<br>";
    } else {
        resultString += "<font color='red'>Ошибка: Точка максимума не найдена</font><br>";
    }

    ui->textBrowser_6->setText(resultString);
    makePlotNewton(a, b, evaluateMathExpression(expression, minimumX), evaluateMathExpression(expression, maximumX), expression);
}



