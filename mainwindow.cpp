#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QRegularExpression>
#include <QValidator>
#include <QKeyEvent>
#include <cmath>
#include <QWidget>
#include "muParser.h"
#include <limits>
#include <algorithm>
#include <random>
#include <QElapsedTimer>
#include <chrono>

std::vector<int> array;

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

void MainWindow::makePlotArray(std::vector<int> arr) {
    QCustomPlot* plotContainer = ui->sortedPlot;

    plotContainer->clearPlottables();

    QCPBars *bars = new QCPBars(plotContainer->xAxis, plotContainer->yAxis);
    bars->setWidth(1);

    QVector<double> x(arr.size());
    QVector<double> y(arr.size());
    double maxY = 0.0;

    for (int i = 0; i < arr.size(); ++i) {
        x[i] = i;
        y[i] = arr[i];

        if (arr[i] > maxY) {
            maxY = arr[i];
        }
    }

    bars->setData(x, y);

    QPen pen;
    pen.setColor(QColor(255, 0, 0));
    bars->setPen(pen);
    bars->setBrush(QColor(255, 0, 0, 50));
    if (maxY == 0 && arr.size() == 0) {
        plotContainer->xAxis->setRange(0, 5);
        plotContainer->yAxis->setRange(0, 5);

        plotContainer->replot();

        return;
    }

    plotContainer->xAxis->setRange(0, arr.size() - 1);
    plotContainer->yAxis->setRange(0, maxY);

    plotContainer->replot();
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

    if (std::abs(leftValue) < 0.01)
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

double goldenSectionMin(double a, double b, double epsilon, const std::string& expression) {
    const double phi = (1 + std::sqrt(5)) / 2;
    int maxIterations = 1000;
    int iterations = 0;
    double x1 = a + (1 - 1 / phi) * (b - a);
    double x2 = a + 1 / phi * (b - a);
    double f_x1 = evaluateMathExpression(expression, x1);
    double f_x2 = evaluateMathExpression(expression, x2);

    while ((b - a) > epsilon || iterations <= maxIterations) {
        if (f_x1 < f_x2) {
            b = x2;
            x2 = x1;
            x1 = a + (1 - 1 / phi) * (b - a);
            f_x2 = f_x1;
            f_x1 = evaluateMathExpression(expression, x1);
        } else {
            a = x1;
            x1 = x2;
            x2 = a + 1 / phi * (b - a);
            f_x1 = f_x2;
            f_x2 = evaluateMathExpression(expression, x2);
        }
        ++iterations;
    }

    return (a + b) / 2;
}

double goldenSectionMax(double a, double b, double epsilon, const std::string& expression) {
    const double phi = (1 + std::sqrt(5)) / 2;
    int maxIterations = 1000;
    int iterations = 0;
    double x1 = a + (1 - 1 / phi) * (b - a);
    double x2 = a + 1 / phi * (b - a);
    double f_x1 = evaluateMathExpression(expression, x1);
    double f_x2 = evaluateMathExpression(expression, x2);

    while ((b - a) > epsilon || iterations <= maxIterations) {
        if (f_x1 < f_x2) {
            a = x1;
            x1 = x2;
            x2 = a + 1 / phi * (b - a);
            f_x1 = f_x2;
            f_x2 = evaluateMathExpression(expression, x2);
        } else {
            b = x2;
            x2 = x1;
            x1 = a + (1 - 1 / phi) * (b - a);
            f_x2 = f_x1;
            f_x1 = evaluateMathExpression(expression, x1);
        }
        ++iterations;
    }

    return (a + b) / 2;
}

int bubbleSort(std::vector<int>& arr) {
    int iterations = 0;
    int n = arr.size();
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                std::swap(arr[j], arr[j + 1]);
            }
            iterations++;
        }
    }
    return iterations;
}

int insertionSort(std::vector<int>& arr) {
    int n = arr.size();
    int iterations = 0;
    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
        iterations++;
    }
    return iterations;
}

int cocktailSort(std::vector<int>& arr) {
    int iterations = 0;
    int n = arr.size();
    bool swapped = true;
    int start = 0;
    int end = n - 1;
    while (swapped) {
        swapped = false;
        for (int i = start; i < end; i++) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }
        if (!swapped) break;
        swapped = false;
        end--;
        for (int i = end - 1; i >= start; i--) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }
        start++;
        iterations++;
    }
    return iterations;
}

int quickSortIterations = 0;

int partition(std::vector<int>& arr, int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);
    for (int j = low; j < high; j++) {
        if (arr[j] <= pivot) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    return (i + 1);
}

void quickSort(std::vector<int>& arr, int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
        quickSortIterations++;
    }
}

bool isSorted(const std::vector<int>& arr) {
    int n = arr.size();
    for (int i = 1; i < n; i++) {
        if (arr[i] < arr[i - 1]) {
            return false;
        }
    }
    return true;
}

double bogoSort(std::vector<int>& arr, int maxIterations) {
    std::random_device rd;
    std::default_random_engine rng(rd());
    int n = arr.size();
    int iterations = 0;

    while (!isSorted(arr)) {
        if (iterations >= maxIterations) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        std::shuffle(arr.begin(), arr.end(), rng);
        iterations++;
    }

    return static_cast<double>(iterations);
}



double findMinCoordinate(double x, double h, const std::string& expression, double epsilon, int maxIterations) {
    double f_prev = evaluateMathExpression(expression, x);
    double f_current;

    int iterations = 0;
    while (iterations < maxIterations) {
        double x_prev = x;

        x -= h;

        f_current = evaluateMathExpression(expression, x);

        if (f_current > f_prev) {
            x = x_prev;
            if (h > epsilon) {
                h /= 2.0; // Уменьшаем шаг вдвое
            } else {
                break; // Минимальный шаг достигнут
            }
        } else {
            f_prev = f_current;
            h *= 1.1; // Увеличиваем шаг на 10%
        }

        iterations++;
    }

    return x;
}

double findMaxCoordinate(double x, double h, const std::string& expression, double epsilon, int maxIterations) {
    double f_prev = evaluateMathExpression(expression, x);
    double f_current;

    int iterations = 0;
    while (iterations < maxIterations) {
        double x_prev = x;

        x += h;

        f_current = evaluateMathExpression(expression, x);

        if (f_current > f_prev) {
            x = x_prev;
            if (h > epsilon) {
                h /= 2.0;
            } else {
                break;
            }
        } else {
            f_prev = f_current;
            h *= 1.1;
        }

        iterations++;
    }

    return x;
}

std::pair<double, double> findExtrema(double a, double b, double epsilon, const std::string& expression) {
    double initialStep = epsilon * 10.0;

    double x_min = findMinCoordinate(a, initialStep, expression, epsilon, 100);
    double min_result = evaluateMathExpression(expression, x_min);

    double x_max = findMaxCoordinate(a, initialStep, expression, epsilon, 100);
    double max_result = evaluateMathExpression(expression, x_max);

    std::cout << "Initial values: x_min = " << x_min << ", min_result = " << min_result << ", x_max = " << x_max << ", max_result = " << max_result << std::endl;

    for (double x = a; x <= b - epsilon; x += initialStep) {
        double current_value = evaluateMathExpression(expression, x);

        std::cout << "Iteration: x = " << x << ", current_value = " << current_value << std::endl;

        if (current_value < min_result) {
            x_min = findMinCoordinate(x, initialStep, expression, epsilon, 100);
            min_result = evaluateMathExpression(expression, x_min);
            if (initialStep <= epsilon) {
                break;
            }
        }

        if (current_value > max_result) {
            x_max = findMaxCoordinate(x, initialStep, expression, epsilon, 100);
            max_result = evaluateMathExpression(expression, x_max);
            if (initialStep <= epsilon) {
                break;
            }
        }

        initialStep *= 1.1;
    }

    std::cout << "Final values: x_min = " << x_min << ", min_result = " << min_result << ", x_max = " << x_max << ", max_result = " << max_result << std::endl;

    return std::make_pair(x_min, x_max);
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

void MainWindow::makePlotCoordinate(double a, double b, double minY, double maxY, const std::string& expression) {
    QCustomPlot* plotContainer = ui->plotCoordinate;

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

void MainWindow::setButtonState(bool dichotomy, bool newton, bool sorting, bool coordinate)
{
    ui->dichotomy->setEnabled(dichotomy);
    ui->newton->setEnabled(newton);
    ui->sorting->setEnabled(sorting);
    ui->coordinate->setEnabled(coordinate);
}

void MainWindow::on_dichotomy_clicked()
{
    ui->stackedWidget->setCurrentIndex(0);
    setButtonState(false, true, true, true);
}

void MainWindow::on_newton_clicked()
{
    ui->stackedWidget->setCurrentIndex(1);
    setButtonState(true, false, true, true);
}

void MainWindow::on_sorting_clicked()
{
    ui->stackedWidget->setCurrentIndex(2);
    setButtonState(true, true, false, true);
}

void MainWindow::on_coordinate_clicked()
{
    ui->stackedWidget->setCurrentIndex(3);
    setButtonState(true, true, true, false);
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

std::string makeTextFromArray(std::vector<int> array) {
    std::string resultString = "";
    for (int i = 0; i < array.size(); ++i) {
        resultString += std::to_string(array[i]);
        if (i != array.size() - 1) {
            resultString += ", ";
        }
    }
    return resultString;
}

void MainWindow::on_addNewArrayButton_clicked()
{
    const int newArrayNumber = ui->newArrayNumber->value();
    array.push_back(newArrayNumber);
    const std::string arrayString = makeTextFromArray(array);
    ui->notSortedArrayTextBox->setText("Исходный массив: " + QString::fromStdString(arrayString));
    makePlotArray(array);
}

void MainWindow::on_sortButton_clicked()
{
    if (array.size() == 0) {
        QMessageBox::warning(this, "Предупреждение", "Массив пуст");
        return;
    }
    if (!ui->bubbleCheck->isChecked() && !ui->insertionCheck->isChecked() && !ui->cocktailCheck->isChecked() && !ui->fastCheck->isChecked() && !ui->bogoCheck->isChecked()) {
        QMessageBox::warning(this, "Предупреждение", "Вы не выбрали ни один способ сортировки");
        return;
    }

    QString resultTimeText = "";
    QString resultText = "";

    if (ui->bubbleCheck->isChecked()) {
        std::vector<int> arr = array;
        auto start = std::chrono::high_resolution_clock::now();
        int iterations = bubbleSort(arr);
        auto end = std::chrono::high_resolution_clock::now();
        qint64 elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        resultTimeText += "Время выполнения пузырьковой сортировки: " + QString::number(elapsedTime) + " нс\n";
                              resultTimeText += "Количество итераций: " + QString::number(iterations) + "\n";

        resultText += "Пузырьковая сортировка: " + QString::fromStdString(makeTextFromArray(arr)) + "<br>";
    }

    if (ui->insertionCheck->isChecked()) {
        std::vector<int> arr = array;
        auto start = std::chrono::high_resolution_clock::now();
        int iterations = insertionSort(arr);
        auto end = std::chrono::high_resolution_clock::now();
        qint64 elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        resultTimeText += "Время выполнения сортировки вставками: " + QString::number(elapsedTime) + " нс\n";
                              resultTimeText += "Количество итераций: " + QString::number(iterations) + "\n";
        resultText += "Сортировка вставками: " + QString::fromStdString(makeTextFromArray(arr)) + "<br>";
    }

    if (ui->cocktailCheck->isChecked()) {
        std::vector<int> arr = array;
        auto start = std::chrono::high_resolution_clock::now();
        int iterations = cocktailSort(arr);
        auto end = std::chrono::high_resolution_clock::now();
        qint64 elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        resultTimeText += "Время выполнения шейкерной сортировки: " + QString::number(elapsedTime) + " нс\n";
                              resultTimeText += "Количество итераций: " + QString::number(iterations) + "\n";
        resultText += "Шейкерная сортировка: " + QString::fromStdString(makeTextFromArray(arr)) + "<br>";
    }

    if (ui->fastCheck->isChecked()) {
        std::vector<int> arr = array;
        auto start = std::chrono::high_resolution_clock::now();
        quickSort(arr, 0, arr.size() - 1);
        auto end = std::chrono::high_resolution_clock::now();
        qint64 elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        resultTimeText += "Время выполнения быстрой сортировки: " + QString::number(elapsedTime) + " нс\n";
                              resultTimeText += "Количество итераций: " + QString::number(quickSortIterations) + "\n";
        quickSortIterations = 0;
        resultText += "Быстрая сортировка: " + QString::fromStdString(makeTextFromArray(arr)) + "<br>";
    }

    if (ui->bogoCheck->isChecked()) {
        std::vector<int> arr = array;
        auto start = std::chrono::high_resolution_clock::now();
        double iterations = bogoSort(arr, 10000);
        auto end = std::chrono::high_resolution_clock::now();
        qint64 elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        if (!std::isnan(iterations)) {
            resultTimeText += "Время выполнения бого-сортировки: " + QString::number(elapsedTime) + " нс\n";
            resultText += "Бого-сортировка: " + QString::fromStdString(makeTextFromArray(arr)) + "<br>";
resultTimeText += "Количество итераций: " + QString::number(iterations) + "\n";
        } else {
            resultTimeText += "Время выполнения бого-сортировки: Превышено максимальное количество итераций (10000)\n";
            resultText += "Бого-сортировка не выполнена\n";
        }

    }
    bubbleSort(array);
    makePlotArray(array);
    ui->sortedArrayTextBox->setText("Итог сортировки: <br>" + resultText);
    QMessageBox::information(this, "Успешно", resultTimeText);
}



void MainWindow::on_clearArray_clicked()
{
    array.clear();
    const std::string arrayString = makeTextFromArray(array);
    ui->notSortedArrayTextBox->setText("Исходный массив: " + QString::fromStdString(arrayString));
    makePlotArray(array);
}



void MainWindow::on_randomArray_clicked()
{
//    const int arraySize = ui->newArrayNumber->value();
    const int arraySize = 5000;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 99);

    array.clear();

    for (int i = 0; i < arraySize; i++) {
        int randomValue = dist(gen);
        array.push_back(randomValue);
    }

    const std::string arrayString = makeTextFromArray(array);
    ui->notSortedArrayTextBox->setText("Исходный массив: " + QString::fromStdString(arrayString));
    makePlotArray(array);
}

void MainWindow::on_goldRatioResult_clicked()
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

    double minimumX = goldenSectionMin(a, b, e, expression);
    double maximumX = goldenSectionMax(a, b, e, expression);

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




void MainWindow::on_coordinateBtn_clicked()
{
    QString resultString = "<b>Результат:<br>";

    double a = ui->coordinateStart->text().toDouble();
    double b = ui->coordinateEnd->text().toDouble();
    double e = ui->coordinateEpsilon->text().toDouble();
    const std::string expression = ui->coordinateExpression->text().toStdString();

    e = pow(10, -e);

    if (ui->lineEdit_6->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите точность</font>";
        ui->coordinateResult->setText(resultString);
        return;
    }

    if (ui->lineEdit_7->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите начальную точку</font>";
        ui->coordinateResult->setText(resultString);
        return;
    }

    if (ui->lineEdit_8->text() == "") {
        resultString += "<font color='red'>Ошибка: Введите конечную точку</font>";
        ui->coordinateResult->setText(resultString);
        return;
    }

    if (a >= b) {
        resultString += "<font color='red'>Ошибка: Начальная точка равна или правее конечной.</font>";
        ui->coordinateResult->setText(resultString);
        return;
    }

    if (std::isnan(evaluateMathExpression(expression, (a+b) / 2.0))) {
        resultString += "<font color='red'>Ошибка: Функция не может быть вычислена.</font>";
        ui->coordinateResult->setText(resultString);
        return;
    }

    std::pair<double, double> extrema = findExtrema(a, b, e, expression);

    double minimumX = extrema.first;
    double maximumX = extrema.second;

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

    ui->coordinateResult->setText(resultString);
    makePlotCoordinate(a, b, evaluateMathExpression(expression, minimumX), evaluateMathExpression(expression, maximumX), expression);
}

