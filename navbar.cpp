#include "navbar.h"
#include "ui_navbar.h"
#include "mainwindow.h"
#include "secondwindow.h"

Navbar::Navbar(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Navbar)
{
    ui->setupUi(this);

    connect(this, &Navbar::showMainWindow, parent, &MainWindow::show);
    connect(this, &Navbar::showSecondWindow, parent, &SecondWindow::show);
}

void Navbar::on_mainButton_clicked()
{
    emit showMainWindow();
}

void Navbar::on_secondButton_clicked()
{
    emit showSecondWindow();
}

Navbar::~Navbar()
{
    delete ui;
}


void Navbar::on_dichotomy_clicked()
{
    MainWindow* mainWindow = qobject_cast<MainWindow*>(parentWidget());
    if (mainWindow) {
        mainWindow->getStackedWidget()->setCurrentIndex(0);
    }
}

void Navbar::on_newton_clicked() {
    MainWindow* mainWindow = qobject_cast<MainWindow*>(parentWidget());
    if (mainWindow) {
        mainWindow->getStackedWidget()->setCurrentIndex(1);
    }
}
