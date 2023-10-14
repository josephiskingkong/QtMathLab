#ifndef NAVBAR_H
#define NAVBAR_H

#include <QWidget>

namespace Ui {
class Navbar;
}

class Navbar : public QWidget
{
    Q_OBJECT

public:
    Navbar(QWidget *parent = nullptr);
    ~Navbar();

signals:
    void showMainWindow();
    void showSecondWindow();

private slots:
    void on_mainButton_clicked();
    void on_secondButton_clicked();

    void on_dichotomy_clicked();

    void on_newton_clicked();
private:
    Ui::Navbar *ui;
};

#endif // NAVBAR_H
