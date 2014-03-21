#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:

    void on_pushButton_Run_clicked();

    void on_button_input_clicked();

    void on_pushButton_clearOutput_clicked();

    void on_pushButton_clearCommandLine_clicked();

    void on_menu_input_triggered();

    void on_actionQuit_triggered();

private:
    Ui::MainWindow *ui;
    QString m_inputFilename;
};

#endif // MAINWINDOW_H
