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
    void on_actionInput_triggered();

    void on_actionSample_triggered();

    void on_actionOutput_triggered();

    void on_toolButton_input_clicked();

    void on_toolButton_output_clicked();

    void on_toolButton_sample_clicked();

    void on_pushButton_run_clicked();

    void on_toolButton_createTable_clicked();

    void on_pushButton_restore_clicked();


private:
    void setClassTable(const QStringList &labels);
    void setDefaults();

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
