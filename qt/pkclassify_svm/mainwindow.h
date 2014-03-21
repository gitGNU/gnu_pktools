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

    void on_actionTraining_triggered();

    void on_actionMask_triggered();

    void on_actionOutput_triggered();

    void on_toolButton_input_clicked();

    void on_toolButton_mask_clicked();

    void on_toolButton_output_clicked();

    void on_toolButton_training_clicked();

    void on_training_returnPressed();

    void on_pushButton_run_clicked();

    void on_lineEdit_2_returnPressed();

private:
    void setClassTable(const QStringList &labels);

    Ui::MainWindow *ui;
    QString m_input;
    QString m_training;
    QString m_mask;
    QString m_output;
};

#endif // MAINWINDOW_H
