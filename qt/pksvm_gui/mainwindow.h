/**********************************************************************
mainwindow.h: GUI for pktools
Copyright (C) 2008-2014 Pieter Kempeneers

This file is part of pktools

pktools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pktools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pktools.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
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

    void on_toolButton_createTable_clicked();

    void on_pushButton_restore_clicked();

    void on_commandLinkButtonPrepareTable_clicked();

    void on_actionColor_table_triggered();

    void on_toolButton_ct_clicked();

private:
    void setClassTable(const QStringList &labels);
    void setDefaults();

    Ui::MainWindow *ui;
//    QString m_input;
//    QString m_training;
//    QString m_mask;
//    QString m_output;
};

#endif // MAINWINDOW_H
