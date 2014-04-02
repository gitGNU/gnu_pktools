/**********************************************************************
mainwindow.h
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
    void on_toolButton_input_clicked();

    void on_toolButton_defaults_clicked();

    void on_actionInput_image_triggered();

    void on_actionOutput_image_triggered();

    void on_actionQuit_triggered();

    void on_toolButton_output_clicked();

    void on_toolButton_Run_clicked();

private:
    Ui::MainWindow *ui;
    void setDefaults();
};

#endif // MAINWINDOW_H
