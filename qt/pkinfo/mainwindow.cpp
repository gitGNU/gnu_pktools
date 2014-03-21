/**********************************************************************
mainwindow.cc
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
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QProcess>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_button_input_clicked()
{ 
    MainWindow::on_menu_input_triggered();
}

void MainWindow::on_pushButton_Run_clicked()
{
    ui->outputEdit->clear();
    ui->commandLineEdit->clear();

    QObject *parent;

    try{
        QString program = "pkinfo";
        if(m_inputFilename.isEmpty())
            MainWindow::on_menu_input_triggered();

        if(m_inputFilename.isEmpty()){
            QString qsError="No input image selected";
            throw(qsError);
        }
        program+=" --input ";
        program+=m_inputFilename;

        QList<QCheckBox*> qcheckBoxList = this->findChildren<QCheckBox *>();

        for(QList<QCheckBox*>::ConstIterator qcbit=qcheckBoxList.begin();qcbit!=qcheckBoxList.end();++qcbit){
            if((*qcbit)->isChecked()){
                QString qsOption;
                qsOption+=" --";
                qsOption+=(*qcbit)->objectName();
                program+=qsOption;
            }
        }

        QList<QComboBox*> qcomboBoxList = this->findChildren<QComboBox *>();

        for(QList<QComboBox*>::ConstIterator qcbit=qcomboBoxList.begin();qcbit!=qcomboBoxList.end();++qcbit){
            QString qsOption;
            qsOption+=" --";
            qsOption+=(*qcbit)->objectName();
            program+=qsOption;
            program+=" ";
            program+=QString::number((*qcbit)->currentIndex());
        }

        QList<QLineEdit*> qlineEditList = this->findChildren<QLineEdit *>();

        for(QList<QLineEdit*>::ConstIterator qlbit=qlineEditList.begin();qlbit!=qlineEditList.end();++qlbit){
            if(!((*qlbit)->text().isEmpty())){
                QString qsOption;
                qsOption+=" --";
                qsOption+=(*qlbit)->objectName();
                qsOption+=" ";
                qsOption+=(*qlbit)->text();
                program+=qsOption;
            }
        }

        ui->commandLineEdit->insertPlainText(program);

//        QProcess *myProcess = new QProcess(parent);
        QProcess *myProcess = new QProcess(this);
        myProcess->start(program);
        myProcess->waitForFinished(-1);
        QString p_stdout = myProcess->readAll();
        ui->outputEdit->appendPlainText(p_stdout);
        delete myProcess;
    }
    catch(QString qsError){
        QMessageBox msgBox;
        msgBox.setText(qsError);
        msgBox.exec();
    }
}

void MainWindow::on_pushButton_clearOutput_clicked()
{
    ui->outputEdit->clear();
}

void MainWindow::on_pushButton_clearCommandLine_clicked()
{
    ui->commandLineEdit->clear();
}

void MainWindow::on_menu_input_triggered()
{
    QObject *parent;
//    QProcess *myProcess = new QProcess(parent);
    QProcess *myProcess = new QProcess(this);
    QString program;
    m_inputFilename = QFileDialog::getOpenFileName(this, "Input image");
    if ( m_inputFilename.isNull() == false ){
        //fill in combobox with number of bands
        program="pkinfo -nb -i ";
        program+=m_inputFilename;
        myProcess->start(program);
        myProcess->waitForFinished(-1);
        QString p_stdout = myProcess->readAll();
        qDebug() << p_stdout;
        int nband=p_stdout.section(' ',1).toInt();
        QStringList list;
        for(int iband=0;iband<nband;++iband){
            QString qsband="band";
            qsband+=QString::number(iband);
            list.append(qsband);
        }
        ui->band->addItems(list);
        myProcess->close();
    }
    delete myProcess;
}

void MainWindow::on_actionQuit_triggered()
{
    close();
}
