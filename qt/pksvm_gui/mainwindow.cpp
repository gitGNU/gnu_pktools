/**********************************************************************
mainwindow.cpp: GUI for pktools
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
#include <QStandardItemModel>
#include <QMessageBox>
#include <QProcess>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QStringList svmlist;
    svmlist << "C_SVC" << "nu_SVC" << "one_class" << "epsilon_SVR" << "nu_SVR";
    ui->svmtype->addItems(svmlist);
    QStringList kernellist;
    kernellist << "radial" << "linear" << "polynomial" << "sigmoid";
    ui->kerneltype->addItems(kernellist);
    setDefaults();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setDefaults()
{
    //tab training
    //m_training="d:\\osgeo\\course\\openstreetmap\\training2.sqlite";
//    ui->training->setText(m_training);
    ui->label->setText("label");
    ui->cv->setText("0");
    //tab input/output
    ui->msknodata->setText("0");
    ui->nodata->setText("0");
    //tab classifier
    ui->coef0->setText("0");
    ui->nu->setText("0.5");
    ui->kd->setText("3");
    ui->label->setText("label");
    ui->cv->setText("0");
    ui->gamma->setText("0");
    ui->ccost->setText("1");
}

void MainWindow::on_actionTraining_triggered()
{
    m_training = QFileDialog::getOpenFileName(this, "Training");
    ui->training->setText(m_training);
    this->on_training_returnPressed();
}

void MainWindow::on_actionMask_triggered()
{
    m_mask = QFileDialog::getOpenFileName(this, "Mask");
    ui->mask->setText(m_mask);
}

void MainWindow::on_actionOutput_triggered()
{
    m_output = QFileDialog::getOpenFileName(this, "Output");
    ui->output->setText(m_output);
}

void MainWindow::on_actionInput_triggered()
{
    m_input = QFileDialog::getOpenFileName(this, "Input");
    ui->input->setText(m_input);
}

void MainWindow::on_toolButton_input_clicked()
{
    on_actionInput_triggered();
}

void MainWindow::on_toolButton_mask_clicked()
{
    on_actionMask_triggered();
}

void MainWindow::on_toolButton_output_clicked()
{
    on_actionOutput_triggered();
}

void MainWindow::on_toolButton_training_clicked()
{
    on_actionTraining_triggered();
}

void MainWindow::on_training_returnPressed()
{
    //eventually read classes from vector file to fill in table...
}

void MainWindow::setClassTable(const QStringList &labels)
{
    QStandardItemModel *model = new QStandardItemModel(labels.size(),3,this); //nlabel rows and 3 columns
    model->setHorizontalHeaderItem(0, new QStandardItem(QString("label name")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("class nr")));
    model->setHorizontalHeaderItem(2, new QStandardItem(QString("prior prob")));
    for(int ilabel=0;ilabel<labels.size();++ilabel){
        QStandardItem *firstCol = new QStandardItem(QString(labels[ilabel]));
        model->setItem(ilabel,0,firstCol);
        QStandardItem *secondCol = new QStandardItem(QString::number(ilabel+1));
        model->setItem(ilabel,1,secondCol);
        QStandardItem *thirdCol = new QStandardItem(QString::number(1.0/labels.size()));
        model->setItem(ilabel,2,thirdCol);
    }
    ui->tableView_labels->setModel(model);
}

void MainWindow::on_pushButton_run_clicked()
{
    try{
        QString program = "pksvm";
        if(m_training.isEmpty())
            MainWindow::on_actionTraining_triggered();

        if(m_training.isEmpty()){
            QString qsError="No training vector file selected";
            throw(qsError);
        }

//        QList<QCheckBox*> qcheckBoxList = this->findChildren<QCheckBox *>();

//        for(QList<QCheckBox*>::ConstIterator qcbit=qcheckBoxList.begin();qcbit!=qcheckBoxList.end();++qcbit){
//            if((*qcbit)->isChecked()){
//                QString qsOption;
//                qsOption+=" --";
//                qsOption+=(*qcbit)->objectName();
//                program+=qsOption;
//            }
//        }

        QList<QComboBox*> qcomboBoxList = this->findChildren<QComboBox *>();

        for(QList<QComboBox*>::ConstIterator qcbit=qcomboBoxList.begin();qcbit!=qcomboBoxList.end();++qcbit){
            QString qsOption;
            qsOption+=" --";
            qsOption+=(*qcbit)->objectName();
            program+=qsOption;
            program+=" ";
            program+=(*qcbit)->currentText();
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

        ui->commandLineEdit->insert(program);

//        QProcess *myProcess = new QProcess(parent);
        QProcess *myProcess = new QProcess(this);
        myProcess->start(program);
        myProcess->waitForFinished(-1);
        QString p_stdout = myProcess->readAll();
        ui->consoleEdit->clear();
        ui->consoleEdit->insertPlainText(p_stdout);
        delete myProcess;
    }
    catch(QString qsError){
        QMessageBox msgBox;
        msgBox.setText(qsError);
        msgBox.exec();
    }
}

void MainWindow::on_toolButton_createTable_clicked()
{
    int nclass=ui->nclass->text().toInt();
    QStringList labels;
    for(int iclass=1;iclass<=nclass;++iclass){
        QString lstring="label";
        lstring+=QString::number(iclass);
        labels << lstring;
    }
    setClassTable(labels);
}

void MainWindow::on_pushButton_restore_clicked()
{
    setDefaults();
}
