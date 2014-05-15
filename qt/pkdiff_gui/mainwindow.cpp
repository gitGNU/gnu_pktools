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
    QStringList formatlist;
    formatlist << "SQLite" << "ESRI Shapefile";
    ui->f->addItems(formatlist);
    setDefaults();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setDefaults()
{
    //tab input/output
    ui->input->clear();
    ui->reference->clear();
    ui->msknodata->setText("0");
    ui->output->clear();
    ui->confusion->setChecked(false);
}

void MainWindow::on_actionReference_triggered()
{
    QString qsreference= QFileDialog::getOpenFileName(this, "Reference");
    ui->reference->setText(qsreference);
}

void MainWindow::on_actionMask_triggered()
{
    QString qsmask = QFileDialog::getOpenFileName(this, "Mask");
    ui->mask->setText(qsmask);
}

void MainWindow::on_actionOutput_triggered()
{
    QString qsoutput = QFileDialog::getSaveFileName(this,"Output image","","*.*");
    ui->output->setText(qsoutput);
}

void MainWindow::on_actionInput_triggered()
{
    QString qsinput = QFileDialog::getOpenFileName(this, "Input");
    ui->input->setText(qsinput);
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

void MainWindow::on_toolButton_reference_clicked()
{
    on_actionReference_triggered();
}

void MainWindow::setClassTable(const QStringList &labels)
{
    QStandardItemModel *model = new QStandardItemModel(labels.size(),2,this); //nlabel rows and 2 columns
    model->setHorizontalHeaderItem(0, new QStandardItem(QString("label name")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("class nr")));
    for(int ilabel=0;ilabel<labels.size();++ilabel){
        QStandardItem *firstCol = new QStandardItem(QString(labels[ilabel]));
        model->setItem(ilabel,0,firstCol);
        QStandardItem *secondCol = new QStandardItem(QString::number(ilabel+1));
        model->setItem(ilabel,1,secondCol);
    }
    ui->tableView_labels->setModel(model);
}

void MainWindow::on_pushButton_run_clicked()
{
    try{
        ui->commandLineEdit->clear();
        ui->consoleEdit->clear();
        QString program = "pkdiff";
        if(ui->input->text().isEmpty()){
            MainWindow::on_actionInput_triggered();
            if(ui->input->text().isEmpty()){
                QString qsError="No input raster dataset selected";
                throw(qsError);
            }
        }
        if(ui->reference->text().isEmpty()){
            MainWindow::on_actionReference_triggered();
            if(ui->reference->text().isEmpty()){
                QString qsError="No reference vector file selected";
                throw(qsError);
            }
        }
        for(int irow=0;irow<ui->tableView_labels->model()->rowCount();++irow){
            QString qsOption;
            qsOption+=" --class ";
            qsOption+=ui->tableView_labels->model()->data(ui->tableView_labels->model()->index(irow,0)).toString();
            qsOption+=" --reclass ";
            qsOption+=ui->tableView_labels->model()->data(ui->tableView_labels->model()->index(irow,1)).toString();
            program+=qsOption;
         }

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

        if(ui->confusion->isChecked())
            program+=" --confusion";

        ui->commandLineEdit->setText(program);

//        QProcess *myProcess = new QProcess(parent);
        QProcess *myProcess = new QProcess(this);
        myProcess->start(program);
        myProcess->setProcessChannelMode(QProcess::MergedChannels);
        this->setCursor(Qt::WaitCursor);
        myProcess->waitForFinished(-1);
        this->setCursor(Qt::ArrowCursor);
        QMessageBox msgBox;
        QString p_stderr = myProcess->readAllStandardError();
        if(!p_stderr.isEmpty()){
            msgBox.setText(p_stderr);
            msgBox.exec();
        }
        QString p_stdout = myProcess->readAll();
        ui->consoleEdit->insertPlainText(p_stdout);
        delete myProcess;
    }
    catch(QString qsError){
        QMessageBox msgBox;
        msgBox.setText(qsError);
        msgBox.exec();
    }
}

void MainWindow::on_pushButton_restore_clicked()
{
    setDefaults();
}

void MainWindow::on_commandLinkButtonPrepareTable_clicked()
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
