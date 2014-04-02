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
#include <QProcess>
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QStringList resamplelist;
    resamplelist << "near" << "bilinear";
    ui->resample->addItems(resamplelist);
    QStringList crulelist;
    crulelist << "overwrite" << "maxndvi" << "maxband" <<"minband" << "mean" << "mode" << "median" << "sum";
    ui->crule->addItems(crulelist);
    QStringList interleavedlist;
    interleavedlist << "BAND" << "LINE" << "PIXEL" <<"BSQ";
    ui->interleaved->addItems(interleavedlist);
    QStringList compressedlist;
    compressedlist << "NONE" << "LZW" << "PACKBITS" <<"DEFLATE";
    ui->compressed->addItems(compressedlist);
    QStringList otypelist;
    otypelist << "Byte" << "Int16" << "UInt16" << "UInt32" << "Int32" << "Float32" << "Float64" << "CInt16" << "CInt32" << "CFloat32" << "CFloat64";
    ui->otype->addItems(otypelist);
    QStringList oformatlist;
    oformatlist << "GTiff" << "HFA" << "ENVI";
    ui->oformat->addItems(oformatlist);
    setDefaults();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setDefaults()
{
    //input
    ui->listWidget_input->clear();
    ui->ulx->clear();
    ui->uly->clear();
    ui->lrx->clear();
    ui->lry->clear();
    //composit
    ui->resample->setCurrentIndex(0);
    ui->crule->setCurrentIndex(0);
    ui->rband->setText("0");
    ui->bndnodata->setText("0");
    ui->srcnodata->setText("0");
    ui->min->clear();
    ui->max->clear();
    //output
    ui->output->clear();
    ui->a_srs->clear();
    ui->otype->setCurrentIndex(0);
    ui->oformat->setCurrentIndex(0);
    ui->ct->clear();
    ui->description->clear();
    ui->dx->clear();
    ui->dy->clear();
    ui->interleaved->setCurrentIndex(0);
    ui->tiled->setChecked(false);
    ui->compressed->setCurrentIndex(0);
    ui->dstnodata->clear();
    ui->file->clear();
}

void MainWindow::on_toolButton_input_clicked()
{
    on_actionInput_image_triggered();
}

void MainWindow::on_toolButton_output_clicked()
{
    on_actionOutput_image_triggered();
}

void MainWindow::on_toolButton_defaults_clicked()
{
    setDefaults();
}

void MainWindow::on_actionInput_image_triggered()
{
    QFileDialog dialog(this);
    dialog.setDirectory(QDir::homePath());
    dialog.setFileMode(QFileDialog::ExistingFiles);
    QStringList fileNames;
    if (dialog.exec())
        fileNames = dialog.selectedFiles();
    ui->listWidget_input->addItems(fileNames);
}

void MainWindow::on_actionOutput_image_triggered()
{
    QString outputfilename=QFileDialog::getSaveFileName(this,"Output image","","*.*");
    ui->output->setText(outputfilename);
}

void MainWindow::on_actionQuit_triggered()
{
    close();
}


void MainWindow::on_toolButton_Run_clicked()
{
    try{
        ui->commandLineEdit->clear();
        ui->consoleEdit->clear();

        QString program = "pkcomposite";

        if(ui->listWidget_input->count()<1)
            MainWindow::on_actionInput_image_triggered();
        if(ui->listWidget_input->count()<1){
            QString qsError="No input image file selected";
            throw(qsError);
        }

        for(int i = 0; i < ui->listWidget_input->count(); ++i)
        {
            QListWidgetItem* item = ui->listWidget_input->item(i);
            program+=" --input "+item->text();
        }

        if(ui->output->text().isEmpty())
            MainWindow::on_actionOutput_image_triggered();
        if(ui->output->text().isEmpty()){
            QString qsError="No output image file selected";
            throw(qsError);
        }

        program+=" --resample "+ui->resample->currentText();
        program+=" --crule "+ui->crule->currentText();
        program+=" --otype "+ui->otype->currentText();
        program+=" --oformat "+ui->oformat->currentText();
        program+=" -co COMPRESS="+ui->compressed->currentText();
        program+=" -co INTERLEAVED="+ui->interleaved->currentText();
        if(ui->tiled->isChecked())
            program+=" -co TILED=YES";

//        QList<QCheckBox*> qcheckBoxList = this->findChildren<QCheckBox *>();

//        for(QList<QCheckBox*>::ConstIterator qcbit=qcheckBoxList.begin();qcbit!=qcheckBoxList.end();++qcbit){
//            if((*qcbit)->isChecked()){
//                QString qsOption;
//                qsOption+=" --";
//                qsOption+=(*qcbit)->objectName();
//                program+=qsOption;
//            }
//        }

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

//        myProcess->start(program);
        myProcess->waitForFinished(-1);
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
