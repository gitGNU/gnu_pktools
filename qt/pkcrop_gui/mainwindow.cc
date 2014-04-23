/**********************************************************************
mainwindow.cc: GUI for pktools
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
    QStringList interleavedlist;
    interleavedlist << "BAND" << "LINE" << "PIXEL" <<"BSQ";
    ui->interleaved->addItems(interleavedlist);
    QStringList compressedlist;
    compressedlist << "NONE" << "LZW" << "PACKBITS" <<"DEFLATE";
    ui->compressed->addItems(compressedlist);
    QStringList otypelist;
    otypelist << "" << "Byte" << "Int16" << "UInt16" << "UInt32" << "Int32" << "Float32" << "Float64" << "CInt16" << "CInt32" << "CFloat32" << "CFloat64";
    ui->otype->addItems(otypelist);
    QStringList oformatlist;
    oformatlist << "" << "GTiff" << "HFA" << "ENVI";
    ui->oformat->addItems(oformatlist);

    setDefaults();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setDefaults()
{
    m_as=false;
    m_manual=false;
    //input
    ui->listWidget_input->clear();
    ui->ulx->clear();
    ui->uly->clear();
    ui->lrx->clear();
    ui->lry->clear();
    ui->extent->clear();
    //scaling
    ui->resample->setCurrentIndex(0);
    ui->as_from->clear();
    ui->as_to->clear();
    ui->scale->clear();
    ui->offset->clear();
    //output
    ui->output->clear();
    ui->a_srs->clear();
    ui->otype->setCurrentIndex(0);
    ui->oformat->setCurrentIndex(0);
    ui->ct->clear();
    ui->dx->clear();
    ui->dy->clear();
    ui->interleaved->setCurrentIndex(0);
    ui->tiled->setChecked(false);
    ui->compressed->setCurrentIndex(0);
    ui->nodata->clear();
}

void MainWindow::on_toolButton_input_clicked()
{
    on_actionInput_triggered();
}


void MainWindow::on_toolButton_extent_clicked()
{
    on_actionExtent_triggered();
}

void MainWindow::on_toolButton_output_clicked()
{
    on_actionOutput_triggered();
}

void MainWindow::on_toolButton_defaults_clicked()
{
    setDefaults();
}

void MainWindow::on_toolButton_ct_clicked()
{
    QString qsct = QFileDialog::getOpenFileName(this, "Color table ASCII");
    ui->ct->setText(qsct);
}

void MainWindow::on_actionInput_triggered()
{
    QFileDialog dialog(this);
    dialog.setDirectory(QDir::homePath());
    dialog.setFileMode(QFileDialog::ExistingFiles);
    QStringList fileNames;
    if (dialog.exec())
        fileNames = dialog.selectedFiles();
    ui->listWidget_input->addItems(fileNames);
    //fill in band list
    QProcess *myProcess = new QProcess(this);
    QString program="pkinfo -nb -i ";
    //todo: loop over all filenames and get the minimum number of bands?
    program+=fileNames[0];
    myProcess->start(program);
    myProcess->waitForFinished(-1);
    QString p_stdout=myProcess->readAll();
    int nband=p_stdout.section(' ',1).toInt();
    QStringList bandlist;
    for(int iband=0;iband<nband;++iband){
        QString qsband=QString::number(iband);
        bandlist << qsband;
    }
    ui->listWidget_band->addItems(bandlist);
    ui->listWidget_band->setSelectionMode(QAbstractItemView::ExtendedSelection);
    ui->listWidget_band->selectAll();
}

void MainWindow::on_actionExtent_triggered()
{
    QString qsextent = QFileDialog::getOpenFileName(this, "extent");
    ui->extent->setText(qsextent);
}

void MainWindow::on_actionOutput_triggered()
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
            MainWindow::on_actionInput_triggered();
        if(ui->listWidget_input->count()<1){
            QString qsError="No input image file selected";
            throw(qsError);
        }

        for(int i = 0; i < ui->listWidget_input->count(); ++i)
        {
            QListWidgetItem* item = ui->listWidget_input->item(i);
            program+=" --input "+item->text();
        }

        for(int i = 0; i < ui->listWidget_band->count(); ++i)
        {
            QListWidgetItem* item = ui->listWidget_band->item(i);
            program+=" --band "+item->text();
        }

        if(ui->output->text().isEmpty())
            MainWindow::on_actionOutput_triggered();
        if(ui->output->text().isEmpty()){
            QString qsError="No output image file selected";
            throw(qsError);
        }

        program+=" --resample "+ui->resample->currentText();
        if(!ui->otype->currentText().isEmpty())
            program+=" --otype "+ui->otype->currentText();
        if(!ui->oformat->currentText().isEmpty())
            program+=" --oformat "+ui->oformat->currentText();
        program+=" -co COMPRESS="+ui->compressed->currentText();
        program+=" -co INTERLEAVE="+ui->interleaved->currentText();
        if(ui->tiled->isChecked())
            program+=" -co TILED=YES";

        //todo: radiobuttons on scaling
        if(m_as){
            program+=" --autoscale ";
            program+=ui->as_from->text();
            program+=" --autoscale ";
            program+=ui->as_to->text();
        }
        else if(m_manual){
            program+=" --scale ";
            program+=ui->scale->text();
            program+=" ---offset ";
            program+=ui->offset->text();
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
        myProcess->setProcessChannelMode(QProcess::MergedChannels);
        myProcess->waitForFinished(-1);
        QString p_stderr = myProcess->readAllStandardError();
        if(!p_stderr.isEmpty()){
            QMessageBox msgBox;
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

void MainWindow::on_autoscale_clicked()
{
    m_as=true;
    m_manual=false;
}

void MainWindow::on_manual_clicked()
{
    m_as=false;
    m_manual=true;
}

void MainWindow::on_noscale_clicked()
{
    m_as=false;
    m_manual=false;
}
