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
    QStringList rulelist;
    rulelist << "point" << "centroid" << "mean" << "proportion" << "minimum" << "minimum" << "maximum" << "maximum voting" << "sum";
    ui->rule->addItems(rulelist);
    QStringList formatlist;
    formatlist << "ESRI Shapefile" << "SQLite";
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
    ui->msknodata->setText("0");
    ui->polygon->setChecked(false);
    ui->f->setCurrentIndex(0);
    //tab extract
    ui->bname->setText("B");
    ui->rule->setCurrentIndex(0);

}

void MainWindow::on_actionInput_triggered()
{
    m_input = QFileDialog::getOpenFileName(this, "Input");
    ui->input->setText(m_input);
}

void MainWindow::on_actionSample_triggered()
{
    m_sample = QFileDialog::getOpenFileName(this, "Sample");
    ui->sample->setText(m_sample);
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

void MainWindow::on_toolButton_input_clicked()
{
    on_actionInput_triggered();
}

void MainWindow::on_toolButton_sample_clicked()
{
    on_actionSample_triggered();
}

void MainWindow::on_toolButton_mask_clicked()
{
    on_actionMask_triggered();
}

void MainWindow::on_toolButton_output_clicked()
{
    on_actionOutput_triggered();
}

void MainWindow::setClassTable(const QStringList &labels)
{
    QStandardItemModel *model = new QStandardItemModel(labels.size(),2,this); //2 Rows and 3 Columns
    model->setHorizontalHeaderItem(0, new QStandardItem(QString("label name")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("select(%)")));
    for(int ilabel=0;ilabel<labels.size();++ilabel){
        QStandardItem *firstCol = new QStandardItem(QString(labels[ilabel]));
        model->setItem(ilabel,0,firstCol);
        QStandardItem *secondCol = new QStandardItem(QString::number(100));
        model->setItem(ilabel,1,secondCol);
    }
    ui->tableView_labels->setModel(model);
}

void MainWindow::on_pushButton_run_clicked()
{
    try{
        ui->input->clear();
        ui->consoleEdit->clear();

        QString program = "pkextract";
        if(m_sample.isEmpty())
            MainWindow::on_actionSample_triggered();

        if(m_sample.isEmpty()){
            QString qsError="No sample file selected";
            throw(qsError);
        }

        if(m_input.isEmpty())
            MainWindow::on_actionInput_triggered();

        if(m_input.isEmpty()){
            QString qsError="No input file selected";
            throw(qsError);
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

        ui->commandLineEdit->insert(program);

//        QProcess *myProcess = new QProcess(parent);
        QProcess *myProcess = new QProcess(this);

        myProcess->start(program);
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

