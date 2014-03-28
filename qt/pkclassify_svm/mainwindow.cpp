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
    ui->coef0->setText("0");
    ui->ccost->setText("1");
    ui->gamma->setText("0");
    ui->nu->setText("0.5");
    //test
    m_training="d:\\osgeo\\course\\openstreetmap\\training2.sqlite";
    ui->training->setText(m_training);
}

MainWindow::~MainWindow()
{
    delete ui;
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
    m_training=ui->training->text();
    QStringList labels;
    labels << "forest" << "non-forest";
    setClassTable(labels);
}

void MainWindow::setClassTable(const QStringList &labels)
{
    QStandardItemModel *model = new QStandardItemModel(labels.size(),2,this); //2 Rows and 3 Columns
    model->setHorizontalHeaderItem(0, new QStandardItem(QString("label name")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("class  number")));
    for(int ilabel=0;ilabel<labels.size();++ilabel){
        QStandardItem *firstRow = new QStandardItem(QString(labels[ilabel]));
        model->setItem(ilabel,0,firstRow);
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

void MainWindow::on_lineEdit_2_returnPressed()
{
    //run
}
