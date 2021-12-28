import qdarkstyle
import pyqtgraph as pg
from pyqtgraph.widgets.TableWidget import TableWidget
from Frr import Ui_MainWindow
import sys
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import numpy as np
from qt_material import apply_stylesheet
from PyQt5.QtWebEngineWidgets import *

from codeFormulas.Zha import *
from codeFormulas.Han import *
from codeFormulas.Eurocode import *
from codeFormulas.Japan import *
from codeFormulas.America import *
import xgboost as xgb
import pandas as pd
import numpy as np

model = xgb.Booster(model_file="hybridModel\\FrrModel_1030")

class MainWindow(QMainWindow,Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        self.setWindowFlags(Qt.FramelessWindowHint)
        self.setWindowTitle('Spook')
        #self.setWindowIcon(QIcon('cat.jpg'))
        
    #Zha's
    @pyqtSlot()
    def on_pushButton_zha_clicked(self):
        D=float(self.lineEdit_diameter_zha.text())
        t=float(self.lineEdit_thickness_zha.text())
        L=float(self.lineEdit_length_zha.text())
        fy=float(self.lineEdit_tube_strength_zha.text())
        fc=float(self.lineEdit_concrete_strength_zha.text())
        FRR=float(self.lineEdit_Frr_zha.text())
        is_circular=self.radioButton_Circular_zha.isChecked()
        if self.radioButton_FF_zha.isChecked():
            bc=1
        elif self.radioButton_FP_zha.isChecked():
            bc=0
        else:
            bc=-1

        z=ZhaFrr(D,t,L,fy,fc,FRR,is_circular,bc)
        z.get_NuT()
        load=z.NuT
        self.lineEdit_Load_zha.setText(str(int(load)))

    #Han
    @pyqtSlot()
    def on_pushButton_han_clicked(self):
        D=float(self.lineEdit_diameter_han.text())
        t=float(self.lineEdit_thickness_han.text())
        L=float(self.lineEdit_length_han.text())
        fy=float(self.lineEdit_tube_strength_han.text())
        fc=float(self.lineEdit_concrete_strength_han.text())
        FRR=float(self.lineEdit_Frr_han.text())
        is_circular=self.radioButton_Circular_han.isChecked()
        if self.radioButton_FF_han.isChecked():
            bc=1
        elif self.radioButton_FP_han.isChecked():
            bc=0
        else:
            bc=-1

        h=HanFrr(D,t,L,fy,fc,FRR,is_circular,bc)
        h.get_NfiRd()
        load=z.NfiRd
        self.lineEdit_Load_han.setText(str(int(load)))

    #Eurocode
    @pyqtSlot()
    def on_pushButton_ec_clicked(self):
        D=float(self.lineEdit_diameter_ec.text())/1000
        t=float(self.lineEdit_thickness_ec.text())/1000
        L=float(self.lineEdit_length_ec.text())/1000
        fy=float(self.lineEdit_tube_strength_ec.text())
        fc=float(self.lineEdit_concrete_strength_ec.text())
        FRR=float(self.lineEdit_Frr_ec.text())
        is_circular=self.radioButton_Circular_ec.isChecked()
        
        if self.radioButton_FF_ec.isChecked():
            buckling=0.5
        elif self.radioButton_FP_ec.isChecked():
            buckling=0.7
        else:
            buckling=1
        if is_circular:
            c=CircularCFST(D,t,fc,fy,L,0,0,0,0,0,0,time=FRR,buckling=buckling)
            load=c.EcNuT(10)
        else:
            s=SquareCFST(D,t,fc,fy,L,0,0,0,0,0,0,time=FRR,buckling=buckling)
            load=s.EcNuT(10)
        
        self.lineEdit_Load_ec.setText(str(int(load)))
        
    #Japanese
    @pyqtSlot()
    def on_pushButton_japan_clicked(self):
        
        D=float(self.lineEdit_diameter_japan.text())
        t=float(self.lineEdit_thickness_japan.text())
        D-=2*t
        fc=float(self.lineEdit_concrete_strength_japan.text())
        FRR=float(self.lineEdit_Frr_japan.text())
        is_circular=self.radioButton_Circular_japan.isChecked()

        j=JapanFrr(D,fc,FRR,is_circular)
        j.get_NfiRd()
        load=j.NfiRd
        self.lineEdit_Load_japan.setText(str(int(load)))

    #America
    @pyqtSlot()
    def on_pushButton_am_clicked(self):
        D=float(self.lineEdit_diameter_am.text())
        fc=float(self.lineEdit_concrete_strength_am.text())
        L=float(self.lineEdit_length_am.text())
        load=float(self.lineEdit_Load_am.text())
        is_circular=self.radioButton_Circular_am.isChecked()
        is_siliceous=self.radioButton_siliceous_am.isChecked()
        if self.radioButton_FF_am.isChecked():
            K=0.65
        elif self.radioButton_FP_am.isChecked():
            K=0.8
        else:
            K=1
        FRR=KodurFrr(fc,D,L,load,K,is_circular,is_siliceous)
        self.lineEdit_Frr_am.setText(str(int(FRR)))

    #XGBOOST
    @pyqtSlot()
    def on_pushButton_xgboost_clicked(self):
        D=float(self.lineEdit_diameter_xgboost.text())
        t=float(self.lineEdit_thickness_xgboost.text())
        L=float(self.lineEdit_length_xgboost.text())
        fy=float(self.lineEdit_tube_strength_xgboost.text())
        fc=float(self.lineEdit_concrete_strength_xgboost.text())
        load_ratio=float(self.lineEdit_Load_xgboost.text())
        eccentricity_ratio=float(self.lineEdit_eccentricity_xgboost.text())
        is_circular=self.radioButton_Circular_xgboost.isChecked()
        if is_circular:
            A=D**2*3.14/4
        else:
            A=D**2
        x=pd.DataFrame(np.array([[A,t,L,fy,fc,load_ratio,eccentricity_ratio]]))
        
        x.columns=['Area', 'Tube Thickness', 'Tube Height', 'Steel Strength',
       'Concrete Strength', 'μ', 'e']
        x=xgb.DMatrix(data=x)
        FRR=model.predict(x)
        self.lineEdit_Frr_xgboost.setText(str(int(FRR)))

app =QApplication(sys.argv)
window = MainWindow()

# setup stylesheet
apply_stylesheet(app, theme='dark_teal.xml')
#app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
# or in new API
#app.setStyleSheet(qdarkstyle.load_stylesheet(qt_api='pyqt5'))
window.show()
