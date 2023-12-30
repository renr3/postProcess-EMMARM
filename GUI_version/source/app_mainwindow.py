import modalAnalysis_mainwindow as modalAnalysisWindow
import eModulus_mainwindow as modulusWindow

import sys, os

#Modules for GUI
from PyQt6 import QtWidgets, uic
from PyQt6.QtCore import Qt, QRunnable, QThreadPool, pyqtSlot, pyqtSignal, QObject
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QApplication, QWidget, QCheckBox, QFileDialog, QDialog, QSplashScreen, QVBoxLayout, QLabel
from PyQt6.QtGui import QDoubleValidator, QIntValidator, QPixmap, QPainter

#Variable to make all addresses relative
basedir = os.path.dirname(__file__)

#Make windows recognize this app in the correct group when it is distributed in an installer
try:
    from ctypes import windll  # Only exists on Windows.
    myappid = 'com.EMMARMpostprocess.v.0.1'
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
except ImportError:
    pass

class MyCentralWidget(QWidget):
    def __init__(self, background_image_path, original_widget):
        super().__init__()
        self.background_image_path=background_image_path
        # Create a QVBoxLayout for the central widget
        layout = QVBoxLayout(self)
        layout.addWidget(original_widget)

    def paintEvent(self, event):
        # Create a QPainter and set it to this widget
        painter = QPainter(self)
        # Load the background image using QPixmap
        pixmap = QPixmap(self.background_image_path)  # Replace with your image file path
        # Draw the background image
        painter.drawPixmap(0, 0, self.width(), self.height(), pixmap)

class startupWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(startupWindow,self).__init__()
        uic.loadUi(os.path.join(basedir, 'UI/startup_mainwindow.ui'), self)

        self.analysisToRun = 'exit'

        original_widget = self.centralwidget
        central_widget = MyCentralWidget(os.path.join(basedir, 'UI/icon/background-400x400.png'), original_widget)
        self.setCentralWidget(central_widget)
        
        #Connect signals
        self.toolButton_modalAnalysis.clicked.connect(self.setModalAnalysis)
        self.toolButton_eModulusEstimation.clicked.connect(self.setModulusEstimationAnalysis)

    def setModalAnalysis(self):
        self.analysisToRun = 'modal'
        self.close()
    def setModulusEstimationAnalysis(self):
        self.analysisToRun = 'modulus'
        self.close()

if __name__ == '__main__':
    app_startup = QApplication(sys.argv)
    startup = startupWindow()
    startup.show()
    app_startup.exec()
    analysisToRun = startup.analysisToRun
    exitFlag = True
    #sys.exit(app_startup.exec())
    if startup.analysisToRun == 'exit':
        sys.exit()
    app_modulusEstimationAnalysis = None
    app_modalAnalysis = None
    while exitFlag:
        if analysisToRun == 'modal':
            if app_modalAnalysis is None:
                app_modalAnalysis = QApplication(sys.argv)
            modalAnalysis = modalAnalysisWindow.modalAnalysis_mainWindow()
            modalAnalysis.show()
            app_modalAnalysis.exec()
            if modalAnalysis.openModulusApp:
                analysisToRun = 'modulus'
            else:
                exitFlag=False
        elif analysisToRun == 'modulus':
            if app_modulusEstimationAnalysis is None:
                app_modulusEstimationAnalysis = QApplication(sys.argv)
            modulusEstimationAnalysis = modulusWindow.eModulus_mainWindow()
            modulusEstimationAnalysis.show()
            app_modulusEstimationAnalysis.exec()
            if modulusEstimationAnalysis.openModalApp:
                analysisToRun = 'modal'
            else:
                exitFlag=False
    else:
        sys.exit()