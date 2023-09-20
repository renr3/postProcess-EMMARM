import modalAnalysis_mainwindow as modalAnalysisWindow
import eModulus_mainwindow as modulusWindow

import sys

#Modules for GUI
from PyQt6 import QtWidgets, uic
from PyQt6.QtCore import Qt, QRunnable, QThreadPool, pyqtSlot, pyqtSignal, QObject
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QApplication, QWidget, QCheckBox, QFileDialog, QDialog, QSplashScreen, QVBoxLayout, QLabel
from PyQt6.QtGui import QDoubleValidator, QIntValidator, QPixmap, QPainter

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

class SplashScreen(QtWidgets.QMainWindow):
    def __init__(self):
        super(SplashScreen,self).__init__()
        uic.loadUi(r'GUI/startup_mainwindow.ui', self)

        original_widget = self.centralwidget
        central_widget = MyCentralWidget(r'GUI/icon/background-400x400.png', original_widget)
        self.setCentralWidget(central_widget)
        
        #Connect signals
        self.toolButton_modalAnalysis
if __name__ == '__main__':
    app = QApplication(sys.argv)
    splash = SplashScreen()
    splash.show()

    sys.exit(app.exec())

'''
app = QtWidgets.QApplication(sys.argv)
main = modulusWindow.eModulus_mainWindow()
main.show()

app2 = QtWidgets.QApplication(sys.argv)
main2 = modalAnalysisWindow.modalAnalysis_mainWindow()
main2.show()
sys.exit(app.exec())
'''