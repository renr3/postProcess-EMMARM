# ------------------------------------------------------
# -------------------- mplwidget.py --------------------
# ------------------------------------------------------
from PyQt6.QtWidgets import*

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from matplotlib.figure import Figure

from PyQt6.QtGui import QFont, QFontInfo

import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('Qt5Agg')
    
class MplWidget(QWidget):
    
    def __init__(self, parent = None):
        super(QWidget, self).__init__()
        #QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvasQTAgg(Figure())
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.canvas.draw()
        
        #Get the default font style from the system and apply to the canvas
        self.resetFontStyle()

        self.setLayout(vertical_layout)

    
    def update_figure(self, new_figure):
        widget_width,widget_height=self.canvas.figure.get_size_inches()
        original_dpi = self.canvas.figure.get_dpi()
        # Clear the existing figure content
        self.canvas.figure.clf()
        self.canvas.axes.cla()

        # Copy the contents of the new figure onto the canvas
        self.canvas.figure = new_figure

        # Set the size of the new figure based on the widget's size
        self.canvas.figure.set_size_inches(widget_width, widget_height)
        self.canvas.figure.set_dpi(original_dpi)
        self.canvas.figure.tight_layout()
        self.resetFontStyle()
        #self.canvas.setMinimumSize(widget_width,widget_height)  # Adjust the scale factor as needed

        # Redraw the canvas
        self.canvas.draw()

    def resetFontStyle(self):
        default_font = QFont()
        default_font_info = QFontInfo(default_font)
        system_font_family = default_font_info.family()
        system_font_size = default_font_info.pointSize()
        plt.rcParams.update({'font.family': system_font_family})
        plt.rcParams.update({'font.size': system_font_size})