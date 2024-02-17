"""
Authors:
Randy Heiland (heiland@iu.edu)
Dr. Paul Macklin (macklinp@iu.edu)
Rf. Credits.md
"""

import sys
import logging
import os
from pathlib import Path
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QFrame,QApplication,QWidget,QTabWidget,QLineEdit,QHBoxLayout,QVBoxLayout,QRadioButton,QPushButton, QLabel,QCheckBox,QComboBox,QScrollArea,QGridLayout, QFileDialog    # , QMessageBox

# class QCheckBox_custom(QCheckBox):  # it's insane to have to do this!
#     def __init__(self,name):
#         super(QCheckBox, self).__init__(name)

#         checkbox_style = """
#                 QCheckBox::indicator:checked {
#                     background-color: rgb(255,255,255);
#                     border: 1px solid #5A5A5A;
#                     width : 15px;
#                     height : 15px;
#                     border-radius : 3px;
#                     image: url(images:checkmark.png);
#                 }
#                 QCheckBox::indicator:unchecked
#                 {
#                     background-color: rgb(255,255,255);
#                     border: 1px solid #5A5A5A;
#                     width : 15px;
#                     height : 15px;
#                     border-radius : 3px;
#                 }
#                 """
#         self.setStyleSheet(checkbox_style)

class BioinfImport(QWidget):
    def __init__(self):
        super().__init__()
        self.default_time_units = "min"

        qlineedit_style = """
        QLineEdit: disabled {
            background-color:#ff0000;
        }

        QLineEdit: enabled {
            background-color:#ffffff;
        }
        """
        self.setStyleSheet(qlineedit_style)

        self.combobox_stylesheet = """ 
            QComboBox{
                color: #000000;
                background-color: #FFFFFF; 
            }
            """

        self.scroll = QScrollArea()  # might contain centralWidget

        self.bioinf_import_params = QWidget()

        self.bioinf_import_tab_layout = QGridLayout()
        # self.config_tab_layout.addWidget(self.tab_widget, 0,0,1,1) # w, row, column, rowspan, colspan

        vbox = QVBoxLayout()
        vbox.addLayout(self.bioinf_import_tab_layout)
        self.bioinf_import_params.setLayout(vbox)

        self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setWidgetResizable(True)

        self.scroll.setWidget(self.bioinf_import_params) # self.config_params = QWidget()

        self.layout = QVBoxLayout(self)  # leave this!
        self.layout.addWidget(self.scroll)
