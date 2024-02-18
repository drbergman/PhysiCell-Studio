"""
Authors:
Randy Heiland (heiland@iu.edu)
Dr. Paul Macklin (macklinp@iu.edu)
Rf. Credits.md
"""

import sys
import logging
import os
import anndata
import copy
import numpy as np
import time

from pathlib import Path
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QFrame,QApplication,QWidget,QTabWidget,QLineEdit,QHBoxLayout,QVBoxLayout,QRadioButton,QPushButton, QLabel,QCheckBox,QComboBox,QScrollArea,QGridLayout, QFileDialog, QButtonGroup, QToolButton   # , QMessageBox
from PyQt5.QtGui import QIcon

class BioinfImportWindow(QWidget):
    def __init__(self):
        super().__init__()
        # self.layout = QVBoxLayout()
        self.setWindowTitle("Bioinformatics Import Walkthrough")

        # self.setLayout(self.layout)

class QCheckBox_custom(QCheckBox):  # it's insane to have to do this!
    def __init__(self,name):
        super(QCheckBox, self).__init__(name)

        checkbox_style = """
                QCheckBox::indicator:checked {
                    background-color: rgb(255,255,255);
                    border: 1px solid #5A5A5A;
                    width : 15px;
                    height : 15px;
                    border-radius : 3px;
                    image: url(images:checkmark.png);
                }
                QCheckBox::indicator:unchecked
                {
                    background-color: rgb(255,255,255);
                    border: 1px solid #5A5A5A;
                    width : 15px;
                    height : 15px;
                    border-radius : 3px;
                }
                """
        self.setStyleSheet(checkbox_style)

class BioinfImport(QWidget):
    def __init__(self, config_tab, ics_tab):
        super().__init__()
        self.config_tab = config_tab
        self.ics_tab = ics_tab
        self.default_time_units = "min"

        self.qlineedit_style = {}
        self.qlineedit_style[False] = "QLineEdit {background-color: rgb(200,200,200); color: black;}"
        self.qlineedit_style[True] = "QLineEdit {background-color: white; color: black;}"
        # self.combobox_stylesheet = """ 
        #     QComboBox{
        #         color: #000000;
        #         background-color: #FFFFFF; 
        #     }
        #     """

        # self.scroll = QScrollArea()  # might contain centralWidget

        # self.bioinf_import_params = QWidget()

        # self.bioinf_import_tab_layout = QGridLayout()
        # # self.config_tab_layout.addWidget(self.tab_widget, 0,0,1,1) # w, row, column, rowspan, colspan

        # vbox = QVBoxLayout()
        # vbox.addLayout(self.bioinf_import_tab_layout)
        # self.bioinf_import_params.setLayout(vbox)

        # self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        # self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        # self.scroll.setWidgetResizable(True)

        # self.scroll.setWidget(self.bioinf_import_params) # self.config_params = QWidget()

        vbox = QVBoxLayout()
        hbox = QHBoxLayout()
        self.import_button = QPushButton("Import")
        self.import_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        self.import_button.clicked.connect(self.import_cb)
        hbox.addWidget(self.import_button)

        label = QLabel("Column name with cell type (prompt will follow if left blank): ")
        hbox.addWidget(label)

        self.column_line_edit = QLineEdit()
        # self.column_line_edit.setStyleSheet(self.stylesheet)
        # fixed_width_value = 80
        # self.num_cells.setFixedWidth(fixed_width_value)
        # self.num_cells.setValidator(QtGui.QIntValidator(1,100000))
        self.column_line_edit.setEnabled(True)
        self.column_line_edit.setText('leiden')
        hbox.addWidget(self.column_line_edit)

        vbox.addLayout(hbox)

        base_widget = QWidget()
        base_widget.setLayout(vbox)
        self.layout = QVBoxLayout(self)  # leave this!
        self.layout.addWidget(base_widget)

    def import_cb(self):
        file_path = "/Users/danielbergman/seq-to-ic-test/data/pbmc3k_clustered.h5ad"
        self.adata = anndata.read_h5ad(file_path)

        print("------------anndata object loaded-------------")
        print(self.adata)

        self.editing_cell_type_names = True
        col_names = list(self.adata.obs.columns)
        if self.column_line_edit.text() in col_names:
            self.current_column = self.column_line_edit.text()
            self.continue_to_cell_type_names_cb()
            return
        
        vbox = QVBoxLayout()
        label = QLabel("Select column that contains cell type info:")
        vbox.addWidget(label)

        self.column_combobox = QComboBox()
        self.column_combobox.currentIndexChanged.connect(self.column_combobox_changed_cb)
        for col_name in col_names:
            self.column_combobox.addItem(col_name)
        vbox.addWidget(self.column_combobox)

        hbox = QHBoxLayout()
        continue_button = QPushButton("Continue")
        continue_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        continue_button.clicked.connect(self.continue_to_cell_type_names_cb)
        hbox.addWidget(continue_button)

        vbox.addLayout(hbox)

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

    def column_combobox_changed_cb(self, idx):
        self.current_column = self.column_combobox.currentText()
        # print(self.current_column)

    def continue_to_cell_type_names_cb(self):
        self.cell_types_original = self.adata.obs[self.current_column]
        self.cell_types_list_original = self.cell_types_original.unique().tolist()
        self.cell_types_list_original.sort()
        self.remaining_cell_types_list_original = copy.deepcopy(self.cell_types_list_original)

        self.cell_types_list_final = []
        self.check_cell_type_names()

    def check_cell_type_names(self):
        label = self.list_current_cell_type_names() 

        self.cell_type_dict = {}
        vbox = QVBoxLayout()
        vbox.addWidget(label)

        label = QLabel("Accept these as is or edit (merge and delete, and rename)?")
        vbox.addWidget(label)

        hbox = QHBoxLayout()
        accept_button = QPushButton("Accept")
        accept_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        accept_button.clicked.connect(self.accept_as_is_cb)
        hbox.addWidget(accept_button)

        edit_button = QPushButton("Edit")
        edit_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        edit_button.clicked.connect(self.edit_cell_types_cb)
        hbox.addWidget(edit_button)

        vbox.addLayout(hbox)

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

    def continue_to_cell_counts(self):

        names_width = 100
        counts_width = 120
        props_width = 120
        manual_width = 120
        
        vbox = QVBoxLayout()

        hbox = QHBoxLayout()
        label = QLabel("Cell Type")
        label.setFixedWidth(names_width)
        hbox.addWidget(label)

        self.counts_button_group = QButtonGroup()
        self.counts_button_group.idToggled.connect(self.counts_button_cb)

        self.use_counts_as_is_radio_button = QRadioButton("Use counts")
        self.use_counts_as_is_radio_button.setFixedWidth(counts_width)
        self.use_counts_as_is_radio_button.setChecked(True)
        self.counts_button_group.addButton(self.use_counts_as_is_radio_button,0)
        hbox.addWidget(self.use_counts_as_is_radio_button)

        self.use_props_radio_button = QRadioButton("Use proportions")
        self.use_props_radio_button.setFixedWidth(props_width)
        self.use_props_radio_button.setChecked(False)
        self.counts_button_group.addButton(self.use_props_radio_button,1)
        hbox.addWidget(self.use_props_radio_button)

        self.use_manual_radio_button = QRadioButton("Set manually")
        self.use_manual_radio_button.setFixedWidth(manual_width)
        self.use_manual_radio_button.setChecked(False)
        self.counts_button_group.addButton(self.use_manual_radio_button,2)
        hbox.addWidget(self.use_manual_radio_button)

        vbox.addLayout(hbox)

        self.cell_counts = {}
        for cell_type in self.cell_types_list_final:
            self.cell_counts[cell_type] = 0
        for ctn in self.cell_types_final:
            self.cell_counts[ctn] += 1

        self.cell_type_props = [self.cell_counts[cell_type]/len(self.cell_types_final) for cell_type in self.cell_types_list_final]

        self.type_prop = {}
        self.type_manual = {}

        num_validator = QtGui.QIntValidator()
        num_validator.setBottom(0)

        for idx, cell_type in enumerate(self.cell_types_list_final):
            hbox = QHBoxLayout()
            label = QLabel(cell_type)
            label.setFixedWidth(names_width)
            hbox.addWidget(label)

            type_count = QLineEdit()
            type_count.setEnabled(False)
            type_count.setText(str(self.cell_counts[cell_type]))
            type_count.setFixedWidth(counts_width)
            type_count.setStyleSheet(self.qlineedit_style[False])
            hbox.addWidget(type_count)

            self.type_prop[cell_type] = QLineEdit()
            self.type_prop[cell_type].setEnabled(False)
            self.type_prop[cell_type].setText(str(self.cell_counts[cell_type]))
            self.type_prop[cell_type].setFixedWidth(props_width)
            self.type_prop[cell_type].setStyleSheet(self.qlineedit_style[False])
            self.type_prop[cell_type].setValidator(num_validator)
            self.type_prop[cell_type].setObjectName(str(idx))
            self.type_prop[cell_type].textChanged.connect(self.prop_box_changed_cb)
            hbox.addWidget(self.type_prop[cell_type])

            self.type_manual[cell_type] = QLineEdit()
            self.type_manual[cell_type].setEnabled(False)
            self.type_manual[cell_type].setText(str(self.cell_counts[cell_type]))
            self.type_manual[cell_type].setFixedWidth(manual_width)
            self.type_manual[cell_type].setStyleSheet(self.qlineedit_style[False])
            self.type_manual[cell_type].setValidator(num_validator)
            self.type_manual[cell_type].setObjectName(str(idx))
            hbox.addWidget(self.type_manual[cell_type])
            
            vbox.addLayout(hbox)
        
        hbox = QHBoxLayout()
        label = QLabel("Total")
        label.setFixedWidth(names_width)
        hbox.addWidget(label)

        type_count = QLineEdit()
        type_count.setEnabled(False)
        type_count.setText(str(len(self.cell_types_final)))
        type_count.setFixedWidth(counts_width)
        type_count.setStyleSheet(self.qlineedit_style[False])
        hbox.addWidget(type_count)

        self.total_prop = QLineEdit()
        self.total_prop.setEnabled(False)
        self.total_prop.setText(str(len(self.cell_types_final)))
        self.total_prop.setFixedWidth(props_width)
        self.total_prop.setStyleSheet(self.qlineedit_style[False])
        self.total_prop.setValidator(num_validator)
        self.total_prop.textChanged.connect(self.prop_box_changed_cb)
        self.total_prop.setObjectName("total_prop")
        hbox.addWidget(self.total_prop)

        self.total_manual = QLineEdit()
        self.total_manual.setEnabled(False)
        self.total_manual.setText(str(len(self.cell_types_final)))
        self.total_manual.setFixedWidth(manual_width)
        self.total_manual.setStyleSheet(self.qlineedit_style[False])
        self.total_manual.setValidator(num_validator)
        hbox.addWidget(self.total_manual)

        vbox.addLayout(hbox)

        self.continue_to_cell_pos = QPushButton("Continue")
        self.continue_to_cell_pos.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        self.continue_to_cell_pos.clicked.connect(self.continue_to_cell_pos_cb)
        vbox.addWidget(self.continue_to_cell_pos)
        
        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

        self.prop_box_callback_paused = False
    
    def prop_box_changed_cb(self, text):
        # print(f"self.prop_box_callback_paused = {self.prop_box_callback_paused}")
        if self.prop_box_callback_paused:
            return
        type_prop_sender = self.sender()
        if type_prop_sender.hasAcceptableInput() is False:
            return
        # text = type_prop_sender.text()
        # print(f"text = {text}")
        current_idx = type_prop_sender.objectName()
        # print(f"type_prop_sender.objectName() = {current_idx}")
        self.prop_box_callback_paused = True
        if type_prop_sender.objectName()=="total_prop":
            mult = int(text)
        else:
            # print(f"int(text) = {int(text)}")
            # print(f"self.cell_type_props[int(current_idx)] = {self.cell_type_props[int(current_idx)]}")
            mult = int(text) / self.cell_type_props[int(current_idx)]
            self.total_prop.setText(str(round(mult)))

        for idx, cell_type in enumerate(self.cell_types_list_final):
            # print(f"idx = {idx}, cell_type = {cell_type}")
            if current_idx==str(idx):
                continue
            self.type_prop[cell_type].setText(str(round(mult * self.cell_type_props[idx])))
        self.prop_box_callback_paused = False

    def counts_button_cb(self):
        print(f"self.counts_button_group.id = {self.counts_button_group.checkedId()}")
        enable_props = self.counts_button_group.checkedId()==1
        for k in self.type_prop.keys():
            self.type_prop[k].setEnabled(enable_props)
            self.type_prop[k].setStyleSheet(self.qlineedit_style[enable_props])
        self.total_prop.setEnabled(enable_props)
        self.total_prop.setStyleSheet(self.qlineedit_style[enable_props])
        enable_manual = self.counts_button_group.checkedId()==2
        for k in self.type_manual.keys():
            self.type_manual[k].setEnabled(enable_manual)
            self.type_manual[k].setStyleSheet(self.qlineedit_style[enable_manual])
        self.total_manual.setEnabled(enable_manual)
        self.total_manual.setStyleSheet(self.qlineedit_style[enable_manual])

    def continue_to_cell_pos_cb(self):
        if self.counts_button_group.checkedId()==0: # use counts found in data file
            pass
        elif self.counts_button_group.checkedId()==1: # use counts found in data file
            for idx, cell_type in enumerate(self.cell_types_list_final):
                # print(f"idx = {idx}, cell_type = {cell_type}")
                self.cell_counts[idx] = int(self.type_prop[cell_type].text())
        elif self.counts_button_group.checkedId()==2:
            for idx, cell_type in enumerate(self.cell_types_list_final):
                # print(f"idx = {idx}, cell_type = {cell_type}")
                self.cell_counts[idx] = int(self.type_manual[cell_type].text())

        self.cell_pos_button_group = QButtonGroup()
        self.cell_pos_button_group.setExclusive(True)
        
        button_width = 150
        button_height = 150
        icon_width = round(0.8 * button_width)
        icon_height = round(0.8 * button_height)
        master_vbox = QVBoxLayout()

        tool_button_style_sheet = """
            QPushButton {
                background-color : lightblue
            }
            QPushButton::unchecked {
                background-color : lightblue
            }
            QPushButton::checked {
                background-color : black
            }
            """

        self.qpushbutton_style = {}
        self.qpushbutton_style[False] = tool_button_style_sheet
        self.qpushbutton_style[True] = "QPushButton {background-color: black; color: black;}"

        label = QLabel("Select an option for how to start placing cells:")
        label.setAlignment(QtCore.Qt.AlignLeft)

        master_vbox.addWidget(label)

        hbox = QHBoxLayout()
            
        vbox = QVBoxLayout()
        full_rectangle_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/scatter_square.svg"))
        full_rectangle_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        full_rectangle_button.setIconSize(size)
        full_rectangle_button.setCheckable(True)
        full_rectangle_button.setStyleSheet(tool_button_style_sheet) 
        self.cell_pos_button_group.addButton(full_rectangle_button,0)
        vbox.addWidget(full_rectangle_button)

        label = QLabel("Whole Environment")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        hbox.addLayout(vbox)
            
        vbox = QVBoxLayout()
        partial_rectangle_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/rectangle.svg"))
        partial_rectangle_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        partial_rectangle_button.setIconSize(size)
        partial_rectangle_button.setCheckable(True)
        partial_rectangle_button.setStyleSheet(tool_button_style_sheet) 
        self.cell_pos_button_group.addButton(partial_rectangle_button,1)
        vbox.addWidget(partial_rectangle_button)

        label = QLabel("Rectangle")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        hbox.addLayout(vbox)

        vbox = QVBoxLayout()
        disc_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/disc.svg"))
        disc_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        disc_button.setIconSize(size)
        disc_button.setCheckable(True)
        disc_button.setStyleSheet(tool_button_style_sheet) 
        self.cell_pos_button_group.addButton(disc_button,2)
        vbox.addWidget(disc_button)

        label = QLabel("Disc")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        hbox.addLayout(vbox)


        master_vbox.addLayout(hbox)

        self.window = BioinfImportWindow()
        self.window.setLayout(master_vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

        
        print(f"self.cell_counts = {self.cell_counts}")
            
        # self.save_ics()

    def save_ics(self):
        # if len(self.csv_array) == 0:
        #     msg = "No cells created. You must Plot first."
        #     print(msg)
        #     msgBox = QMessageBox()
        #     msgBox.setIcon(QMessageBox.Information)
        #     msgBox.setText(msg)
        #     msgBox.setStandardButtons(QMessageBox.Ok)
        #     returnValue = msgBox.exec()
        #     return

        # print("\n------- ics_tab.py: save_cb() -------")
        # x = y = z = np.arange(0.0,5.0,1.0)
        # np.savetxt('cells.csv', (x,y,z), delimiter=',')
        # print(self.csv_array)
        dir_name = "./config/"
        # print(f'dir_name={dir_name}<end>')
        if len(dir_name) > 0 and not os.path.isdir(dir_name):
            os.makedirs(dir_name)
            time.sleep(1)
        full_fname = os.path.join(dir_name,"cells.csv")
        print("save_cb(): full_fname=",full_fname)
        print("----- Writing v2 (with cell names) .csv file for cells")
        print("----- full_fname=",full_fname)
        
        xmin = float(self.config_tab.xmin.text())
        xmax = float(self.config_tab.xmax.text())
        dx = xmax-xmin
        ymin = float(self.config_tab.ymin.text())
        ymax = float(self.config_tab.ymax.text())
        dy = ymax-ymin
        with open(full_fname, 'w') as f:
            f.write('x,y,z,type\n')  # PhysiCell checks for "x" or "X"
            for ctn in self.cell_types_final:
                x = xmin + dx * np.random.uniform()
                y = ymin + dy * np.random.uniform()
                z = 0
                f.write(f'{x},{y},{z},{ctn}\n')

        # else:
        #     print("----- Writing v1 (with cell indices) .csv file for cells")
        #     print("----- full_fname=",full_fname)
        #     np.savetxt(full_fname, self.csv_array, delimiter=',')

    def accept_as_is_cb(self):
        self.cell_types_final = self.cell_types_original
        self.cell_types_list_final = self.cell_types_list_original
        return self.continue_to_cell_counts()

    def edit_cell_types_cb(self):
        label = self.list_current_cell_type_names()
        
        vbox = QVBoxLayout()
        vbox.addWidget(label)

        label = QLabel("How would you like to proceed? <html><ul><li>Keep: select some types to keep as is</li><li>Delete: select some types to delete</li><li>Merge: select some types to merge</li><li>Accept: accept the remaining cell types as is, i.e. Keep them all</li></ul></html>")
        vbox.addWidget(label)

        hbox = QHBoxLayout()
        keep_button = QPushButton("Keep")
        keep_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        keep_button.clicked.connect(self.keep_cb)
        hbox.addWidget(keep_button)

        delete_button = QPushButton("Delete")
        delete_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        delete_button.clicked.connect(self.delete_cb)
        hbox.addWidget(delete_button)

        merge_button = QPushButton("Merge")
        merge_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        merge_button.clicked.connect(self.merge_cb)
        hbox.addWidget(merge_button)

        accept_button = QPushButton("Accept")
        accept_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        accept_button.clicked.connect(self.accept_cb)
        hbox.addWidget(accept_button)

        vbox.addLayout(hbox)

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

    def keep_cb(self):
        continue_button = self.list_remaining_with_checkboxes("Keep")
        continue_button.clicked.connect(self.continue_on_keep_cb)

    def delete_cb(self):
        continue_button = self.list_remaining_with_checkboxes("Delete")
        continue_button.clicked.connect(self.continue_on_delete_cb)

    def merge_cb(self):
        continue_button = self.list_remaining_with_checkboxes("Merge")
        continue_button.clicked.connect(self.continue_on_merge_cb)

    def accept_cb(self):
        for ctn in self.remaining_cell_types_list_original:
            self.cell_type_dict[ctn] = ctn
        self.continue_to_rename()

    def continue_to_rename(self):
        vbox = QVBoxLayout()
        label = QLabel("Rename your chosen cell types if you like:")
        vbox.addWidget(label)
        self.final_types = []
        self.final_type_pre_image = {}
        all_keys = list(self.cell_type_dict.keys())
        all_keys.sort()
        for ctn in all_keys:
            final_type = self.cell_type_dict[ctn]
            if final_type is not None:
                if final_type not in self.final_types:
                    self.final_types.append(self.cell_type_dict[ctn])
                    self.final_type_pre_image[final_type] = [ctn]
                else:
                    self.final_type_pre_image[final_type].append(ctn)
        labels = {}
        self.new_name_line_edit = {}
        for final_type in self.final_types:
            hbox = QHBoxLayout()
            label_text = ", ".join(self.final_type_pre_image[final_type])
            label_text += " --> "
            labels[final_type] = QLabel(label_text)
            self.new_name_line_edit[final_type] = QLineEdit()
            self.new_name_line_edit[final_type].setText(self.final_type_pre_image[final_type][0])
            hbox.addWidget(labels[final_type])
            hbox.addWidget(self.new_name_line_edit[final_type])
            vbox.addLayout(hbox)

        continue_button = QPushButton("Continue")
        continue_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        continue_button.clicked.connect(self.continue_on_rename_cb)
        vbox.addWidget(continue_button)

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

    def continue_on_rename_cb(self):
        # keep editing here
        for final_type in self.final_types:
            self.cell_types_list_final.append(self.new_name_line_edit[final_type].text())
            for ctn in self.final_type_pre_image[final_type]:
                self.cell_type_dict[ctn] = self.new_name_line_edit[final_type].text()

        self.cell_types_final = [self.cell_type_dict[ctn] for ctn in self.cell_types_original if self.cell_type_dict[ctn] is not None]
        return self.continue_to_cell_counts()

    def list_remaining_with_checkboxes(self, s):
        vbox = QVBoxLayout()
        label = QLabel(f"Select cell types to {s}:")
        vbox.addWidget(label)
        self.checkbox_dict = {}
        for ctn in self.remaining_cell_types_list_original:
            self.checkbox_dict[ctn] = QCheckBox_custom(ctn)
            self.checkbox_dict[ctn].setChecked(False)
            self.checkbox_dict[ctn].setEnabled(True)
            vbox.addWidget(self.checkbox_dict[ctn])

        hbox = QHBoxLayout()
        continue_button = QPushButton("Continue")
        continue_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        hbox.addWidget(continue_button)

        go_back = QPushButton("Go back")
        go_back.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        go_back.clicked.connect(self.go_back_cb)
        hbox.addWidget(go_back)

        vbox.addLayout(hbox)

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()    

        return continue_button

    def continue_on_keep_cb(self):
        for ctn in self.checkbox_dict.keys():
            if self.checkbox_dict[ctn].isChecked():
                self.cell_type_dict[ctn] = ctn
        self.remaining_cell_types_list_original = [ctn for ctn in self.remaining_cell_types_list_original if not self.checkbox_dict[ctn].isChecked()]
        if self.remaining_cell_types_list_original:
            self.edit_cell_types_cb()
        else:
            self.continue_to_rename()

    def continue_on_delete_cb(self):
        for ctn in self.checkbox_dict.keys():
            if self.checkbox_dict[ctn].isChecked():
                self.cell_type_dict[ctn] = None
        self.remaining_cell_types_list_original = [ctn for ctn in self.remaining_cell_types_list_original if not self.checkbox_dict[ctn].isChecked()]
        if self.remaining_cell_types_list_original:
            self.edit_cell_types_cb()
        else:
            self.continue_to_rename()

    def continue_on_merge_cb(self):
        first_name = None
        for ctn in self.checkbox_dict.keys():
            if self.checkbox_dict[ctn].isChecked():
                if first_name is None:
                    first_name = ctn
                self.cell_type_dict[ctn] = first_name
        self.remaining_cell_types_list_original = [ctn for ctn in self.remaining_cell_types_list_original if not self.checkbox_dict[ctn].isChecked()]
        if self.remaining_cell_types_list_original:
            self.edit_cell_types_cb()
        else:
            self.continue_to_rename()

    def go_back_cb(self):
        self.edit_cell_types_cb()
    
    def list_current_cell_type_names(self):
        if self.editing_cell_type_names:
            label_text = "The following cell types remain for editing"
            cell_type_list = self.remaining_cell_types_list_original
        else:
            label_text = "The following cell types were found:"
            cell_type_list = self.cell_types_list_original
        label_text +=  "<html><ul>"
        for ctn in cell_type_list:
            label_text += f"<li>{ctn}</li>"
        label_text += "</ul></html>"
        return QLabel(label_text)
    