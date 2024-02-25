"""
Authors:
Randy Heiland (heiland@iu.edu)
Dr. Paul Macklin (macklinp@iu.edu)
Rf. Credits.md
"""

import sys
import logging
import os
try:
    import anndata
    HAVE_ANNDATA = True
except:
    HAVE_ANNDATA = False

import copy
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch, Rectangle, Annulus, Wedge
from matplotlib.collections import PatchCollection
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

from pathlib import Path
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QFrame,QApplication,QWidget,QTabWidget,QLineEdit,QHBoxLayout,QVBoxLayout,QRadioButton,QPushButton, QLabel,QCheckBox,QComboBox,QScrollArea,QGridLayout, QFileDialog, QButtonGroup, QToolButton, QSplitter  # , QMessageBox
from PyQt5.QtGui import QIcon

from studio_classes import QHLine, QVLine

class GoBackButton(QPushButton):
    def __init__(self, parent):
        super().__init__()
        self.setText("\u2190 Go back")
        self.setStyleSheet(f"QPushButton {{background-color: lightgreen; color: black;}}")
        self.clicked.connect(parent.go_back_to_prev_window)

class BioinfImportWindow(QWidget):
    def __init__(self):
        super().__init__()
        # self.layout = QVBoxLayout()
        self.setWindowTitle("Bioinformatics Import Walkthrough")

class BioinfImportPlotWindow(QWidget):
    def __init__(self, main_window, config_tab):
        super().__init__()
        self.main_window = main_window
        self.config_tab = config_tab
        self.preview_patch = None
        self.csv_array = {}
        for cell_type in self.main_window.cell_types_list_final:
            self.csv_array[cell_type] = np.empty((0,3))

        self.alpha_value = 1.0
        value = ['gray','red','yellow','green','blue','magenta','orange','lime','cyan','hotpink','peachpuff','darkseagreen','lightskyblue']
        value = value[0:len(self.main_window.cell_types_list_final)]
        self.color_by_celltype = dict(zip(self.main_window.cell_types_list_final, value))

        self.plot_xmin = float(self.config_tab.xmin.text())
        self.plot_xmax = float(self.config_tab.xmax.text())
        self.plot_dx = self.plot_xmax - self.plot_xmin
        self.plot_ymin = float(self.config_tab.ymin.text())
        self.plot_ymax = float(self.config_tab.ymax.text())
        self.plot_dy = self.plot_ymax - self.plot_ymin

        vbox = QVBoxLayout()
       
        hbox = self.create_par_area()
        vbox.addLayout(hbox)

        hbox = QHBoxLayout()
        label = QLabel("""
                       Use your mouse and keyboard to draw your selected region:\
                       <html><ul>\
                       <li>Click: set (x0,y0)</li>\
                       <li>Shift-click: set (w,h), r, or r1</li>\
                       <li>Click-and-drag: set (x0,y0) + {(w,h), r, or r1}</li>\
                       <li>Shift-click-and-drag: set (r0,r1)</li>\
                       <li>Ctrl-click-and-drag: set (r0,r1)</li>\
                       <li>Ctrl-click: set r0</li>\
                       <li>Alt-click-and-drag: set (\u03b81,\u03b82)</li>\
                       <li>Alt+ctrl-click-and-drag: set \u03b81</li>\
                       <li>Alt+shift-click-and-drag: set \u03b82</li>\
                       </ul></html>
                       """
                       )
        hbox.addWidget(label)
        self.create_figure()
        hbox.addWidget(self.canvas)
        
        vbox.addLayout(hbox)

        hbox = QHBoxLayout()
        self.plot_cells_button = QPushButton("Plot", enabled=True)
        self.plot_cells_button.setStyleSheet(self.main_window.qpushbutton_style_sheet)
        self.plot_cells_button.clicked.connect(self.plot_cell_pos)
        hbox.addWidget(self.plot_cells_button)
        
        self.sync_par_area()

        self.cancel_button = QPushButton("Hide")
        self.cancel_button.setStyleSheet("QPushButton {background-color: white; color: black;}")
        self.cancel_button.clicked.connect(self.cancel_cb)
        hbox.addWidget(self.cancel_button)

        vbox.addLayout(hbox)

        self.finish_write_button = QPushButton("Overwrite",enabled=False)
        self.finish_write_button.setStyleSheet(self.main_window.qpushbutton_style_sheet)
        self.finish_write_button.clicked.connect(self.finish_write_button_cb)

        self.finish_append_button = QPushButton("Append",enabled=False)
        self.finish_append_button.setStyleSheet(self.main_window.qpushbutton_style_sheet)
        self.finish_append_button.clicked.connect(self.finish_append_button_cb)

        hbox = QHBoxLayout()
        hbox.addWidget(self.finish_write_button)
        hbox.addWidget(self.finish_append_button)
        vbox.addLayout(hbox)

        self.setLayout(vbox)
        self.hide()
        self.show()

    def create_par_area(self):
        par_label_width = 50
        par_text_width = 75
        self.par_label = []
        self.par_text = []
        hbox = QHBoxLayout()
        for i in range(6):
            self.par_label.append(QLabel())
            self.par_label[i].setAlignment(QtCore.Qt.AlignRight)
            self.par_label[i].setFixedWidth(par_label_width)
            self.par_text.append(QLineEdit())
            self.par_text[i].setFixedWidth(par_text_width)
            self.par_text[i].setStyleSheet(self.main_window.qlineedit_style_sheet)
            hbox.addWidget(self.par_label[i])
            hbox.addWidget(self.par_text[i])

        coord_validator = QtGui.QDoubleValidator()
        self.par_text[0].setValidator(coord_validator)
        self.par_text[1].setValidator(coord_validator)
        self.par_text[4].setValidator(coord_validator) # theta 1
        self.par_text[5].setValidator(coord_validator) # theta 2
        pos_par_validator = QtGui.QDoubleValidator()
        pos_par_validator.setBottom(0)
        for i in range(2,4):
            self.par_text[i].setValidator(pos_par_validator)

        return hbox
    
    def sync_par_area(self):
        if self.preview_patch is not None:
            self.preview_patch.remove()
            self.canvas.update()
            self.canvas.draw()
            self.preview_patch = None

        for cid in self.mpl_cid:
            self.canvas.mpl_disconnect(cid) # might throw error if none exist...
        self.mpl_cid = []

        if self.main_window.cell_pos_button_group.checkedId()==0:
            for pt in self.par_text:
                pt.setEnabled(False)
            for pl in self.par_label:
                pl.setText("")
            self.everywhere_plotter()
        
        else: # set callbacks common to all
            self.mpl_cid.append(self.canvas.mpl_connect("button_release_event", self.mouse_released_cb)) 

            if self.main_window.cell_pos_button_group.checkedId()==1:
                npars = 4
                self.par_label[0].setText("x0")
                self.par_label[1].setText("y0")
                self.par_label[2].setText("width")
                self.par_label[3].setText("height")
                for idx, pt in enumerate(self.par_text):
                    pt.setEnabled(idx < npars)
                    try:
                        pt.textChanged.disconnect()
                    except:
                        pass
                    pt.textChanged.connect(self.rectangle_plotter)
                for i in range(npars,len(self.par_label)):
                    self.par_label[i].setText("")

                print("Switching to rectangle_plotter")
                self.mpl_cid.append(self.canvas.mpl_connect("button_press_event", self.rectangle_mouse_press))
                self.mpl_cid.append(self.canvas.mpl_connect("motion_notify_event", self.rectangle_mouse_motion))
                self.rectangle_plotter()

            elif self.main_window.cell_pos_button_group.checkedId()==2:
                npars = 3
                self.par_label[0].setText("x0")
                self.par_label[1].setText("y0")
                self.par_label[2].setText("r")
                for idx, pt in enumerate(self.par_text):
                    pt.setEnabled(idx < npars)
                    try:
                        pt.textChanged.disconnect()
                    except:
                        pass
                    pt.textChanged.connect(self.disc_plotter)
                for i in range(npars,len(self.par_label)):
                    self.par_label[i].setText("")
                print("Switching to disc_plotter")
                self.mpl_cid.append(self.canvas.mpl_connect("button_press_event", self.disc_mouse_press))
                self.mpl_cid.append(self.canvas.mpl_connect("motion_notify_event", self.disc_mouse_motion))
                self.disc_plotter()

            elif self.main_window.cell_pos_button_group.checkedId()==3:
                npars = 4
                self.par_label[0].setText("x0")
                self.par_label[1].setText("y0")
                self.par_label[2].setText("r0")
                self.par_label[3].setText("r1")
                for idx, pt in enumerate(self.par_text):
                    pt.setEnabled(idx < npars)
                    try:
                        pt.textChanged.disconnect()
                    except:
                        pass
                    pt.textChanged.connect(self.annulus_plotter)
                for i in range(npars,len(self.par_label)):
                    self.par_label[i].setText("")
                print("Switching to annulus_plotter")
                self.mpl_cid.append(self.canvas.mpl_connect("button_press_event", self.annulus_mouse_press))
                self.mpl_cid.append(self.canvas.mpl_connect("motion_notify_event", self.annulus_mouse_motion))
                self.annulus_plotter()

            elif self.main_window.cell_pos_button_group.checkedId()==4:
                npars = 6
                self.par_label[0].setText("x0")
                self.par_label[1].setText("y0")
                self.par_label[2].setText("r0")
                self.par_label[3].setText("r1")
                self.par_label[4].setText("\u03b81 (\u00b0)")
                self.par_label[5].setText("\u03b82 (\u00b0)")
                for idx, pt in enumerate(self.par_text):
                    pt.setEnabled(idx < npars)
                    try:
                        pt.textChanged.disconnect()
                    except:
                        pass
                    pt.textChanged.connect(self.wedge_plotter)
                for i in range(npars,len(self.par_label)):
                    self.par_label[i].setText("")
                print("Switching to wedge_plotter")
                self.mpl_cid.append(self.canvas.mpl_connect("button_press_event", self.wedge_mouse_press))
                self.mpl_cid.append(self.canvas.mpl_connect("motion_notify_event", self.wedge_mouse_motion))
                self.wedge_plotter()

            elif self.main_window.cell_pos_button_group.checkedId()==5:
                npars = 6
                self.par_label[0].setText("x0")
                self.par_label[1].setText("y0")
                self.par_label[2].setText("r0")
                self.par_label[3].setText("r1")
                self.par_label[4].setText("\u03b81 (\u00b0)")
                self.par_label[5].setText("\u03b82 (\u00b0)")
                for idx, pt in enumerate(self.par_text):
                    pt.setEnabled(idx > 1 and idx < npars)
                    try:
                        pt.textChanged.disconnect()
                    except:
                        pass
                    pt.textChanged.connect(self.wedge_plotter)
                for i in range(npars,len(self.par_label)):
                    self.par_label[i].setText("")
                self.wedge_plotter()

        for pt in self.par_text:
            if pt.isEnabled():
                pt.editingFinished.connect(self.move_cursor_to_start)
    
    def move_cursor_to_start(self):
        self.sender().setCursorPosition(0)

    def everywhere_plotter(self):
        if self.preview_patch is None:
            self.preview_patch = self.ax0.add_patch(Rectangle((self.plot_xmin,self.plot_ymin),self.plot_dx,self.plot_dy,alpha=0.2))
        else:
            self.preview_patch.set_bounds(self.plot_xmin,self.plot_ymin,self.plot_dx,self.plot_dy)

        self.plot_cells_button.setEnabled(self.main_window.is_any_cell_type_button_group_checked())

        self.canvas.update()
        self.canvas.draw()

    def get_x0y0(self):
        return float(self.par_text[0].text()), float(self.par_text[1].text())
        
    def set_x0y0(self, event):
        self.x0y0 = (event.xdata,event.ydata)
        self.assign_par(self.x0y0[0],0)
        self.assign_par(self.x0y0[1],1)

    def mouse_released_cb(self, event):
        self.mouse_pressed = False
        self.x0y0 = None

    def time_check(self):
        bval = self.need_to_check_time and (time.time() < self.time_of_press + self.single_click_grace_period)
        if bval is False:
            self.need_to_check_time = False
        return bval
    
    def assign_par(self, value, idx):
        self.par_text[idx].setText(str(value))
        self.par_text[idx].setCursorPosition(0)

    def rectangle_mouse_press(self, event):
        if event.inaxes is None:
            self.x0y0 = None
            self.mouse_pressed = False
            return # then mouse is not over axes, move on
        self.mouse_pressed = True
        self.current_key_is_shift = event.key=="shift"
        self.need_to_check_time = True
        self.time_of_press = time.time()
        if self.current_key_is_shift:
            self.x0y0 = None
            self.rectangle_helper(event, *self.get_x0y0)
        else:
            self.set_x0y0(event)
        if self.shape_started is False:
            self.default_rectangle_pars()

    def rectangle_mouse_motion(self, event):
        if (event.inaxes is None) or (self.mouse_pressed is False) or (self.current_key_is_shift) or self.time_check(): # make sure we are in the axes, that the mouse is pressed, and we're not holding shift
            return
        self.rectangle_helper(event, *self.x0y0)
 
    def rectangle_helper(self, event, xL, yL):
        xR = event.xdata
        yR = event.ydata
        self.assign_par(np.min([xL,xR]),0)
        self.assign_par(np.min([yL,yR]),1)
        self.assign_par(abs(xR-xL),2)
        self.assign_par(abs(yR-yL),3)

    def default_rectangle_pars(self):
        self.default_x0y0()
        self.default_wh()
        self.shape_started = True

    def default_x0y0(self):
        if self.par_text[0].text() == "":
            self.assign_par(0.5*(self.plot_xmin+self.plot_xmax), 0)
        if self.par_text[1].text() == "":
            self.assign_par(0.5*(self.plot_ymin+self.plot_ymax), 1)
            
    def default_wh(self, factor = 0.5):
        x0, y0 = self.get_x0y0()
        if self.par_text[2].text() == "":
            dL = abs(self.plot_xmin-x0)
            dR = abs(self.plot_xmax-x0)
            if dL > dR:
                w = factor * dL
                x0 -= w
                self.assign_par(x0, 0)
            else:
                w = factor * dR
            self.assign_par(w, 2)
        if self.par_text[3].text() == "":
            dL = abs(self.plot_ymin-y0)
            dR = abs(self.plot_ymax-y0)
            if dL > dR:
                h = factor * dL
                y0 -= h
                self.assign_par(y0, 1)
            else:
                h = factor * dR
            self.assign_par(h, 3)

    def disc_mouse_press(self, event):
        if event.inaxes is None:
            self.x0y0 = None
            self.mouse_pressed = False
            return # then mouse is not over axes, move on
        self.mouse_pressed = True
        self.current_key_is_shift = event.key=="shift"
        self.need_to_check_time = True
        self.time_of_press = time.time()
        if self.current_key_is_shift:
            self.x0y0 = None
            self.set_radius_helper(event, *self.get_x0y0(), 2)
        else:
            self.set_x0y0(event)
        if self.shape_started is False:
            self.default_disc_pars()

    def disc_mouse_motion(self, event):
        if (event.inaxes is None) or (self.mouse_pressed is False) or (self.current_key_is_shift) or self.time_check(): # make sure we are in the axes, that the mouse is pressed, and we're not holding shift
            return
        self.set_radius_helper(event, *self.x0y0, 2)
    
    def default_disc_pars(self):
        self.default_x0y0()
        self.default_radius(2)
        self.shape_started = True
       
    def default_radius(self, idx, factor = 0.9):
        if self.par_text[idx].text() == "":
            x0, y0 = self.get_x0y0()
            r = factor * np.min(np.abs([self.plot_xmax-x0,self.plot_xmin-x0,self.plot_ymax-y0,self.plot_ymin-y0]))
            self.assign_par(r, idx)
            
    def annulus_mouse_press(self, event):
        if event.inaxes is None:
            self.x0y0 = None
            self.mouse_pressed = False
            return # then mouse is not over axes, move on
        self.mouse_pressed = True
        self.current_key = event.key
        self.need_to_check_time = True
        self.time_of_press = time.time()
        if self.current_key is None:
            self.set_x0y0(event)
        elif self.current_key == "shift":
            self.x0y0 = self.get_x0y0()
            self.r_start = self.compute_radius(event, *self.x0y0)
            self.set_radius_helper(event, *self.x0y0, 3)
        elif self.current_key == "control":
            self.x0y0 = self.get_x0y0()
            self.r_start = self.compute_radius(event, *self.x0y0)
            self.set_radius_helper(event, *self.x0y0, 2)
        if self.shape_started is False:
            self.default_annulus_pars()

    def annulus_mouse_motion(self, event):
        if (event.inaxes is None) or (self.mouse_pressed is False) or self.time_check() or (self.current_key not in [None,"shift","control"]) : # make sure we are in the axes, that the mouse is pressed, and we're not holding shift or control
            return
        self.radii_motion_helper(event)

    def default_annulus_pars(self):
        self.default_x0y0()
        self.default_radius(2, factor = 0) # set inner radius to half the distance to the nearest boundary
        self.default_radius(3)
        self.shape_started = True

    def compute_radius(self,event, x0, y0):
        x1 = event.xdata
        y1 = event.ydata
        return np.sqrt((x1-x0)**2 + (y1-y0)**2)

    def radii_motion_helper(self, event):
        if self.current_key is None:
            self.set_radius_helper(event, *self.x0y0, 3)
            return
        # if we're here, then we clicked (and continue to click) with the shift or control key pressed for at least a little time; we are setting the inner and outer radii by the start and current mouse position
        self.r_end = self.compute_radius(event, *self.x0y0)
        r0, r1 = [self.r_start,self.r_end] if self.r_start < self.r_end else [self.r_end,self.r_start]
        self.assign_par(r0, 2)
        self.assign_par(r1, 3)

    def set_radius_helper(self, event, x0, y0, idx):
        r = self.compute_radius(event, x0, y0)
        if r is None: # then left axes or something to cause this
            return
        self.assign_par(r, idx)
        return r

    def wedge_mouse_press(self, event):
        if event.inaxes is None:
            self.x0y0 = None
            self.mouse_pressed = False
            return # then mouse is not over axes, move on
        self.mouse_pressed = True
        self.current_key = event.key
        self.need_to_check_time = True
        self.time_of_press = time.time()
        if self.current_key is None:
            self.set_x0y0(event)
        elif self.current_key == "shift":
            self.x0y0 = self.get_x0y0()
            self.r_start = self.compute_radius(event, *self.x0y0)
            self.set_radius_helper(event, *self.x0y0, 3)
        elif self.current_key == "control":
            self.x0y0 = self.get_x0y0()
            self.r_start = self.compute_radius(event, *self.x0y0)
            self.set_radius_helper(event, *self.x0y0, 2)
        elif "alt" in self.current_key: # then work with thetas somehow
            self.x0y0 = self.get_x0y0()
            if "shift" in self.current_key:
                self.assign_par(self.get_angle(event, *self.x0y0), 5) # set theta 2 because shift is above control on my keyboard
            elif "control" in self.current_key or "ctrl" in self.current_key: # I cannot fathom why someone decided that control should show up in two different ways...
                self.assign_par(self.get_angle(event, *self.x0y0), 4) # set theta 2 because shift is above control on my keyboard
            elif self.current_key == "alt":
                theta = self.get_angle(event, *self.x0y0)
                self.assign_par(theta, 4)
                self.assign_par(theta, 5)
        if self.shape_started is False:
            self.default_wedge_pars()

    def wedge_mouse_motion(self, event):
        if (event.inaxes is None) or (self.mouse_pressed is False) or self.time_check() or (self.current_key not in [None,"shift","control","alt","alt+control","ctrl+alt","alt+shift","shift+alt"]) : # make sure we are in the axes, that the mouse is pressed, and we're not holding shift or control or alt
            return
        if self.current_key == "alt" or self.current_key == "alt+shift" or self.current_key == "shift+alt":
            self.assign_par(self.get_angle(event, *self.x0y0), 5)
            return
        if self.current_key == "alt+control" or self.current_key == "ctrl+alt":
            self.assign_par(self.get_angle(event, *self.x0y0), 4)
            return
        self.radii_motion_helper(event)

    def default_wedge_pars(self):
        self.default_annulus_pars()
        if self.par_text[4].text() == "":
            self.assign_par("0",4)
        if self.par_text[5].text() == "":
            self.assign_par("360",5)
        
    def get_angle(self, event, x0, y0):
        x1 = event.xdata
        y1 = event.ydata
        return 57.295779513082323 * np.arctan2(y1-y0,x1-x0) # convert to degrees

    def rectangle_plotter(self):
        pars = []
        for idx, pt in enumerate(self.par_text):
            if idx > 3:
                break
            if pt.hasAcceptableInput() is False:
                self.plot_cells_button.setEnabled(False)
                return # do not update unless all are ready
            pars.append(float(pt.text()))
        # self.plot_cells_button.setEnabled(True)
        x0, y0, width, height = pars
        if self.preview_patch is None:
            self.preview_patch = self.ax0.add_patch(Rectangle((x0,y0),width,height,alpha=0.2))
        else:
            self.preview_patch.set_bounds(x0,y0,width,height)

        # check left edge of rect is left of right edge of domain, right edge of rect is right of left edge of domain (similar in y direction)
        bval = (x0 < self.plot_xmax) and (x0+width > self.plot_xmin) and (y0 < self.plot_ymax) and (y0+height > self.plot_ymin) # make sure the rectangle intersects the domain with positive area

        # self.plot_cells_button.setEnabled(True)
        self.plot_cells_button.setEnabled(bval and self.main_window.is_any_cell_type_button_group_checked())

        self.canvas.update()
        self.canvas.draw()

    def disc_plotter(self):
        pars = []
        for idx, pt in enumerate(self.par_text):
            if idx > 2:
                break
            if pt.hasAcceptableInput() is False:
                self.plot_cells_button.setEnabled(False)
                return # do not update unless all are ready
            pars.append(float(pt.text()))
        # self.plot_cells_button.setEnabled(True)
        x0, y0, r = pars
        if self.preview_patch is None:
            self.preview_patch = self.ax0.add_patch(Circle((x0,y0),r,alpha=0.2))
        else:
            self.preview_patch.set(center=(x0,y0),radius=r)

        # check the disc intersects the domain in non-trivial manner
        r2 = self.get_distance2_to_domain(x0, y0)[0]
        bval = r2 < r*r # make sure the distance from center of Circle to domain is less than radius of circle
        
        self.plot_cells_button.setEnabled(bval and self.main_window.is_any_cell_type_button_group_checked())
        self.canvas.update()
        self.canvas.draw()

    def get_distance2_to_domain(self, x0, y0):
        if x0 < self.plot_xmin:
            dx = x0 - self.plot_xmin # negative
        elif x0 <= self.plot_xmax:
            dx = 0
        else:
            dx = x0 - self.plot_xmax # positive

        if y0 < self.plot_ymin:
            dy = y0 - self.plot_ymin # negative
        elif y0 <= self.plot_ymax:
            dy = 0
        else:
            dy = y0 - self.plot_ymax # positive
        return dx*dx + dy*dy, dx, dy
    
    def annulus_plotter(self):
        pars = []
        for idx, pt in enumerate(self.par_text):
            if idx > 3:
                break
            if pt.hasAcceptableInput() is False:
                self.plot_cells_button.setEnabled(False)
                return # do not update unless all are ready
            pars.append(float(pt.text()))
        x0, y0, r0, r1 = pars
        if r1==0 or (r1 < r0):
            if self.preview_patch: # probably a way to impose this using validators, but that would require dynamically updating the validators...
                self.preview_patch.remove()
                self.canvas.update()
                self.canvas.draw()
                self.preview_patch = None
            self.plot_cells_button.setEnabled(False)
            return
        
        r2 = self.get_distance2_to_domain(x0,y0)[0]
        cr2 = self.get_circumscribing_radius(x0, y0)
        # outer_radius_reaches_domain = r2 < r1*r1
        # inner_radius_does_not_contain_entire_domain = cr2 > r0*r0
        bval = (r2 < r1*r1) and (cr2 > r0*r0)

        self.plot_cells_button.setEnabled(bval and self.main_window.is_any_cell_type_button_group_checked())

        # print(f"r0 = {r0}\t r1 = {r1}\twidth = {r1-r0}\twidth <= r1 = {r1-r0 <= r1}")
        width = r1-r0
        try:
            self.annulus_setter(x0,y0,r1,width)
        except:
            # my PR to matplotlib should resolve the need for this check!
            # if width==r1: # hack to address the bug in matplotlib.patches.Annulus which checks if width < r0 rather than <= r0 as the error message suggests it does
            #     width *= 1 - np.finfo(width).eps # reduce width by the littlest bit possible to make sure this bug doesn't hit
            print("\tBioinfImport WARNING: You likely can use an update to matplotlib to fix a bug in their Annulus plots.\n\tWe'll take care of it for now.")
            self.annulus_setter(x0,y0,r1,width*(1-np.finfo(width).eps))

        self.canvas.update()
        self.canvas.draw()

    def annulus_setter(self, x0, y0, r1, width):
        if self.preview_patch is None:
            self.preview_patch = self.ax0.add_patch(Annulus((x0,y0),r1,width,alpha=0.2))
        else:
            self.preview_patch.set(center=(x0,y0),radii=r1,width=width)

    def get_circumscribing_radius(self,x0,y0):
        if 2*x0 < self.plot_xmin + self.plot_xmax: # if left of midpoint
            dx = self.plot_xmax - x0
        else:
            dx = x0 - self.plot_xmin
        if 2*y0 < self.plot_ymin + self.plot_ymax: # if left of midpoint
            dy = self.plot_ymax - y0
        else:
            dy = y0 - self.plot_ymin
        return dx*dx + dy*dy

    def wedge_plotter(self):
        pars = []
        for idx, pt in enumerate(self.par_text):
            if idx >= 6:
                break
            if pt.hasAcceptableInput() is False:
                self.plot_cells_button.setEnabled(False)
                return # do not update unless all are ready
            pars.append(float(pt.text()))
        x0, y0, r0, r1, th1, th2 = pars
        if r1 < r0:
            # print(f"r1 = {r1}, r0 = {r0}")
            if self.preview_patch: # probably a way to impose this using validators, but that would require dynamically updating the validators...
                self.preview_patch.remove()
                self.canvas.update()
                self.canvas.draw()
                self.preview_patch = None
            self.plot_cells_button.setEnabled(False)
            return
        
        r2, dx, dy = self.get_distance2_to_domain(x0,y0)
        cr2 = self.get_circumscribing_radius(x0, y0)
        # outer_radius_reaches_domain = r2 < r1*r1
        # inner_radius_does_not_contain_entire_domain = cr2 > r0*r0
        bval = (r2 < r1*r1) and (cr2 > r0*r0)

        bval = bval and self.wedge_in_domain(x0,y0,r0,r1,th1,th2,dx,dy,r2)

        self.plot_cells_button.setEnabled(bval and self.main_window.is_any_cell_type_button_group_checked())

        if self.preview_patch is None:
            self.preview_patch = self.ax0.add_patch(Wedge((x0,y0),r1,th1,th2,width=r1-r0,alpha=0.2))
        else:
            self.preview_patch.set(center=(x0,y0),radius=r1,theta1=th1,theta2=th2,width=r1-r0)
        
        self.canvas.update()
        self.canvas.draw()

    def wedge_in_domain(self,x0,y0,r0,r1,th1,th2,dx,dy,r2):
        th1, th2 = normalize_thetas(th1,th2)
        if r2==0: # then (x0,y0) is in domain
            # first find shortest distances to edge of domain
            r_th1 = self.distance_to_domain_from_within(x0,y0,th1)
            if r_th1 > r0:
                return True
            r_th2 = self.distance_to_domain_from_within(x0,y0,th2)
            if r_th2 > r0:
                return True

            # If the above don't work, then hopefully checking these easy distances may work
            starting_theta_step = 1 + (th1 // 90) # first 90 deg angle that could go to side of domain for shortest distance purposes
            end_theta_step = 1 + (th2 // 90) # first 90 deg angle that could go to side of domain for shortest distance purposes
            mid_thetas_step = np.arange(starting_theta_step,end_theta_step)
            for th in mid_thetas_step:
                if th % 4 == 0: # right
                    d = self.plot_xmax - x0
                elif th % 4 == 1: # up
                    d = self.plot_ymax - y0
                elif th % 4 == 2: # left
                    d = x0 - self.plot_xmin
                elif th % 4 == 3: # down
                   d = y0 - self.plot_ymin
                if d > r0:
                    return True
            
            # If those easy-to-calculate distances fail, then we try the corners, which are our last hope
            r02 = r0*r0
            th1_rad = th1*0.017453292519943
            th2_rad = th2*0.017453292519943
            th = np.arctan2(self.plot_ymax-y0,self.plot_xmax-x0)
            if th > th1_rad and th < th2_rad:
                if  (self.plot_xmax - x0)**2 + (self.plot_ymax-y0)**2 > r02:
                    return True
            th = np.arctan2(self.plot_ymax-y0,self.plot_xmin-x0)
            if th > th1_rad and th < th2_rad:
                if  (self.plot_xmin - x0)**2 + (self.plot_ymax-y0)**2 > r02:
                    return True
            th = np.arctan2(self.plot_ymin-y0,self.plot_xmin-x0)
            if th > th1_rad and th < th2_rad:
                if  (self.plot_xmin - x0)**2 + (self.plot_ymin-y0)**2 > r02:
                    return True
            th = np.arctan2(self.plot_ymin-y0,self.plot_xmax-x0)
            if th > th1_rad and th < th2_rad:
                if  (self.plot_xmax - x0)**2 + (self.plot_ymin-y0)**2 > r02:
                    return True
            return False
        else: # then (x0,y0) is not in the domain
            return self.wedge_in_domain_center_out(x0,y0,r0,r1,th1,th2,dx,dy)

    def wedge_in_domain_center_out(self,x0,y0,r0,r1,th1,th2,dx,dy):
        # check to see if the wedge is at all in the domain when the center is outside the domain
        xL, xR, yL, yR = [self.plot_xmin, self.plot_xmax, self.plot_ymin, self.plot_ymax]
        # WLOG set center on left or bottom-left of domain
        if dx > 0: 
            if dy == 0: # reflect so it is on left
                xL, xR = [-xR,-xL]
                x0 *= -1
                th1, th2 = [180-th2,180-th1] # reflections flip orientation
                dx *= -1
            elif dy > 0: # rotate 180
                xL, xR, yL, yR = [-xR, -xL, -yR, -yL]
                x0 *= -1
                y0 *= -1
                th1 += 180
                th2 += 180
                dx *= -1
                dy *= -1
            else: # dy < 0 rotate 270
                xL, xR, yL, yR = [yL, yR, -xR, -xL]
                x0, y0 = [y0, -x0]
                th1 += 270
                th2 += 270
                dx, dy = [dy, -dx]
        elif dx == 0:
            if dy < 0: # rotate 270
                xL, xR, yL, yR = [yL, yR, -xR, -xL]
                x0, y0 = [y0, -x0]
                th1 += 270
                th2 += 270
                dx, dy = [dy, 0]
            else: # dy > 0 rotate 90
                xL, xR, yL, yR = [-yR, -yL, xL, xR]
                x0, y0 = [-y0, x0]
                th1 += 90
                th2 += 90
                dx, dy = [-dy, 0]
        else: # dx < 0
            if dy > 0: # reflect in y axis so on the bottom
                yL, yR = [-yR,-yL]
                y0 *= -1
                th1, th2 = [-th2,-th1] # reflections flip orientation
                dy *= -1
        th1, th2 = normalize_thetas(th1,th2)

        # now affine shift so (x0,y0) at (0,0)
        xL -= x0
        xR -= x0
        yL -= y0
        yR -= y0
        x0 = 0
        y0 = 0

        # print(f"Transformed coords: (x0,y0,xL,xR,yL,yR,th1,th2,dx,dy) = {(x0,y0,xL,xR,yL,yR,th1,th2,dx,dy)}")

        # Now I can proceed as if the center is left or bottom-left of domain, i.e. dx<0 and dy<=0

        th1_rad = th1*0.017453292519943
        th2_rad = th2*0.017453292519943

        bounding_th_R = np.arctan2(yR,xL) # top-left corner is always the upper (Right) theta bound in this reference frame
        bounding_d_R2 = xL**2 + yR**2 + np.zeros((2,1))
        inner_th_1 = np.arctan2(yR,xR) # top-right is always interior in this frame of reference
        if dy == 0:
            inner_1_d2 = np.array([(xL/np.cos(inner_th_1))**2,xR**2 + yR**2]).reshape((2,1))
            bounding_th_L = np.arctan2(yL,xL)
            bounding_d_L2 = xL**2 + yL**2 + np.zeros((2,1))
            inner_th_2 = np.arctan2(yL,xR)
            inner_2_d2 = np.array([(xL/np.cos(inner_th_2))**2,xR**2 + yL**2]).reshape((2,1))
            th_right_d2 = np.array([xL**2,xR**2]).reshape((2,1))
            TH = np.array([bounding_th_L,inner_th_2,0,inner_th_1,bounding_th_R])
            D = np.concatenate((bounding_d_L2,inner_2_d2,th_right_d2,inner_1_d2,bounding_d_R2),axis=1)
        else:
            temp_dx = xL/np.cos(inner_th_1)
            temp_dy = yL/np.sin(inner_th_1)
            inner_1_d2 = np.array([(np.max([temp_dx,temp_dy]))**2,xR**2 + yR**2]).reshape((2,1))
            bounding_th_L = np.arctan2(yL,xR)
            bounding_d_L2 = xR**2 + yL**2 + np.zeros((2,1))
            inner_th_2 = np.arctan2(yL,xL)
            temp_dx =  xR/np.cos(inner_th_1)
            temp_dy = yR/np.sin(inner_th_1)
            inner_2_d2 = np.array([xL**2 + yL**2,(np.min([temp_dx,temp_dy]))**2]).reshape((2,1))
            if inner_th_1 > inner_th_2:
                inner_th_1, inner_th_2 = [inner_th_2, inner_th_1]
                inner_1_d2, inner_2_d2 = [inner_2_d2, inner_1_d2]
            TH = np.array([bounding_th_L,inner_th_1,inner_th_2,bounding_th_R])
            D = np.concatenate((bounding_d_L2,inner_1_d2,inner_2_d2,bounding_d_R2),axis=1)
        # print(f"TH = {TH * 180/np.pi}")
        # print(f"D = {D}")
        th1_inbounds = th1_rad > bounding_th_L and th1_rad < bounding_th_R
        if th1_inbounds and (th1_rad not in TH):
            TH, D = compute_theta_intersection_distances(TH,D,th1_rad,x0,y0,xL,xR,yL,yR,dy)
        th2_inbounds = th2_rad > bounding_th_L and th2_rad < bounding_th_R
        th2_inbounds = th2_inbounds or ((th2_rad-2*np.pi) > bounding_th_L and (th2_rad-2*np.pi) < bounding_th_R) # it's possible that th2 being on [th1,th1+360] might not lie between these theta value, but intersect nonetheless
        if th2_inbounds and (th2_rad not in TH):
            TH, D = compute_theta_intersection_distances(TH,D,th2_rad,x0,y0,xL,xR,yL,yR,dy)
        # print(f"TH = {TH * 180/np.pi}")
        # print(f"D = {D}")
       
        # print(f"th1_rad = {th1_rad}, th2_rad = {th2_rad}")
        TH_rel = (TH-th1_rad) % (2*np.pi)
        # print(f"TH_rel = {TH_rel * 180/np.pi}")
        # print(f"TH_rel = {TH_rel}")
        d_old2 = None

        A = np.concatenate((TH_rel.reshape(1,len(TH_rel)),TH.reshape(1,len(TH)),D),axis=0)
        ord = np.argsort(A[0])
        A = A[:,ord]

        # print(f"A[0] = {A[0] * 180/np.pi}")
        # print(f"A[1] = {A[1] * 180/np.pi}")
        # print(f"A[2:] = {A[2:]}")

        r02 = r0*r0
        r12 = r1*r1
        # print(f"r02 = {r02}, r12 = {r12}")
        # print(f"th1_rad = {th1_rad}, th2_rad = {th2_rad}")
        for idx, abs_th in enumerate(A[1]):
            if A[0,idx] > th2_rad - th1_rad: # then we've finished our rotations
                # print("1. break")
                break
            d_new2 = A[2:4,idx]
            # print(f"\tidx = {idx}\n\tabs_th = {abs_th}\n\td_new2 = {d_new2}\n\td_old2 = {d_old2}")
            if d_new2[0] < r12 and d_new2[1] > r02:
                # print("2. return true")
                return True
            if d_old2 is not None:
                if (d_old2[0] >= r12 and d_new2[1] <= r02) or (d_old2[1] <= r02 and d_new2[0] >= r12):
                    # then the interval passed through the region
                    # print("3. return true")
                    return True
            if abs_th == bounding_th_R: # the sweep of the ray is leaving the domain
                d_old2 = None
            else:
                d_old2 = copy.deepcopy(d_new2)
        # if control passes here, then we have found that the wedge does not intersect the domain
        # print("4. return false")
        return False

    def distance_to_domain_from_within(self, x0, y0, th):
        v1_x = np.cos(th)
        v1_y = np.sin(th)
        if v1_x > 0:
            r_x = (self.plot_xmax - x0) / v1_x
        elif v1_x < 0:
            r_x = (self.plot_xmin - x0) / v1_x
        else:
            r_x = np.inf
        if v1_y > 0:
            r_y = (self.plot_ymax - y0) / v1_y
        elif v1_y < 0:
            r_y = (self.plot_ymin - y0) / v1_y
        else:
            r_y = np.inf
        return r_x if r_x <= r_y else r_y

    def create_figure(self):
        self.figure = plt.figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setStyleSheet("background-color:transparent;")

        self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.canvas.setFocus()

        self.ax0 = self.figure.add_subplot(111, adjustable='box')

        self.format_axis()

        self.mpl_cid = []

        self.mouse_pressed = False
        self.x0y0 = None
        self.current_key_is_shift = False
        self.single_click_grace_period = 0.1 # seconds for a single-click grace period, don't allow movement within this window to cause a redraw
        self.shape_started = False

        self.canvas.update()
        self.canvas.draw()

    def format_axis(self):
        self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
        self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)
        self.ax0.set_aspect(1.0)

    def plot_cell_pos(self):
        self.constrain_preview_to_axes = False
        # print(f"self.main_window.checkbox_dict.keys() = {self.main_window.checkbox_dict.keys()}")
        for ctn in self.main_window.checkbox_dict.keys():
            if self.main_window.checkbox_dict[ctn].isChecked():
                # print(f"writing for {ctn}")
                self.plot_cell_pos_single(ctn)
        self.canvas.update()
        self.canvas.draw()

        self.plot_cells_button.setEnabled(False)

        for b in self.main_window.checkbox_dict.values():
            if b.isEnabled() is True:
                return
        # If control passes here, then all the buttons are disabled and the plotting is done
        self.finish_write_button.setEnabled(True)
        self.finish_append_button.setEnabled(True)

    def max_dist_to_domain(self,x0,y0):
        xL, xR, yL, yR = [self.plot_xmin, self.plot_xmax, self.plot_ymin, self.plot_ymax]
        dx = xL-x0 if (2*x0 > xL+xR) else x0-xR # distance to furtherst vertical edge
        dy = yL-y0 if (2*y0 > yL+yR) else y0-yR # distance to furtherst horizontal edge
        return np.sqrt(dx*dx + dy*dy)
    
    def plot_cell_pos_single(self, cell_type):
        # print(f"cell_type = {cell_type}")
        # print(f"self.main_window.cell_counts.keys() = {self.main_window.cell_counts.keys()}")
        N = self.main_window.cell_counts[cell_type]
        if type(self.preview_patch) is Rectangle:
            # first make sure the rectangle is all in bounds
            if self.constrain_preview_to_axes is False:
                corners = self.preview_patch.get_corners()
                corners = np.array([[min(max(x,self.plot_xmin),self.plot_xmax),min(max(y,self.plot_ymin),self.plot_ymax)] for x,y in corners[[0,2]]])
                self.preview_patch.set_bounds(corners[0,0],corners[0,1],corners[1,0]-corners[0,0],corners[1,1]-corners[0,1])
                self.constrain_preview_to_axes = True
            x0, y0 = self.preview_patch.get_xy()
            width = self.preview_patch.get_width()
            height = self.preview_patch.get_height()
            x = x0 + width * np.random.uniform(size=(N,1))
            y = y0 + height * np.random.uniform(size=(N,1))
            z = np.zeros((N,1))
            self.new_pos = np.concatenate((x,y,z),axis=1)
        elif type(self.preview_patch) is Circle:
            x0, y0 = self.preview_patch.get_center()
            r = self.preview_patch.get_radius()
            r = np.min([r,self.max_dist_to_domain(x0,y0)])
            self.wedge_sample(N, x0, y0, r)
        elif type(self.preview_patch) is Annulus:
            x0, y0 = self.preview_patch.get_center()
            r1 = self.preview_patch.get_radii()[0] # annulus is technically an ellipse, get_radii returns (semi-major,semi-minor) axis lengths, since I'm using circles, these will be the same
            r1 = np.min([r1,self.max_dist_to_domain(x0,y0)])
            width = self.preview_patch.get_width()
            r0 = r1 - width
            r0 = np.max([r0,np.sqrt(self.get_distance2_to_domain(x0,y0)[0])])
            self.wedge_sample(N, x0, y0, r1, r0=r0)
        elif type(self.preview_patch) is Wedge:
            x0, y0 = self.preview_patch.center
            r1 = self.preview_patch.r
            r1 = np.min([r1,self.max_dist_to_domain(x0,y0)])
            th1 = self.preview_patch.theta1  
            th2 = self.preview_patch.theta2
            if th2 == th1:
                pass
            elif ((th2-th1) % 360) == 0:
                # print(f"th1 = {th1}, th2 = {th2}, th2-th1 = {th2-th1}, (th2-th11) % 360 = {(th2-th1) % 360}")
                th2 = th1 + 360
            else:
                th2 -= 360 * ((th2-th1) // 360) # I promise this works if dth=th2-th1 < 0, 0<dth<360, and dth>360. 
            width = self.preview_patch.width
            r0 = r1 - width
            r0 = np.max([r0,np.sqrt(self.get_distance2_to_domain(x0,y0)[0])])
            self.wedge_sample(N, x0, y0, r1, r0=r0, th_lim=(th1*0.017453292519943,th2*0.017453292519943))
        else:
            print("unknown patch")
        self.csv_array[cell_type] = np.append(self.csv_array[cell_type],self.new_pos,axis=0)

        self.circles(self.new_pos, s=8., color=self.color_by_celltype[cell_type], edgecolor='black', linewidth=0.5, alpha=self.alpha_value)

        self.main_window.checkbox_dict[cell_type].setEnabled(False)
        self.main_window.checkbox_dict[cell_type].setChecked(False)
        self.main_window.undo_button[cell_type].setEnabled(True)

    def wedge_sample(self,N,x0,y0,r1, r0=0, th_lim=(0,2*np.pi)):
        i_start = 0
        self.new_pos = np.empty((N,3))
        while i_start < N:
            if r0 == 0:
                d = r1*np.sqrt(np.random.uniform(size=N-i_start))
            else:
                d = np.sqrt(r0*r0 + (r1*r1-r0*r0)*np.random.uniform(size=N-i_start))
            th = th_lim[0] + (th_lim[1]-th_lim[0]) * np.random.uniform(size=N-i_start)
            x = x0 + d * np.cos(th)
            y = y0 + d * np.sin(th)
            xy = np.array([[a,b] for a,b in zip(x,y) if a>=self.plot_xmin and a<=self.plot_xmax and b>=self.plot_ymin and b<=self.plot_ymax])
            if len(xy)==0:
                continue
            # z = np.zeros((len(xy),1))
            self.new_pos[i_start:(i_start+xy.shape[0]),0:2] = xy
            # self.new_pos[range(i_start,i_start+xy.shape[0]),2] = z
            i_start += xy.shape[0]

    def circles(self, pos, s, c='b', vmin=None, vmax=None, **kwargs):
        """
        See https://gist.github.com/syrte/592a062c562cd2a98a83 

        Make a scatter plot of circles. 
        Similar to plt.scatter, but the size of circles are in data scale.
        Parameters
        ----------
        x, y : scalar or array_like, shape (n, )
            Input data
        s : scalar or array_like, shape (n, ) 
            Radius of circles.
        c : color or sequence of color, optional, default : 'b'
            `c` can be a single color format string, or a sequence of color
            specifications of length `N`, or a sequence of `N` numbers to be
            mapped to colors using the `cmap` and `norm` specified via kwargs.
            Note that `c` should not be a single numeric RGB or RGBA sequence 
            because that is indistinguishable from an array of values
            to be colormapped. (If you insist, use `color` instead.)  
            `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
        vmin, vmax : scalar, optional, default: None
            `vmin` and `vmax` are used in conjunction with `norm` to normalize
            luminance data.  If either are `None`, the min and max of the
            color array is used.
        kwargs : `~matplotlib.collections.Collection` properties
            Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
            norm, cmap, transform, etc.
        Returns
        -------
        paths : `~matplotlib.collections.PathCollection`
        Examples
        --------
        a = np.arange(11)
        circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
        plt.colorbar()
        License
        --------
        This code is under [The BSD 3-Clause License]
        (http://opensource.org/licenses/BSD-3-Clause)
        """
        x, y, _ = pos.T
        if np.isscalar(c):
            kwargs.setdefault('color', c)
            c = None

        if 'fc' in kwargs:
            kwargs.setdefault('facecolor', kwargs.pop('fc'))
        if 'ec' in kwargs:
            kwargs.setdefault('edgecolor', kwargs.pop('ec'))
        if 'ls' in kwargs:
            kwargs.setdefault('linestyle', kwargs.pop('ls'))
        if 'lw' in kwargs:
            kwargs.setdefault('linewidth', kwargs.pop('lw'))
        # You can set `facecolor` with an array for each patch,
        # while you can only set `facecolors` with a value for all.

        zipped = np.broadcast(x, y, s)
        patches = [Circle((x_, y_), s_)
                for x_, y_, s_ in zipped]
        collection = PatchCollection(patches, **kwargs)
        if c is not None:
            c = np.broadcast_to(c, zipped.shape).ravel()
            collection.set_array(c)
            collection.set_clim(vmin, vmax)

        # ax = plt.gca()
        # ax.add_collection(collection)
        # ax.autoscale_view()
        self.ax0.add_collection(collection)
        self.ax0.autoscale_view()
        # plt.draw_if_interactive()
        if c is not None:
            # plt.sci(collection)
            self.ax0.sci(collection)
        # return collection

    def cancel_cb(self):
        self.hide() # this will work for now, but maybe a better way to handle closing the window?
        pass

    def finish_write_button_cb(self):
        self.check_for_new_celldefs()
        self.set_file_name()
        with open(self.full_fname, 'w') as f:
            f.write('x,y,z,type\n')
        self.add_cell_positions_to_file()

    def finish_append_button_cb(self):
        self.check_for_new_celldefs()
        self.set_file_name()
        with open(self.full_fname, 'r') as f:
            first_line = f.readline()
            if "x,y,z,type" not in first_line:
                print(f"{self.full_fname} is not properly formatted for appending.\nIt needs to start with 'x,y,z,type,...'")
                return
        self.add_cell_positions_to_file()

    def check_for_new_celldefs(self):
        for cell_type in self.main_window.cell_types_list_final:
            if cell_type in self.main_window.celldef_tab.celltypes_list:
                print(f"BioinfImportPlotWindow: {cell_type} found in current list of cell types. Not appending this...")
            else:
                self.main_window.celldef_tab.new_cell_def_named(cell_type)

    def set_file_name(self):
        dir_name = self.main_window.csv_folder.text()
        # print(f'dir_name={dir_name}<end>')
        if len(dir_name) > 0 and not os.path.isdir(dir_name):
            os.makedirs(dir_name)
            time.sleep(1)
        self.full_fname = os.path.join(dir_name,self.main_window.csv_file.text())
        Path(self.full_fname).touch() # make sure the file exists
        print("save_cb(): self.full_fname=",self.full_fname)

    def add_cell_positions_to_file(self):
        with open(self.full_fname, 'a') as f:
            for cell_type in self.csv_array.keys():
                for pos in self.csv_array[cell_type]:
                    f.write(f'{pos[0]},{pos[1]},{pos[2]},{cell_type}\n')
        self.main_window.ics_tab.clear_cb()
        self.main_window.ics_tab.import_from_file(self.full_fname)
        self.main_window.ics_tab.tab_widget.setCurrentIndex(self.main_window.ics_tab.base_tab_id)
        self.main_window.ics_tab.csv_folder.setText(self.main_window.csv_folder.text())
        self.main_window.ics_tab.output_file.setText(self.main_window.csv_file.text())
        self.close()
        self.main_window.window.close()
        print("BioinfImportWindow: Colors will likely change in the ICs tab due to previous cell types being present.")

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
                QCheckBox:disabled
                {
                    background-color:lightgray;
                }
                """
        self.setStyleSheet(checkbox_style)

class BioinfImport(QWidget):
    def __init__(self, config_tab, celldef_tab, ics_tab, bioinf_import_test):
        super().__init__()
        if HAVE_ANNDATA is False:
            vbox = QVBoxLayout()
            s = "This tab allows the import of an anndata object to generate cell initial conditions."
            s += "\nThis tab will read the adata.obs column selected for cell types and walk you through how to place them."
            s += "\nHowever, you do not have anndata installed. You need to install that in your environment:\n\t1. pip install anndata\n---or---\n\t2. conda install anndata -c conda-forge"
            s += "\n\nAfter installing, restart studio."
            label = QLabel(s)
            label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            vbox.addWidget(label)
            base_widget = QWidget()
            base_widget.setLayout(vbox)
            self.layout = QVBoxLayout(self)  # leave this!
            self.layout.addWidget(base_widget)
            return

        self.config_tab = config_tab
        self.celldef_tab = celldef_tab
        self.ics_tab = ics_tab
        self.default_time_units = "min"

        # self.full_fname = f"{self.csv_folder.text()}/{self.csv_file.text()}"
        self.qlineedit_style_sheet = """
            QLineEdit:disabled {
                background-color: rgb(200,200,200);
                color: black;
            }
            QLineEdit:enabled {
                background-color: white;
                color: black;
            }
            """
        self.qpushbutton_style_sheet = """
            QPushButton:enabled {
                background-color : lightgreen;
            }
            QPushButton:disabled {
                background-color : grey;
            }
            """

        vbox = QVBoxLayout()
        vbox.addStretch(1)
        vbox.addWidget(QLabel("Importing",styleSheet="QLabel {background-color : orange;}",alignment=QtCore.Qt.AlignCenter,maximumHeight=20))
        hbox = QHBoxLayout()
        self.import_button = QPushButton("Import from AnnData")
        self.import_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        self.import_button.clicked.connect(self.import_cb)
        hbox.addWidget(self.import_button)

        label = QLabel("Column name with cell type (prompt will follow if left blank): ")
        hbox.addWidget(label)

        self.column_line_edit = QLineEdit()
        self.column_line_edit.setEnabled(True)
        self.column_line_edit.setText('leiden')
        hbox.addWidget(self.column_line_edit)

        vbox.addLayout(hbox)

        vbox.addWidget(QHLine())

        vbox.addWidget(QLabel("Cell initial conditions file",styleSheet="QLabel {background-color : orange;}",alignment=QtCore.Qt.AlignCenter,maximumHeight=20))

        hbox = QHBoxLayout()
        label = QLabel("folder")
        hbox.addWidget(label)

        self.csv_folder = QLineEdit("./config")
        hbox.addWidget(self.csv_folder)
        label = QLabel("file")
        hbox.addWidget(label)

        self.csv_file = QLineEdit("cells.csv")
        hbox.addWidget(self.csv_file)

        vbox.addLayout(hbox)

        vbox.addStretch(1)

        base_widget = QWidget()
        base_widget.setLayout(vbox)
        self.layout = QVBoxLayout(self)  # leave this!
        self.layout.addWidget(base_widget)

        if bioinf_import_test:
            self.import_file("./data/pbmc3k_clustered.h5ad")
            self.continue_to_rename()
            self.continue_on_rename_cb()
            self.continue_to_cell_pos_cb()


    def select_all_button_cb(self):
        for cbd in self.checkbox_dict.values():
            if cbd.isEnabled():
                cbd.setChecked(True)
        bval = self.is_any_cell_type_button_group_checked()
        if self.ics_plot_area:
            self.ics_plot_area.plot_cells_button.setEnabled(bval)
    
    def deselect_all_button_cb(self):
        for cbd in self.checkbox_dict.values():
            if cbd.isEnabled():
                cbd.setChecked(False)
        self.show_plot_button.setEnabled(False)
        if self.ics_plot_area:
            self.ics_plot_area.plot_cells_button.setEnabled(False)

    def import_cb(self):
        full_file_path = QFileDialog.getOpenFileName(self,'',".")
        file_path = full_file_path[0]
        # file_path = "./data/pbmc3k_clustered.h5ad"

        self.import_file(file_path)

    def import_file(self,file_path):
        try:
            self.adata = anndata.read_h5ad(file_path)
        except:
            print(f"Import failed...")
            return

        print("------------anndata object loaded-------------")
        # print(self.adata)

        self.editing_cell_type_names = True
        col_names = list(self.adata.obs.columns)
        auto_continue = False
        if self.column_line_edit.text() in col_names:
            s = "Select column that contains cell type info:"
            auto_continue = True
        elif self.column_line_edit.text() != "":
            s = f"{self.column_line_edit.text()} not found in the columns of the obs matrix.\nSelect from the following:"
        else:
            s = "Select column that contains cell type info:"
        
        vbox = QVBoxLayout()
        label = QLabel(s)
        vbox.addWidget(label)

        self.column_combobox = QComboBox()
        self.column_combobox.currentIndexChanged.connect(self.column_combobox_changed_cb)
        for col_name in col_names:
            self.column_combobox.addItem(col_name)
        vbox.addWidget(self.column_combobox)

        hbox = QHBoxLayout()
        continue_button = QPushButton("Continue \u2192")
        continue_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        continue_button.clicked.connect(self.continue_to_cell_type_names_cb)
        hbox.addWidget(continue_button)

        vbox.addLayout(hbox)

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        if auto_continue:
            self.current_column = self.column_line_edit.text()
            self.continue_to_cell_type_names_cb()
            return
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
        # continue_button = self.list_remaining_with_checkboxes("Keep")
        # label = self.list_current_cell_type_names() 

        keep_color = "lightgreen"
        merge_color = "yellow"
        delete_color = "#FFCCCB"
        # label = QLabel(f"""The following cell types were found. \
        #                 Choose which to <b style="color: {keep_color}; background-color: black">KEEP</b>, \
        #                 <b style="color: {merge_color}; background-color: black">MERGE</b>, and \
        #                 <b style="color: {delete_color}; background-color: black">DELETE</b>.
        #                 By default, all are kept.\
        #                 """)
        label = QLabel(f"""The following cell types were found.<br>\
                        Choose which to <b style="background-color: {keep_color};">KEEP</b>, \
                        <b style="background-color: {merge_color};">MERGE</b>, and \
                        <b style="background-color: {delete_color};">DELETE</b>.<br>\
                        By default, all are kept.\
                        """)
        # label.setWordWrap(True)
        # label.setM
        checkbox_style_template = lambda x : f"""
                QCheckBox {{
                    color : black;
                    font-weight: bold;
                    background-color: {x};
                }}
                QCheckBox::indicator:checked {{
                    background-color: rgb(255,255,255);
                    border: 1px solid #5A5A5A;
                    width : 15px;
                    height : 15px;
                    border-radius : 3px;
                    image: url(images:checkmark.png);
                }}
                QCheckBox::indicator:unchecked
                {{
                    background-color: {"rgb(255,255,255)" if x==keep_color else "rgb(0,0,0)"};
                    border: 1px solid #5A5A5A;
                    width : 15px;
                    height : 15px;
                    border-radius : 3px;
                }}
                """
        self.checkbox_style = {}
        self.checkbox_style["keep"] = checkbox_style_template(keep_color)
        self.checkbox_style["merge"] = checkbox_style_template(merge_color)
        self.checkbox_style["delete"] = checkbox_style_template(delete_color)
        self.cell_type_dict = {}

        vbox = QVBoxLayout()
        # label = QLabel(f"Select cell types :")
        vbox.addWidget(label)
        self.checkbox_dict = {}
        self.keep_button = {}
        self.checkbox_group = QButtonGroup(exclusive=False)
        self.checkbox_group.buttonToggled.connect(self.cell_type_toggled_cb)
        self.merge_id = 0
        rename_button_style_sheet_template = lambda x : f"""
            QPushButton {{
                color : black;
                font-weight : bold;
            }}
            QPushButton:enabled {{
                background-color: {x};
            }}
            QPushButton:disabled {{
                background-color : grey;
            }}
            """
        self.rename_button_style_sheet = {}
        self.rename_button_style_sheet["keep"] = rename_button_style_sheet_template(keep_color)
        self.rename_button_style_sheet["merge"] = rename_button_style_sheet_template(merge_color)
        self.rename_button_style_sheet["delete"] = rename_button_style_sheet_template(delete_color)
        row_height = 30
        print(f"self.cell_types_list_original = {self.cell_types_list_original}")
        for cell_type in self.cell_types_list_original:
            hbox = QHBoxLayout()
            self.checkbox_dict[cell_type] = QCheckBox(cell_type,styleSheet=self.checkbox_style["keep"])
            self.checkbox_dict[cell_type].setChecked(False)
            self.checkbox_dict[cell_type].setEnabled(True)
            self.checkbox_dict[cell_type].setFixedHeight(row_height)
            self.checkbox_group.addButton(self.checkbox_dict[cell_type])
            
            self.keep_button[cell_type] = QPushButton("Keep",enabled=False,objectName=cell_type)
            self.keep_button[cell_type].setStyleSheet(rename_button_style_sheet_template(keep_color))
            self.keep_button[cell_type].setFixedHeight(row_height)
            self.keep_button[cell_type].clicked.connect(self.continue_on_keep_cb)

            hbox.addWidget(self.checkbox_dict[cell_type])
            hbox.addStretch(1)
            hbox.addWidget(self.keep_button[cell_type])
            vbox.addLayout(hbox)
            self.cell_type_dict[cell_type] = cell_type

        hbox = QHBoxLayout()

        self.merge_button = QPushButton("Merge",enabled=False,styleSheet=rename_button_style_sheet_template(merge_color))
        self.merge_button.clicked.connect(self.continue_on_merge_cb)
        hbox.addWidget(self.merge_button)

        self.delete_button = QPushButton("Delete",enabled=False,styleSheet=rename_button_style_sheet_template(delete_color))
        self.delete_button.clicked.connect(self.continue_on_delete_cb)
        hbox.addWidget(self.delete_button)

        vbox.addLayout(hbox)
        vbox.addWidget(QHLine())

        hbox = QHBoxLayout()

        go_back_button = GoBackButton(self)
        # go_back_button = QPushButton("\u2190 Go back")
        # go_back_button.setStyleSheet(f"QPushButton {{background-color: lightgreen; color: black;}}")
        # go_back_button.clicked.connect(self.go_back_to_prev_window)

        accept_button = QPushButton("Finish \u2192")
        accept_button.setStyleSheet(f"QPushButton {{background-color: lightgreen; color: black;}}")
        accept_button.clicked.connect(self.continue_to_rename)

        hbox.addWidget(go_back_button)
        hbox.addWidget(accept_button)
        vbox.addLayout(hbox)
        # self.checkbox_dict = create_checkboxes_for_cell_types(vbox, self.remaining_cell_types_list_original)

        # hbox = QHBoxLayout()
        # continue_button = QPushButton("Continue")
        # continue_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        # hbox.addWidget(continue_button)

        # go_back = QPushButton("Go back")
        # go_back.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        # go_back.clicked.connect(self.go_back_cb)
        # hbox.addWidget(go_back)

        # vbox.addLayout(hbox)


        
        # vbox = QVBoxLayout()
        # vbox.addWidget(label)

        # label = QLabel("Accept these as is or edit (merge and delete, and rename)?")
        # vbox.addWidget(label)

        # hbox = QHBoxLayout()
        # accept_button = QPushButton("Accept")
        # accept_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        # accept_button.clicked.connect(self.accept_as_is_cb)
        # hbox.addWidget(accept_button)

        # edit_button = QPushButton("Edit")
        # edit_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        # edit_button.clicked.connect(self.edit_cell_types_cb)
        # hbox.addWidget(edit_button)

        # vbox.addLayout(hbox)

        self.previous_windows = [self.window]
        self.previous_windows[0].hide()

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

    def cell_type_toggled_cb(self):
        checked_count = 0
        enable_delete_button = False
        enable_merge_button = False
        for checkbox in self.checkbox_dict.values():
            checked_count += checkbox.isChecked()
            if checked_count == 1:
                enable_delete_button = True
            if checked_count == 2:
                enable_merge_button = True
                break
        self.delete_button.setEnabled(enable_delete_button)
        self.merge_button.setEnabled(enable_merge_button)

    def continue_to_cell_counts(self):

        names_width = 100
        counts_width = 120
        props_width = 120
        manual_width = 120
        confluence_width = 150
        
        vbox = QVBoxLayout()

        vbox.addWidget(QLabel("Set how many of each cell type to place. You may use the counts found in your data, or use proportions or confluence to scale these counts.\nYou may also set the counts manually; these values will track the others, allowing you to adjust from those baselines."))

        vbox_cols = [QVBoxLayout() for i in range(5)]

        # hbox = QHBoxLayout()
        label = QLabel("Cell Type")
        label.setFixedWidth(names_width)
        vbox_cols[0].addWidget(label)

        self.counts_button_group = QButtonGroup()
        self.counts_button_group.idToggled.connect(self.counts_button_cb)
        counts_button_group_next_id = 0

        self.use_counts_as_is_radio_button = QRadioButton("Use counts")
        self.use_counts_as_is_radio_button.setFixedWidth(counts_width)
        self.use_counts_as_is_radio_button.setChecked(True)
        self.counts_button_group.addButton(self.use_counts_as_is_radio_button,counts_button_group_next_id)
        counts_button_group_next_id += 1

        self.use_props_radio_button = QRadioButton("Use proportions")
        self.use_props_radio_button.setFixedWidth(props_width)
        self.use_props_radio_button.setChecked(False)
        self.counts_button_group.addButton(self.use_props_radio_button,counts_button_group_next_id)
        counts_button_group_next_id += 1

        self.use_confluence_radio_button = QRadioButton("Set confluence (%)")
        self.use_confluence_radio_button.setFixedWidth(confluence_width)
        self.use_confluence_radio_button.setChecked(False)
        self.counts_button_group.addButton(self.use_confluence_radio_button,counts_button_group_next_id)
        counts_button_group_next_id += 1

        self.use_manual_radio_button = QRadioButton("Set manually")
        self.use_manual_radio_button.setFixedWidth(manual_width)
        self.use_manual_radio_button.setChecked(False)
        self.counts_button_group.addButton(self.use_manual_radio_button,counts_button_group_next_id)
        counts_button_group_next_id += 1

        vbox_cols[1].addWidget(self.use_counts_as_is_radio_button)
        vbox_cols[2].addWidget(self.use_props_radio_button)
        vbox_cols[3].addWidget(self.use_confluence_radio_button)
        vbox_cols[4].addWidget(self.use_manual_radio_button)

        # vbox.addLayout(hbox)

        self.cell_counts = {}
        for cell_type in self.cell_types_list_final:
            self.cell_counts[cell_type] = 0
        for ctn in self.cell_types_final:
            self.cell_counts[ctn] += 1

        self.cell_type_props = [self.cell_counts[cell_type]/len(self.cell_types_final) for cell_type in self.cell_types_list_final]

        self.type_prop = {}
        self.type_manual = {}
        self.type_confluence = {}

        self.prop_box_callback_paused = False
        self.conf_box_callback_paused = False

        num_validator = QtGui.QIntValidator()
        num_validator.setBottom(0)

        self.setup_confluence_info()

        for idx, cell_type in enumerate(self.cell_types_list_final):
            hbox = QHBoxLayout()
            label = QLabel(cell_type)
            label.setFixedWidth(names_width)

            type_count = QLineEdit(enabled=False)
            type_count.setText(str(self.cell_counts[cell_type]))
            type_count.setFixedWidth(counts_width)
            type_count.setStyleSheet(self.qlineedit_style_sheet)

            self.type_prop[cell_type] = QLineEdit(enabled=False)
            self.type_prop[cell_type].setText(str(self.cell_counts[cell_type]))
            self.type_prop[cell_type].setFixedWidth(props_width)
            self.type_prop[cell_type].setStyleSheet(self.qlineedit_style_sheet)
            self.type_prop[cell_type].setValidator(num_validator)
            self.type_prop[cell_type].setObjectName(str(idx))
            self.type_prop[cell_type].textChanged.connect(self.prop_box_changed_cb)

            self.type_confluence[cell_type] = QLineEdit(enabled=False)
            self.type_confluence[cell_type].setFixedWidth(confluence_width)
            self.type_confluence[cell_type].setStyleSheet(self.qlineedit_style_sheet)
            self.type_confluence[cell_type].setValidator(QtGui.QDoubleValidator(bottom=0))
            self.type_confluence[cell_type].setObjectName(str(idx))
            self.type_confluence[cell_type].textChanged.connect(self.confluence_box_changed_cb)

            self.type_manual[cell_type] = QLineEdit(enabled=False)
            self.type_manual[cell_type].setText(str(self.cell_counts[cell_type]))
            self.type_manual[cell_type].setFixedWidth(manual_width)
            self.type_manual[cell_type].setStyleSheet(self.qlineedit_style_sheet)
            self.type_manual[cell_type].setValidator(num_validator)
            self.type_manual[cell_type].setObjectName(str(idx))
            
            vbox_cols[0].addWidget(label)
            vbox_cols[1].addWidget(type_count)
            vbox_cols[2].addWidget(self.type_prop[cell_type])
            vbox_cols[3].addWidget(self.type_confluence[cell_type])
            vbox_cols[4].addWidget(self.type_manual[cell_type])

            # vbox.addLayout(hbox)
        
        # hbox = QHBoxLayout()
        label = QLabel("Total")
        label.setFixedWidth(names_width)

        type_count = QLineEdit(enabled=False)
        type_count.setText(str(len(self.cell_types_final)))
        type_count.setFixedWidth(counts_width)
        type_count.setStyleSheet(self.qlineedit_style_sheet)

        self.total_prop = QLineEdit(enabled=False)
        self.total_prop.setText(str(len(self.cell_types_final)))
        self.total_prop.setFixedWidth(props_width)
        self.total_prop.setStyleSheet(self.qlineedit_style_sheet)
        self.total_prop.setValidator(num_validator)
        self.total_prop.textChanged.connect(self.prop_box_changed_cb)
        self.total_prop.setObjectName("total_prop")

        self.total_manual = QLineEdit(enabled=False)
        self.total_manual.setText(str(len(self.cell_types_final)))
        self.total_manual.setFixedWidth(manual_width)
        self.total_manual.setStyleSheet(self.qlineedit_style_sheet)
        self.total_manual.setValidator(num_validator)

        self.total_conf = QLineEdit(enabled=False)
        self.total_conf.setFixedWidth(confluence_width)
        self.total_conf.setStyleSheet(self.qlineedit_style_sheet)
        self.total_conf.setValidator(QtGui.QDoubleValidator(bottom=0))
        self.total_conf.textChanged.connect(self.confluence_box_changed_cb)
        self.total_conf.setObjectName("total_conf")
        self.total_conf.setText("100")

        vbox_cols[0].addWidget(label)
        vbox_cols[1].addWidget(type_count)
        vbox_cols[2].addWidget(self.total_prop)
        vbox_cols[3].addWidget(self.total_conf)
        vbox_cols[4].addWidget(self.total_manual)


        hbox = QHBoxLayout()
        for idx, vb in enumerate(vbox_cols):
            hbox.addLayout(vb)
            if idx != (len(vbox_cols)-1):
                hbox.addWidget(QVLine())
        vbox.addLayout(hbox)

        hbox = QHBoxLayout()

        self.go_back_to_rename = GoBackButton(self)
        # self.go_back_to_rename = QPushButton("\u2190 Go back")
        # self.go_back_to_rename.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        # self.go_back_to_rename.clicked.connect(self.go_back_to_prev_window)

        self.continue_to_cell_pos = QPushButton("Continue \u2192")
        self.continue_to_cell_pos.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        self.continue_to_cell_pos.clicked.connect(self.continue_to_cell_pos_cb)

        hbox.addWidget(self.go_back_to_rename)
        hbox.addWidget(self.continue_to_cell_pos)
        
        vbox.addLayout(hbox)

        self.previous_windows.append(self.window)
        self.previous_windows[-1].hide()

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

        print(f"self.counts_button_group.checkedButton().text() = {self.counts_button_group.checkedButton().text()}")
    
    # def go_back_to_rename(self):

    def prop_box_changed_cb(self, text):
        # print(f"self.prop_box_callback_paused = {self.prop_box_callback_paused}")
        if self.prop_box_callback_paused:
            return
        type_prop_sender = self.sender()
        if type_prop_sender.hasAcceptableInput() is False:
            return
        # text = type_prop_sender.text()
        # print(f"text = {text}")
        current_name = type_prop_sender.objectName()
        # print(f"type_prop_sender.objectName() = {current_name}")
        self.prop_box_callback_paused = True
        if current_name=="total_prop":
            mult = int(text)
            self.total_manual.setText(text)
        else:
            # print(f"int(text) = {int(text)}")
            # print(f"self.cell_type_props[int(current_name)] = {self.cell_type_props[int(current_name)]}")
            mult = int(text) / self.cell_type_props[int(current_name)]
            self.total_prop.setText(str(round(mult)))
            self.total_manual.setText(str(round(mult)))

        for idx, cell_type in enumerate(self.cell_types_list_final):
            # print(f"idx = {idx}, cell_type = {cell_type}")
            if current_name==str(idx):
                continue
            self.type_prop[cell_type].setText(str(round(mult * self.cell_type_props[idx])))
            self.type_manual[cell_type].setText(str(round(mult * self.cell_type_props[idx])))
        self.prop_box_callback_paused = False

    def confluence_box_changed_cb(self, text):
        # print(f"self.conf_box_callback_paused = {self.conf_box_callback_paused}")
        if self.conf_box_callback_paused:
            return
        type_conf_sender = self.sender()
        if type_conf_sender.hasAcceptableInput() is False:
            return
        # text = type_conf_sender.text()
        # print(f"text = {text}")
        current_name = type_conf_sender.objectName()
        # print(f"type_conf_sender.objectName() = {current_name}")
        self.conf_box_callback_paused = True
        current_conf = float(text)
        if current_name=="total_conf":
            mult = current_conf
            mult /= self.prop_dot_ratios
        else:
            # print(f"int(text) = {int(text)}")
            # print(f"self.cell_type_props[int(current_name)] = {self.cell_type_props[int(current_name)]}")
            current_idx = int(current_name)
            mult = current_conf
            mult /= self.prop_total_area_one[self.cell_types_list_final[current_idx]]
            mult /= self.cell_type_props[current_idx]

        total_conf = 0
        total_n = 0
        for idx, cell_type in enumerate(self.cell_types_list_final):
            # print(f"idx = {idx}, cell_type = {cell_type}")
            if current_name==str(idx):
                n = round(0.01*current_conf/self.prop_total_area_one[cell_type])
                self.type_manual[cell_type].setText(str(n))
                total_conf += current_conf
                total_n += n
                continue
            new_conf = mult * self.cell_type_props[idx] * self.prop_total_area_one[cell_type]
            n = round(0.01*new_conf/self.prop_total_area_one[cell_type])
            total_conf += new_conf
            total_n += n
            self.type_confluence[cell_type].setText(str(new_conf))
            self.type_manual[cell_type].setText(str(n))
        if current_name!="total_conf":
            self.total_conf.setText(str(total_conf))
        self.total_manual.setText(str(total_n))
        
        if total_conf > 100:
            self.total_conf.setStyleSheet("QLineEdit {background-color : red; color : black;}")
        else:
            self.total_conf.setStyleSheet(self.qlineedit_style_sheet)
            
        self.conf_box_callback_paused = False

    def counts_button_cb(self):
        enable_props = self.counts_button_group.checkedId()==1
        enable_confluence = self.counts_button_group.checkedId()==2
        enable_manual = self.counts_button_group.checkedId()==3

        for k in self.type_prop.keys():
            self.type_prop[k].setEnabled(enable_props)
        self.total_prop.setEnabled(enable_props)
        
        for k in self.type_confluence.keys():
            self.type_confluence[k].setEnabled(enable_confluence)
        self.total_conf.setEnabled(enable_confluence)
        if enable_confluence and float(self.total_conf.text()) > 100:
            self.total_conf.setStyleSheet("QLineEdit {background-color : red; color : black;}")
        else:
            self.total_conf.setStyleSheet(self.qlineedit_style_sheet)
        
        for k in self.type_manual.keys():
            self.type_manual[k].setEnabled(enable_manual)
        self.total_manual.setEnabled(enable_manual)

    def continue_to_cell_pos_cb(self):
        if self.counts_button_group.checkedId()==0: # use counts found in data file
            pass
        elif self.counts_button_group.checkedId()==1: # use props found in data file
            for cell_type in self.cell_types_list_final:
                self.cell_counts[cell_type] = int(self.type_prop[cell_type].text())
        elif self.counts_button_group.checkedId()==2: # set by confluence
            for cell_type in self.cell_types_list_final:
                self.cell_counts[cell_type] = round(0.01 * float(self.type_confluence[cell_type].text()) / self.prop_total_area_one[cell_type])
        elif self.counts_button_group.checkedId()==3: # manually set
            for cell_type in self.cell_types_list_final:
                self.cell_counts[cell_type] = int(self.type_manual[cell_type].text())

        # print(f"self.cell_counts = {self.cell_counts}")

        self.cell_types_to_place = self.cell_types_list_final
        self.ics_plot_area = None
        self.create_cell_type_scroll_area()
        self.create_pos_scroll_area()

        splitter = QSplitter()
        splitter.addWidget(self.cell_type_scroll_area)
        splitter.addWidget(self.pos_scroll_area)

        vbox = QVBoxLayout()
        vbox.addWidget(splitter)

        hbox = QHBoxLayout()
        go_back_button = GoBackButton(self)
        # go_back_button = QPushButton("\u2190 Go back")
        # go_back_button.setStyleSheet(self.qpushbutton_style_sheet)
        # go_back_button.clicked.connect(self.go_back_to_prev_window)

        self.show_plot_button = QPushButton("Show plot window \u2193",enabled=False)
        self.show_plot_button.setStyleSheet(self.qpushbutton_style_sheet)
        self.show_plot_button.clicked.connect(self.show_plot_button_cb)

        hbox.addWidget(go_back_button)
        hbox.addWidget(self.show_plot_button)

        vbox.addLayout(hbox)

        self.previous_windows.append(self.window)
        self.previous_windows[-1].hide()

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()
        
        # print(f"self.cell_counts = {self.cell_counts}")
            
    def show_plot_button_cb(self):
        if self.ics_plot_area is None:
            self.create_ics_plot_area()
        self.ics_plot_area.hide()
        self.ics_plot_area.show()

    def create_cell_type_scroll_area(self):
        vbox_main = QVBoxLayout()
        label = QLabel("Select cell type(s) to place.\nGreyed out cell types have already been placed.")
        vbox_main.addWidget(label)

        hbox_mid = QHBoxLayout()
        vbox_mid_checkboxes = QVBoxLayout()
        self.cell_type_button_group = QButtonGroup(exclusive=False)
        # self.cell_type_button_group.setExclusive(False)
        self.cell_type_button_group.buttonClicked.connect(self.cell_type_button_group_cb)
        self.checkbox_dict = create_checkboxes_for_cell_types(vbox_mid_checkboxes, self.cell_types_to_place)
        self.undo_button = {}
        for cbd in self.checkbox_dict.values():
            self.cell_type_button_group.addButton(cbd)
        vbox_mid_undos = QVBoxLayout()
        for cell_type in self.cell_types_to_place:
            self.undo_button[cell_type] = QPushButton("Undo",enabled=False,styleSheet=self.qpushbutton_style_sheet,objectName=cell_type)
            self.undo_button[cell_type].clicked.connect(self.undo_button_cb)
            vbox_mid_undos.addWidget(self.undo_button[cell_type])
        
        hbox_mid.addLayout(vbox_mid_checkboxes)
        hbox_mid.addLayout(vbox_mid_undos)

        vbox_main.addLayout(hbox_mid)

        self.select_all_button = QPushButton("Select remaining",styleSheet=self.qpushbutton_style_sheet)
        self.select_all_button.clicked.connect(self.select_all_button_cb)
        
        self.deselect_all_button = QPushButton("Unselect all",styleSheet=self.qpushbutton_style_sheet)
        self.deselect_all_button.clicked.connect(self.deselect_all_button_cb)
        hbox = QHBoxLayout()
        hbox.addWidget(self.select_all_button)
        hbox.addWidget(self.deselect_all_button)
        
        vbox_main.addLayout(hbox)

        cell_type_scroll_area_widget = QWidget()
        cell_type_scroll_area_widget.setLayout(vbox_main)

        self.cell_type_scroll_area = QScrollArea()
        self.cell_type_scroll_area.setWidget(cell_type_scroll_area_widget)

    def undo_button_cb(self):
        undone_cell_type = self.sender().objectName()
        self.ics_plot_area.csv_array[undone_cell_type] = np.empty((0,3))
        self.ics_plot_area.ax0.cla()
        self.ics_plot_area.format_axis()
        for cell_type in self.ics_plot_area.csv_array.keys():
            self.ics_plot_area.circles(self.ics_plot_area.csv_array[cell_type], s=8., color=self.ics_plot_area.color_by_celltype[cell_type], edgecolor='black', linewidth=0.5, alpha=self.ics_plot_area.alpha_value)

        self.ics_plot_area.sync_par_area() # easy way to redraw the patch for current plotting
        
        self.checkbox_dict[undone_cell_type].setEnabled(True)
        self.checkbox_dict[undone_cell_type].setChecked(False)
        self.undo_button[undone_cell_type].setEnabled(False)
        self.ics_plot_area.finish_write_button.setEnabled(False)
        self.ics_plot_area.finish_append_button.setEnabled(False)

    def cell_type_button_group_cb(self):
        bval = self.is_any_cell_type_button_group_checked()
        if self.ics_plot_area is None:
            return
        if bval:
            # The plot is created and at least one is checked. See if the parameters are ready
            self.ics_plot_area.sync_par_area() # this call is overkill, I just want to see if the parameters call for the Plot button being enabled
        else:
            self.ics_plot_area.plot_cells_button.setEnabled(False)

    def is_any_cell_type_button_group_checked(self):
        bval = self.cell_type_button_group.checkedButton() is not None
        self.show_plot_button.setEnabled(bval)
        return bval

    def cell_pos_button_group_cb(self):
        if self.ics_plot_area:
            # print("syncing...")
            self.ics_plot_area.sync_par_area()
            self.ics_plot_area.shape_started = False # for helping with default parameters when using the mouse
        
    def create_pos_scroll_area(self):
        self.cell_pos_button_group = QButtonGroup()
        self.cell_pos_button_group.setExclusive(True)
        self.cell_pos_button_group.idToggled.connect(self.cell_pos_button_group_cb)
        next_button_id = 0
        
        button_width = 250
        button_height = 250
        icon_width = round(0.8 * button_width)
        icon_height = round(0.8 * button_height)
        master_vbox = QVBoxLayout()

        qpushbutton_style_sheet = """
            QPushButton {
                background-color : lightblue;
                color : black;
            }
            QPushButton::unchecked {
                background-color : lightblue;
            }
            QPushButton::checked {
                background-color : black;
            }
            """

        label = QLabel("Select an option for how to place selected cell types:")
        label.setAlignment(QtCore.Qt.AlignLeft)

        master_vbox.addWidget(label)

        # grid_layout = QHBoxLayout()
        grid_layout = QGridLayout()
        rI = 0
        cI = 0
        cmax = 1
            
        vbox = QVBoxLayout()
        full_rectangle_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/scatter_square.svg"))
        full_rectangle_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        full_rectangle_button.setIconSize(size)
        full_rectangle_button.setCheckable(True)
        full_rectangle_button.setChecked(True)
        full_rectangle_button.setStyleSheet(qpushbutton_style_sheet) 
        self.cell_pos_button_group.addButton(full_rectangle_button,next_button_id)
        next_button_id += 1
        vbox.addWidget(full_rectangle_button)

        label = QLabel("Everywhere")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        grid_layout.addLayout(vbox,rI,cI,1,1)
        rI, cI = [rI,cI+1] if cI < cmax else [rI+1,0]
            
        vbox = QVBoxLayout()
        partial_rectangle_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/rectangle.svg"))
        partial_rectangle_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        partial_rectangle_button.setIconSize(size)
        partial_rectangle_button.setCheckable(True)
        partial_rectangle_button.setStyleSheet(qpushbutton_style_sheet) 
        self.cell_pos_button_group.addButton(partial_rectangle_button,next_button_id)
        next_button_id += 1
        vbox.addWidget(partial_rectangle_button)

        label = QLabel("Rectangle")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        grid_layout.addLayout(vbox,rI,cI,1,1)
        rI, cI = [rI,cI+1] if cI < cmax else [rI+1,0]

        vbox = QVBoxLayout()
        disc_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/disc.svg"))
        disc_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        disc_button.setIconSize(size)
        disc_button.setCheckable(True)
        disc_button.setStyleSheet(qpushbutton_style_sheet) 
        self.cell_pos_button_group.addButton(disc_button,next_button_id)
        next_button_id += 1
        vbox.addWidget(disc_button)

        label = QLabel("Disc")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        grid_layout.addLayout(vbox,rI,cI,1,1)
        rI, cI = [rI,cI+1] if cI < cmax else [rI+1,0]

        vbox = QVBoxLayout()
        annulus_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/annulus.svg"))
        annulus_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        annulus_button.setIconSize(size)
        annulus_button.setCheckable(True)
        annulus_button.setStyleSheet(qpushbutton_style_sheet) 
        self.cell_pos_button_group.addButton(annulus_button,next_button_id)
        next_button_id += 1
        vbox.addWidget(annulus_button)

        label = QLabel("Annulus")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        grid_layout.addLayout(vbox,rI,cI,1,1)
        rI, cI = [rI,cI+1] if cI < cmax else [rI+1,0]

        vbox = QVBoxLayout()
        wedge_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/wedge.svg"))
        wedge_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        wedge_button.setIconSize(size)
        wedge_button.setCheckable(True)
        wedge_button.setStyleSheet(qpushbutton_style_sheet) 
        self.cell_pos_button_group.addButton(wedge_button,next_button_id)
        next_button_id += 1
        vbox.addWidget(wedge_button)

        label = QLabel("Wedge")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        grid_layout.addLayout(vbox,rI,cI,1,1)
        rI, cI = [rI,cI+1] if cI < cmax else [rI+1,0]

        vbox = QVBoxLayout()
        rainbow_button = QPushButton(icon=QIcon(sys.path[0] + "/icon/rainbow.svg"),enabled=False) # not ready for this yet
        rainbow_button.setFixedSize(button_width,button_height)
        size = QtCore.QSize(icon_width, icon_height) 
        rainbow_button.setIconSize(size)
        rainbow_button.setCheckable(True)
        rainbow_button.setStyleSheet(qpushbutton_style_sheet) 
        self.cell_pos_button_group.addButton(rainbow_button,next_button_id)
        next_button_id += 1
        vbox.addWidget(rainbow_button)

        label = QLabel("Rainbow\n(why not?!)\nWork in progess...")
        label.setFixedWidth(button_width)
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(label)

        grid_layout.addLayout(vbox,rI,cI,1,1)
        rI, cI = [rI,cI+1] if cI < cmax else [rI+1,0]

        master_vbox.addLayout(grid_layout)

        pos_scroll_area_widget = QWidget()
        pos_scroll_area_widget.setLayout(master_vbox)

        self.pos_scroll_area = QScrollArea()
        self.pos_scroll_area.setWidget(pos_scroll_area_widget)

    def create_ics_plot_area(self):
        self.ics_plot_area = BioinfImportPlotWindow(self, self.config_tab)

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
            label_text += " \u21d2 "
            labels[final_type] = QLabel(label_text)
            self.new_name_line_edit[final_type] = QLineEdit()
            self.new_name_line_edit[final_type].setText(self.final_type_pre_image[final_type][0])
            hbox.addWidget(labels[final_type])
            hbox.addWidget(self.new_name_line_edit[final_type])
            vbox.addLayout(hbox)

        hbox = QHBoxLayout()
        go_back_button = GoBackButton(self)
        # go_back_button = QPushButton("\u2190 Go back")
        # go_back_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        # go_back_button.clicked.connect(self.go_back_to_prev_window)

        continue_button = QPushButton("Continue \u2192")
        continue_button.setStyleSheet("QPushButton {background-color: lightgreen; color: black;}")
        continue_button.clicked.connect(self.continue_on_rename_cb)

        hbox.addWidget(go_back_button)
        hbox.addWidget(continue_button)

        vbox.addLayout(hbox)

        self.previous_windows.append(self.window)
        self.previous_windows[-1].hide()

        self.window = BioinfImportWindow()
        self.window.setLayout(vbox)

        # hack to bring to foreground
        self.window.hide()
        self.window.show()

    def go_back_to_prev_window(self):
        self.window = self.previous_windows.pop()
        self.window.hide()
        self.window.show()

    def continue_on_rename_cb(self):
        # keep editing here
        self.cell_types_list_final = []
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
        self.checkbox_dict = create_checkboxes_for_cell_types(vbox, self.remaining_cell_types_list_original)

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
        cell_type = self.sender().objectName()
        # for cell_type in self.checkbox_dict.keys():
        # if self.checkbox_dict[cell_type].isChecked():
        self.set_cell_type_to_keep(cell_type)

    def set_cell_type_to_keep(self, cell_type, check_merge_gp=True):
        self.cell_type_dict[cell_type] = cell_type
        self.checkbox_dict[cell_type].setEnabled(True)
        self.checkbox_dict[cell_type].setStyleSheet(self.checkbox_style["keep"])
        if  check_merge_gp and ("\u21d2 Merge Gp. #" in self.checkbox_dict[cell_type].text()):
            # check for merge group that is no longer merging
            gp_preimage = [ctn for ctn, new_name in self.cell_type_dict.items() if (ctn != cell_type and new_name == self.cell_type_dict[cell_type])] # get cell types that map into this merge group
            if len(gp_preimage)==1: # then delete this merge gp:
                # gp_preimage = [ctn for ctn in gp_preimage if ctn != cell_type]
                self.set_cell_type_to_keep(gp_preimage[0], check_merge_gp=False)
            elif cell_type == self.cell_type_dict[cell_type]: # make sure this current cell type did not set the merge group name
                first_name = None
                for ctn in gp_preimage:
                    if first_name is None:
                        first_name = ctn
                    self.cell_type_dict[ctn] = first_name

            # new_name = self.cell_type_dict[cell_type]
            # for ctn,  in self.cell_type_dict.items():
            #     if self.cell_type_dict[ctn]
        self.cell_type_dict[cell_type] = cell_type
        self.checkbox_dict[cell_type].setText(f"{cell_type}")
        self.keep_button[cell_type].setEnabled(False)
        # self.remaining_cell_types_list_original = [cell_type for cell_type in self.remaining_cell_types_list_original if not self.checkbox_dict[cell_type].isChecked()]
        # if self.remaining_cell_types_list_original:
        #     self.edit_cell_types_cb()
        # else:
        #     self.continue_to_rename()

    def continue_on_delete_cb(self):
        for cell_type in self.checkbox_dict.keys():
            if self.checkbox_dict[cell_type].isChecked():
                self.cell_type_dict[cell_type] = None
                self.checkbox_dict[cell_type].setChecked(False)
                self.checkbox_dict[cell_type].setEnabled(False)
                self.checkbox_dict[cell_type].setStyleSheet(self.checkbox_style["delete"])

                self.keep_button[cell_type].setEnabled(True)
        # self.remaining_cell_types_list_original = [cell_type for cell_type in self.remaining_cell_types_list_original if not self.checkbox_dict[cell_type].isChecked()]
        # if self.remaining_cell_types_list_original:
        #     self.edit_cell_types_cb()
        # else:
        #     self.continue_to_rename()

    def continue_on_merge_cb(self):
        self.merge_id += 1
        first_name = None
        for cell_type in self.checkbox_dict.keys():
            if self.checkbox_dict[cell_type].isChecked():
                if first_name is None:
                    first_name = cell_type
                self.cell_type_dict[cell_type] = first_name
                self.checkbox_dict[cell_type].setChecked(False)
                self.checkbox_dict[cell_type].setEnabled(False)
                self.checkbox_dict[cell_type].setStyleSheet(self.checkbox_style["merge"])
                self.checkbox_dict[cell_type].setText(f"{cell_type} \u21d2 Merge Gp. #{self.merge_id}")

                self.keep_button[cell_type].setEnabled(True)
        # self.remaining_cell_types_list_original = [cell_type for cell_type in self.remaining_cell_types_list_original if not self.checkbox_dict[cell_type].isChecked()]
        # if self.remaining_cell_types_list_original:
        #     self.edit_cell_types_cb()
        # else:
        #     self.continue_to_rename()

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
    
    def setup_confluence_info(self):
        volume_env = (float(self.config_tab.xmax.text()) - float(self.config_tab.xmin.text())) * (float(self.config_tab.ymax.text()) - float(self.config_tab.ymin.text()))
        self.prop_total_area_one = {}
        self.prop_dot_ratios = 0
        for idx, cell_type in enumerate(self.cell_types_list_final):
            if cell_type in self.celldef_tab.param_d.keys():
                cell_volume = float(self.celldef_tab.param_d[cell_type]['volume_total'])
            else:
                cell_volume = 2494 # use PhysiCell default
            self.prop_total_area_one[cell_type] = (((9*np.pi*cell_volume**2) / 16) ** (1./3)) / volume_env
            self.prop_dot_ratios += (self.cell_type_props[idx] * self.prop_total_area_one[cell_type])
        pass

    def fill_gui(self):
        self.csv_folder.setText(self.config_tab.csv_folder.text())
        self.csv_file.setText(self.config_tab.csv_file.text())

def create_checkboxes_for_cell_types(vbox, cell_types):
    checkbox_dict = {}
    for cell_type in cell_types:
        checkbox_dict[cell_type] = QCheckBox_custom(cell_type)
        checkbox_dict[cell_type].setChecked(False)
        checkbox_dict[cell_type].setEnabled(True)
        vbox.addWidget(checkbox_dict[cell_type])

    return checkbox_dict

def normalize_thetas(th1,th2):
    th1 = th1 % 360
    th1 = th1 - 360 if th1 >= 180 else th1 # get th1 in [-180,180]
    th2 -= 360 * ((th2-th1) // 360) # get th2 in interval (th1,th1+360)
    return th1, th2

def compute_theta_intersection_distances(TH,D,th,x0,y0,xL,xR,yL,yR,dy):
    if dy == 0:
        ds = (np.array([xL,xR])-x0)/np.cos(th)
        if th > 0:
            ds = np.append(ds, (yR-y0)/np.sin(th))
        else:
            ds = np.append(ds, (yL-y0)/np.sin(th))
        d2 = np.array([(xL-x0)/np.cos(th),np.min(ds)]).reshape((2,1))
    else:
        d2 = np.max([(xL-x0)/np.cos(th),(yL-y0)/np.sin(th)])
        d2 = np.append(d2,np.min([(xR-x0)/np.cos(th),(yR-y0)/np.sin(th)]))
    TH = np.append(TH,th)
    D = np.append(D,d2,axis=1)
    return TH, D