#!/home/pcholine/anaconda3/bin/python

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

from functools import partial
import os
import sys

from filter_msp import filter_msp
from lipidDB import ConfigDB
from main_routine import main_routine_func

class Dummy(): pass

class ThreadSearch(QtCore.QThread):
    def __init__(self,ms1_error, ms2_error, ms_tol_ppm, plot_XIC, plot_ms2, aggregate_results, search_negative_ions, search_positive_ions, 
                    score_ms1_error,score_ms2_error,score_hypergeom,score_intensity_expl,score_pred_tr,min_intensity_explained,min_hypergeom,
                    rt_file,mzml_files_loc,class_fragments,fa_fragments,db_file,gui,parent = None):
        super(ThreadSearch,self).__init__(parent)
        self.ms1_error = ms1_error
        self.ms2_error = ms2_error
        self.ms_tol_ppm = ms_tol_ppm
        self.plot_XIC = plot_XIC
        self.plot_ms2 = plot_ms2
        self.aggregate_results = aggregate_results
        self.search_negative_ions = search_negative_ions
        self.search_positive_ions = search_positive_ions
        self.score_ms1_error = score_ms1_error
        self.score_ms2_error = score_ms2_error
        self.score_hypergeom = score_hypergeom
        self.score_intensity_expl = score_intensity_expl
        self.score_pred_tr = score_pred_tr
        self.min_intensity_explained = min_intensity_explained
        self.min_hypergeom = min_hypergeom
        self.rt_file = rt_file
        self.mzml_files_loc = mzml_files_loc
        self.class_fragments = class_fragments
        self.fa_fragments = fa_fragments
        self.db_file = db_file
        self.gui = gui
    
    def run(self):
        main_routine_func(ms1_error = self.ms1_error,
            ms2_error = self.ms2_error,
            ms_tol_ppm = self.ms_tol_ppm,
            plot_XIC = self.plot_XIC,
            plot_ms2 = self.plot_ms2,
            aggregate_results = self.aggregate_results,
            search_negative_ions = self.search_negative_ions,
            search_positive_ions = self.search_positive_ions,
            score_ms1_error = self.score_ms1_error,
            score_ms2_error = self.score_ms2_error,
            score_hypergeom = self.score_hypergeom,
            score_intensity_expl = self.score_intensity_expl,
            score_pred_tr = self.score_pred_tr,
            min_intensity_explained = self.min_intensity_explained,
            min_hypergeom = self.min_hypergeom,
            rt_file = self.rt_file,
            mzml_files_loc = self.mzml_files_loc,
            filter_head_spec=self.class_fragments,
            filter_fa_spec=self.fa_fragments,
            db_file=self.db_file,
            gui = self.gui)

class OutLog:
    def __init__(self, edit, out=None, color=None):
        """(edit, out=None, color=None) -> can write stdout, stderr to a
        QTextEdit.
        edit = QTextEdit
        out = alternate stream ( can be the original sys.stdout )
        color = alternate color (i.e. color stderr a different color)
        """
        self.edit = edit
        self.out = None
        self.color = color

    def write(self, m):
        if self.color:
            tc = self.edit.textColor()
            self.edit.setTextColor(self.color)

        self.edit.moveCursor(QtGui.QTextCursor.End)
        self.edit.insertPlainText( m )

        if self.color:
            self.edit.setTextColor(tc)

        if self.out:
            self.out.write(m)

class Ui_Dialog(object):
    def showErrorDialog(self,where_error):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)

        msg.setText("Please fill in all the required information")
        msg.setInformativeText(where_error)
        msg.setWindowTitle("Input error")
        #msg.setDetailedText("The details are as follows:")
        msg.setStandardButtons(QMessageBox.Ok) #| QMessageBox.Cancel
        #msg.buttonClicked.connect(msgbtn)

        retval = msg.exec_()
    
    def browse_for_file(self,set_obj,use_dir=False):
        if use_dir:  
            dlg = QFileDialog() #QFileDialog.getExistingDirectory(QFileDialog.Directory)
            dlg.setFileMode(QFileDialog.Directory)
            dlg.setOption(QFileDialog.ShowDirsOnly, True)
        else: dlg = QFileDialog()
        #dlg.setFileMode(QFileDialog.AnyFile)
        #dlg.setFilter(["Text files (*.txt)"])
        #filenames = QStringList()
        
        if dlg.exec_(): filenames = dlg.selectedFiles()
        try: set_obj.setText(filenames[0])
        except: pass
       
    def run_search(self):
        mzml = self.inputfield_mzml.displayText()
        tr_pred = self.inputfield_tr_pred.displayText()
        native_lipids = self.inputfield_native_lipids.displayText()
        output_dir = self.inputfield_output_dir.displayText()
        db_file = self.inputfield_db_file.displayText()
        ms1_tolerance = float(self.inputfield_ms1_tolerance.displayText())
        ms2_tolerance = float(self.inputfield_ms2_tolerance.displayText())
        threshold_intensity = float(self.inputfield_threshold_intensity.displayText())
        threshold_hypergeo_score = float(self.inputfield_threshold_hypergeo_score.displayText())

        if self.radio_ppm.isChecked() == False: 
            ppm = False
            dalton = True
        else:
            ppm = True
            dalton = False
        
        if self.check_negative.checkState() == 0: negative = False
        else: negative = True
        if self.check_positive.checkState() == 0: positive = False
        else: positive = True
        if self.check_ms1_tolerance.checkState() == 0: score_ms1_tolerance = False
        else: score_ms1_tolerance = True
        if self.check_ms2_tolerance.checkState() == 0: score_ms2_tolerance = False
        else: score_ms2_tolerance = True
        if self.check_hypergeo_score.checkState() == 0: hypergeo_score = False
        else: hypergeo_score = True
        if self.check_intensity.checkState() == 0: intensity = False
        else: intensity = True
        if self.check_tr_pred.checkState() == 0: score_tr_pred = False
        else: score_tr_pred = True
        if self.check_ox_pred.checkState() == 0: ox_pred = False
        else: ox_pred = True
        if self.check_ms2_plot.checkState() == 0: ms2_plot = False
        else: ms2_plot = True
        if self.check_xic_plot.checkState() == 0: xic_plot = False
        else: xic_plot = True
        if self.check_class_fragments.checkState() == 0: class_fragments = False
        else: class_fragments = True
        if self.check_fa_fragments.checkState() == 0: fa_fragments = False
        else: fa_fragments = True
        if self.check_aggregate_results.checkState() == 0: aggregate_results = False
        else: aggregate_results = True
        
        self.threadsearch = ThreadSearch(ms1_tolerance, ms2_tolerance, ppm, xic_plot, ms2_plot, aggregate_results,
                                            negative, positive, score_ms1_tolerance, score_ms2_tolerance, hypergeo_score,
                                            intensity, score_tr_pred, threshold_intensity, threshold_hypergeo_score,
                                            tr_pred,mzml,class_fragments,fa_fragments,db_file,self)
        self.threadsearch.start()
        
    
    def run_db_creation(self):
        input_config_file = self.inputfield_input_config.displayText()
        output_file = self.inputfield_output_db.displayText()
        
        db_in = ConfigDB(config_file=input_config_file)
        db_in.create_db(outfile_name=output_file)
        
        input_config_file = self.inputfield_input_config.displayText()
        output_file = self.inputfield_output_db.displayText()
               
        max_db = int(self.inputfield_max_db.displayText())
        max_c = int(self.inputfield_max_c.displayText())
        min_c = int(self.inputfield_min_c.displayText())
        max_mods = int(self.inputfield_max_sum_mods.displayText())
        
        if self.check_even_c.checkState() == 0: even = False
        else: even = True
        
        if self.check_ignore_mods.checkState() == 0: ignore_mods = False
        else: ignore_mods = True
        
        filter_class = []
        
        filter_msp(output_file,output_file,output_file+"_toformula.csv",max_db=max_db,max_c=max_c,min_c=min_c,max_mods=max_mods,even=even,ignore_mods=ignore_mods,filter_class=False)
    
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(438, 650)
        self.tabWidget = QtWidgets.QTabWidget(Dialog)
        self.tabWidget.setGeometry(QtCore.QRect(0, 0, 441, 651))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayoutWidget_5 = QtWidgets.QWidget(self.tab)
        self.gridLayoutWidget_5.setGeometry(QtCore.QRect(20, 20, 401, 236))
        self.gridLayoutWidget_5.setObjectName("gridLayoutWidget_5")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.gridLayoutWidget_5)
        self.gridLayout_5.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.label_16 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_16.setObjectName("label_16")
        self.gridLayout_5.addWidget(self.label_16, 10, 0, 1, 2)
        self.browse_mzml = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.browse_mzml.setObjectName("browse_mzml")
        self.gridLayout_5.addWidget(self.browse_mzml, 1, 0, 1, 1)
        self.browse_native_lipids = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.browse_native_lipids.setObjectName("browse_native_lipids")
        self.gridLayout_5.addWidget(self.browse_native_lipids, 9, 0, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_7.setObjectName("label_7")
        self.gridLayout_5.addWidget(self.label_7, 0, 0, 1, 2)
        self.browse_tr_pred = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.browse_tr_pred.setObjectName("browse_tr_pred")
        self.gridLayout_5.addWidget(self.browse_tr_pred, 5, 0, 1, 1)
        self.inputfield_native_lipids = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.inputfield_native_lipids.setObjectName("inputfield_native_lipids")
        self.gridLayout_5.addWidget(self.inputfield_native_lipids, 9, 1, 1, 1)
        self.inputfield_mzml = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.inputfield_mzml.setObjectName("inputfield_mzml")
        self.gridLayout_5.addWidget(self.inputfield_mzml, 1, 1, 1, 1)
        self.inputfield_tr_pred = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.inputfield_tr_pred.setObjectName("inputfield_tr_pred")
        self.gridLayout_5.addWidget(self.inputfield_tr_pred, 5, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_6.setObjectName("label_6")
        self.gridLayout_5.addWidget(self.label_6, 7, 0, 1, 2)
        self.label_10 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_10.setObjectName("label_10")
        self.gridLayout_5.addWidget(self.label_10, 4, 0, 1, 2)
        self.browse_output_dir = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.browse_output_dir.setObjectName("browse_output_dir")
        self.gridLayout_5.addWidget(self.browse_output_dir, 11, 0, 1, 1)
        self.inputfield_output_dir = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.inputfield_output_dir.setObjectName("inputfield_output_dir")
        self.gridLayout_5.addWidget(self.inputfield_output_dir, 11, 1, 1, 1)
        self.browse_db_file = QtWidgets.QPushButton(self.gridLayoutWidget_5)
        self.browse_db_file.setObjectName("browse_db_file")
        self.gridLayout_5.addWidget(self.browse_db_file, 3, 0, 1, 1)
        self.inputfield_db_file = QtWidgets.QLineEdit(self.gridLayoutWidget_5)
        self.inputfield_db_file.setObjectName("inputfield_db_file")
        self.gridLayout_5.addWidget(self.inputfield_db_file, 3, 1, 1, 1)
        self.label_20 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_20.setObjectName("label_20")
        self.gridLayout_5.addWidget(self.label_20, 2, 0, 1, 2)
        self.gridLayoutWidget_7 = QtWidgets.QWidget(self.tab)
        self.gridLayoutWidget_7.setGeometry(QtCore.QRect(20, 460, 401, 117))
        self.gridLayoutWidget_7.setObjectName("gridLayoutWidget_7")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.gridLayoutWidget_7)
        self.gridLayout_7.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.check_xic_plot = QtWidgets.QCheckBox(self.gridLayoutWidget_7)
        self.check_xic_plot.setObjectName("check_xic_plot")
        self.gridLayout_7.addWidget(self.check_xic_plot, 0, 1, 1, 1)
        self.check_ms2_plot = QtWidgets.QCheckBox(self.gridLayoutWidget_7)
        self.check_ms2_plot.setObjectName("check_ms2_plot")
        self.gridLayout_7.addWidget(self.check_ms2_plot, 0, 0, 1, 1)
        self.check_aggregate_results = QtWidgets.QCheckBox(self.gridLayoutWidget_7)
        self.check_aggregate_results.setObjectName("check_aggregate_results")
        self.gridLayout_7.addWidget(self.check_aggregate_results, 2, 0, 1, 2)
        self.inputfield_threshold_intensity = QtWidgets.QLineEdit(self.gridLayoutWidget_7)
        self.inputfield_threshold_intensity.setObjectName("inputfield_threshold_intensity")
        self.gridLayout_7.addWidget(self.inputfield_threshold_intensity, 4, 1, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.label_14.setObjectName("label_14")
        self.gridLayout_7.addWidget(self.label_14, 4, 0, 1, 1)
        self.inputfield_threshold_hypergeo_score = QtWidgets.QLineEdit(self.gridLayoutWidget_7)
        self.inputfield_threshold_hypergeo_score.setObjectName("inputfield_threshold_hypergeo_score")
        self.gridLayout_7.addWidget(self.inputfield_threshold_hypergeo_score, 5, 1, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.label_15.setObjectName("label_15")
        self.gridLayout_7.addWidget(self.label_15, 5, 0, 1, 1)
        self.check_class_fragments = QtWidgets.QCheckBox(self.gridLayoutWidget_7)
        self.check_class_fragments.setObjectName("check_class_fragments")
        self.gridLayout_7.addWidget(self.check_class_fragments, 1, 0, 1, 1)
        self.check_fa_fragments = QtWidgets.QCheckBox(self.gridLayoutWidget_7)
        self.check_fa_fragments.setObjectName("check_fa_fragments")
        self.gridLayout_7.addWidget(self.check_fa_fragments, 1, 1, 1, 1)
        self.gridLayoutWidget_9 = QtWidgets.QWidget(self.tab)
        self.gridLayoutWidget_9.setGeometry(QtCore.QRect(20, 360, 401, 74))
        self.gridLayoutWidget_9.setObjectName("gridLayoutWidget_9")
        self.gridLayout_9 = QtWidgets.QGridLayout(self.gridLayoutWidget_9)
        self.gridLayout_9.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.check_hypergeo_score = QtWidgets.QCheckBox(self.gridLayoutWidget_9)
        self.check_hypergeo_score.setObjectName("check_hypergeo_score")
        self.gridLayout_9.addWidget(self.check_hypergeo_score, 1, 0, 1, 1)
        self.check_tr_pred = QtWidgets.QCheckBox(self.gridLayoutWidget_9)
        self.check_tr_pred.setObjectName("check_tr_pred")
        self.gridLayout_9.addWidget(self.check_tr_pred, 2, 0, 1, 1)
        self.check_intensity = QtWidgets.QCheckBox(self.gridLayoutWidget_9)
        self.check_intensity.setObjectName("check_intensity")
        self.gridLayout_9.addWidget(self.check_intensity, 1, 1, 1, 1)
        self.check_ox_pred = QtWidgets.QCheckBox(self.gridLayoutWidget_9)
        self.check_ox_pred.setObjectName("check_ox_pred")
        self.gridLayout_9.addWidget(self.check_ox_pred, 2, 1, 1, 1)
        self.check_ms1_tolerance = QtWidgets.QCheckBox(self.gridLayoutWidget_9)
        self.check_ms1_tolerance.setObjectName("check_ms1_tolerance")
        self.gridLayout_9.addWidget(self.check_ms1_tolerance, 0, 0, 1, 1)
        self.check_ms2_tolerance = QtWidgets.QCheckBox(self.gridLayoutWidget_9)
        self.check_ms2_tolerance.setObjectName("check_ms2_tolerance")
        self.gridLayout_9.addWidget(self.check_ms2_tolerance, 0, 1, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.tab)
        self.label_11.setGeometry(QtCore.QRect(20, 340, 161, 16))
        self.label_11.setObjectName("label_11")
        self.label_13 = QtWidgets.QLabel(self.tab)
        self.label_13.setGeometry(QtCore.QRect(20, 440, 161, 16))
        self.label_13.setObjectName("label_13")
        self.label_12 = QtWidgets.QLabel(self.tab)
        self.label_12.setGeometry(QtCore.QRect(20, 0, 161, 16))
        self.label_12.setObjectName("label_12")
        self.line = QtWidgets.QFrame(self.tab)
        self.line.setGeometry(QtCore.QRect(20, 350, 401, 16))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.line_2 = QtWidgets.QFrame(self.tab)
        self.line_2.setGeometry(QtCore.QRect(20, 450, 401, 20))
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.line_3 = QtWidgets.QFrame(self.tab)
        self.line_3.setGeometry(QtCore.QRect(20, 10, 401, 20))
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.run_button_search = QtWidgets.QPushButton(self.tab)
        self.run_button_search.setGeometry(QtCore.QRect(20, 600, 401, 23))
        self.run_button_search.setObjectName("run_button_search")
        self.progress_search = QtWidgets.QProgressBar(self.tab)
        self.progress_search.setGeometry(QtCore.QRect(20, 580, 401, 16))
        self.progress_search.setProperty("value", 0)
        self.progress_search.setObjectName("progress_search")
        self.gridLayoutWidget_6 = QtWidgets.QWidget(self.tab)
        self.gridLayoutWidget_6.setGeometry(QtCore.QRect(20, 260, 401, 77))
        self.gridLayoutWidget_6.setObjectName("gridLayoutWidget_6")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.gridLayoutWidget_6)
        self.gridLayout_6.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.label_9 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_9.setObjectName("label_9")
        self.gridLayout_6.addWidget(self.label_9, 6, 0, 2, 1)
        self.label_8 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_8.setObjectName("label_8")
        self.gridLayout_6.addWidget(self.label_8, 1, 0, 2, 1)
        self.inputfield_ms1_tolerance = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.inputfield_ms1_tolerance.setObjectName("inputfield_ms1_tolerance")
        self.gridLayout_6.addWidget(self.inputfield_ms1_tolerance, 1, 1, 2, 1)
        self.inputfield_ms2_tolerance = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.inputfield_ms2_tolerance.setObjectName("inputfield_ms2_tolerance")
        self.gridLayout_6.addWidget(self.inputfield_ms2_tolerance, 6, 1, 2, 1)
        self.radio_dalton = QtWidgets.QRadioButton(self.gridLayoutWidget_6)
        self.radio_dalton.setObjectName("radio_dalton")
        self.gridLayout_6.addWidget(self.radio_dalton, 6, 2, 1, 1)
        self.radio_ppm = QtWidgets.QRadioButton(self.gridLayoutWidget_6)
        self.radio_ppm.setObjectName("radio_ppm")
        self.gridLayout_6.addWidget(self.radio_ppm, 2, 2, 1, 1)
        self.check_positive = QtWidgets.QCheckBox(self.gridLayoutWidget_6)
        self.check_positive.setObjectName("check_positive")
        self.gridLayout_6.addWidget(self.check_positive, 8, 1, 1, 1)
        self.check_negative = QtWidgets.QCheckBox(self.gridLayoutWidget_6)
        self.check_negative.setObjectName("check_negative")
        self.gridLayout_6.addWidget(self.check_negative, 8, 0, 1, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.gridLayoutWidget = QtWidgets.QWidget(self.tab_2)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(20, 20, 401, 118))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.inputfield_input_config = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.inputfield_input_config.setObjectName("inputfield_input_config")
        self.gridLayout.addWidget(self.inputfield_input_config, 1, 1, 1, 1)
        self.browse_input_config = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.browse_input_config.setObjectName("browse_input_config")
        self.gridLayout.addWidget(self.browse_input_config, 1, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 2)
        self.label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.browse_output_db = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.browse_output_db.setObjectName("browse_output_db")
        self.gridLayout.addWidget(self.browse_output_db, 4, 0, 1, 1)
        self.inputfield_output_db = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.inputfield_output_db.setObjectName("inputfield_output_db")
        self.gridLayout.addWidget(self.inputfield_output_db, 4, 1, 1, 1)
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.tab_2)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(20, 170, 401, 88))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.check_PA = QtWidgets.QCheckBox(self.gridLayoutWidget_2)
        self.check_PA.setObjectName("check_PA")
        self.gridLayout_2.addWidget(self.check_PA, 0, 0, 1, 1)
        self.check_PG = QtWidgets.QCheckBox(self.gridLayoutWidget_2)
        self.check_PG.setObjectName("check_PG")
        self.gridLayout_2.addWidget(self.check_PG, 1, 0, 1, 1)
        self.check_PC = QtWidgets.QCheckBox(self.gridLayoutWidget_2)
        self.check_PC.setObjectName("check_PC")
        self.gridLayout_2.addWidget(self.check_PC, 0, 1, 1, 1)
        self.check_PS = QtWidgets.QCheckBox(self.gridLayoutWidget_2)
        self.check_PS.setObjectName("check_PS")
        self.gridLayout_2.addWidget(self.check_PS, 1, 1, 1, 1)
        self.check_no_filter = QtWidgets.QCheckBox(self.gridLayoutWidget_2)
        self.check_no_filter.setChecked(True)
        self.check_no_filter.setObjectName("check_no_filter")
        self.gridLayout_2.addWidget(self.check_no_filter, 3, 1, 1, 1)
        self.check_PI = QtWidgets.QCheckBox(self.gridLayoutWidget_2)
        self.check_PI.setObjectName("check_PI")
        self.gridLayout_2.addWidget(self.check_PI, 2, 1, 1, 1)
        self.check_PE = QtWidgets.QCheckBox(self.gridLayoutWidget_2)
        self.check_PE.setObjectName("check_PE")
        self.gridLayout_2.addWidget(self.check_PE, 2, 0, 1, 1)
        self.check_TG = QtWidgets.QCheckBox(self.gridLayoutWidget_2)
        self.check_TG.setObjectName("check_TG")
        self.gridLayout_2.addWidget(self.check_TG, 3, 0, 1, 1)
        self.gridLayoutWidget_3 = QtWidgets.QWidget(self.tab_2)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(20, 270, 401, 100))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.gridLayout_3.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_3 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_3.setObjectName("label_3")
        self.gridLayout_3.addWidget(self.label_3, 0, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_4.setObjectName("label_4")
        self.gridLayout_3.addWidget(self.label_4, 1, 0, 1, 1)
        self.inputfield_min_c = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.inputfield_min_c.setObjectName("inputfield_min_c")
        self.gridLayout_3.addWidget(self.inputfield_min_c, 0, 1, 1, 1)
        self.inputfield_max_c = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.inputfield_max_c.setObjectName("inputfield_max_c")
        self.gridLayout_3.addWidget(self.inputfield_max_c, 1, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_5.setObjectName("label_5")
        self.gridLayout_3.addWidget(self.label_5, 3, 0, 1, 1)
        self.inputfield_max_sum_mods = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.inputfield_max_sum_mods.setObjectName("inputfield_max_sum_mods")
        self.gridLayout_3.addWidget(self.inputfield_max_sum_mods, 3, 1, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_17.setObjectName("label_17")
        self.gridLayout_3.addWidget(self.label_17, 2, 0, 1, 1)
        self.inputfield_max_db = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.inputfield_max_db.setObjectName("inputfield_max_db")
        self.gridLayout_3.addWidget(self.inputfield_max_db, 2, 1, 1, 1)
        self.gridLayoutWidget_4 = QtWidgets.QWidget(self.tab_2)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(20, 380, 401, 80))
        self.gridLayoutWidget_4.setObjectName("gridLayoutWidget_4")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.gridLayoutWidget_4)
        self.gridLayout_4.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.check_even_c = QtWidgets.QCheckBox(self.gridLayoutWidget_4)
        self.check_even_c.setObjectName("check_even_c")
        self.gridLayout_4.addWidget(self.check_even_c, 0, 0, 1, 1)
        self.check_ignore_mods = QtWidgets.QCheckBox(self.gridLayoutWidget_4)
        self.check_ignore_mods.setChecked(True)
        self.check_ignore_mods.setObjectName("check_ignore_mods")
        self.gridLayout_4.addWidget(self.check_ignore_mods, 1, 0, 1, 1)
        self.run_button_db = QtWidgets.QPushButton(self.tab_2)
        self.run_button_db.setGeometry(QtCore.QRect(20, 580, 401, 23))
        self.run_button_db.setObjectName("run_button_db")
        self.line_4 = QtWidgets.QFrame(self.tab_2)
        self.line_4.setGeometry(QtCore.QRect(20, 10, 401, 20))
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.label_18 = QtWidgets.QLabel(self.tab_2)
        self.label_18.setGeometry(QtCore.QRect(20, 0, 161, 16))
        self.label_18.setObjectName("label_18")
        self.label_19 = QtWidgets.QLabel(self.tab_2)
        self.label_19.setGeometry(QtCore.QRect(20, 150, 161, 16))
        self.label_19.setObjectName("label_19")
        self.line_5 = QtWidgets.QFrame(self.tab_2)
        self.line_5.setGeometry(QtCore.QRect(20, 160, 401, 20))
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.outputfield_text = QtWidgets.QTextBrowser(self.tab_2)
        self.outputfield_text.setGeometry(QtCore.QRect(20, 470, 401, 101))
        self.outputfield_text.setAcceptDrops(False)
        self.outputfield_text.setAcceptRichText(False)
        self.outputfield_text.setObjectName("outputfield_text")
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.tabWidget.addTab(self.tab_3, "")
        
        sys.stdout = OutLog(self.outputfield_text, sys.stdout)
        
        self.run_button_search.clicked.connect(self.run_search)
        self.run_button_db.clicked.connect(self.run_db_creation)
        self.browse_input_config.clicked.connect(partial(self.browse_for_file,self.inputfield_input_config))
        self.browse_db_file.clicked.connect(partial(self.browse_for_file,self.inputfield_db_file))
        self.browse_output_db.clicked.connect(partial(self.browse_for_file,self.inputfield_output_db))
        self.browse_mzml.clicked.connect(partial(self.browse_for_file,self.inputfield_mzml))
        self.browse_tr_pred.clicked.connect(partial(self.browse_for_file,self.browse_tr_pred))
        self.browse_native_lipids.clicked.connect(partial(self.browse_for_file,self.browse_native_lipids))
        self.browse_output_dir.clicked.connect(partial(self.browse_for_file,self.inputfield_output_dir,use_dir=True))
                
        self.retranslateUi(Dialog)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("LIT search engine", "LIT search engine"))
        self.label_16.setText(_translate("Dialog", "Output directory"))
        self.browse_mzml.setText(_translate("Dialog", "Browse"))
        self.browse_native_lipids.setText(_translate("Dialog", "Browse"))
        self.label_7.setText(_translate("Dialog", "Input mzml/mgf file or directory"))
        self.browse_tr_pred.setText(_translate("Dialog", "Browse"))
        self.label_6.setText(_translate("Dialog", "Input native lipid identifications"))
        self.label_10.setText(_translate("Dialog", "Input retention time prediction file"))
        self.browse_output_dir.setText(_translate("Dialog", "Browse"))
        self.browse_db_file.setText(_translate("Dialog", "Browse"))
        self.label_20.setText(_translate("Dialog", "Input database file"))
        self.check_xic_plot.setText(_translate("Dialog", "Plot XIC"))
        self.check_ms2_plot.setText(_translate("Dialog", "Plot MS2 alignments"))
        self.check_aggregate_results.setText(_translate("Dialog", "Aggregate results from multiple files into one"))
        self.label_14.setText(_translate("Dialog", "Intensity explained threshold"))
        self.label_15.setText(_translate("Dialog", "Hypergeometric score threshold (-log10(p-value))"))
        self.check_class_fragments.setText(_translate("Dialog", "Filter class fragment peaks"))
        self.check_fa_fragments.setText(_translate("Dialog", "Filter FA fragment peaks"))
        self.check_hypergeo_score.setText(_translate("Dialog", "Hypergeometric score"))
        self.check_tr_pred.setText(_translate("Dialog", "Predicted retention time"))
        self.check_intensity.setText(_translate("Dialog", "Intensity explained (%)"))
        self.check_ox_pred.setText(_translate("Dialog", "Predicted oxidative modification"))
        self.check_ms1_tolerance.setText(_translate("Dialog", "MS1 mass error"))
        self.check_ms2_tolerance.setText(_translate("Dialog", "MS2 mass error"))
        self.label_11.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:10pt; font-weight:600;\">Scoring parameters</span></p></body></html>"))
        self.label_13.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:10pt; font-weight:600;\">Output parameters</span></p></body></html>"))
        self.label_12.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:10pt; font-weight:600;\">Search parameters</span></p></body></html>"))
        self.run_button_search.setText(_translate("Dialog", "Run"))
        self.label_9.setText(_translate("Dialog", "MS2 mass tolerance                                  "))
        self.label_8.setText(_translate("Dialog", "MS1 mass tolerance                                  "))
        self.radio_dalton.setText(_translate("Dialog", "Dalton"))
        self.radio_ppm.setText(_translate("Dialog", "PPM"))
        self.check_positive.setText(_translate("Dialog", "Positive ions"))
        self.check_negative.setText(_translate("Dialog", "Negative ions"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("Dialog", "Search engine"))
        self.browse_input_config.setText(_translate("Dialog", "Browse"))
        self.label_2.setText(_translate("Dialog", "Output file"))
        self.label.setText(_translate("Dialog", "Input database configuration file"))
        self.browse_output_db.setText(_translate("Dialog", "Browse"))
        self.check_PA.setText(_translate("Dialog", "PA"))
        self.check_PG.setText(_translate("Dialog", "PG"))
        self.check_PC.setText(_translate("Dialog", "PC"))
        self.check_PS.setText(_translate("Dialog", "PS"))
        self.check_no_filter.setText(_translate("Dialog", "No filtering"))
        self.check_PI.setText(_translate("Dialog", "PI"))
        self.check_PE.setText(_translate("Dialog", "PE"))
        self.check_TG.setText(_translate("Dialog", "TG"))
        self.label_3.setText(_translate("Dialog", "Minimal length in carbon atoms of FA"))
        self.label_4.setText(_translate("Dialog", "Maximum length in carbon atoms of FA"))
        self.label_5.setText(_translate("Dialog", "Maximum amount of (summed) modifications"))
        self.label_17.setText(_translate("Dialog", "Maximal double bonds in carbon atoms of FA"))
        self.check_even_c.setText(_translate("Dialog", "Even number of carbon atoms in FA"))
        self.check_ignore_mods.setText(_translate("Dialog", "Ignore previous stated rules if lipid contains modification"))
        self.run_button_db.setText(_translate("Dialog", "Run"))
        self.label_18.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:10pt; font-weight:600;\">Input and output</span></p></body></html>"))
        self.label_19.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:10pt; font-weight:600;\">Filter settings database</span></p></body></html>"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("Dialog", "Create DB"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("Dialog", "About"))

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
