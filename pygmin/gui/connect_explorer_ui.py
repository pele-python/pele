# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'connect_explorer_ui.ui'
#
# Created: Thu Jan 24 10:24:34 2013
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(928, 718)
        self.horizontalLayout = QtGui.QHBoxLayout(Form)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setContentsMargins(10, -1, -1, -1)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.wgt_neb = NEBWidget(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.wgt_neb.sizePolicy().hasHeightForWidth())
        self.wgt_neb.setSizePolicy(sizePolicy)
        self.wgt_neb.setObjectName(_fromUtf8("wgt_neb"))
        self.verticalLayout_2.addWidget(self.wgt_neb)
        self.wgt_ogl_slider = Show3DWithSlider(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.wgt_ogl_slider.sizePolicy().hasHeightForWidth())
        self.wgt_ogl_slider.setSizePolicy(sizePolicy)
        self.wgt_ogl_slider.setObjectName(_fromUtf8("wgt_ogl_slider"))
        self.verticalLayout_2.addWidget(self.wgt_ogl_slider)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setContentsMargins(10, -1, -1, -1)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.formLayout = QtGui.QFormLayout()
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setContentsMargins(10, -1, -1, -1)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.btn_refineTS = QtGui.QPushButton(Form)
        self.btn_refineTS.setObjectName(_fromUtf8("btn_refineTS"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.btn_refineTS)
        self.btn_show_neb_path = QtGui.QPushButton(Form)
        self.btn_show_neb_path.setObjectName(_fromUtf8("btn_show_neb_path"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.btn_show_neb_path)
        self.btn_show_ts_path = QtGui.QPushButton(Form)
        self.btn_show_ts_path.setObjectName(_fromUtf8("btn_show_ts_path"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.btn_show_ts_path)
        self.btn_show_pushoff1 = QtGui.QPushButton(Form)
        self.btn_show_pushoff1.setObjectName(_fromUtf8("btn_show_pushoff1"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.btn_show_pushoff1)
        self.verticalLayout_3.addLayout(self.formLayout)
        self.horizontalLayout.addLayout(self.verticalLayout_3)

        self.retranslateUi(Form)
        QtCore.QObject.connect(self.btn_show_ts_path, QtCore.SIGNAL(_fromUtf8("clicked()")), Form.show_TS_path)
        QtCore.QObject.connect(self.btn_show_neb_path, QtCore.SIGNAL(_fromUtf8("clicked()")), Form.show_neb_path)
        QtCore.QObject.connect(self.btn_show_pushoff1, QtCore.SIGNAL(_fromUtf8("clicked()")), Form.show_pushoff_path)
        QtCore.QObject.connect(self.btn_refineTS, QtCore.SIGNAL(_fromUtf8("clicked()")), Form.refine_transition_state)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_refineTS.setText(QtGui.QApplication.translate("Form", "refine transition state", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_show_neb_path.setText(QtGui.QApplication.translate("Form", "show NEB path", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_show_ts_path.setText(QtGui.QApplication.translate("Form", "show TS refinement trajectory", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_show_pushoff1.setText(QtGui.QApplication.translate("Form", "show left pushoff path", None, QtGui.QApplication.UnicodeUTF8))

from nebdlg import NEBWidget
from show3d import Show3DWithSlider
