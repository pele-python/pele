# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'graph_view_ui.ui'
#
# Created: Mon Jun 10 12:16:02 2013
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
        Form.resize(699, 575)
        self.verticalLayout_3 = QtGui.QVBoxLayout(Form)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.label_status = QtGui.QLabel(Form)
        self.label_status.setText(_fromUtf8(""))
        self.label_status.setObjectName(_fromUtf8("label_status"))
        self.verticalLayout_3.addWidget(self.label_status)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.canvas = MPLWidgetWithToolbar(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.canvas.sizePolicy().hasHeightForWidth())
        self.canvas.setSizePolicy(sizePolicy)
        self.canvas.setObjectName(_fromUtf8("canvas"))
        self.horizontalLayout.addWidget(self.canvas)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setContentsMargins(1, -1, -1, -1)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.formLayout = QtGui.QFormLayout()
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.btn_show_all = QtGui.QPushButton(Form)
        self.btn_show_all.setObjectName(_fromUtf8("btn_show_all"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.btn_show_all)
        self.checkBox_zoom = QtGui.QCheckBox(Form)
        self.checkBox_zoom.setObjectName(_fromUtf8("checkBox_zoom"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.checkBox_zoom)
        self.verticalLayout_2.addLayout(self.formLayout)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout_3.addLayout(self.horizontalLayout)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_show_all.setText(QtGui.QApplication.translate("Form", "show all", None, QtGui.QApplication.UnicodeUTF8))
        self.checkBox_zoom.setText(QtGui.QApplication.translate("Form", "zoom on click", None, QtGui.QApplication.UnicodeUTF8))

from pele.gui.ui.mplwidget import MPLWidgetWithToolbar
