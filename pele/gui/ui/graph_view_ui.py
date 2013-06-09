# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'graph_view_ui.ui'
#
# Created: Sun Jun  9 14:29:39 2013
#      by: PyQt4 UI code generator 4.9.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

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
        self.verticalLayout_2.addLayout(self.formLayout)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout_3.addLayout(self.horizontalLayout)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.btn_show_all.setText(_translate("Form", "show all", None))

from pele.gui.ui.mplwidget import MPLWidgetWithToolbar
