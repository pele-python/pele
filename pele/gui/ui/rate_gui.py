# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'rate_gui.ui'
#
# Created: Fri Jan 17 16:24:05 2014
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
        Form.resize(618, 513)
        self.verticalLayout_2 = QtGui.QVBoxLayout(Form)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_status = QtGui.QLabel(Form)
        self.label_status.setText(_fromUtf8(""))
        self.label_status.setObjectName(_fromUtf8("label_status"))
        self.verticalLayout.addWidget(self.label_status)
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setContentsMargins(-1, 10, -1, -1)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.btn_compute = QtGui.QPushButton(Form)
        self.btn_compute.setObjectName(_fromUtf8("btn_compute"))
        self.gridLayout.addWidget(self.btn_compute, 3, 0, 1, 1)
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.label = QtGui.QLabel(Form)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.lineEdit_B = QtGui.QLineEdit(Form)
        self.lineEdit_B.setObjectName(_fromUtf8("lineEdit_B"))
        self.gridLayout.addWidget(self.lineEdit_B, 1, 1, 1, 1)
        self.lineEdit_A = QtGui.QLineEdit(Form)
        self.lineEdit_A.setObjectName(_fromUtf8("lineEdit_A"))
        self.gridLayout.addWidget(self.lineEdit_A, 0, 1, 1, 1)
        self.lineEdit_T = QtGui.QLineEdit(Form)
        self.lineEdit_T.setObjectName(_fromUtf8("lineEdit_T"))
        self.gridLayout.addWidget(self.lineEdit_T, 2, 1, 1, 1)
        self.label_3 = QtGui.QLabel(Form)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 2, 0, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.textBrowser = QtGui.QTextBrowser(Form)
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.verticalLayout.addWidget(self.textBrowser)
        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_compute.setText(QtGui.QApplication.translate("Form", "Compute Rates", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("Form", "Product State", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Form", "Reactant State", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("Form", "Temperature", None, QtGui.QApplication.UnicodeUTF8))

