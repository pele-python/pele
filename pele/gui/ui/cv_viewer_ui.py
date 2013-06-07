# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cv_viewer_ui.ui'
#
# Created: Fri Jun  7 16:36:48 2013
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
        Form.resize(675, 522)
        self.verticalLayoutWidget = QtGui.QWidget(Form)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(30, 20, 631, 491))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_status = QtGui.QLabel(self.verticalLayoutWidget)
        self.label_status.setText(_fromUtf8(""))
        self.label_status.setObjectName(_fromUtf8("label_status"))
        self.verticalLayout.addWidget(self.label_status)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.mplwidget = MPLWidgetWithToolbar(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mplwidget.sizePolicy().hasHeightForWidth())
        self.mplwidget.setSizePolicy(sizePolicy)
        self.mplwidget.setObjectName(_fromUtf8("mplwidget"))
        self.horizontalLayout.addWidget(self.mplwidget)
        self.formLayout = QtGui.QFormLayout()
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.btn_recalculate = QtGui.QPushButton(self.verticalLayoutWidget)
        self.btn_recalculate.setObjectName(_fromUtf8("btn_recalculate"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.btn_recalculate)
        self.label = QtGui.QLabel(self.verticalLayoutWidget)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.label)
        self.lineEdit_nmin_max = QtGui.QLineEdit(self.verticalLayoutWidget)
        self.lineEdit_nmin_max.setObjectName(_fromUtf8("lineEdit_nmin_max"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.lineEdit_nmin_max)
        self.label_2 = QtGui.QLabel(self.verticalLayoutWidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.label_2)
        self.label_3 = QtGui.QLabel(self.verticalLayoutWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.LabelRole, self.label_3)
        self.label_4 = QtGui.QLabel(self.verticalLayoutWidget)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.LabelRole, self.label_4)
        self.lineEdit_Tmin = QtGui.QLineEdit(self.verticalLayoutWidget)
        self.lineEdit_Tmin.setObjectName(_fromUtf8("lineEdit_Tmin"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.lineEdit_Tmin)
        self.lineEdit_Tmax = QtGui.QLineEdit(self.verticalLayoutWidget)
        self.lineEdit_Tmax.setObjectName(_fromUtf8("lineEdit_Tmax"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.FieldRole, self.lineEdit_Tmax)
        self.lineEdit_nT = QtGui.QLineEdit(self.verticalLayoutWidget)
        self.lineEdit_nT.setObjectName(_fromUtf8("lineEdit_nT"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.FieldRole, self.lineEdit_nT)
        self.horizontalLayout.addLayout(self.formLayout)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_recalculate.setText(QtGui.QApplication.translate("Form", "recalculate", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Form", "max # minima", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("Form", "Tmin", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("Form", "Tmax", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("Form", "# points", None, QtGui.QApplication.UnicodeUTF8))

from pele.gui.ui.mplwidget import MPLWidgetWithToolbar
