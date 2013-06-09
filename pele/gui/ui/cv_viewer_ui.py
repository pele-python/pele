# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cv_viewer_ui.ui'
#
# Created: Sun Jun  9 14:03:34 2013
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
        Form.resize(785, 580)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(Form)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_status = QtGui.QLabel(Form)
        self.label_status.setText(_fromUtf8(""))
        self.label_status.setObjectName(_fromUtf8("label_status"))
        self.verticalLayout.addWidget(self.label_status)
        self.splitter = QtGui.QSplitter(Form)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName(_fromUtf8("splitter"))
        self.mplwidget = MPLWidgetWithToolbar(self.splitter)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mplwidget.sizePolicy().hasHeightForWidth())
        self.mplwidget.setSizePolicy(sizePolicy)
        self.mplwidget.setObjectName(_fromUtf8("mplwidget"))
        self.widget = QtGui.QWidget(self.splitter)
        self.widget.setObjectName(_fromUtf8("widget"))
        self.formLayout = QtGui.QFormLayout(self.widget)
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setMargin(0)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.btn_recalculate = QtGui.QPushButton(self.widget)
        self.btn_recalculate.setObjectName(_fromUtf8("btn_recalculate"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.btn_recalculate)
        self.label = QtGui.QLabel(self.widget)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.label)
        self.lineEdit_nmin_max = QtGui.QLineEdit(self.widget)
        self.lineEdit_nmin_max.setObjectName(_fromUtf8("lineEdit_nmin_max"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.lineEdit_nmin_max)
        self.label_2 = QtGui.QLabel(self.widget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.label_2)
        self.label_3 = QtGui.QLabel(self.widget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.LabelRole, self.label_3)
        self.label_4 = QtGui.QLabel(self.widget)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.LabelRole, self.label_4)
        self.lineEdit_Tmin = QtGui.QLineEdit(self.widget)
        self.lineEdit_Tmin.setObjectName(_fromUtf8("lineEdit_Tmin"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.lineEdit_Tmin)
        self.lineEdit_Tmax = QtGui.QLineEdit(self.widget)
        self.lineEdit_Tmax.setObjectName(_fromUtf8("lineEdit_Tmax"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.FieldRole, self.lineEdit_Tmax)
        self.lineEdit_nT = QtGui.QLineEdit(self.widget)
        self.lineEdit_nT.setObjectName(_fromUtf8("lineEdit_nT"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.FieldRole, self.lineEdit_nT)
        self.verticalLayout.addWidget(self.splitter)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.btn_recalculate.setText(_translate("Form", "recalculate", None))
        self.label.setText(_translate("Form", "max # minima", None))
        self.label_2.setText(_translate("Form", "Tmin", None))
        self.label_3.setText(_translate("Form", "Tmax", None))
        self.label_4.setText(_translate("Form", "# points", None))

from pele.gui.ui.mplwidget import MPLWidgetWithToolbar
