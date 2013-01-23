# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'local_connect_browser.ui'
#
# Created: Wed Jan 23 09:49:55 2013
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
        Form.resize(635, 718)
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
        self.formLayout.setContentsMargins(10, -1, -1, -1)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.pushButton = QtGui.QPushButton(Form)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.pushButton)
        self.verticalLayout_3.addLayout(self.formLayout)
        self.horizontalLayout.addLayout(self.verticalLayout_3)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("Form", "PushButton", None, QtGui.QApplication.UnicodeUTF8))

from nebdlg import NEBWidget
from show3d import Show3DWithSlider
