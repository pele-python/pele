# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '_blj_dialog.ui'
#
# Created: Fri Nov  8 12:16:44 2013
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_DialogLJSetup(object):
    def setupUi(self, DialogLJSetup):
        DialogLJSetup.setObjectName(_fromUtf8("DialogLJSetup"))
        DialogLJSetup.resize(351, 348)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DialogLJSetup.sizePolicy().hasHeightForWidth())
        DialogLJSetup.setSizePolicy(sizePolicy)
        DialogLJSetup.setModal(True)
        self.buttonBox = QtGui.QDialogButtonBox(DialogLJSetup)
        self.buttonBox.setGeometry(QtCore.QRect(30, 280, 301, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayoutWidget = QtGui.QWidget(DialogLJSetup)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(20, 20, 301, 241))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.lineEdit_sigBB = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_sigBB.setObjectName(_fromUtf8("lineEdit_sigBB"))
        self.gridLayout.addWidget(self.lineEdit_sigBB, 5, 1, 1, 1)
        self.label_2 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.lineEdit_natoms = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_natoms.setInputMask(_fromUtf8(""))
        self.lineEdit_natoms.setObjectName(_fromUtf8("lineEdit_natoms"))
        self.gridLayout.addWidget(self.lineEdit_natoms, 1, 1, 1, 1)
        self.lineEdit_ntypeA = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_ntypeA.setObjectName(_fromUtf8("lineEdit_ntypeA"))
        self.gridLayout.addWidget(self.lineEdit_ntypeA, 2, 1, 1, 1)
        self.label = QtGui.QLabel(self.gridLayoutWidget)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.lineEdit_epsAB = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_epsAB.setObjectName(_fromUtf8("lineEdit_epsAB"))
        self.gridLayout.addWidget(self.lineEdit_epsAB, 4, 1, 1, 1)
        self.lineEdit_sigAB = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_sigAB.setObjectName(_fromUtf8("lineEdit_sigAB"))
        self.gridLayout.addWidget(self.lineEdit_sigAB, 3, 1, 1, 1)
        self.lineEdit_epsBB = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_epsBB.setObjectName(_fromUtf8("lineEdit_epsBB"))
        self.gridLayout.addWidget(self.lineEdit_epsBB, 6, 1, 1, 1)
        self.label_3 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 3, 0, 1, 1)
        self.label_4 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 4, 0, 1, 1)
        self.label_5 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout.addWidget(self.label_5, 5, 0, 1, 1)
        self.label_6 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout.addWidget(self.label_6, 6, 0, 1, 1)

        self.retranslateUi(DialogLJSetup)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), DialogLJSetup.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), DialogLJSetup.reject)
        QtCore.QMetaObject.connectSlotsByName(DialogLJSetup)

    def retranslateUi(self, DialogLJSetup):
        DialogLJSetup.setWindowTitle(QtGui.QApplication.translate("DialogLJSetup", "Create binary Lennard-Jones system", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_sigBB.setText(QtGui.QApplication.translate("DialogLJSetup", "0.88", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("DialogLJSetup", "number type A", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_natoms.setText(QtGui.QApplication.translate("DialogLJSetup", "13", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_ntypeA.setInputMask(QtGui.QApplication.translate("DialogLJSetup", "999; ", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_ntypeA.setText(QtGui.QApplication.translate("DialogLJSetup", "10", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("DialogLJSetup", "Number of particles", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_epsAB.setText(QtGui.QApplication.translate("DialogLJSetup", "1.5", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_sigAB.setText(QtGui.QApplication.translate("DialogLJSetup", "0.8", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_epsBB.setText(QtGui.QApplication.translate("DialogLJSetup", "0.5", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("DialogLJSetup", "sigma AB", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("DialogLJSetup", "eps AB", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("DialogLJSetup", "sigma BB", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("DialogLJSetup", "eps BB", None, QtGui.QApplication.UnicodeUTF8))

