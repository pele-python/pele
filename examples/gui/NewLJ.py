# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'NewLJ.ui'
#
# Created: Thu May 10 03:10:06 2012
#      by: PyQt4 UI code generator 4.8.6
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
        DialogLJSetup.resize(349, 144)
        DialogLJSetup.setWindowTitle(QtGui.QApplication.translate("DialogLJSetup", "Create new Lennard-Jones system", None, QtGui.QApplication.UnicodeUTF8))
        DialogLJSetup.setModal(True)
        self.buttonBox = QtGui.QDialogButtonBox(DialogLJSetup)
        self.buttonBox.setGeometry(QtCore.QRect(20, 100, 301, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayoutWidget = QtGui.QWidget(DialogLJSetup)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(20, 20, 301, 61))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_2 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_2.setText(QtGui.QApplication.translate("DialogLJSetup", "Number of minima to save", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.lineNatoms = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineNatoms.setInputMask(_fromUtf8(""))
        self.lineNatoms.setText(QtGui.QApplication.translate("DialogLJSetup", "13", None, QtGui.QApplication.UnicodeUTF8))
        self.lineNatoms.setObjectName(_fromUtf8("lineNatoms"))
        self.gridLayout.addWidget(self.lineNatoms, 1, 1, 1, 1)
        self.lineNsave = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineNsave.setInputMask(QtGui.QApplication.translate("DialogLJSetup", "999; ", None, QtGui.QApplication.UnicodeUTF8))
        self.lineNsave.setText(QtGui.QApplication.translate("DialogLJSetup", "50", None, QtGui.QApplication.UnicodeUTF8))
        self.lineNsave.setObjectName(_fromUtf8("lineNsave"))
        self.gridLayout.addWidget(self.lineNsave, 2, 1, 1, 1)
        self.label = QtGui.QLabel(self.gridLayoutWidget)
        self.label.setText(QtGui.QApplication.translate("DialogLJSetup", "Number of particles", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)

        self.retranslateUi(DialogLJSetup)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), DialogLJSetup.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), DialogLJSetup.reject)
        QtCore.QMetaObject.connectSlotsByName(DialogLJSetup)

    def retranslateUi(self, DialogLJSetup):
        pass

