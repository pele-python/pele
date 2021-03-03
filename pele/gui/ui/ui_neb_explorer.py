# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui/ui_neb_explorer.ui'
#
# Created: Thu Feb  7 00:13:53 2013
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from __future__ import absolute_import
from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(800, 600)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 29))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.toolBar = QtGui.QToolBar(MainWindow)
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionRun = QtGui.QAction(MainWindow)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/run.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionRun.setIcon(icon)
        self.actionRun.setObjectName(_fromUtf8("actionRun"))
        self.actionE = QtGui.QAction(MainWindow)
        self.actionE.setCheckable(True)
        self.actionE.setChecked(True)
        self.actionE.setObjectName(_fromUtf8("actionE"))
        self.actionS = QtGui.QAction(MainWindow)
        self.actionS.setCheckable(True)
        self.actionS.setChecked(True)
        self.actionS.setObjectName(_fromUtf8("actionS"))
        self.actionK = QtGui.QAction(MainWindow)
        self.actionK.setCheckable(True)
        self.actionK.setChecked(True)
        self.actionK.setObjectName(_fromUtf8("actionK"))
        self.actionRms = QtGui.QAction(MainWindow)
        self.actionRms.setCheckable(True)
        self.actionRms.setChecked(True)
        self.actionRms.setObjectName(_fromUtf8("actionRms"))
        self.actionNimages = QtGui.QAction(MainWindow)
        self.actionNimages.setCheckable(True)
        self.actionNimages.setChecked(True)
        self.actionNimages.setObjectName(_fromUtf8("actionNimages"))
        self.action3D = QtGui.QAction(MainWindow)
        self.action3D.setCheckable(True)
        self.action3D.setObjectName(_fromUtf8("action3D"))
        self.actionParams = QtGui.QAction(MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/preferences.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionParams.setIcon(icon1)
        self.actionParams.setObjectName(_fromUtf8("actionParams"))
        self.actionReset = QtGui.QAction(MainWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/restart.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionReset.setIcon(icon2)
        self.actionReset.setObjectName(_fromUtf8("actionReset"))
        self.actionSave = QtGui.QAction(MainWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/save.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionSave.setIcon(icon3)
        self.actionSave.setObjectName(_fromUtf8("actionSave"))
        self.actionLoad = QtGui.QAction(MainWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/open.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionLoad.setIcon(icon4)
        self.actionLoad.setObjectName(_fromUtf8("actionLoad"))
        self.actionTS = QtGui.QAction(MainWindow)
        self.actionTS.setObjectName(_fromUtf8("actionTS"))
        self.toolBar.addAction(self.actionRun)
        self.toolBar.addAction(self.actionReset)
        self.toolBar.addAction(self.actionParams)
        self.toolBar.addAction(self.actionSave)
        self.toolBar.addAction(self.actionLoad)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionE)
        self.toolBar.addAction(self.actionS)
        self.toolBar.addAction(self.actionK)
        self.toolBar.addAction(self.actionNimages)
        self.toolBar.addAction(self.actionRms)
        self.toolBar.addAction(self.action3D)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionTS)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "NEB explorer", None, QtGui.QApplication.UnicodeUTF8))
        self.toolBar.setWindowTitle(QtGui.QApplication.translate("MainWindow", "toolBar", None, QtGui.QApplication.UnicodeUTF8))
        self.actionRun.setText(QtGui.QApplication.translate("MainWindow", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.actionRun.setToolTip(QtGui.QApplication.translate("MainWindow", "run the neb (continue)", None, QtGui.QApplication.UnicodeUTF8))
        self.actionE.setText(QtGui.QApplication.translate("MainWindow", "E", None, QtGui.QApplication.UnicodeUTF8))
        self.actionE.setToolTip(QtGui.QApplication.translate("MainWindow", "show energies", None, QtGui.QApplication.UnicodeUTF8))
        self.actionS.setText(QtGui.QApplication.translate("MainWindow", "s", None, QtGui.QApplication.UnicodeUTF8))
        self.actionS.setToolTip(QtGui.QApplication.translate("MainWindow", "path length", None, QtGui.QApplication.UnicodeUTF8))
        self.actionK.setText(QtGui.QApplication.translate("MainWindow", "k", None, QtGui.QApplication.UnicodeUTF8))
        self.actionRms.setText(QtGui.QApplication.translate("MainWindow", "rms", None, QtGui.QApplication.UnicodeUTF8))
        self.actionNimages.setText(QtGui.QApplication.translate("MainWindow", "nimages", None, QtGui.QApplication.UnicodeUTF8))
        self.action3D.setText(QtGui.QApplication.translate("MainWindow", "3D", None, QtGui.QApplication.UnicodeUTF8))
        self.actionParams.setText(QtGui.QApplication.translate("MainWindow", "params", None, QtGui.QApplication.UnicodeUTF8))
        self.actionParams.setToolTip(QtGui.QApplication.translate("MainWindow", "edit parameters", None, QtGui.QApplication.UnicodeUTF8))
        self.actionReset.setText(QtGui.QApplication.translate("MainWindow", "reset", None, QtGui.QApplication.UnicodeUTF8))
        self.actionReset.setToolTip(QtGui.QApplication.translate("MainWindow", "start from fresh interpolation", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave.setText(QtGui.QApplication.translate("MainWindow", "save", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave.setToolTip(QtGui.QApplication.translate("MainWindow", "save path to file", None, QtGui.QApplication.UnicodeUTF8))
        self.actionLoad.setText(QtGui.QApplication.translate("MainWindow", "load", None, QtGui.QApplication.UnicodeUTF8))
        self.actionLoad.setToolTip(QtGui.QApplication.translate("MainWindow", "load path from file", None, QtGui.QApplication.UnicodeUTF8))
        self.actionTS.setText(QtGui.QApplication.translate("MainWindow", "TS refinement", None, QtGui.QApplication.UnicodeUTF8))
        self.actionTS.setToolTip(QtGui.QApplication.translate("MainWindow", "switch to ts refinement", None, QtGui.QApplication.UnicodeUTF8))

from . import resources_rc

