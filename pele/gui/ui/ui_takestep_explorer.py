# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_takestep_explorer.ui'
#
# Created: Wed Apr 24 11:11:57 2013
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(854, 611)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setContentsMargins(5, -1, -1, -1)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout_2.addWidget(self.label)
        self.show3d = Show3DWithSlider(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.show3d.sizePolicy().hasHeightForWidth())
        self.show3d.setSizePolicy(sizePolicy)
        self.show3d.setObjectName(_fromUtf8("show3d"))
        self.verticalLayout_2.addWidget(self.show3d)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setContentsMargins(0, -1, -1, -1)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setMaximumSize(QtCore.QSize(200, 16777215))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout.addWidget(self.label_2)
        self.listMinima = QtGui.QListWidget(self.centralwidget)
        self.listMinima.setMaximumSize(QtCore.QSize(200, 16777215))
        self.listMinima.setObjectName(_fromUtf8("listMinima"))
        self.verticalLayout.addWidget(self.listMinima)
        self.horizontalLayout.addLayout(self.verticalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 854, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.toolBar = QtGui.QToolBar(MainWindow)
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionDisplace = QtGui.QAction(MainWindow)
        self.actionDisplace.setObjectName(_fromUtf8("actionDisplace"))
        self.actionQuench = QtGui.QAction(MainWindow)
        self.actionQuench.setObjectName(_fromUtf8("actionQuench"))
        self.actionShow_path = QtGui.QAction(MainWindow)
        self.actionShow_path.setCheckable(True)
        self.actionShow_path.setChecked(False)
        self.actionShow_path.setObjectName(_fromUtf8("actionShow_path"))
        self.toolBar.addAction(self.actionDisplace)
        self.toolBar.addAction(self.actionQuench)
        self.toolBar.addAction(self.actionShow_path)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "Energy (id)", None, QtGui.QApplication.UnicodeUTF8))
        self.toolBar.setWindowTitle(QtGui.QApplication.translate("MainWindow", "toolBar", None, QtGui.QApplication.UnicodeUTF8))
        self.actionDisplace.setText(QtGui.QApplication.translate("MainWindow", "displace", None, QtGui.QApplication.UnicodeUTF8))
        self.actionQuench.setText(QtGui.QApplication.translate("MainWindow", "quench", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_path.setText(QtGui.QApplication.translate("MainWindow", "show path", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_path.setToolTip(QtGui.QApplication.translate("MainWindow", "display the quench trajectory", None, QtGui.QApplication.UnicodeUTF8))

from pele.gui.show3d_with_slider import Show3DWithSlider
