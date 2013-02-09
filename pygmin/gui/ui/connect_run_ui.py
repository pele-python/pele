# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'connect_run_ui.ui'
#
# Created: Sat Feb  9 16:57:57 2013
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

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(764, 829)
        MainWindow.setDockNestingEnabled(False)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setEnabled(True)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 764, 22))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.view_Log = QtGui.QDockWidget(MainWindow)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.view_Log.sizePolicy().hasHeightForWidth())
        self.view_Log.setSizePolicy(sizePolicy)
        self.view_Log.setFeatures(QtGui.QDockWidget.DockWidgetFloatable|QtGui.QDockWidget.DockWidgetMovable)
        self.view_Log.setAllowedAreas(QtCore.Qt.AllDockWidgetAreas)
        self.view_Log.setObjectName(_fromUtf8("view_Log"))
        self.dockWidgetContents_5 = QtGui.QWidget()
        self.dockWidgetContents_5.setObjectName(_fromUtf8("dockWidgetContents_5"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.dockWidgetContents_5)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.textEdit = QtGui.QTextEdit(self.dockWidgetContents_5)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.textEdit.sizePolicy().hasHeightForWidth())
        self.textEdit.setSizePolicy(sizePolicy)
        self.textEdit.setMinimumSize(QtCore.QSize(0, 0))
        self.textEdit.setBaseSize(QtCore.QSize(0, 500))
        self.textEdit.setObjectName(_fromUtf8("textEdit"))
        self.horizontalLayout.addWidget(self.textEdit)
        self.view_Log.setWidget(self.dockWidgetContents_5)
        MainWindow.addDockWidget(QtCore.Qt.DockWidgetArea(8), self.view_Log)
        self.view_ogl = QtGui.QDockWidget(MainWindow)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.view_ogl.sizePolicy().hasHeightForWidth())
        self.view_ogl.setSizePolicy(sizePolicy)
        self.view_ogl.setFeatures(QtGui.QDockWidget.DockWidgetFloatable|QtGui.QDockWidget.DockWidgetMovable)
        self.view_ogl.setAllowedAreas(QtCore.Qt.AllDockWidgetAreas)
        self.view_ogl.setObjectName(_fromUtf8("view_ogl"))
        self.dockWidgetContents_3 = QtGui.QWidget()
        self.dockWidgetContents_3.setObjectName(_fromUtf8("dockWidgetContents_3"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.dockWidgetContents_3)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.ogl = Show3DWithSlider(self.dockWidgetContents_3)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ogl.sizePolicy().hasHeightForWidth())
        self.ogl.setSizePolicy(sizePolicy)
        self.ogl.setObjectName(_fromUtf8("ogl"))
        self.horizontalLayout_2.addWidget(self.ogl)
        self.view_ogl.setWidget(self.dockWidgetContents_3)
        MainWindow.addDockWidget(QtCore.Qt.DockWidgetArea(4), self.view_ogl)
        self.toolBar = QtGui.QToolBar(MainWindow)
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionEnergy = QtGui.QAction(MainWindow)
        self.actionEnergy.setCheckable(True)
        self.actionEnergy.setChecked(False)
        self.actionEnergy.setObjectName(_fromUtf8("actionEnergy"))
        self.action3D = QtGui.QAction(MainWindow)
        self.action3D.setCheckable(True)
        self.action3D.setObjectName(_fromUtf8("action3D"))
        self.actionGraph = QtGui.QAction(MainWindow)
        self.actionGraph.setCheckable(True)
        self.actionGraph.setObjectName(_fromUtf8("actionGraph"))
        self.actionStop = QtGui.QAction(MainWindow)
        self.actionStop.setCheckable(True)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/pause.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/pause.png")), QtGui.QIcon.Normal, QtGui.QIcon.On)
        self.actionStop.setIcon(icon)
        self.actionStop.setObjectName(_fromUtf8("actionStop"))
        self.actionD_Graph = QtGui.QAction(MainWindow)
        self.actionD_Graph.setCheckable(True)
        self.actionD_Graph.setObjectName(_fromUtf8("actionD_Graph"))
        self.actionSummary = QtGui.QAction(MainWindow)
        self.actionSummary.setCheckable(True)
        self.actionSummary.setObjectName(_fromUtf8("actionSummary"))
        self.actionLog = QtGui.QAction(MainWindow)
        self.actionLog.setCheckable(True)
        self.actionLog.setObjectName(_fromUtf8("actionLog"))
        self.toolBar.addAction(self.actionStop)
        self.toolBar.addAction(self.actionLog)
        self.toolBar.addAction(self.action3D)
        self.toolBar.addAction(self.actionEnergy)
        self.toolBar.addAction(self.actionGraph)
        self.toolBar.addAction(self.actionD_Graph)
        self.toolBar.addAction(self.actionSummary)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.view_Log.setWindowTitle(_translate("MainWindow", "Log", None))
        self.view_ogl.setWindowTitle(_translate("MainWindow", "3D view", None))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar", None))
        self.actionEnergy.setText(_translate("MainWindow", "Energy", None))
        self.actionEnergy.setToolTip(_translate("MainWindow", "toggle energy window", None))
        self.action3D.setText(_translate("MainWindow", "3D", None))
        self.action3D.setToolTip(_translate("MainWindow", "toggle 3D viewer", None))
        self.actionGraph.setText(_translate("MainWindow", "Graph", None))
        self.actionGraph.setToolTip(_translate("MainWindow", "toggle graph view", None))
        self.actionStop.setText(_translate("MainWindow", "stop", None))
        self.actionD_Graph.setText(_translate("MainWindow", "D-Graph", None))
        self.actionD_Graph.setToolTip(_translate("MainWindow", "disconnectivity graph", None))
        self.actionSummary.setText(_translate("MainWindow", "Summary", None))
        self.actionSummary.setToolTip(_translate("MainWindow", "display summary information", None))
        self.actionLog.setText(_translate("MainWindow", "Log", None))
        self.actionLog.setToolTip(_translate("MainWindow", "display log information", None))

from pygmin.gui.show3d import Show3DWithSlider
import resources_rc
