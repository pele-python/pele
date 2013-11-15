# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_normalmode_explorer.ui'
#
# Created: Wed Nov 13 11:48:24 2013
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
        MainWindow.resize(863, 613)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.splitter = QtGui.QSplitter(self.centralwidget)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName(_fromUtf8("splitter"))
        self.verticalLayoutWidget_2 = QtGui.QWidget(self.splitter)
        self.verticalLayoutWidget_2.setObjectName(_fromUtf8("verticalLayoutWidget_2"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setMargin(0)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.mplwidget = MPLWidget(self.verticalLayoutWidget_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mplwidget.sizePolicy().hasHeightForWidth())
        self.mplwidget.setSizePolicy(sizePolicy)
        self.mplwidget.setMinimumSize(QtCore.QSize(0, 0))
        self.mplwidget.setMaximumSize(QtCore.QSize(16777215, 300))
        self.mplwidget.setObjectName(_fromUtf8("mplwidget"))
        self.verticalLayout_2.addWidget(self.mplwidget)
        self.view3D = Show3DWithSlider(self.verticalLayoutWidget_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.view3D.sizePolicy().hasHeightForWidth())
        self.view3D.setSizePolicy(sizePolicy)
        self.view3D.setMinimumSize(QtCore.QSize(200, 200))
        self.view3D.setObjectName(_fromUtf8("view3D"))
        self.verticalLayout_2.addWidget(self.view3D)
        self.verticalLayoutWidget = QtGui.QWidget(self.splitter)
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(self.verticalLayoutWidget)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.listNormalmodes = QtGui.QListWidget(self.verticalLayoutWidget)
        self.listNormalmodes.setObjectName(_fromUtf8("listNormalmodes"))
        self.verticalLayout.addWidget(self.listNormalmodes)
        self.pushClose = QtGui.QPushButton(self.verticalLayoutWidget)
        self.pushClose.setObjectName(_fromUtf8("pushClose"))
        self.verticalLayout.addWidget(self.pushClose)
        self.horizontalLayout.addWidget(self.splitter)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 863, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.toolBar = QtGui.QToolBar(MainWindow)
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionRun = QtGui.QAction(MainWindow)
        self.actionRun.setCheckable(True)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/run.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/run.png")), QtGui.QIcon.Normal, QtGui.QIcon.On)
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/run.png")), QtGui.QIcon.Disabled, QtGui.QIcon.Off)
        self.actionRun.setIcon(icon)
        self.actionRun.setObjectName(_fromUtf8("actionRun"))
        self.actionSave = QtGui.QAction(MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/save.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionSave.setIcon(icon1)
        self.actionSave.setObjectName(_fromUtf8("actionSave"))
        self.actionParameters = QtGui.QAction(MainWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(_fromUtf8(":/icons/icons/preferences.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionParameters.setIcon(icon2)
        self.actionParameters.setObjectName(_fromUtf8("actionParameters"))
        self.actionShow_energies = QtGui.QAction(MainWindow)
        self.actionShow_energies.setCheckable(True)
        self.actionShow_energies.setObjectName(_fromUtf8("actionShow_energies"))
        self.actionHessian_eigs = QtGui.QAction(MainWindow)
        self.actionHessian_eigs.setCheckable(True)
        self.actionHessian_eigs.setObjectName(_fromUtf8("actionHessian_eigs"))
        self.toolBar.addAction(self.actionRun)
        self.toolBar.addAction(self.actionSave)
        self.toolBar.addAction(self.actionParameters)
        self.toolBar.addAction(self.actionShow_energies)
        self.toolBar.addAction(self.actionHessian_eigs)

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.pushClose, QtCore.SIGNAL(_fromUtf8("clicked()")), MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Normalmode Browser", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Normalmodes", None, QtGui.QApplication.UnicodeUTF8))
        self.pushClose.setText(QtGui.QApplication.translate("MainWindow", "Close", None, QtGui.QApplication.UnicodeUTF8))
        self.toolBar.setWindowTitle(QtGui.QApplication.translate("MainWindow", "toolBar", None, QtGui.QApplication.UnicodeUTF8))
        self.actionRun.setText(QtGui.QApplication.translate("MainWindow", "Play", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave.setText(QtGui.QApplication.translate("MainWindow", "Save", None, QtGui.QApplication.UnicodeUTF8))
        self.actionParameters.setText(QtGui.QApplication.translate("MainWindow", "Parameters", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_energies.setText(QtGui.QApplication.translate("MainWindow", "show energies", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_energies.setToolTip(QtGui.QApplication.translate("MainWindow", "<html><head/><body><p>Show a plot of the energies.  Note the energy will only agree with the harmonic approximation if the metric tensor is the identity.  Click the button &quot;Hessian eigs&quot; to make the energies agree.</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.actionHessian_eigs.setText(QtGui.QApplication.translate("MainWindow", "Hessian eigs", None, QtGui.QApplication.UnicodeUTF8))
        self.actionHessian_eigs.setToolTip(QtGui.QApplication.translate("MainWindow", "<html><head/><body><p>Show the Hessian eigenvalues and eigenvectors instead of the normal modes. These will be different if the metric tensor is not the identity (e.g. if with curvilinear coordinates like angle-axis).</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))

from pele.gui.show3d_with_slider import Show3DWithSlider
from pele.gui.ui.mplwidget import MPLWidget
import resources_rc
