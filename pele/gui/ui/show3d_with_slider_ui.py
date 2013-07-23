# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'show3d_with_slider_ui.ui'
#
# Created: Wed Apr 24 11:36:52 2013
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_show3d_with_slider(object):
    def setupUi(self, show3d_with_slider):
        show3d_with_slider.setObjectName(_fromUtf8("show3d_with_slider"))
        show3d_with_slider.resize(658, 587)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(show3d_with_slider.sizePolicy().hasHeightForWidth())
        show3d_with_slider.setSizePolicy(sizePolicy)
        show3d_with_slider.setMinimumSize(QtCore.QSize(200, 200))
        self.gridLayout = QtGui.QGridLayout(show3d_with_slider)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.btn_animate = QtGui.QPushButton(show3d_with_slider)
        self.btn_animate.setObjectName(_fromUtf8("btn_animate"))
        self.horizontalLayout.addWidget(self.btn_animate)
        self.label = QtGui.QLabel(show3d_with_slider)
        self.label.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.oglwgt = Show3D(show3d_with_slider)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.oglwgt.sizePolicy().hasHeightForWidth())
        self.oglwgt.setSizePolicy(sizePolicy)
        self.oglwgt.setObjectName(_fromUtf8("oglwgt"))
        self.verticalLayout.addWidget(self.oglwgt)
        self.slider = QtGui.QSlider(show3d_with_slider)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.slider.sizePolicy().hasHeightForWidth())
        self.slider.setSizePolicy(sizePolicy)
        self.slider.setOrientation(QtCore.Qt.Horizontal)
        self.slider.setObjectName(_fromUtf8("slider"))
        self.verticalLayout.addWidget(self.slider)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.retranslateUi(show3d_with_slider)
        QtCore.QMetaObject.connectSlotsByName(show3d_with_slider)

    def retranslateUi(self, show3d_with_slider):
        show3d_with_slider.setWindowTitle(QtGui.QApplication.translate("show3d_with_slider", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_animate.setText(QtGui.QApplication.translate("show3d_with_slider", "animate", None, QtGui.QApplication.UnicodeUTF8))

from pele.gui.show3d import Show3D
