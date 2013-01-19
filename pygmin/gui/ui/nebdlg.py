import numpy as np
from PyQt4.QtGui import QDialog, QApplication
import nebbrowser

import sys

app = QApplication(sys.argv)
    
wnd = QDialog()

ui = nebbrowser.Ui_Form()
ui.setupUi(wnd)

wnd.show()

ui.widget.axes.plot(np.sin(np.arange(0, 6, 0.01)))

sys.exit(app.exec_()) 
    