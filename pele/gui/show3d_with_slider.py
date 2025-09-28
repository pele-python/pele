from PyQt5 import QtCore
from PyQt5.Qt import QWidget
from PyQt5.QtCore import pyqtSlot

from pele.gui.ui.show3d_with_slider_ui import Ui_show3d_with_slider
from pele.utils.events import Signal


class Show3DWithSlider(QWidget):
    """
    The main OpenGL viewer for pele
    
    This viewer can display a single structure, two structures
    overlaid, or a pathway of structures (with a slider) 
    """
    def __init__(self, *args, **kwargs):
        super(Show3DWithSlider, self).__init__(*args, **kwargs)
        self.setMinimumSize(200, 200)

        
        self.ui = Ui_show3d_with_slider()
        self.ui.setupUi(self)
         
        self.label = self.ui.label
        self.label.hide()
        
        self.ui.btn_animate.hide()
        self.ui.btn_animate.setCheckable(True)
        
        self.oglwgt = self.ui.oglwgt
        

        self.slider = self.ui.slider
        
        self.on_frame_updated = Signal()

        self.animate = False
        self._animate_dir = 1
    
    def setSystem(self, system):
        """
        set the system which has information about how to draw the structures.
        
        it should have a function:
        
            system.draw(coords)
        """
        self.oglwgt.setSystem(system)
        
    def setCoords(self, coords, index=1):
        """
        display a structure (not a path) in the ogl viewer
        
        Parameters
        ----------
        coords : 1d array
            the structure to display
        index : 1 or 2
            the viewer can display two structures overlayed.  This specifies
            whether this structure should be in slot 1 or slot 2
        """
        if index not in (1, 2):
            raise ValueError("index must be either 1 or 2")
        self.animate = False
        self.messages = None
        self.coordspath = None
        self.slider.hide()
        self.oglwgt.setCoords(coords, index=index)
        self.label.hide()
        self.ui.btn_animate.hide()

    def setCoordsPath(self, coordspath, frame=None, labels=None):
        """
        show a path in the viewer
        
        Parameters
        ----------
        coordspath : 2d numpy array
            coordspath[i,j] is the coordinates of the j'th degree of freedom of the i'th structure
        frame : int
            frame number to show initially
        labels : list of strings
            labels for the structures that will be shown above the ogl viewer
        """
        self.oglwgt.setCoords(None, index=2)
        if frame is None:
            frame = self.slider.value()
        
        self.ui.btn_animate.show()
        
        if labels is None:
            self.label.hide()
        else:
            self.label.show()
        
        if frame < 0:
            frame = coordspath.shape[0]-1
        else:
            frame = min(frame, coordspath.shape[0]-1)
        
        self.coordspath = coordspath
        self.messages = labels
        self.slider.show()
        self.slider.setRange(0, coordspath.shape[0]-1)
        self.showFrame(frame)
    
    @pyqtSlot(int)
    def on_slider_valueChanged(self, i):
        return self._showFrame(i)

    def _showFrame(self, i):
        self.oglwgt.setCoords(self.coordspath[i,:], index=1)
        if self.messages is not None:
            self.label.setText(self.messages[i])
        self.on_frame_updated(i, sender=self)

    def showFrame(self, i):
        if i == -1:
            i = self.coordspath.shape[0] - 1
        self.slider.setValue(i)
        self._showFrame(i)

    def get_slider_index(self):
        return self.slider.value()

    def on_btn_animate_clicked(self, checked=None):
        if checked is None: return
        if checked:
            self.start_animation()
        else:
            self.stop_animation()

    def start_animation(self):
        self.animate=True
        self._animate_dir = 1
        QtCore.QTimer.singleShot(0., self._next_frame)
    
    def stop_animation(self):
        self.animate = False
    
    def _next_frame(self):
        if not self.animate: return
        cur = self.slider.value()
        if cur == self.slider.maximum():
            self._animate_dir = -1
        elif cur ==  self.slider.minimum():
            self._animate_dir = 1
        cur += self._animate_dir
#        self.slider.setValue(cur)
        self.showFrame(cur)
        
        if self.animate:
            frames_per_second = 10.
            QtCore.QTimer.singleShot(1000. / frames_per_second, self._next_frame)

    def sizeHint(self):
        w, h = 500,500 #self.get_width_height()
        return QtCore.QSize(w, h)

    def minimumSizeHint(self):
        return QtCore.QSize(10, 10)

