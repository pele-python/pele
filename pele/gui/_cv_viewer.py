import numpy as np

from PyQt5 import QtWidgets, QtCore

from pele.gui.ui.cv_viewer_ui import Ui_Form
from pele.thermodynamics import GetThermodynamicInfoParallel, minima_to_cv
from pele.utils.events import Signal


class GetThermodynamicInfoParallelQT(GetThermodynamicInfoParallel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.on_finish = Signal()
        self.on_error = Signal()
        self.running = False

    def poll(self):
        if not self.running:
            return
        if self.njobs == 0:
            self.refresh_timer.stop()
            self.finish()
            return
        if not self.done_queue.empty():
            self.njobs -= 1
            ret = self.done_queue.get()
            try:
                self._process_return_value(ret)
            except Exception as error:
                self.cancel()
                self.on_error(error)
        elif all(not worker.is_alive() for worker in self.workers):
            # a worker's last result can arrive between the empty() check
            # above and its exit; once it has exited, the queue holds it
            if not self.done_queue.empty():
                return
            self.cancel()
            self.on_error(RuntimeError(
                "Thermodynamic worker exited before returning all results."
            ))

    def finish(self):
        self.refresh_timer.stop()
        super().finish()
        self.running = False
        self._close_queues()
        self.on_finish()

    def _close_queues(self):
        for queue in (self.send_queue, self.done_queue):
            queue.cancel_join_thread()
            queue.close()

    def cancel(self):
        if not self.running:
            return
        self.refresh_timer.stop()
        self._kill_workers()
        for worker in self.workers:
            worker.close()
        self.running = False
        self._close_queues()

    def start(self):
        # Fork before put() starts the queue feeder thread.
        for worker in self.workers:
            worker.start()
        self._populate_queue()
        self.running = True

        self.refresh_timer = QtCore.QTimer()
        self.refresh_timer.timeout.connect(self.poll)
        self.refresh_timer.start(50)  # time in msec


class HeatCapacityWidget(QtWidgets.QWidget):
    def __init__(self, system, database, parent=None):
        super().__init__(parent=parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self.system = system
        self.database = database

        self.canvas = self.ui.mplwidget.canvas
        self.axes = self.canvas.axes

    def rebuild_cv_plot(self):
        try:
            self._get_T_range()
            self._get_nmin_max()
            self._compute_thermodynamic_info(on_finish=self.make_cv_plot)
        except ValueError as error:
            self._show_error(error)

    def make_cv_plot(self):
        try:
            self._compute_cv()
            self._plot_cv()
        except ValueError as error:
            self._show_error(error)

    def _show_error(self, error):
        self.ui.label_status.setText(str(error))

    def closeEvent(self, event):
        if hasattr(self, "worker"):
            self.worker.cancel()
        super().closeEvent(event)

    def _get_ndof(self):
        return self.system.get_ndof()

    def _get_nmin_max(self):
        txt = self.ui.lineEdit_nmin_max.text()
        if len(txt) > 0:
            nmin_max = int(txt)
            if nmin_max <= 0:
                raise ValueError("Maximum number of minima must be positive.")
        else:
            nmin_max = None
        return nmin_max

    def _compute_thermodynamic_info(
        self, nproc=2, on_finish=None, verbose=False
    ):
        nmin = self._get_nmin_max()
        if nmin is not None and nmin < self.database.number_of_minima():
            self.minima = self.database.minima()[:nmin]
        else:
            self.minima = self.database.minima()
        if not self.minima:
            raise ValueError("No minima are available for heat capacity.")

        if hasattr(self, "worker"):
            self.worker.cancel()

        self.worker = GetThermodynamicInfoParallelQT(
            self.system,
            self.database,
            npar=nproc,
            verbose=verbose,
            only_minima=True,
        )
        if on_finish is not None:
            self.worker.on_finish.connect(on_finish)
        self.worker.on_error.connect(self._show_error)
        self.worker.start()

        njobs = self.worker.njobs
        self.ui.label_status.setText(
            "computing thermodynamic information for %d minima" % njobs
        )

    def _compute_cv(self):
        Tlist = self._get_T_range()
        if not any(not minimum.invalid for minimum in self.minima):
            raise ValueError("No valid minima are available for heat capacity.")
        lZ, U, U2, Cv = minima_to_cv(self.minima, Tlist, self._get_ndof())
        self.Tlist = Tlist
        self.Cv = Cv

    def _get_Tmin(self):
        res = 0.01
        txt = self.ui.lineEdit_Tmin.text()
        if len(txt) > 0:
            res = float(txt)
        return res

    def _get_Tmax(self):
        res = 1.0
        txt = self.ui.lineEdit_Tmax.text()
        if len(txt) > 0:
            res = float(txt)
        return res

    def _get_nT(self):
        res = 100
        txt = self.ui.lineEdit_nT.text()
        if len(txt) > 0:
            res = int(txt)
        return res

    def _get_T_range(self):
        Tmin = self._get_Tmin()
        Tmax = self._get_Tmax()
        nT = self._get_nT()

        if not np.isfinite([Tmin, Tmax]).all() or not 0 < Tmin < Tmax:
            raise ValueError("Temperatures must satisfy 0 < Tmin < Tmax.")
        if nT <= 0:
            raise ValueError("Number of points must be a positive integer.")
        return np.linspace(Tmin, Tmax, nT, endpoint=False)

    def _plot_cv(self):
        self.ui.label_status.setText(
            "showing heat capacity calculated with %d minima" % len(self.minima)
        )
        axes = self.axes
        axes.clear()
        axes.plot(self.Tlist, self.Cv)
        self.canvas.draw()

    def on_btn_recalculate_clicked(self, clicked=None):
        if clicked is None:
            return
        self.rebuild_cv_plot()


class HeatCapacityViewer(QtWidgets.QMainWindow):
    def __init__(self, system, database, parent=None, app=None):
        super().__init__(parent=parent)
        self.cv_widget = HeatCapacityWidget(system, database, parent=self)
        self.setCentralWidget(self.cv_widget)
        self.setWindowTitle("Harmonic Superposition Heat Capacity")

    def rebuild_cv_plot(self):
        self.cv_widget.rebuild_cv_plot()

    def closeEvent(self, event):
        if hasattr(self.cv_widget, "worker"):
            self.cv_widget.worker.cancel()
        super().closeEvent(event)


def test():
    import sys
    from pele.systems import LJCluster

    app = QtWidgets.QApplication(sys.argv)
    system = LJCluster(13)

    db = system.create_database()
    bh = system.get_basinhopping(db, outstream=None)
    bh.run(200)

    obj = HeatCapacityViewer(system, db)
    obj.show()

    def test_start():
        obj.rebuild_cv_plot()

    QtCore.QTimer.singleShot(10, test_start)
    sys.exit(app.exec_())


if __name__ == "__main__":
    test()
