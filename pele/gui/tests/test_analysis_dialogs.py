"""Headless regressions for analysis dialog actions and Python 3 data handling."""
import pickle
import os
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pytest

pytest.importorskip("PyQt5")
pytest.importorskip("OpenGL")
from PyQt5 import QtCore, QtWidgets

from pele.gui.dlg_params import DlgParams, EditParamsWidget
from pele.gui.neb_explorer import NEBExplorer, NEBRunner
from pele.gui.normalmode_browser import NormalmodeBrowser
from pele.gui.takestep_explorer import TakestepExplorer
from pele.gui._cv_viewer import HeatCapacityWidget, HeatCapacityViewer, GetThermodynamicInfoParallelQT
from pele.gui._rate_gui import RateWidget, RateViewer
from pele.storage import Database
from pele.thermodynamics import GetThermodynamicInfoParallel


class AnalysisDialogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self):
        self.system = Mock()
        self.system.get_orthogonalize_to_zero_eigenvectors.return_value = None
        self.system.get_potential.return_value.getEnergy.return_value = 0.0
        self.system.get_potential.return_value.getEnergyGradient.return_value = (
            0.0, np.ones(3)
        )
        self.database = Database()
        self.addCleanup(self.database.engine.dispose)
        self.addCleanup(self.database.connection.close)
        self.addCleanup(self.database.session.close)

    def keep(self, widget):
        self.addCleanup(widget.close)
        return widget

    def test_parameter_edit_updates_original_and_nested_values(self):
        for cls in (DlgParams, EditParamsWidget):
            with self.subTest(dialog=cls.__name__):
                params = {"count": 2, "nested": {"enabled": True}, "optional": None}
                widget = self.keep(cls(params=params))
                count = widget.model.findItems("count")[0]
                value = widget.model.item(count.row(), 1)
                widget.model.blockSignals(True)
                value.setText("7")
                widget.item_changed(value)
                self.assertEqual(params["count"], 7)
                value.setText("invalid")
                widget.item_changed(value)
                self.assertEqual(value.text(), "7")
                nested = widget.model.findItems("nested")[0]
                self.assertEqual(nested.rowCount(), 1)
                enabled = nested.child(0, 1)
                enabled.setCheckState(QtCore.Qt.Unchecked)
                widget.item_changed(enabled)
                self.assertFalse(params["nested"]["enabled"])
                optional = widget.model.findItems("optional")[0]
                value = widget.model.item(optional.row(), 1)
                value.setText("invalid")
                widget.item_changed(value)
                self.assertIsNone(params["optional"])
                self.assertEqual(value.text(), "None")

    def test_neb_saved_path_can_be_loaded(self):
        widget = self.keep(NEBExplorer(system=self.system, app=self.app))
        path = np.arange(9.0).reshape(3, 3)
        widget.coords1, widget.coords2 = path[0], path[-1]
        widget.nebrunner.run = Mock()
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "path.pickle"
            with filename.open("wb") as output:
                pickle.dump(path, output)
            with patch("pele.gui.neb_explorer.QtWidgets.QFileDialog") as dialog:
                dialog.return_value.exec_.return_value = True
                dialog.return_value.selectedFiles.return_value = [str(filename)]
                widget.on_actionLoad_triggered(False)
        np.testing.assert_array_equal(widget.initial_path, path)
        self.assertFalse(widget.nebrunner.run.call_args.kwargs["run"])

    def test_unimplemented_parameter_menu_actions_are_disabled(self):
        for cls in (DlgParams, EditParamsWidget):
            widget = self.keep(cls(params={"count": 1, "nested": {}}))
            view = widget.ui.treeParams if cls is DlgParams else widget.view
            for key in ("count", "nested"):
                with self.subTest(dialog=cls.__name__, parameter=key):
                    item = widget.model.findItems(key)[0]
                    view.setCurrentIndex(widget.model.indexFromItem(item))
                    menu = QtWidgets.QMenu(widget)
                    with patch.object(menu, "exec_"), patch(
                        "pele.gui.dlg_params.QtWidgets.QMenu", return_value=menu
                    ):
                        widget.open_context_menu(QtCore.QPoint())
                    self.assertEqual(len(menu.actions()), 1)
                    self.assertFalse(menu.actions()[0].isEnabled())

    def test_neb_zero_update_frequency_keeps_endpoints(self):
        runner = NEBRunner(self.app, self.system, freq=0)
        runner.step_shift = 0
        for attr in ("stepnum", "k", "rms", "nimages", "energies", "distances"):
            setattr(runner, attr, [])
        args = dict(energies=np.zeros(2), distances=np.ones(1),
                    path=np.zeros((2, 3)), stepnum=0, rms=0, k=1)
        runner._neb_update(event="initial", **args)
        runner._neb_update(event="update", **args)
        self.assertEqual(len(runner.energies), 1)

    def test_neb_energy_pick_selects_path_frame(self):
        widget = self.keep(NEBExplorer(system=self.system, app=self.app))
        path = np.arange(9.0).reshape(3, 3)
        widget.show3d.setCoordsPath(path)
        widget.energies.on_pick(SimpleNamespace(ind=[2]))
        self.assertEqual(widget.show3d.get_slider_index(), 2)

    def test_neb_failure_restores_run_actions(self):
        runner = NEBRunner(self.app, self.system)
        neb = Mock()
        neb.run.side_effect = RuntimeError("optimizer failed")
        runner.create_neb = Mock(return_value=neb)
        finished = Mock()
        runner.on_run_finished.connect(finished)
        with self.assertRaisesRegex(RuntimeError, "optimizer failed"):
            runner.run(np.zeros(3), np.ones(3))
        finished.assert_called_once()

    def test_normalmode_export_requires_a_selection(self):
        widget = self.keep(NormalmodeBrowser(system=self.system, app=self.app))
        self.assertFalse(widget.ui.actionSave.isEnabled())
        widget.set_coords(np.zeros(3), normalmodes=[(2.0, np.ones(3))])
        widget.ui.listNormalmodes.setCurrentRow(0)
        self.assertTrue(widget.ui.actionSave.isEnabled())
        widget.on_listNormalmodes_currentItemChanged(None)
        self.assertIsNone(widget.current_selection)
        self.assertFalse(widget.ui.actionSave.isEnabled())

    def test_normalmode_reloading_retains_selected_animation(self):
        widget = self.keep(NormalmodeBrowser(system=self.system, app=self.app))
        modes = [(2.0, np.ones(3))]
        widget.set_coords(np.zeros(3), normalmodes=modes)
        widget.ui.listNormalmodes.setCurrentRow(0)
        widget.set_coords(np.ones(3), normalmodes=modes)
        self.assertIsNotNone(widget.ui.view3D.coordspath)
        self.assertEqual(len(widget.ui.view3D.coordspath), 30)

    def test_takestep_uses_supplied_database_and_handles_cleared_selection(self):
        del self.system.database
        self.database.addMinimum(0.0, np.zeros(3))
        widget = self.keep(TakestepExplorer(system=self.system, app=self.app,
                                           database=self.database))
        self.assertEqual(widget.ui.listMinima.count(), 1)
        widget.ui.listMinima.setCurrentRow(0)
        widget.on_listMinima_currentItemChanged(None, None)
        self.assertIsNone(widget.quenched)

    def test_thermo_worker_starts_before_populating_queues(self):
        worker = GetThermodynamicInfoParallelQT(self.system, self.database, npar=0)
        calls = []
        worker.workers = [SimpleNamespace(start=lambda: calls.append("start"))]
        worker._populate_queue = lambda: calls.append("populate")
        with patch("pele.gui._cv_viewer.QtCore.QTimer"):
            worker.start()
        self.assertEqual(calls, ["start", "populate"])
        worker.send_queue.close()
        worker.done_queue.close()

    def test_thermo_worker_completes_and_reports_errors(self):
        self.database.addMinimum(0.0, np.zeros(3))
        self.system.get_pgorder.return_value = 1
        self.system.get_log_product_normalmode_freq.return_value = 0.0
        for failure in (None, RuntimeError("cannot compute modes")):
            with self.subTest(failure=failure):
                self.system.get_pgorder.side_effect = failure
                worker = GetThermodynamicInfoParallelQT(
                    self.system, self.database, npar=1, recalculate=True
                )
                finished, failed = Mock(), Mock()
                worker.on_finish.connect(finished)
                worker.on_error.connect(failed)
                worker.start()
                deadline = time.monotonic() + 5
                while worker.running and time.monotonic() < deadline:
                    self.app.processEvents()
                    time.sleep(0.01)
                still_running = worker.running
                worker.cancel()
                self.assertFalse(still_running, "thermodynamic worker did not finish")
                if failure is None:
                    finished.assert_called_once()
                    failed.assert_not_called()
                else:
                    finished.assert_not_called()
                    failed.assert_called_once()
                    self.assertEqual(str(failed.call_args.args[0]), str(failure))
                self.assertFalse(worker.refresh_timer.isActive())
                self.assertTrue(worker.send_queue._closed)
                self.assertTrue(worker.done_queue._closed)

    def test_heat_capacity_reports_an_empty_database(self):
        widget = self.keep(HeatCapacityWidget(self.system, self.database))
        widget.rebuild_cv_plot()
        self.assertIn("No minima", widget.ui.label_status.text())
        self.assertFalse(hasattr(widget, "worker"))

    def test_minima_only_thermodynamics_skips_transition_states(self):
        first = self.database.addMinimum(0.0, np.zeros(3))
        second = self.database.addMinimum(1.0, np.ones(3))
        self.database.addTransitionState(2.0, np.full(3, 2.0), first, second)
        worker = GetThermodynamicInfoParallel(
            self.system, self.database, npar=0, only_minima=True
        )
        worker._populate_queue()
        jobs = [worker.send_queue.get(timeout=1) for _ in range(worker.njobs)]
        worker.send_queue.close()
        worker.done_queue.close()
        self.assertEqual([job[0] for job in jobs], ["m", "m"])

    def test_thermo_cancel_stops_processes_and_timer(self):
        worker = GetThermodynamicInfoParallelQT(self.system, self.database, npar=1)
        worker.start()
        try:
            worker.cancel()
            self.assertFalse(worker.refresh_timer.isActive())
            self.assertTrue(all(process._closed for process in worker.workers))
            self.assertTrue(worker.send_queue._closed)
            self.assertTrue(worker.done_queue._closed)
            worker.cancel()  # closing a finished viewer is harmless
        finally:
            for process in worker.workers:
                if not process._closed:
                    process.terminate()
                    process.join(timeout=2)
                    process.close()
            worker.refresh_timer.stop()

    def test_thermo_reports_worker_exit_without_result(self):
        self.database.addMinimum(0.0, np.zeros(3))
        self.system.get_pgorder.side_effect = lambda coords: os._exit(3)
        worker = GetThermodynamicInfoParallelQT(self.system, self.database, npar=1)
        failed = Mock()
        worker.on_error.connect(failed)
        worker.start()
        try:
            worker.workers[0].join(timeout=2)
            worker.poll()
            self.assertFalse(worker.running)
            failed.assert_called_once()
            self.assertIn("exited", str(failed.call_args.args[0]))
        finally:
            worker.cancel()

    def test_thermo_keeps_result_that_arrives_as_worker_exits(self):
        # the result lands between poll()'s empty() check and its liveness check
        self.database.addMinimum(0.0, np.zeros(3))
        self.system.get_pgorder.return_value = 1
        self.system.get_log_product_normalmode_freq.return_value = 0.0
        worker = GetThermodynamicInfoParallelQT(self.system, self.database, npar=1)
        failed, finished = Mock(), Mock()
        worker.on_error.connect(failed)
        worker.on_finish.connect(finished)
        worker.start()
        try:
            worker.workers[0].join(timeout=5)
            empty = worker.done_queue.empty
            worker.done_queue.empty = Mock(side_effect=[True, empty(), empty()])
            worker.poll()  # sees an empty queue, then a dead worker
            worker.poll()  # takes the result
            worker.poll()  # finishes
            failed.assert_not_called()
            finished.assert_called_once()
        finally:
            worker.cancel()

    def test_closing_analysis_viewers_cancels_workers(self):
        for cls, attr in ((HeatCapacityViewer, "cv_widget"),
                          (RateViewer, "rate_widget")):
            viewer = self.keep(cls(self.system, self.database))
            worker = Mock()
            getattr(viewer, attr).worker = worker
            viewer.close()
            worker.cancel.assert_called_once()
            self.assertFalse(getattr(viewer, attr).isHidden())

    def test_heat_capacity_valid_default_plot(self):
        widget = self.keep(HeatCapacityWidget(self.system, self.database))
        self.system.get_ndof.return_value = 3
        minimum = self.database.addMinimum(0.0, np.zeros(3))
        minimum.fvib, minimum.pgorder = 0.0, 1
        widget.minima = [minimum]
        widget.make_cv_plot()
        self.assertEqual(len(widget.Tlist), 100)
        np.testing.assert_allclose(widget.Cv, 3)

    def test_heat_capacity_resize_gives_graph_the_extra_height(self):
        widget = self.keep(HeatCapacityWidget(self.system, self.database))
        widget.ui.label_status.setText("showing heat capacity calculated with 5 minima")
        widget.resize(480, 1050)
        widget.show()
        self.app.processEvents()
        self.assertFalse(widget.grab().isNull())
        self.assertLess(widget.ui.label_status.height(), widget.height() / 10)
        self.assertGreater(widget.ui.splitter.height(), widget.height() * 0.8)

    def test_heat_capacity_rejects_invalid_input_before_starting_workers(self):
        widget = self.keep(HeatCapacityWidget(self.system, self.database))
        cases = [("lineEdit_nT", "0"), ("lineEdit_nT", "2.5"),
                 ("lineEdit_Tmin", "-1"), ("lineEdit_Tmax", "nan"),
                 ("lineEdit_Tmax", "0.01"), ("lineEdit_nmin_max", "0")]
        with patch.object(widget, "_compute_thermodynamic_info") as start:
            for name, value in cases:
                with self.subTest(field=name, value=value):
                    field = getattr(widget.ui, name)
                    field.setText(value)
                    widget.rebuild_cv_plot()
                    start.assert_not_called()
                    self.assertTrue(widget.ui.label_status.text())
                    field.clear()

    def test_rates_reject_missing_or_identical_states_before_workers(self):
        widget = self.keep(RateWidget(self.system, self.database))
        with patch.object(widget, "_compute_thermodynamic_info") as start:
            widget.compute_rates()
            start.assert_not_called()
            self.assertTrue(widget.ui.label_status.text())
            minimum = self.database.addMinimum(0.0, np.zeros(3))
            widget.update_A(minimum)
            widget.update_B(minimum)
            widget.compute_rates()
            start.assert_not_called()

    def test_rates_clear_selections_and_enable_compute_only_for_distinct_states(self):
        widget = self.keep(RateWidget(self.system, self.database))
        self.assertFalse(widget.ui.btn_compute.isEnabled())
        first = self.database.addMinimum(0.0, np.zeros(3))
        second = self.database.addMinimum(1.0, np.ones(3))
        widget.update_A(first)
        widget.update_B(second)
        self.assertTrue(widget.ui.btn_compute.isEnabled())
        widget.update_B(first)
        self.assertFalse(widget.ui.btn_compute.isEnabled())
        widget.update_A(None)
        widget.update_B(None)
        self.assertEqual(widget.ui.lineEdit_A.text(), "")
        self.assertEqual(widget.ui.lineEdit_B.text(), "")
        self.assertFalse(widget.ui.btn_compute.isEnabled())

    def test_rates_reject_invalid_temperature_before_workers(self):
        widget = self.keep(RateWidget(self.system, self.database))
        widget.update_A(self.database.addMinimum(0.0, np.zeros(3)))
        widget.update_B(self.database.addMinimum(1.0, np.ones(3)))
        with patch.object(widget, "_compute_thermodynamic_info") as start:
            for temperature in ("bad", "0", "-1", "nan", "inf"):
                with self.subTest(temperature=temperature):
                    widget.ui.lineEdit_T.setText(temperature)
                    widget.compute_rates()
                    start.assert_not_called()
                    self.assertTrue(widget.ui.label_status.text())

    def test_rates_compute_both_directions_and_report_disconnected_network(self):
        widget = self.keep(RateWidget(self.system, self.database))
        first = self.database.addMinimum(0.0, np.zeros(3))
        second = self.database.addMinimum(1.0, np.ones(3))
        ts = self.database.addTransitionState(2.0, np.full(3, 2.0), first, second)
        for point in (first, second, ts):
            point.fvib, point.pgorder = 0.0, 1
        widget.update_A(first)
        widget.update_B(second)
        widget.transition_states = [ts]
        widget._compute_rates()
        result = widget.ui.textBrowser.toPlainText()
        self.assertIn("rate [1] -> [2]", result)
        self.assertIn("rate [2] -> [1]", result)
        self.assertEqual(widget.ui.label_status.text(), "")
        widget.transition_states = []
        widget._compute_rates()
        self.assertTrue(widget.ui.label_status.text())


if __name__ == "__main__":
    unittest.main()
