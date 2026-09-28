"""Background-worker regressions, using real processes, pipes and Qt timers."""

import io
import multiprocessing as mp
import os
import sys
import time
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("PyQt5")
pytest.importorskip("OpenGL")

from pele.gui.bhrunner import BHManager, BHRunner
from pele.gui.connect_all import ConnectAllDialog
from pele.gui.connect_explorer_dlg import ConnectExplorerDialog, _TSListItem
from pele.gui.connect_run_dlg import ConnectViewer
from pele.gui.double_ended_connect_runner import DECRunner
from pele.storage import Database


class WorkerSystem:
    def __init__(self, delay=0, fail=False):
        self.params = SimpleNamespace(gui=SimpleNamespace(_sort_lists=True))
        self.delay = delay
        self.fail = fail
        self.ready = mp.Event()

    def create_database(self):
        return Database()

    def get_basinhopping(self, database, outstream):
        self.database = database
        self.step = 0
        return self

    def run(self, nsteps):
        self.ready.set()
        time.sleep(self.delay)
        self.step += 1
        self.database.addMinimum(float(self.step), np.array([self.step]))

    def get_double_ended_connect(self, min1, min2, db, fresh_connect):
        self.ready.set()
        if self.fail:
            raise RuntimeError("intentional worker failure")
        self.min1, self.min2 = min1, min2
        self.database = db
        self.graph = SimpleNamespace(areConnected=lambda a, b: True)
        return self

    def connect(self):
        time.sleep(self.delay)
        self.database.addTransitionState(
            2., np.array([.5]), self.min1, self.min2,
            eigenval=-1., eigenvec=np.array([1.]),
        )

    def returnPath(self):
        return [self.min1, self.min2], [0., 1.], [0., 1.]

    def smooth_path(self, path):
        raise AssertionError("smoothing was explicitly disabled")


def wait_for(qapp, predicate, timeout=5):
    deadline = time.monotonic() + timeout
    while not predicate() and time.monotonic() < deadline:
        qapp.processEvents()
        time.sleep(.005)
    assert predicate(), "worker did not finish before the deadline"


def cleanup(runner):
    process = getattr(runner, "bhprocess", None) or getattr(runner, "decprocess", None)
    if process is not None:
        if process.is_alive():
            process.kill()
        process.join()
    if hasattr(runner, "refresh_timer"):
        runner.refresh_timer.stop()


def make_connect_runner(system, **kwargs):
    db = Database()
    min1 = db.addMinimum(0., np.array([0.]))
    min2 = db.addMinimum(1., np.array([1.]))
    return DECRunner(system, db, min1, min2, **kwargs)


def test_basin_hopping_drains_results_before_manager_removes_worker(qapp):
    db = Database()
    manager = BHManager(WorkerSystem(), db)
    manager.start_worker(nsteps=8)
    worker = manager.workers[0]
    try:
        worker.bhprocess.join(5)
        assert not worker.is_alive()
        manager._remove_dead()
        assert db.number_of_minima() == 8 or worker in manager.workers
        wait_for(qapp, lambda: db.number_of_minima() == 8)
        assert manager.number_of_workers() == 0
        assert not worker.refresh_timer.isActive()
    finally:
        manager.refresh_timer.stop()
        cleanup(worker)


def test_basin_hopping_stop_does_not_block_the_gui(qapp):
    system = WorkerSystem(delay=.5)
    runner = BHRunner(system, Database(), nsteps=10)
    runner.start()
    try:
        assert system.ready.wait(5)
        started = time.monotonic()
        runner.kill()
        assert time.monotonic() - started < .2
        wait_for(qapp, lambda: not runner.refresh_timer.isActive())
        assert not runner.is_alive()
    finally:
        cleanup(runner)


def test_connect_worker_crash_finishes_once_and_reports_error(qapp):
    log = io.StringIO()
    runner = make_connect_runner(WorkerSystem(fail=True), outstream=log)
    finished = []
    callback = lambda: finished.append(runner.is_running)
    runner.on_finished.connect(callback)
    runner.start()
    try:
        wait_for(qapp, lambda: bool(finished))
        assert finished == [False]
        assert not runner.success
        assert not runner.refresh_timer.isActive()
        assert "intentional worker failure" in log.getvalue()
        runner.poll()
        assert finished == [False]
    finally:
        cleanup(runner)


def test_connect_worker_respects_disabled_smoothing_and_transfers_results(qapp):
    log = io.StringIO()
    runner = make_connect_runner(
        WorkerSystem(), outstream=log, return_smoothed_path=False,
    )
    runner.start()
    try:
        wait_for(qapp, lambda: not runner.is_running)
        assert runner.success
        assert runner.database.number_of_transition_states() == 1
        assert len(runner.newminima) == 2
        assert len(runner.newtransition_states) == 1
        assert "smoothing was explicitly disabled" not in log.getvalue()
    finally:
        cleanup(runner)


def test_unstarted_workers_can_be_stopped(qapp):
    bh = BHRunner(WorkerSystem(), Database())
    assert not bh.is_alive()
    bh.kill()
    dec = make_connect_runner(WorkerSystem())
    dec.terminate_early()
    assert not dec.is_running


def test_close_waits_for_children_and_drains_results(qapp):
    manager = BHManager(WorkerSystem(), Database())
    manager.start_worker(nsteps=100)
    bh = manager.workers[0]
    dec = make_connect_runner(WorkerSystem(delay=.1), return_smoothed_path=False)
    dec.start()
    try:
        assert dec.system.ready.wait(5)
        manager.kill_all_workers(wait=True)
        dec.terminate_early(wait=True)
        assert not bh.is_alive()
        assert not dec.decprocess.is_alive()
        assert not dec.is_running
        assert not bh.refresh_timer.isActive()
        assert not dec.refresh_timer.isActive()
    finally:
        manager.refresh_timer.stop()
        cleanup(bh)
        cleanup(dec)


@pytest.mark.parametrize("dialog_class", [ConnectViewer, ConnectAllDialog])
def test_connect_dialog_stop_and_close_before_first_job(qapp, dialog_class):
    dialog = dialog_class(WorkerSystem(), Database(), app=qapp)
    try:
        dialog.on_actionKill_triggered(True)
        dialog.close()
    finally:
        dialog.hide()


def test_connect_all_views_are_safe_before_first_result(qapp, monkeypatch):
    errors = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: errors.append(exc[1]))
    dialog = ConnectAllDialog(WorkerSystem(), Database(), app=qapp)
    try:
        dialog.show()
        dialog.ui.actionEnergy.setChecked(True)
        dialog.ui.action3D.setChecked(True)
        assert errors == []
    finally:
        dialog.hide()


def test_connect_explorer_empty_pushoff_and_cleared_selection(qapp, monkeypatch):
    errors = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: errors.append(exc[1]))
    dialog = ConnectExplorerDialog(WorkerSystem(), qapp)
    try:
        dialog.show_pushoff_path()
        dialog.show_TS_path()
        path = np.array([[0., 0., 0.]])
        item = _TSListItem(0, path, ["ts"], path, ["pushoff"])
        dialog.ts_list.addItem(item)
        # Highlighting needs a NEB plot; selection-cleared is independent of it.
        dialog.ts_list.blockSignals(True)
        dialog.ts_list.setCurrentItem(item)
        dialog.ts_list.blockSignals(False)
        dialog.reset()
        assert errors == []
    finally:
        dialog.close()


class SmoothedWorkerSystem(WorkerSystem):
    def smooth_path(self, path):
        return path


class ExitingWorkerSystem(WorkerSystem):
    def connect(self):
        os._exit(17)


class SlowSmoothingSystem(WorkerSystem):
    def __init__(self):
        super().__init__()
        self.smoothing = mp.Event()

    def smooth_path(self, path):
        self.smoothing.set()
        time.sleep(5)
        return path


def test_stopping_during_smoothing_does_not_publish_incomplete_success(qapp):
    runner = make_connect_runner(SlowSmoothingSystem())
    runner.start()
    try:
        assert runner.system.smoothing.wait(5)
        runner.terminate_early(wait=True)
        assert not runner.success
        assert not runner.is_running
        assert runner.database.number_of_transition_states() == 1
    finally:
        cleanup(runner)


def test_connect_worker_abnormal_exit_reports_failure(qapp):
    log = io.StringIO()
    runner = make_connect_runner(ExitingWorkerSystem(), outstream=log)
    runner.start()
    try:
        wait_for(qapp, lambda: not runner.is_running)
        assert not runner.success
        assert "17" in log.getvalue()
    finally:
        cleanup(runner)


def test_connect_all_summary_measures_elapsed_wall_time(qapp):
    runner = make_connect_runner(SmoothedWorkerSystem(delay=.2))
    dialog = ConnectAllDialog(runner.system, runner.database, app=qapp)
    dialog.combobox.setCurrentText("global min")
    dialog.ui.actionPause.setChecked(True)
    dialog.start()
    try:
        wait_for(qapp, lambda: bool(dialog.connect_summary.attempts))
        attempt = dialog.connect_summary.attempts[0]
        assert attempt.success
        assert attempt.time >= .2
    finally:
        dialog.close()
        cleanup(dialog.decrunner)


def test_connect_all_restart_waits_for_previous_child(qapp):
    runner = make_connect_runner(SmoothedWorkerSystem(delay=1))
    dialog = ConnectAllDialog(runner.system, runner.database, app=qapp)
    dialog.combobox.setCurrentText("global min")
    dialog.start()
    previous = dialog.decrunner
    try:
        assert runner.system.ready.wait(5)
        dialog.on_actionKill_triggered(True)
        dialog.start()
        assert dialog.decrunner is previous
    finally:
        dialog.close()
        cleanup(previous)
        cleanup(dialog.decrunner)


def test_connect_restart_clears_previous_results(qapp):
    runner = make_connect_runner(SmoothedWorkerSystem())
    runner.start()
    try:
        wait_for(qapp, lambda: not runner.is_running)
        assert runner.success
        np.testing.assert_array_equal(runner.smoothed_path, [[0.], [1.]])
        runner.system.fail = True
        runner.start()
        assert not runner.success
        assert not runner.newminima
        assert not runner.newtransition_states
        wait_for(qapp, lambda: not runner.is_running)
        assert not runner.success
    finally:
        cleanup(runner)


@pytest.mark.parametrize("nmin", [0, 1, 2])
def test_connect_all_stops_when_no_pairs_remain(qapp, nmin):
    db = Database()
    for i in range(nmin):
        db.addMinimum(float(i), np.array([i]))
    dialog = ConnectAllDialog(WorkerSystem(), db, app=qapp)
    dialog.combobox.setCurrentText("global min")
    if nmin == 2:
        dialog.connect_manager.get_connect_job(strategy="gmin")
    try:
        dialog.start()
        assert not dialog.is_running
        assert dialog.decrunner is None
        assert dialog.ui.actionPause.isChecked()
        assert dialog.textEdit_summary.toPlainText()
    finally:
        dialog.close()
