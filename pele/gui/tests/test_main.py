from unittest.mock import patch

import numpy as np
import pytest

pytest.importorskip("PyQt5")
pytest.importorskip("OpenGL")
from PyQt5 import QtCore, QtWidgets

from pele.gui._list_views import MinimumStandardItemModel
from pele.gui.run import MainGUI
from pele.storage import Database, Minimum
from pele.systems import LJCluster


@pytest.fixture
def window(qapp):
    with patch("OpenGL.GLUT.glutInit"):
        window = MainGUI(qapp, LJCluster(3))
    yield window
    window.list_manager._sort_timer.stop()
    window.close()


def add_landscape(db):
    coords = np.array([0., 0., 0., 1., 0., 0., 0., 1., 0.])
    m1 = db.addMinimum(-3., coords)
    m2 = db.addMinimum(-2., coords * 1.1)
    ts = db.addTransitionState(-1., coords * 1.05, m1, m2)
    return m1, m2, ts


def test_database_removals_clear_rows_and_selection(window):
    db = window.system.database
    m1, m2, ts = add_landscape(db)
    window.SelectMinimum(m1)
    window._SelectMinimum1(m1)
    window._SelectMinimum2(m2)
    model = window.list_manager.ts_list_model
    window.ui.list_TS.setCurrentIndex(model.index(0, 0))
    db.removeMinimum(m1)
    assert db.number_of_transition_states() == 0
    assert model.rowCount() == 0
    assert window.list_manager.minima_list_model.rowCount() == 1
    assert window.list_manager.ts_selected is None
    # Qt may select the next row after removal; the deleted object must be gone.
    assert window.minima_selection.minimum1 is not m1
    assert window.ui.ogl_main.minima[1] is not m1


def test_database_switch_detaches_old_database_and_resets_views(window):
    old = window.system.database
    m1, m2, ts = add_landscape(old)
    window.SelectMinimum(m1)
    window._SelectMinimum1(m1)
    window._SelectMinimum2(m2)
    new = Database()
    window.connect_db(new)
    assert window.minima_selection.minimum1 is None
    assert window.minima_selection.minimum2 is None
    assert window.ui.ogl_main.coords[1] is None
    old.addMinimum(-4., np.zeros(9))
    assert window.list_manager.minima_list_model.rowCount() == 0
    new.addMinimum(-5., np.ones(9))
    assert window.list_manager.minima_list_model.rowCount() == 1


def test_failed_database_switch_preserves_current_rows(window, tmp_path):
    old = window.system.database
    add_landscape(old)
    with pytest.raises(Exception):
        window.connect_db(str(tmp_path / "missing" / "database.sqlite"))
    assert window.system.database is old
    assert window.list_manager.minima_list_model.rowCount() == 2


def test_unselected_actions_are_disabled_and_stop_is_safe(window):
    window.on_btn_stop_basinhopping_clicked(False)
    for button in (window.ui.btnAlign, window.ui.btnConnect,
                   window.ui.btnNEB, window.ui.pushNormalmodesMin,
                   window.ui.pushNormalmodesTS):
        assert not button.isEnabled()
    assert not window.ui.action_delete_minimum.isEnabled()


def test_empty_context_menu_is_harmless(window):
    window.list_manager.list_view_on_context(QtCore.QPoint(50, 50))
    window.list_manager.transition_state_on_context(QtCore.QPoint(50, 50))


def test_list_limit_keeps_lowest_energies_and_cleans_lookup(qapp):
    model = MinimumStandardItemModel(nmax=2)
    minima = [Minimum(e, [e]) for e in (1., 2., 3., -1.)]
    for i, minimum in enumerate(minima):
        minimum._id = i + 1
        model.addMinimum(minimum)
    assert sorted(model.item(i).minimum.energy for i in range(model.rowCount())) == [-1., 1.]
    assert set(model._minimum_to_item) == {minima[0], minima[3]}
    model.clear()
    assert not model._minimum_to_item
    assert model.minimum_from_index(QtCore.QModelIndex()) is None
    model.set_nmax(None)
    model.addMinimum(minima[0])
    assert model.rowCount() == 1


def test_list_limit_applies_before_initial_database_load(qapp):
    system = LJCluster(3)
    system.params.gui.list_nmax = 1
    db = Database()
    add_landscape(db)
    with patch("OpenGL.GLUT.glutInit"), patch.object(system, "create_database", return_value=db):
        window = MainGUI(qapp, system)
    try:
        assert window.list_manager.minima_list_model.rowCount() == 1
    finally:
        window.close()


def test_export_coordinates(window, tmp_path):
    from pele.gui._list_views import SaveCoordsAction
    minimum = add_landscape(window.system.database)[0]
    filename = tmp_path / "coords.dat"
    action = SaveCoordsAction(minimum, window)
    with patch.object(QtWidgets.QFileDialog, "getSaveFileName", return_value=(str(filename), "")):
        action(False)
    np.testing.assert_array_equal(np.loadtxt(filename), minimum.coords)


def test_select_minimum_outside_display_limit(window):
    window.list_manager.minima_list_model.set_nmax(1)
    m1, m2, ts = add_landscape(window.system.database)
    window.SelectMinimum(m2)
    window._SelectMinimum1(m1)
    window._SelectMinimum2(m2)
    assert window.ui.ogl_main.minima[1] is m2
    assert window.minima_selection.minimum2 is m2
    assert window.ui.btnConnect.isEnabled()


def test_invalid_basin_hopping_steps_do_not_start_worker(window):
    with patch("pele.gui.run.BHManager") as manager, patch.object(QtWidgets.QMessageBox, "warning"):
        for value in ("abc", "-1", "0"):
            window.ui.lineEdit_bh_nsteps.setText(value)
            window.on_btn_start_basinhopping_clicked(False)
        assert not manager.called


def test_close_windows_closes_every_analysis_dialog(window):
    window.on_btn_rates_clicked(False)
    window.on_action_edit_params_triggered(False)
    rate, params = window.rate_viewer, window.paramsdlg
    window.on_btn_close_all_clicked(False)
    assert not rate.isVisible()
    assert not params.isVisible()
    assert not hasattr(window, "rate_viewer")
    assert not hasattr(window, "paramsdlg")


def test_reopening_parameters_does_not_orphan_an_earlier_dialog(window):
    window.on_action_edit_params_triggered(False)
    first = window.paramsdlg
    window.on_action_edit_params_triggered(False)
    window.on_btn_close_all_clicked(False)
    assert not first.isVisible()


def test_merge_updates_transition_state_endpoint_columns(window):
    m1, m2, ts = add_landscape(window.system.database)
    window.ui.list_TS.setCurrentIndex(window.list_manager.ts_list_model.index(0, 0))
    with patch.object(QtWidgets.QMessageBox, "question", return_value=QtWidgets.QMessageBox.Ok):
        window._merge_minima(m1, m2)
    model = window.list_manager.ts_list_model
    assert model.item(0, 2).n == m1.id()
    assert model.item(0, 3).n == m1.id()
    assert window.list_manager.ts_selected is None
    assert not window.ui.pushNormalmodesTS.isEnabled()
