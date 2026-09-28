"""Viewer regressions; run with QT_QPA_PLATFORM=offscreen when headless."""

from unittest.mock import patch

import matplotlib
import numpy as np
import pytest

pytest.importorskip("PyQt5")
pytest.importorskip("OpenGL")
from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtTest import QTest

from pele.gui.graph_viewer import GraphViewWidget
from pele.gui.show3d import Show3D
from pele.gui.show3d_with_slider import Show3DWithSlider
from pele.gui.ui.dgraph_dlg import DGraphWidget, minimum_energy_path
from pele.storage import Database


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def test_slider_updates_each_frame_once_and_clamps(app):
    viewer = Show3DWithSlider()
    frames = []
    callback = lambda i, sender: frames.append(i)
    viewer.on_frame_updated.connect(callback)
    path = np.arange(9).reshape(3, 3)
    viewer.setCoordsPath(path, frame=1, labels=["first", "middle", "last"])
    assert frames == [1]
    assert viewer.label.text() == "middle"
    np.testing.assert_array_equal(viewer.oglwgt.coords[1], path[1])
    viewer.showFrame(100)
    assert frames == [1, 2]
    assert viewer.slider.value() == 2
    assert viewer.label.text() == "last"
    viewer.setCoordsPath(path[:1], labels=["only"])
    assert frames == [1, 2, 0]
    viewer.close()


def test_single_frame_cannot_animate_and_stopping_resets_button(app):
    viewer = Show3DWithSlider()
    viewer.setCoordsPath(np.ones((1, 3)))
    viewer.start_animation()
    viewer._next_frame()
    assert not viewer.animate
    viewer.setCoordsPath(np.ones((3, 3)))
    viewer.ui.btn_animate.click()
    assert viewer.animate
    viewer.setCoords(np.zeros(3))
    assert not viewer.ui.btn_animate.isChecked()
    assert not viewer.animate
    viewer.close()


def test_empty_path_is_rejected_before_changing_viewer(app):
    viewer = Show3DWithSlider()
    with pytest.raises(ValueError):
        viewer.setCoordsPath(np.empty((0, 3)))
    viewer.close()


def test_resize_tolerates_zero_height(app):
    viewer = Show3D(None)
    # This verifies projection arguments without requiring a display server.
    with patch("pele.gui.show3d.GL"), patch("pele.gui.show3d.GLU") as glu:
        viewer.resizeGL(200, 0)
        assert np.isfinite(glu.gluPerspective.call_args.args[1])
    viewer.close()


def make_database():
    db = Database()
    minima = [db.addMinimum(float(i), np.array([i])) for i in range(3)]
    for m in minima:
        m.fvib, m.pgorder = 0., 1
    for u, v, energy in [(0, 1, 5.), (1, 2, 5.), (0, 2, 10.)]:
        ts = db.addTransitionState(energy, np.zeros(1), minima[u], minima[v])
        ts.fvib, ts.pgorder = 0., 1
    return db, minima


def test_connectivity_graph_handles_empty_database_without_app_argument(app):
    viewer = GraphViewWidget(Database())
    viewer.show_all()
    assert viewer.graph.number_of_nodes() == 0
    db, minima = make_database()
    viewer.database = db
    viewer.show_all()
    assert len(viewer._mimima_layout_list) == 3
    viewer.close()


def test_minimum_energy_path_uses_current_networkx(app):
    db, minima = make_database()
    viewer = GraphViewWidget(db, app=app)
    viewer.make_graph()
    assert minimum_energy_path(viewer.graph, minima[0], minima[2]) == minima
    viewer.close()


def test_zoomed_graph_keeps_boundary_nodes_visible(app):
    db, minima = make_database()
    viewer = GraphViewWidget(db, app=app)
    viewer.make_graph()
    viewer.make_graph_from([minima[0]])
    viewer.show_graph()
    assert len(viewer._boundary_points.get_offsets()) == 2
    assert len(viewer._boundary_points.get_edgecolors()) > 0
    viewer.close()


def test_disconnectivity_redraw_preserves_parameters(app):
    db, minima = make_database()
    viewer = DGraphWidget(db)
    viewer.rebuild_disconnectivity_graph()
    viewer.ui.lineEdit_nlevels.setText("5")
    viewer.redraw_disconnectivity_graph()
    assert viewer.params["nlevels"] == 5
    viewer.close()


def test_graph_committor_uses_all_transition_states(app):
    db, minima = make_database()
    viewer = GraphViewWidget(db, app=app)
    viewer.make_graph()
    viewer._color_by_committor(minima[0], minima[2])
    # Equal barriers from the middle minimum give equal hitting probabilities.
    assert viewer._minima_color_value(minima[1]) == pytest.approx(.5)
    viewer.close()


@pytest.mark.parametrize("method", ["_color_by_committor", "_layout_by_committor"])
def test_disconnectivity_committor_uses_all_transition_states(app, method):
    db, minima = make_database()
    viewer = DGraphWidget(db)
    viewer.rebuild_disconnectivity_graph()
    getattr(viewer, method)(minima[0], minima[2])
    color = viewer.dg.minimum_to_leave[minima[1]].data["colour"]
    assert color == pytest.approx(matplotlib.colormaps["winter"](.5))
    viewer.close()


@pytest.mark.parametrize("index", [0, 2])
def test_disconnectivity_mfpt_is_zero_at_selected_target(app, index):
    db, minima = make_database()
    viewer = DGraphWidget(db)
    viewer.rebuild_disconnectivity_graph()
    viewer._color_by_mfpt(minima[index])
    color = viewer.dg.minimum_to_leave[minima[index]].data["colour"]
    assert color == pytest.approx(matplotlib.colormaps["winter"](0.))
    viewer.close()


def test_animation_start_is_idempotent_and_stop_cancels_ticks(app):
    viewer = Show3DWithSlider()
    viewer.setCoordsPath(np.arange(9).reshape(3, 3))
    frames = []
    callback = lambda i, sender: frames.append(i)
    viewer.on_frame_updated.connect(callback)
    viewer.start_animation()
    viewer.start_animation()
    QTest.qWait(120)
    assert frames == [1]
    viewer.stop_animation()
    QTest.qWait(120)
    assert frames == [1]
    viewer.close()


@pytest.mark.parametrize("close_parent", [False, True])
def test_closing_viewer_or_parent_stops_animation(app, close_parent):
    parent = QWidget()
    viewer = Show3DWithSlider(parent)
    # Exercise widget visibility without requiring an OpenGL display context.
    viewer.oglwgt.hide()
    viewer.setCoordsPath(np.arange(9).reshape(3, 3))
    frames = []
    callback = lambda i, sender: frames.append(i)
    viewer.on_frame_updated.connect(callback)
    parent.show()
    viewer.start_animation()
    try:
        (parent if close_parent else viewer).close()
        QTest.qWait(120)
        assert frames == []
        assert not viewer.animate
        assert not viewer.ui.btn_animate.isChecked()
    finally:
        viewer.stop_animation()
        parent.close()
