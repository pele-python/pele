
from PyQt4 import QtCore, QtGui, Qt


class MinimumStandardItem(Qt.QStandardItem):
    """defines an item to populate the lists of minima in the GUI
    
    These items will be collected in a QStandardItemModel and viewed using
    QListView
    """
    def __init__(self, minimum):
        text="%.4f (%d)"%(minimum.energy, minimum._id)
        super(MinimumStandardItem, self).__init__(text)
        self.minimum = minimum
    def __lt__(self, item2):
        #sort the energies in the list lowest to highest
        return self.minimum.energy < item2.minimum.energy

class TransitionStateStandardItem(Qt.QStandardItem):
    """defines an item to populate the lists of minima in the GUI
    
    These items will be collected in a QStandardItemModel and viewed using
    QListView
    """
    def __init__(self, ts):
        text="%.4f (%d<-%d->%d)"%(ts.energy, ts._minimum1_id, ts._id, ts._minimum2_id)
        super(TransitionStateStandardItem, self).__init__(text)
        self.ts = ts
    def __lt__(self, item2):
        #sort the energies in the list lowest to highest
        return self.ts.energy < item2.ts.energy

    def __getattr__(self, name):
        """return the transition state if the minimum is asked for"""
        if name == "minimum":
            return self.ts
        return super(TransitionStateStandardItem, self).__getattr__(name)

class MinimumStandardItemModel(Qt.QStandardItemModel):
    """a class to manage the list of minima for display in the gui"""
    def __init__(self, nmax=None, **kwargs):
        super(MinimumStandardItemModel, self).__init__(**kwargs)
        self.nmax = nmax # the maximum number of minima
        self.issued_warning = False
        self._minimum_to_item = dict()

    def set_nmax(self, nmax):
        self.nmax = nmax
    
    def item_from_minimum(self, minimum):
        return self._minimum_to_item[minimum]
          
    def appendRow(self, item, *args, **kwargs):
        self._minimum_to_item[item.minimum] = item
        Qt.QStandardItemModel.appendRow(self, item, *args, **kwargs)
        look_back = min(10, self.nmax)
        if self.nmax is not None:
            nrows = self.rowCount()
            if nrows > self.nmax:
                if not self.issued_warning:
                    print "warning: limiting the number of minima displayed in the gui to", self.nmax
                    self.issued_warning = True
                # choose an item to remove from the list.  we can't usume it's totally sorted
                # because it might have been a while since it was last sorted
                candidates = [(self.item(r).minimum.energy, r) for r in xrange(nrows-look_back,nrows)]
                toremove = max(candidates)
                self.takeRow(toremove[1])


class ListViewManager(object):
    def __init__(self, parent):
        self.parent = parent
        self.ui = self.parent.ui
    
        # set up the sorting of minima
        self.need_sorting = False
        self.need_sorting_minima = False
        self.need_sorting_ts = False
        self._sort_timer = QtCore.QTimer()
        self._sort_timer.timeout.connect(self._delayed_sort)
        
        self.ts_selected = None
        
        # keep track of minima using a QStandardItemModel 
        self.minima_list_model = MinimumStandardItemModel()
        # use QListView objects to view the lists of minima in the gui
        self.ui.list_minima_main.setModel(self.minima_list_model)
#        self.ui.list_minima_main.set
        self.ui.listMinima1.setModel(self.minima_list_model)
        self.ui.listMinima2.setModel(self.minima_list_model)

        # connect to the signals of the QSelectionModel objects of the list view's   
        self.ui.list_minima_main.selectionModel().selectionChanged.connect(self.on_list_minima_main_selectionChanged)
        self.ui.listMinima1.selectionModel().selectionChanged.connect(self.on_listMinima1_selectionChanged)
        self.ui.listMinima2.selectionModel().selectionChanged.connect(self.on_listMinima2_selectionChanged)
        
        # set up the list of transition states in the gui
        self.ts_list_model = MinimumStandardItemModel()
        self.ui.list_TS.setModel(self.ts_list_model)
        self.ui.list_TS.selectionModel().selectionChanged.connect(self.on_list_TS_selectionChanged)

        # add actions
        self.ui.list_minima_main.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.list_minima_main.connect(self.ui.list_minima_main, 
                                         QtCore.SIGNAL("customContextMenuRequested(QPoint)"), 
                                                       self.list_view_on_context)

    
    def finish_setup(self):
        """this must be called after NewSystem() is called"""
                # determine the maximum number of minima to keep in the lists.
        # this must be done after NewSystem() is called
        try:
            nmax = self.parent.system.params.gui.list_nmax
            self.minima_list_model.set_nmax(nmax)
            self.ts_list_model.set_nmax(nmax)
        except AttributeError:
            pass
    
    def clear(self):
        self.minima_list_model.clear()
        self.ts_list_model.clear()

    def list_view_on_context(self, point):
        view = self.ui.list_minima_main
        index = view.indexAt(point)
        item = self.minima_list_model.itemFromIndex(index)
        minimum = item.minimum
        
        # create the menu
        menu = QtGui.QMenu("list menu", self.parent)
        
        def save_coords(val):
            filename = QtGui.QFileDialog.getSaveFileName(self.parent, 'Save file name', '.')
            if len(filename) > 0:
                print "saving coords to file", filename
                with open(filename, "w") as fout:
                    fout.write("# id " + str(minimum._id) + " energy " + str(minimum.energy) + "\n")
                    for x in minimum.coords:
                        fout.write( str(x) + "\n")
                
        
        action1 = QtGui.QAction("save coords", self.parent)
        action1.triggered.connect(save_coords)
        menu.addAction(action1)
#        self.ui.list_minima_main.map
#        menu.a
        # show the context menu
        menu.exec_(view.mapToGlobal(point))
    
    def on_list_minima_main_selectionChanged(self, new, old):
        index = new.indexes()[0]
        item = self.minima_list_model.itemFromIndex(index)
        minimum = item.minimum
        self.parent.SelectMinimum(minimum, set_selected=False)

    def on_listMinima1_selectionChanged(self, new, old):
        """when an item in the first list in the connect tab is selected"""
        index = new.indexes()[0]
        item = self.minima_list_model.itemFromIndex(index)
        minimum = item.minimum
        self.parent._SelectMinimum1(minimum, set_selected=False)

    
    def on_listMinima2_selectionChanged(self, new, old):
        """when an item in the second list in the connect tab is selected"""
        index = new.indexes()[0]
        item = self.minima_list_model.itemFromIndex(index)
        minimum = item.minimum
        self.parent._SelectMinimum2(minimum, set_selected=False)

    def _select_main(self, minimum):
        """set the minimum as selected in the basinhopping tab"""
        # I'm surprised we don't need to catch exceptions if the minimum is not in the 
        # model (e.g. if the maximum list length is exceeded)
        item = self.minima_list_model.item_from_minimum(minimum)
        index = self.minima_list_model.indexFromItem(item)
        self.ui.list_minima_main.setCurrentIndex(index)

    def _select1(self, minimum):
        """set the minimum as selected in the first list of the connect tab"""
        item = self.minima_list_model.item_from_minimum(minimum)
        index = self.minima_list_model.indexFromItem(item)
        self.ui.listMinima1.setCurrentIndex(index)

    def _select2(self, minimum):
        """set the minimum as selected in the second list of the connect tab"""
        item = self.minima_list_model.item_from_minimum(minimum)
        index = self.minima_list_model.indexFromItem(item)
        self.ui.listMinima2.setCurrentIndex(index)


    def on_list_TS_selectionChanged(self, new, old):
        index = new.indexes()[0]
        item = self.ts_list_model.item(index.row())
        ts = item.ts
        self.ts_selected = ts
        self.parent.show_TS(ts)


    def _sort_ts(self):
        self.need_sorting_ts = True
        self._sort_lists()
    
    def _sort_minima(self):
        self.need_sorting_minima = True
        self._sort_lists()
    
    def _sort_lists(self):
        """calling this function indicates that the lists need sorting
        
        since sorting can be a *huge* bottleneck for large lists we try to do it
        as infrequently as possible.  Here we wait several seconds too see
        if more minima are added before sorting the lists
        """
        self.need_sorting = True
        if not self._sort_timer.isActive():
            self._sort_timer.start(2000)

    def _delayed_sort(self):
        if not self.need_sorting:
            self.need_sorting_minima = False
            self.need_sorting_ts = False
            print "sorting not needed"
            return
        try:
            # the system params flag flag gui._sort_lists can optionally
            # be set to False to indicate not to sort the lists.  This is
            # useful if many minima are added at once, e.g. at the end of
            # a connect run.  The flag should be set to True when all the
            # minima are added. 
            s = self.parent.system.params.gui._sort_lists
#            print "self.system.params.gui._sort_lists", s
            if not s:
                # wait a few seconds and call this function again
                if not self._sort_timer.isActive():
                    self._sort_timer.start(2000)
                print "delaying sort"
                return
        except AttributeError:
            pass
        
#        print "sorting lists"
        if self.need_sorting_minima:
            self.minima_list_model.sort(0)
        if self.need_sorting_ts:
            self.ts_list_model.sort(0)
        self.need_sorting = False
        self.need_sorting_minima = False
        self.need_sorting_ts = False
        self._sort_timer.stop()
#        print "done sorting lists"

    def NewMinimum(self, minimum, sort_items=True):
        """ add a new minimum to the system """
        self.minima_list_model.appendRow(MinimumStandardItem(minimum))
        if sort_items:
            self._sort_minima()

    def RemoveMinimum(self, minimum):
        """remove a minimum from self.minima_list_model"""
        minid = minimum._id
        items = self.minima_list_model.findItems('*', QtCore.Qt.MatchWildcard)
        for i in items:
#            print "item", i.minimum._id, minid, minimum._id
            if i.minimum._id == minid:
#                print "taking item", i.minimum._id
                self.minima_list_model.takeRow(i.row())
                break

    def NewTS(self, ts, sort=True):
        """add new transition state, or list of transition states"""
        try:
            len(ts)
            is_iterable = True
        except TypeError:
            is_iterable = False
        if is_iterable:   
            tslist = ts
        else: 
            tslist = [ts]
        for ts in tslist:
            tsitem = TransitionStateStandardItem(ts)
            self.ts_list_model.appendRow(tsitem)
        if sort:
            self._sort_ts()