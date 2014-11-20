from PyQt4 import QtCore, QtGui, Qt

class NumberStandardItem(Qt.QStandardItem):
    def __init__(self, n):
        self.n = n
        Qt.QStandardItem.__init__(self, str(self.n))
    
    def __lt__(self, item2):
        return self.n < item2.n


class MinimumStandardItem(Qt.QStandardItem):
    """defines an item to populate the lists of minima in the GUI
    
    These items will be collected in a QStandardItemModel and viewed using
    QListView
    """
    def __init__(self, minimum):
        text="%.4f"% minimum.energy
        super(MinimumStandardItem, self).__init__(text)
        self.minimum = minimum
    def __lt__(self, item2):
        # sort the energies in the list lowest to highest
        return self.minimum.energy < item2.minimum.energy

class TransitionStateStandardItem(Qt.QStandardItem):
    """defines an item to populate the lists of transition states in the GUI
    
    These items will be collected in a QStandardItemModel and viewed using
    QListView
    """
    def __init__(self, ts):
#        text="%.4f (%d<-%d->%d)"%(ts.energy, ts._minimum1_id, ts._id, ts._minimum2_id)
        text="%.4f"% ts.energy
        super(TransitionStateStandardItem, self).__init__(text)
        self.ts = ts
    def __lt__(self, item2):
        # sort the energies in the list lowest to highest
        return self.ts.energy < item2.ts.energy

    def __getattr__(self, name):
        """return the transition state if the minimum is asked for"""
        if name == "minimum":
            return self.ts
        return super(TransitionStateStandardItem, self).__getattr__(name)

class MinimumStandardItemModel(Qt.QStandardItemModel):
    """a class to manage the list of minima for display in the gui"""
    def __init__(self, nmax=None):
        super(MinimumStandardItemModel, self).__init__()
        self.nmax = nmax # the maximum number of minima
        self.issued_warning = False
        self._minimum_to_item = dict()
        
        self.setColumnCount(2)
        
    def headerData(self, section, orientation, role):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                if section == 0:
                    return "Energy"
                elif section == 1:
                    return "ID"
                elif section == 2:
                    return "Min1"
                elif section == 3:
                    return "Min2 "

    def set_nmax(self, nmax):
        self.nmax = nmax
    
    def item_from_minimum(self, minimum):
        return self._minimum_to_item[minimum]
    
    
    def minimum_from_index(self, index):
        item = self.item(index.row())
        return item.minimum
    
    def minimum_from_selection(self, selection):
        try:
            index = selection.indexes()[0]
        except IndexError:
            return None
        return self.minimum_from_index(index)

    def addMinimum(self, m):
        item = MinimumStandardItem(m)
        iditem = NumberStandardItem(item.minimum._id)
        self.appendRow([item, iditem])
    
    def addMinima(self, mlist):
        if self.nmax is not None:
            if len(mlist) > self.nmax:
                mlist.sort(key=lambda m: m.energy)
                mlist = mlist[:self.nmax]
                print "warning: limiting the number of minima displayed in the gui to", self.nmax
        for m in mlist:
            self.addMinimum(m)

    
    def appendRow(self, items, *args, **kwargs):
        Qt.QStandardItemModel.appendRow(self, items)
        mitem = items[0]
#        print "adding minimum", mitem.minimum, mitem.minimum._id
        self._minimum_to_item[mitem.minimum] = mitem
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

class MinimumSortFilterProxyModel(Qt.QSortFilterProxyModel):
    def minimum_from_selection(self, selection):
        source_selection = self.mapSelectionToSource(selection)
        return self.sourceModel().minimum_from_selection(source_selection)
    def minimum_from_index(self, index):
        source_index = self.mapToSource(index)
        return self.sourceModel().minimum_from_index(source_index)

    def lessThan(self, index1, index2):
        item1 = self.sourceModel().itemFromIndex(index1)
        item2 = self.sourceModel().itemFromIndex(index2)
        return item1 < item2
        
    def index_from_minimum(self, minimum):
        source_model = self.sourceModel()
        item = source_model.item_from_minimum(minimum)
        index = source_model.indexFromItem(item)
        return self.mapFromSource(index)

class TransitionStateStandardItemModel(MinimumStandardItemModel):
    def __init__(self, nmax=None):
        MinimumStandardItemModel.__init__(self, nmax=nmax)
    
        self.setColumnCount(4)
    
    def addTS(self, ts):
        item = TransitionStateStandardItem(ts)
        iditem = NumberStandardItem(ts._id)
        m1item = NumberStandardItem(ts.minimum1._id)
        m2item = NumberStandardItem(ts.minimum2._id)
        self.appendRow([item, iditem, m1item, m2item])

    def addTransitionStates(self, tslist):
        if self.nmax is not None:
            if len(tslist) > self.nmax:
                tslist.sort(key=lambda ts: ts.energy)
                tslist = tslist[:self.nmax]
                print "warning: limiting the number of transition states displayed in the gui to", self.nmax
        for ts in tslist:
            self.addTS(ts)

class SaveCoordsAction(QtGui.QAction):
    def __init__(self, minimum, parent=None):
        super(SaveCoordsAction, self).__init__("save coords", parent)
        self.parent = parent
        self.minimum = minimum
        self.triggered.connect(self.__call__)

    def __call__(self, val):
        filename = QtGui.QFileDialog.getSaveFileName(self.parent, 'Save coords to', '.')
        if len(filename) > 0:
            print "saving coords to file", filename
            with open(filename, "w") as fout:
                fout.write("# id " + str(self.minimum._id) + " energy " + str(self.minimum.energy) + "\n")
                for x in self.minimum.coords:
                    fout.write( str(x) + "\n")
    


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
        
        # set up the sort/filter proxies which allow the lists to be
        # sorted separately
        self.mproxy_main = MinimumSortFilterProxyModel()
        self.mproxy_main.setSourceModel(self.minima_list_model)
        self.mproxy_1 = MinimumSortFilterProxyModel()
        self.mproxy_1.setSourceModel(self.minima_list_model)
        self.mproxy_2 = MinimumSortFilterProxyModel()
        self.mproxy_2.setSourceModel(self.minima_list_model)
        
        # use QListView objects to view the lists of minima in the gui
        self.ui.list_minima_main.setModel(self.mproxy_main)
#        self.ui.list_minima_main.set

        self.ui.listMinima1.setModel(self.mproxy_1)
        self.ui.listMinima2.setModel(self.mproxy_2)

        # connect to the signals of the QSelectionModel objects of the list view's   
        self.ui.list_minima_main.selectionModel().selectionChanged.connect(self.on_list_minima_main_selectionChanged)
        self.ui.listMinima1.selectionModel().selectionChanged.connect(self.on_listMinima1_selectionChanged)
        self.ui.listMinima2.selectionModel().selectionChanged.connect(self.on_listMinima2_selectionChanged)
        
        # set up the list of transition states in the gui
        self.ts_list_model = TransitionStateStandardItemModel()
        self.ui.list_TS.setModel(self.ts_list_model)
        self.ui.list_TS.selectionModel().selectionChanged.connect(self.on_list_TS_selectionChanged)

        # add actions
        self.ui.list_minima_main.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.list_minima_main.connect(self.ui.list_minima_main, 
                                         QtCore.SIGNAL("customContextMenuRequested(QPoint)"), 
                                                       self.list_view_on_context)
        self.ui.list_TS.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.list_TS.connect(self.ui.list_TS, 
                                         QtCore.SIGNAL("customContextMenuRequested(QPoint)"), 
                                                       self.transition_state_on_context)
    
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

    def transition_state_on_context(self, point):
        view = self.ui.list_TS
        index = view.indexAt(point)
        ts = self.ts_list_model.minimum_from_index(index)
        
        # create the menu
        menu = QtGui.QMenu("list menu", self.parent)
        
        action1 = SaveCoordsAction(ts, parent=self.parent)
        menu.addAction(action1)
        
        def prepare_in_connect(val):
            print "value", val
            self.parent._SelectMinimum1(ts.minimum1)
            print "selected minimum 1"
            self.parent._SelectMinimum2(ts.minimum2)
        action2 = QtGui.QAction("show in connect tab", self.parent)
        action2.triggered.connect(prepare_in_connect)
        menu.addAction(action2)

        # show the context menu
        menu.exec_(view.mapToGlobal(point))

    def list_view_on_context(self, point):
        view = self.ui.list_minima_main
        index = view.indexAt(point)
        
        minimum = self.mproxy_main.minimum_from_index(index)
        
        # create the menu
        menu = QtGui.QMenu("list menu", self.parent)
        
        action1 = SaveCoordsAction(minimum, parent=self.parent)     
        menu.addAction(action1)

        menu.exec_(view.mapToGlobal(point))
    
    def on_list_minima_main_selectionChanged(self, new, old):
        minimum = self.mproxy_main.minimum_from_selection(new)
        if minimum is not None:
            self.parent.SelectMinimum(minimum, set_selected=False)

    def on_listMinima1_selectionChanged(self, new, old):
        """when an item in the first list in the connect tab is selected"""
        minimum = self.mproxy_1.minimum_from_selection(new)
        if minimum is not None:
            self.parent._SelectMinimum1(minimum, set_selected=False)
    
    def on_listMinima2_selectionChanged(self, new, old):
        """when an item in the second list in the connect tab is selected"""
        minimum = self.mproxy_2.minimum_from_selection(new)
        if minimum is not None:
            self.parent._SelectMinimum2(minimum, set_selected=False)

    def _select_main(self, minimum):
        """set the minimum as selected in the basinhopping tab"""
        # I'm surprised we don't need to catch exceptions if the minimum is not in the 
        # model (e.g. if the maximum list length is exceeded)
        index = self.mproxy_main.index_from_minimum(minimum)
        self.ui.list_minima_main.setCurrentIndex(index)

    def _select1(self, minimum):
        """set the minimum as selected in the first list of the connect tab"""
        index = self.mproxy_1.index_from_minimum(minimum)
        self.ui.listMinima1.setCurrentIndex(index)

    def _select2(self, minimum):
        """set the minimum as selected in the second list of the connect tab"""
        index = self.mproxy_2.index_from_minimum(minimum)
        self.ui.listMinima2.setCurrentIndex(index)


    def on_list_TS_selectionChanged(self, new, old):
#        index = new.indexes()[0]
        ts = self.ts_list_model.minimum_from_selection(new)
        self.ts_selected = ts
        self.parent.show_TS(ts)

    def get_selected_ts(self):
        ts = self.ts_selected
        if ts is None:
            raise Exception("you must select a transition state first")
        return ts
    
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

    def resize_columns_minima(self):
#        print "resizing minima columns"
        self.ui.list_minima_main.resizeColumnsToContents()
        self.ui.listMinima1.resizeColumnsToContents()
        self.ui.listMinima2.resizeColumnsToContents()
    def resize_columns_ts(self):
#        print "resizing ts columns"
        self.ui.list_TS.resizeColumnsToContents()

    def NewMinimum(self, minimum, sort_items=True):
        """ add a new minimum to the system """
        try:
            len(minimum)
            minima = minimum
        except TypeError:
            minima = [minimum]
        self.minima_list_model.addMinima(minima)

        if sort_items:
            self._sort_minima()
        if not hasattr(self, "_minima_columns_resized"):
            self.resize_columns_minima()
            self._minima_columns_resized = True
            

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
        self.ts_list_model.addTransitionStates(tslist)
#        for ts in tslist:
##            tsitem = TransitionStateStandardItem(ts)
#            self.ts_list_model.addTS(ts)
        if sort:
            self._sort_ts()
            
        if not hasattr(self, "_ts_columns_resized"):
            self.resize_columns_ts()
            self._ts_columns_resized = True
