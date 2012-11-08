class GUISystem(object):
    def __init__(self):
        pass
    
    def initialize_database(self):
        pass
    
    def create_basinhopping(self):
        raise BaseException("not implemented")
    
    def draw(self, coordslinear, index):
        raise BaseException("not implemented")
    
    def set_database(self, database):
        self.database = database