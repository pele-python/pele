
class LineSearchError(Exception):
    """
    The exception to return if there is a problem with the line search.  This is also
    used in functions like takeStepNoLineSearch which effectively play the role
    of a line search
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value
