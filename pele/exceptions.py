    
class SignError(ValueError):
    """
    this is an exception to be raised if a parameter is passed to a method
    with the wrong sign (e.g. Ecriterion in DontLeaveBasin must be positive)
    """
    pass
