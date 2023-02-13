"""Functionality to raise sampling errors"""
def raise_invalid_sampling_method(sampling):
    """Raises ValueError due to invalid sample method being specified.
    
    Parameters
    ----------
        sampling : str
            Invalid sampling method
    """
    err_msg = ('Unsupported sampling method, {}. Currently, descmap only '
               'supports "lhs" and "linear".'.format(sampling))
    raise ValueError(err_msg)
