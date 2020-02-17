


def raise_invalid_sampling_method(sampling):
    err_msg = ('Unsupported sampling method, {}. Currently, descmap only '
               'supports "lhs" and "linear".'.format(sampling))
    raise ValueError(err_msg)
