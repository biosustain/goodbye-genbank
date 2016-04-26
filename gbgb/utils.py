def single(qualifiers, name, on_multiple_ignore=False):
    """

    :param qualifiers:
    :param bool on_multiple_ignore: ignore the qualifier if it has multiple values; otherwise returns the first item.
    :return: A single value for the qualifier; ``None`` if the qualifier has no value
    """
    if name not in qualifiers:
        return None
    if isinstance(qualifiers[name], (list, tuple, set)):
        if on_multiple_ignore and len(qualifiers[name]) > 1:
            return None
        if len(qualifiers[name]) == 0:
            return None
        return next(iter(qualifiers[name]))
    return qualifiers[name]


def as_set(qualifiers):
    if qualifiers is None:
        return set()
    if isinstance(qualifiers, set):
        return qualifiers
    if isinstance(qualifiers, (list, tuple)):
        return set(qualifiers)
    return {qualifiers}


def as_list(qualifiers):
    if qualifiers is True or qualifiers is None:
        return []
    if isinstance(qualifiers, list):
        return qualifiers
    if isinstance(qualifiers, (tuple, set)):
        return list(qualifiers)
    return [qualifiers]
