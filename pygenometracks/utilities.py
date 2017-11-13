import sys


def to_string(s):
    """
    This takes care of python2/3 differences
    """
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):
        if sys.version_info[0] == 2:
            return str(s)
        return s.decode('ascii')
    if isinstance(s, list):
        return [to_string(x) for x in s]
    return s


def to_bytes(s):
    """
    Like toString, but for functions requiring bytes in python3
    """
    if sys.version_info[0] == 2:
        return s
    if isinstance(s, bytes):
        return s
    if isinstance(s, str):
        return bytes(s, 'ascii')
    if isinstance(s, list):
        return [to_bytes(x) for x in s]
    return s
