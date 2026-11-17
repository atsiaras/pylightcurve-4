
class PyLCError(BaseException):
    """Base class of all PyLightcurve exceptions.

    It inherits from :class:`BaseException`, not :class:`Exception`, so
    ``except Exception`` does not catch it. Catch ``plc.PyLCError`` instead.

    """
    pass


class PyLCLibraryError(PyLCError):
    """Raised when something is not found in the PyLightcurve data.

    For example, by :func:`locate_planet` when no planet is near the position.

    """
    pass


class PyLCFileError(PyLCError):
    """Raised for problems with files."""
    pass


class PyLCProcessError(PyLCError):
    """Raised when a calculation cannot be completed.

    For example, limb-darkening coefficients for stellar parameters outside the
    model grid.

    """
    pass


class PyLCInputError(PyLCError):
    """Raised when an input is not valid, for example an unknown filter or time format.

    """
    pass
