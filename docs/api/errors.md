# Errors

All PyLightcurve exceptions derive from {class}`~pylightcurve.PyLCError`.

:::{warning}
{class}`~pylightcurve.PyLCError` inherits from `BaseException`, not `Exception`,
so `except Exception:` does not catch PyLightcurve errors. Catch them explicitly:

```python
try:
    planet = plc.get_planet('not a planet')
except plc.PyLCError as error:
    print(error)
```
:::

```{eval-rst}
.. currentmodule:: pylightcurve

.. autoexception:: PyLCError
.. autoexception:: PyLCInputError
.. autoexception:: PyLCLibraryError
.. autoexception:: PyLCProcessError
.. autoexception:: PyLCFileError
```
