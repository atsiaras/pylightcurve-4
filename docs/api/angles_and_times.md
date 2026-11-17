# Angles and times

```{eval-rst}
.. currentmodule:: pylightcurve
```

## Angles

All three angle classes have the same methods. They are listed under
{class}`Degrees`.

```{eval-rst}
.. autoclass:: Degrees
   :members: deg, hours, rad, dms, hms, deg_coord, dms_coord, sin, cos, tan
   :inherited-members:

.. autoclass:: Hours
   :no-members:

.. autoclass:: Rad
   :no-members:

.. autofunction:: arcsin
.. autofunction:: arccos
.. autofunction:: arctan

.. data:: pi

   The angle of 180 degrees, as a :class:`Degrees` object.
```

## Targets and time conversions

```{eval-rst}
.. autoclass:: FixedTarget
   :members: convert_to_bjd_tdb, convert_to_jd, convert_to_mjd, distance_on_sphere
```
