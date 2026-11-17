# API reference

Everything listed here is available directly from the top-level package:

```python
import pylightcurve as plc
```

| Page | Contents |
|:--|:--|
| {doc}`planet` | The `Planet` class, and `get_planet`, `get_all_planets` and `locate_planet` |
| {doc}`models` | Functions for transit and eclipse models, orbits and limb darkening |
| {doc}`fitting` | The general-purpose MCMC class `EmceeFitting` |
| {doc}`angles_and_times` | Angle classes and the `FixedTarget` class for time conversions |
| {doc}`utilities` | Gaussian fitting, distributions, numerical integration, files and other helpers |
| {doc}`errors` | Exceptions raised by PyLightcurve |

```{toctree}
:maxdepth: 1
:hidden:

planet
models
fitting
angles_and_times
utilities
errors
```
