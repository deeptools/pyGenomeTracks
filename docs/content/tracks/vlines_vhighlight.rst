vlines and vhighlight
=====================

Description
-----------

These two types: 'vlines' and 'vhighlight' are not really tracks as they do not occupy a horizontal space.
They go over all tracks plotted in the vertical orientation.
You can define them anywhere in the ini file.
'vlines' can be defined only once and will plot dotted vertical lines corresponding to the beginning of intervals of a bed file.
'vhighlight' can be defined multiple times and will plot vertical rectangle corresponding to the intervals of a bed file.

Parameters
----------

Necessary:
^^^^^^^^^^
- **file**: A bed file.
- **type**: vlines or vhighlight depending what you want to plot.

Optional for vlines:
^^^^^^^^^^^^^^^^^^^^
- **line_width**: `0.5` (default) or float


Optional for vhighlight:
^^^^^^^^^^^^^^^^^^^^^^^^
- **alpha**: `0.5` (default) or float between 0 and 1
- **zorder**: `-100` (default) or float
- **color**: `yellow` (default)
