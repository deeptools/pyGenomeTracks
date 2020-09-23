FAQ
===

.. contents:: 
    :local:

Why the scale of my Hi-C plot suddenly changed
----------------------------------------------
pyGenomeTracks is using `HiCMatrix <https://github.com/deeptools/HiCMatrix>`_ to read the matrix from ``h5`` and ``cool`` format.
From version 12 to version 13, a normalization step when reading ``cool`` file was removed. This normalization was mostly used 
when you were providing ``cool`` file from `cooler balance <https://cooler.readthedocs.io/en/latest/cli.html#cooler-balance>`_.

If you want to keep the old scale you need to downgrade to HiCMatrix version 12 but version 13 also correct some bugs so we advice
to change your ``max_value`` in your parameter file to adjust to the new scale.

No output generated with version 3.5 installed with pip
-------------------------------------------------------
If you used pyGenomeTracks version 3.5 and the last line you get is:

.. code:: bash

    INFO:pygenometracks.tracksClass:initialize x. [xxxxx]

It is highly probable that BEDTools is not installed or not loaded in your environment.
