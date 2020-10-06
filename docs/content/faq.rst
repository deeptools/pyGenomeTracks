FAQ
===

.. contents:: 
    :local:

Why the scale of my Hi-C plot suddenly changed?
-----------------------------------------------
pyGenomeTracks is using `HiCMatrix <https://github.com/deeptools/HiCMatrix>`_ to read the matrix from ``h5`` and ``cool`` format.
From version 12 to version 13, a normalization step when reading ``cool`` file was removed. This normalization was mostly used 
when you were providing ``cool`` file from `cooler balance <https://cooler.readthedocs.io/en/latest/cli.html#cooler-balance>`_.

If you want to keep the old scale you need to downgrade to HiCMatrix version 12 but version 13 also correct some bugs so we advice
to change your ``max_value`` in your parameter file to adjust to the new scale.

My Hi-C plot looks like no correction was applied when using cool matrix
------------------------------------------------------------------------
pyGenomeTracks is using `HiCMatrix <https://github.com/deeptools/HiCMatrix>`_ to read the matrix from ``cool`` format.
Unfortunately, a bug was introduced in version 14 ignoring the correction factors.
This bug was fixed in version 15 so update HiCMatrix to last version should fix it.
