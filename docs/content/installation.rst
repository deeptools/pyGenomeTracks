Installation
============

Remember -- pyGenomeTracks is available for **command line usage** as well as for
**integration into Galaxy servers**!

.. contents::
    :local:

Requirements
-------------

Python dependencies:

* Python >= 3.6
* numpy >= 1.16
* intervaltree >= 2.1.0
* pyBigWig >= 0.3.16
* hicmatrix >= 15
* pysam >= 0.14
* matplotlib >=3.1.1,<=3.3.2
* gffutils >= 0.9
* pybedtools >= 0.8.1
* tqdm >= 4.20

External dependencies:

* BEDTools

Command line installation using ``conda``
-----------------------------------------

We encourage users to use ``conda`` installation to install pygenometracks. All of the requirements for pyGenomeTracks can be installed via Anaconda with:

.. code:: bash

    $ conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks python=3.7

To get a specific version, one can specify it. For example:

.. code:: bash

    $ conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks=3.5 python=3.7

Command line installation using ``pip``
-----------------------------------------

Install pyGenomeTracks using the following command:
::

	$ pip install pyGenomeTracks

All python requirements should be automatically installed.

Since version 3.5, pyGenomeTracks require BEDTools, do not forget to install it or load it into your environment.

If you need to specify a specific path for the installation of the tools, make use of `pip install`'s numerous options:

.. code:: bash

    $ pip install --install-option="--prefix=/MyPath/Tools/pyGenomeTracks" git+https://github.com/deeptools/pyGenomeTracks.git

Command line installation without ``pip``
-------------------------------------------

You are highly recommended to use `conda install` rather than the following complicated steps.

1. Install the requirements listed above in the "requirements" section. This is done automatically by `pip` (except BEDTools).

2. Download source code
::

	$ git clone https://github.com/deeptools/pyGenomeTracks.git

or if you want a particular release, choose one from https://github.com/deeptools/pygenometracks/releases:
::

	$ wget https://github.com/deeptools/pyGenomeTracks/archive/3.1.tar.gz
	$ tar -xzvf

3. install the source code (if you don't have root permission, you can set
a specific folder using the ``--prefix`` option)
::

	$ python setup.py install --prefix /User/Tools/pyGenomeTracks3.1

Galaxy installation
--------------------

pyGenomeTracks can be easily integrated into a local `Galaxy <http://galaxyproject.org>`_.
The wrapper and its dependencies are available in the `Galaxy Tool
Shed <http://toolshed.g2.bx.psu.edu/view/iuc/pygenometracks>`_.

Installation via Galaxy API (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First generate an `API Key <http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key>`_
for your admin user and run the the installation script:
::

	$ python ./scripts/api/install_tool_shed_repositories.py \
		--api YOUR_API_KEY -l http://localhost/ \
		--url http://toolshed.g2.bx.psu.edu/ \
		-o iuc -r <revision> --name pygenometracks \
		--tool-deps --repository-deps --panel-section-name plots

The ``-r`` argument specifies the version of pygenometracks.

You can watch the installation status under: Top Panel --> Admin --> Manage
installed tool shed repositories

Installation via web browser
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  go to the `admin page <http://localhost:8080/admin>`_
-  select *Search and browse tool sheds*
-  Galaxy tool shed --> Visualization --> pygenometracks
-  install pygenometracks
