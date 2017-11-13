# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import re

from setuptools import setup
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.

__version__ = '%s'
"""


def update_version_py():
    if not os.path.isdir(".git"):
        print("This does not appear to be a Git repository.")
        return
    try:
        p = subprocess.Popen(["git", "describe",
                              "--tags", "--always"],
                             stdout=subprocess.PIPE)
    except EnvironmentError:
        print("unable to run git, leaving pygenometracks/_version.py alone")
        return
    stdout = p.communicate()[0]
    if p.returncode != 0:
        print("unable to run git, leaving pygenometracks/_version.py alone")
        return
    ver = stdout.strip()
    f = open("pygenometracks/_version.py", "w")
    f.write(VERSION_PY % ver)
    f.close()
    print("set pygenometracks/_version.py to '%s'" % ver)


def get_version():
    try:
        f = open("pygenometracks/_version.py")
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


class sdist(_sdist):

    def run(self):
        # update_version_py()
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)

# Install class to check for external dependencies from OS environment


class install(_install):

    def run(self):
        # update_version_py()
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return

    def checkProgramIsInstalled(self, program, args, where_to_download,
                                affected_tools):
        try:
            subprocess.Popen([program, args],
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE)
            return True
        except EnvironmentError:
            # handle file not found error.
            # the config file is installed in:
            msg = "\n**{0} not found. This " \
                  "program is needed for the following "\
                  "tools to work properly:\n"\
                  " {1}\n"\
                  "{0} can be downloaded from here:\n " \
                  " {2}\n".format(program, affected_tools,
                                  where_to_download)
            sys.stderr.write(msg)

        except Exception as e:
            sys.stderr.write("Error: {}".format(e))

install_requires_py = ["numpy >= 1.12.1",
                       "matplotlib >= 2.0.0",
                       "intervaltree >= 2.1.0",
                       "pyBigWig >=0.3.4",
                       "future >= 0.16.0",
                       "pytest"
                       ]

if sys.version_info[0] == 2 or (sys.version_info[0] == 3 and sys.version_info[1] == 4):
    install_requires_py.append("configparser >= 3.5.0")

setup(
    name='pyGenomeTracks',
    version=get_version(),
    author='Fidel Ramirez',
    author_email='deeptools@googlegroups.com',
    packages=['pygenometracks'],
    scripts=['bin/make_tracks_file', 'bin/pgt', 'bin/pyGenomeTracks'],
    include_package_data=True,
    package_dir={'pygenometracks': 'pygenometracks'},
    url='http://pygenometracks.readthedocs.io',
    license='LICENSE.txt',
    description='Set of programs to process, analyze and visualize Hi-C data',
    long_description=open('README.md').read(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=install_requires_py,
    zip_safe=False,
    python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    cmdclass={'sdist': sdist, 'install': install}
)
