#!/bin/bash

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh

for TRAVIS_PYTHON_VERSION in 3.6 3.7 3.8; do
    conda create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda -c conda-forge python=$TRAVIS_PYTHON_VERSION --file requirements.txt
    conda activate pgt_test_${TRAVIS_PYTHON_VERSION}
    conda install --yes -c conda-forge -c bioconda pytest ghostscript coverage coverage-badge
    python setup.py install
    coverage run -m py.test
    coverage html
    coverage-badge -f -o docs/coverage.svg
    conda deactivate
done
