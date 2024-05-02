#!/bin/bash

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh

for TRAVIS_PYTHON_VERSION in 3.8 3.9 3.10; do
    for dir in $(find . -name "__pycache__"); do
        rm -r $dir
    done
    if [ -e build ]; then
        rm -r build
    fi
    conda create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda -c conda-forge python=$TRAVIS_PYTHON_VERSION mamba
    conda activate pgt_test_${TRAVIS_PYTHON_VERSION}
    requirement_file=requirements_CI.txt
    mamba install --yes -c conda-forge -c bioconda --file ${requirement_file}
    python setup.py install
    coverage run -m py.test
    coverage html
    coverage-badge -f -o docs/coverage.svg
    conda deactivate
done

TRAVIS_PYTHON_VERSION=3.11
conda create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda -c conda-forge python=$TRAVIS_PYTHON_VERSION
conda activate pgt_test_${TRAVIS_PYTHON_VERSION}
conda install --yes -c conda-forge -c bioconda bedtools
pip install -r requirements_CI.txt
python setup.py install
py.test pygenometracks --doctest-modules
conda deactivate
