#!/bin/bash

# Test with micromamba
eval "$(micromamba shell hook --shell bash)"

# The channel_priority needs to be set to flexible

## Set environment to generate test images
TRAVIS_PYTHON_VERSION=3.12
for dir in $(find . -name "__pycache__"); do
    rm -r $dir
done
if [ -e build ]; then
    rm -r build
fi
micromamba create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda python=$TRAVIS_PYTHON_VERSION
micromamba activate pgt_test_${TRAVIS_PYTHON_VERSION}
requirement_file=requirements_CI.txt
micromamba install --yes -c bioconda --file ${requirement_file}
python -m pip install .
coverage run -m py.test
coverage html
coverage-badge -f -o docs/coverage.svg
micromamba deactivate

## Set environment to test lowMatplotlib
TRAVIS_PYTHON_VERSION=3.10
for dir in $(find . -name "__pycache__"); do
    rm -r $dir
done
if [ -e build ]; then
    rm -r build
fi
micromamba create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda python=$TRAVIS_PYTHON_VERSION
micromamba activate pgt_test_${TRAVIS_PYTHON_VERSION}
awk -v low='3.8.4' '{if($0~/matplotlib/){$0="matplotlib =="low};if($0=="pytest"){$0="pytest<9.2"};print}' requirements_CI.txt > requirements_CI_low.txt
requirement_file=requirements_CI_low.txt
micromamba install --yes -c bioconda --file ${requirement_file}
python -m pip install .
py.test pygenometracks --doctest-modules -n 4

## Set environment to test pip installation for higher python version
TRAVIS_PYTHON_VERSION=3.14
for dir in $(find . -name "__pycache__"); do
    rm -r $dir
done
if [ -e build ]; then
    rm -r build
fi
micromamba create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda python=$TRAVIS_PYTHON_VERSION
micromamba activate pgt_test_${TRAVIS_PYTHON_VERSION}
micromamba install --yes -c bioconda bedtools
requirement_file=requirements_CI.txt
python -m pip install -r ${requirement_file}
python -m pip install .
py.test pygenometracks --doctest-modules
micromamba deactivate
