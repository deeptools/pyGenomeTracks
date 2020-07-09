#!/bin/bash
conda create -n pgt_dev --file requirements.txt -y
conda activate pgt_dev
conda install -y pytest nose ghostscript pathlib coverage coverage-badge
python setup.py install
coverage run -m py.test
coverage html
coverage-badge -f -o docs/coverage.svg
