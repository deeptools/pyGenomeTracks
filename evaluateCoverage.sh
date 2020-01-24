conda create -n pgt_dev --file requirements.txt -y
conda activate pgt_dev
conda install -y pytest nose ghostscript pathlib coverage
python setup.py install
coverage run -m py.test
coverage html
pip install coverage-badge
coverage-badge -o docs/coverage.svg
