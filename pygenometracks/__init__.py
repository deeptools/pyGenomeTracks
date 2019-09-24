import logging
logging.basicConfig(level=logging.INFO)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('hicmatrix').setLevel(logging.WARNING)
logging.getLogger('numexpr').setLevel(logging.WARNING)
