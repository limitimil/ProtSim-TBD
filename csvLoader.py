import csv
_fieldnames = ['chainA','chainB','PDBid','codeA','codeB','path','locA','locB','SegA','SegB']
class _csv_r(object):
    fopen = False
    def __init__(self, stream):
        self.reader = csv.DictReader(stream)
        self.fopen = True
class _csv_w(object):
    fopen = False
    def __init__(self, stream):
        self.writer = csv.DictWriter(stream, _fieldnames)
        self.fopen = False
