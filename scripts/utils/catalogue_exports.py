"""Structural checks for published CODE orbit-table exports."""
from html.parser import HTMLParser


class _Rows(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.cell, self.in_cell = [], [], '', False

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.row = []
        elif tag == 'td':
            self.cell, self.in_cell = '', True

    def handle_data(self, data):
        if self.in_cell:
            self.cell += data

    def handle_endtag(self, tag):
        if tag == 'td':
            self.row.append(self.cell.strip())
            self.in_cell = False
        elif tag == 'tr' and len(self.row) >= 14:
            self.rows.append(self.row)


def code_export_models(body):
    """Require one preferred solution per comet; prevent silent overwrite."""
    parser = _Rows()
    parser.feed(body.decode('utf-8', errors='replace'))
    models = {}
    for row in parser.rows:
        name = row[0]
        if name in models:
            raise ValueError(f'CODE export has multiple solutions for {name}')
        models[name] = row[1]
    if not models:
        raise ValueError('CODE export contains no orbital rows')
    return models


def require_matching_code_models(tables):
    """Shared objects must use the same published solution across legs."""
    common = set.intersection(*(set(t) for t in tables.values()))
    mismatched = [name for name in common
                  if len({table[name] for table in tables.values()}) != 1]
    if mismatched:
        raise ValueError(f'CODE solution models disagree across legs: {mismatched[:8]}')
    return len(common)
