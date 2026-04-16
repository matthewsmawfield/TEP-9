#!/usr/bin/env python
"""Batch query for JPL/Horizons CGI - Adapted for TEP-3I pipeline"""

import re
import sys
import requests
from pathlib import Path

class QueryError(Exception):
    """Base class for query exceptions"""

    def __init__(self, message):
        super().__init__(message)
        self.message = message


class JPLHorizonsQuery:
    """Batch query for JPL/Horizons CGI interface

    This script queries JPL's official trajectory repository to obtain
    real range and range rate data for spacecraft flybys.
    
    Example:
        query = JPLHorizonsQuery()
        data = query.fetch('query.q')
        print(data)

    Processing:
        - Queries JPL Horizons CGI interface
        - Returns raw response text for further processing
        - Used to obtain real trajectory data for TEP analysis
    """

    def __init__(self):
        self._re_spaces = re.compile(r' *= *')
        self._url = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'
        self._pfx = 'batch=1'

    def encode(self, string):
        """Encode parameter string for CGI GET"""
        string = self._re_spaces.sub(r'=', string)
        string = string.replace(r' ', '%20')
        string = string.replace(r';', '%3b')
        string = string.replace(r'&', '%26')
        string = string.replace(r'?', '%3f')
        return string

    def compose(self, fname):
        """Compile CGI GET parameters from file"""
        with open(fname) as pfile:
            query = [self.encode(line.strip()) if '$$' not in line else ''
                     for line in pfile]
        return '&'.join(query)

    def fetch(self, qfname):
        """Fetch data from the JPL/Horizons server, using query parameters from file"""
        query = self.compose(qfname)
        req = requests.get(self._url, params=self._pfx + query)

        if req.status_code != 200:
            raise QueryError("Unexpected response status: {}".format(req.status_code))

        if "Cannot" in req.text:
            issue = ''.join([s if "Cannot" in s else '' for s in req.text.splitlines()])
            raise QueryError("Query failed: {}".format(issue))

        return req.text


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(JPLHorizonsQuery().fetch(sys.argv[1]))
    else:
        print("Usage: python jpl_horizons_query.py <query_file.q>")
