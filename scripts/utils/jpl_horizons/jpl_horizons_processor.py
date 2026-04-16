#!/usr/bin/env python
"""Batch processor for JPL/Horizons CGI response - Adapted for TEP-3I pipeline"""

from io import StringIO
from scipy import constants as c
import sys
import re
import pandas as pd
import json
from pathlib import Path


class JPLHorizonsProcessor:
    """Batch query processor for JPL/Horizons CGI interface

    This script parses JPL Horizons responses and extracts real trajectory
    data for TEP analysis.
    
    Processing:
        - Parses JPL Horizons response format
        - Converts range from km to meters
        - Converts velocity from km/s to m/s
        - Converts light time from minutes to seconds
        - Calculates acceleration, doppler rate
        - Calculates range and velocity lags
    """

    def __init__(self):
        self._re_data = re.compile(r'\$\$SOE.*\$\$EOE', re.DOTALL)
        self._re_angles = re.compile(r'\s*([+-]?\d{2}).(\d{2}).(\d{2}\.\d+)')
        self._speed = 299792458  # speed of light in m/s

    def parse(self, response):
        """Parse the response from the JPL Horizons server"""
        match_result = self._re_data.search(response)
        if not match_result:
            raise ValueError("No data section found in JPL Horizons response")
        
        # Extract data between $$SOE and $$EOE
        data_text = match_result.group()[6:match_result.end()-6]
        lines = data_text.strip().split('\n')
        
        # Month abbreviation to number mapping
        month_map = {
            'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
            'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
            'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'
        }
        
        # Parse data lines manually
        data_lines = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith('$'):
                # Split by comma
                parts = line.split(',')
                if len(parts) >= 8 and parts[0].strip():
                    # Convert date format from "1998-Jan-22 20:00" to "1998-01-22 20:00"
                    date_str = parts[0].strip()
                    for month_abbr, month_num in month_map.items():
                        date_str = date_str.replace(month_abbr, month_num)
                    parts[0] = date_str
                    data_lines.append(parts)
        
        if not data_lines:
            raise ValueError("No valid data lines found in JPL Horizons response")
        
        # Create DataFrame from parsed lines
        data = pd.DataFrame(data_lines, 
                           columns=['timestamp', 'empty1', 'empty2', 'azi_ra', 
                                    'el_dec', 'lighttime', 'range', 'velocity', 'extra'])
        
        # Select only the columns we need and convert types
        data = data[['timestamp', 'azi_ra', 'el_dec', 'lighttime', 'range', 'velocity']].copy()
        data['timestamp'] = pd.to_datetime(data['timestamp'])
        data['lighttime'] = pd.to_numeric(data['lighttime'])
        data['range'] = pd.to_numeric(data['range'])
        data['velocity'] = pd.to_numeric(data['velocity'])
        data = data.set_index('timestamp')
        
        return data

    def hms(self, sh, sm, ss):
        """Convert HMS format to decimal degrees"""
        return float(sh) + float(sm)/60.0 + float(ss)/3600

    def normalize_angle(self, datum):
        """Normalize angle data from HMS format"""
        m = self._re_angles.match(str(datum))
        if m:
            return self.hms(m.group(1), m.group(2), m.group(3))
        return datum

    def normalize(self, data):
        """Normalize data returned by the JPL/Horizons server"""
        # Calculate time differences
        data['dt'] = data.index.to_series().diff().dt.seconds \
                     + data.index.to_series().diff().dt.microseconds.div(1000000, fill_value=0)

        # Convert units
        data['range'] *= 1000  # km to m
        data['velocity'] *= 1000  # km/s to m/s
        data['lighttime'] = data['lighttime'] * 60.0  # minutes to seconds

        # Calculate derived quantities
        data['acceleration'] = data['velocity'].diff()/data['dt']
        data['traveltime'] = data['range'] / self._speed  # one-way light time

        # Calculate range and velocity lags (for flyby anomaly analysis)
        data['rangelag'] = data['velocity'] * data['traveltime']
        data['velocitylag'] = data['acceleration'] * data['traveltime']

        data['ratedoppler'] = data['acceleration'] / c.c

        # Normalize angle data if present
        if 'el/dec' in data.columns:
            data['el/dec'] = data['el/dec'].map(lambda a: self.normalize_angle(a))

        return data

    def load(self, dname):
        """Load server response data file for processing"""
        with open(dname) as dfile:
            resp = dfile.read()
        output = StringIO()
        self.normalize(self.parse(resp)).to_csv(output, sep='\t')
        output.seek(0)
        return output.read()

    def process_to_json(self, response, output_file):
        """Process JPL Horizons response and save as JSON for pipeline use"""
        data = self.normalize(self.parse(response))
        
        # Convert to dict for JSON serialization
        result = {
            'timestamp': data.index.astype(str).tolist(),
            'range_m': data['range'].tolist(),
            'velocity_m_s': data['velocity'].tolist(),
            'lighttime_s': data['lighttime'].tolist(),
            'acceleration_m_s2': data['acceleration'].tolist(),
            'rangelag_m': data['rangelag'].tolist(),
            'velocitylag_m_s': data['velocitylag'].tolist()
        }
        
        # Save to JSON
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)
        
        return result


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(JPLHorizonsProcessor().load(sys.argv[1]))
    else:
        print("Usage: python jpl_horizons_processor.py <response_file.txt>")
