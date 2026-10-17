"""Validated NASA ASCII retrieval with an explicit-TLS FTP fallback."""
import subprocess
import tempfile
from pathlib import Path
from urllib.parse import urlparse
import requests


def validate_daily_ascii(body, columns, min_rows=10000):
    lines = [line.split() for line in body.decode('ascii').splitlines() if line.strip()]
    if len(lines) < min_rows:
        raise ValueError('NASA record is truncated or is not a daily archive')
    for row in lines:
        if len(row) != columns:
            raise ValueError('NASA record has an unexpected column count')
        year, day, hour = map(int, row[:3])
        if not (1963 <= year <= 2100 and 1 <= day <= 366 and hour == 0):
            raise ValueError('NASA record has invalid daily time fields')
        for value in row[3:]:
            float(value)
    return body


def fetch_daily_ascii(url, columns):
    parsed = urlparse(url)
    if parsed.hostname not in {'spdf.gsfc.nasa.gov', 'cdaweb.gsfc.nasa.gov'}:
        raise ValueError('FTP fallback is restricted to the NASA archive')
    https_error = None
    try:
        result = requests.get(url, timeout=20, headers={'User-Agent': 'tep9/1.0'})
        result.raise_for_status()
        return validate_daily_ascii(result.content, columns), url, 'HTTPS'
    except (requests.RequestException, ValueError, UnicodeDecodeError) as exc:
        https_error = str(exc)
    ftp_url = 'ftp://spdf.gsfc.nasa.gov'+parsed.path
    with tempfile.TemporaryDirectory() as directory:
        target = Path(directory)/'archive.dat'
        command = ['curl', '--ssl-reqd', '--fail', '--silent', '--show-error',
                   '--connect-timeout', '15', '--max-time', '120', '--retry', '2',
                   ftp_url, '--output', str(target)]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(f'NASA HTTPS failed ({https_error}); FTP-TLS failed: {result.stderr[-500:]}')
        body = validate_daily_ascii(target.read_bytes(), columns)
    return body, ftp_url, 'FTP with required explicit TLS'
