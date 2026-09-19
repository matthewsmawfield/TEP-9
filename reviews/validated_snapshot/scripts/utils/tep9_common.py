#!/usr/bin/env python3
"""
TEP-9 shared utilities -- geometry, catalogue parsers, constants.

All orbital-geometry helpers and catalogue readers used across the
TEP-9 pipeline live here so every step consumes identical parsing
and axis definitions.
"""
import io
import json
import math
import os
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROC = PROJECT_ROOT / "data" / "processed"
RESULTS = PROJECT_ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
(RESULTS / "figures").mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# Coordinate frames
# ------------------------------------------------------------------

from scripts.utils.coordinates import ECL2GAL, GAL2ECL, angular_separation


def perih_dir(om, Om, inc):
    """Unit vector toward perihelion from osculating elements (radians)."""
    co, so = np.cos(om), np.sin(om)
    cO, sO, ci, si = np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci,
                     so * si])


def sep(a, b):
    """Angular separation of two unit vectors, degrees."""
    return angular_separation(a, b)


def lv(l, b):
    """Ecliptic longitude/latitude (deg) -> unit vector."""
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l),
                     math.cos(b) * math.sin(l),
                     math.sin(b)])


def gv(l, b):
    """Galactic longitude/latitude (deg) -> ecliptic unit vector."""
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l),
                               math.sin(b)])


def lb(v):
    """Unit vector -> (ecliptic longitude, latitude) in degrees."""
    return (math.degrees(math.atan2(v[1], v[0])) % 360,
            math.degrees(math.asin(np.clip(v[2], -1, 1))))


# Reference axes (degrees, ecliptic frame)
AXES = {
    "tno":   lv(49.0, -17.0),          # detached-TNO cluster axis
    "anti":  lv(229.0, 17.0),          # antipode
    "ism":   lv(255.8, 5.16),          # ISM inflow (Bzowski+2015)
    "cmb+":  gv(264.02, 48.25),        # CMB dipole apex (Planck 2018)
    "cmb-":  gv(84.02, -48.25),        # CMB antapex
    "gctr":  gv(0, 0),
    "gpole": gv(0, 90),
    "epole": np.array([0, 0, 1.0]),
}
TNO_AXIS = AXES["tno"]

# ------------------------------------------------------------------
# SBDB catalogue
# ------------------------------------------------------------------

def load_sbdb(name="sbdb_outer_ss.json"):
    """Load a JPL SBDB query dump -> list of dicts with float elements."""
    d = json.load(open(DATA_RAW / "sbdb" / name))
    fields = d["fields"]
    rows = []
    for rec in d["data"]:
        r = dict(zip(fields, rec))
        rows.append(r)
    return rows


def sbdb_float(r, key):
    try:
        return float(r[key])
    except (TypeError, ValueError, KeyError):
        return float("nan")

# ------------------------------------------------------------------
# Warsaw catalogue (fixed-width .dat tables)
# ------------------------------------------------------------------

def parse_warsaw_orbits(path):
    """tableb/c/d format -> list of dicts."""
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(sample=line[0:2].strip(), com=line[3].strip(),
                desig=line[5:17].strip(), tyr=int(line[27:31]),
                q=float(line[42:56]), e=float(line[56:70]),
                w=float(line[70:82]), Om=float(line[82:94]),
                i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows


def parse_warsaw_a1(path):
    """tablea1 -> {desig: [model/datat/qnew dicts]}."""
    out = {}
    for line in open(path):
        if len(line) < 160:
            continue
        d = line[3:15].strip()
        if d:
            out.setdefault(d, []).append(dict(
                model=line[132:140].strip(), datat=line[122:132].strip(),
                qnew=line[156:159].strip()))
    return out


def parse_warsaw_b4(path):
    """tableb4 -> {desig: NG parameter dict}."""
    out = {}
    for line in open(path):
        if len(line) < 90:
            continue
        d = line[3:15].strip()
        try:
            out[d] = dict(A1=float(line[17:29]), A2=float(line[41:52]),
                          A3=float(line[66:77]), eA1=float(line[29:41]),
                          eA2=float(line[52:64]), eA3=float(line[77:86]),
                          model=line[111:115].strip())
        except ValueError:
            continue
    return out


WARSAW_PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
WARSAW_PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}


def warsaw_dedup(rows, pref):
    out = {}
    for r in rows:
        k = r["desig"]
        if k not in out or pref.get(r["com"], 9) < pref.get(out[k]["com"], 9):
            out[k] = r
    return out


def warsaw_designations():
    """Set of Warsaw comet designations (from tablec)."""
    return {l[5:17].strip() for l in open(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
            if len(l) > 115}

# ------------------------------------------------------------------
# CODE catalogue (HTML table export)
# ------------------------------------------------------------------

class _TableParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.cur, self.buf, self.in_td = [], [], "", False
    def handle_starttag(self, t, a):
        if t == "tr":
            self.cur = []
        elif t == "td":
            self.in_td, self.buf = True, ""
    def handle_endtag(self, t):
        if t == "td":
            self.in_td = False
            self.cur.append(self.buf.strip())
        elif t == "tr" and self.cur:
            self.rows.append(self.cur)
    def handle_data(self, d):
        if self.in_td:
            self.buf += d


def parse_code(path):
    """CODE HTML export -> {desig: row dict}."""
    p = _TableParser()
    p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14:
            continue
        try:
            out[r[0].strip()] = dict(
                desig=r[0].strip(), model=r[1].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                tyr=int(r[7].split()[0]) if r[7].split() else 0,
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out


# ------------------------------------------------------------------
# OSSOS VOTable
# ------------------------------------------------------------------

def load_ossos(name="ossos_t3char.vot"):
    """Parse the OSSOS characterized VOTable -> list of dicts."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(DATA_RAW / "ossos" / name)
    root = tree.getroot()
    ns = {"v": "http://www.ivoa.net/xml/VOTable/v1.3"}
    fields = [f.get("name") for f in root.iter("{*}FIELD")]
    rows = []
    for tr in root.iter("{*}TR"):
        vals = [td.text for td in tr.iter("{*}TD")]
        rows.append(dict(zip(fields, vals)))
    return rows

# ------------------------------------------------------------------
# Logging tee: mirror stdout into the step log file
# ------------------------------------------------------------------

class _Tee(io.TextIOBase):
    def __init__(self, *streams):
        self.streams = streams
    def write(self, s):
        for st in self.streams:
            try:
                st.write(s)
            except Exception:
                pass
        return len(s)
    def flush(self):
        for st in self.streams:
            try:
                st.flush()
            except Exception:
                pass


def tee_stdout(logger):
    """Send all print() output to the step log file as well as console."""
    sys.stdout = _Tee(sys.__stdout__, open(logger.get_log_file(), "a"))
    return logger
