"""Recover the persisted-but-dropped leg-fit channels for the b91
independent-refit records: d_in_leg / d_out_leg (independent leg
asymptotes vs the full-arc osculating direction -- the commensurate
analogue of the catalogue orig-vs-osc / fut-vs-osc construction).

Observations and SBDB seeds are already cached under data/raw/mpc/;
this pass is pure compute.  Runs fit_comet on every b91 dual-leg
comet and writes results/step_b91_leg_channels.jsonl.
"""
import json
import multiprocessing as mp
import sys
from pathlib import Path

sys.path.insert(0, ".")
from scripts.utils import mpc_refit
from scripts.utils.mpc_refit import fit_comet
from scripts.utils.tep9_common import DATA_RAW

MPC_DIR = DATA_RAW / "mpc"


def one(des):
    try:
        res, err = fit_comet(des)
        if res is None:
            return {"des": des, "error": err}
        out = {"des": des}
        for k in ("d_in_leg", "d_out_leg", "ddirf", "drot"):
            if k in res:
                out[k] = res[k]
        return out
    except Exception as exc:
        return {"des": des, "error": str(exc)[:120]}


if __name__ == "__main__":
    mpc_refit.configure(
        MPC_DIR / "obs", MPC_DIR / "sbdb_fp",
        MPC_DIR / "obscodes.json", MPC_DIR / "provenance.json",
        DATA_RAW / "spice" / "de440s.bsp",
        DATA_RAW / "naif" / "naif0012.tls",
        "recover_leg_channels")
    recs = [json.loads(l) for l in open("results/step_b91_refit.jsonl")]
    dual = [r["des"] for r in recs if "our_ddirf" in r]
    print("recovering leg channels for %d dual-leg comets" % len(dual),
          flush=True)
    ctx = mp.get_context("fork")
    with ctx.Pool(8, initializer=mpc_refit.refit_worker_init) as pool:
        rows = pool.map(one, dual)
    ok = [r for r in rows if "d_in_leg" in r]
    bad = [r for r in rows if "d_in_leg" not in r]
    print("recovered %d / %d (failures %d)" % (len(ok), len(rows), len(bad)))
    with open("results/step_b91_leg_channels.jsonl", "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
