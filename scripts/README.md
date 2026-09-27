# Analysis and publication tools

`run_all.py` defines 151 registered steps, from source acquisition through orbital, spatial and spacecraft diagnostics. It stops on a failed step and records execution status and source hashes. The exact publication-claims audit runs last because it depends on later result files.

The step names and historical output names use different numbering schemes. Resolve headline manuscript numbers through `../site/claims.json`, not by matching any nearby number in a result file.

Every step constructs `StepLogger("<step_name>")` and calls `tee_stdout(logger)` immediately after; both live in `scripts/utils/`. The logger writes `logs/<step_name>.log` with a UTC-timestamped, step-prefixed record per line: a purpose/environment/script-hash prologue at start and a duration/status/output-file summary (sizes and sha256) at exit. `tee_stdout` routes ordinary `print()` lines through the same handler, so console output is never lost from the log; it is a no-op under nested teeing (e.g. step 134's exec'd step-127 machinery). Register every written artefact with `logger.data_save(path)` so it appears in the summary and provenance. Download steps retry transient HTTP failures with exponential backoff before failing hard.

Useful review commands, from the repository root:

```bash
python3 -m pytest -q tests
python3 scripts/audits/audit_tep9_observables.py
python3 scripts/audits/audit_refit_replay.py
python3 scripts/audits/audit_publication.py
npm run build:markdown --prefix site
python3 scripts/generate_site_pdf.py
```

The nine-object refit replay is a stratified regression check on cached real astrometry, not a full refit of every catalogue member. The observable audit separates direct catalogue-minus-dynamics residuals from reduced-regression population diagnostics. Mathematical and simulated controls are not counted as empirical TEP evidence.

See `../reviews/FULL_REVIEW_20260919.md` for corrected issues, rerun coverage and unresolved physical/statistical requirements. Historical result verdict strings may express stronger interpretations than the revised main paper; they are not independent evidence or a replacement for the audited values.
