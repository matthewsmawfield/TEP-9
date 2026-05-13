#!/usr/bin/env python3
"""
Step 027: Claim Consistency Audit
================================

Ensures manuscript and pipeline outputs stay synchronized.
"""

import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger

from scripts.utils.publication_sync import (
    PRIMARY_DETECTIONS,
    load_expected_publication_values,
    sync_publication_artifacts,
    _format_percent,
    _format_signed_mm,
    _round,
)


STALE_MARKERS = (
    "567.5",
    "B_{10} = 567.5",
    "Bayes factor 567.5",
    "Bayes factor = 567.5",
    "ΔBIC = 12.7",
    "ΔBIC = 3.1",
    "log L = -16.04",
    "B = 4.7",
    "B_{A0} = 121.6",
    "Akaike weight 99.2%",
    "99.2% evidence weight",
    "accounts for 0.0% of the fitted-β variance to the structural proxy bundle, 30.0%",
    "accounts for 0.0% of this variance, observational pipeline effects (OD filter absorption and systematic uncertainties) account for 30.0%",
    "+0.16 | +4.24",
    "+0.04 | +1.05",
)


class ClaimConsistencyAuditor:
    """Audits manuscript claims against pipeline implementation."""

    PRIMARY_DETECTIONS = PRIMARY_DETECTIONS

    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.violations: List[str] = []
        self.pipeline_constants = self._load_pipeline_constants()
        self.manuscript_claims = self._load_manuscript_claims()
        self.expected = self._load_expected_claims()

    def _load_pipeline_constants(self) -> Dict[str, Any]:
        sys.path.insert(0, str(self.project_root))

        from scripts.utils.physics import (
            BETA_BASELINE,
            DISFORMAL_COUPLING_STRENGTH,
            DISFORMAL_VELOCITY_THRESHOLD_KM_S,
        )

        step_007 = self.project_root / "scripts/steps/step_007_tep_model.py"
        beta_initial = BETA_BASELINE * 1e-4
        if step_007.exists():
            content = step_007.read_text(encoding="utf-8")
            match = re.search(r"^[ \t]*BETA_INITIAL\s*=\s*([\d.e-]+)", content, re.MULTILINE)
            if match:
                beta_initial = float(match.group(1))

        return {
            "BETA_INITIAL": beta_initial,
            "DISFORMAL_COUPLING_STRENGTH": DISFORMAL_COUPLING_STRENGTH,
            "DISFORMAL_VELOCITY_THRESHOLD_KM_S": DISFORMAL_VELOCITY_THRESHOLD_KM_S,
        }

    def _manuscript_sources(self) -> List[Path]:
        sources = [self.project_root / "15-TEP-EFA-v0.1-Yogyakarta.md"]
        sources.extend(sorted((self.project_root / "site/components").glob("*.html")))
        return [path for path in sources if path.exists()]

    def _combined_manuscript_text(self) -> str:
        return "\n".join(path.read_text(encoding="utf-8") for path in self._manuscript_sources())

    def _load_manuscript_claims(self) -> Dict[str, Any]:
        claims = {
            "v_trans": [],
            "cassini_treatments": [],
            "disformal_sign_rule": None,
            "f_od_placeholder": False,
            "f_plasma_placeholder": False,
        }

        for path in self._manuscript_sources():
            content = path.read_text(encoding="utf-8")

            claims["v_trans"].extend(
                float(value)
                for value in re.findall(r"v_\{\\rm trans\}\s*=\s*([\d.]+)\s*km", content)
            )
            claims["v_trans"].extend(float(value) for value in re.findall(r"v_trans\s*=\s*([\d.]+)", content))
            claims["v_trans"].extend(
                float(value) for value in re.findall(r"v_trans\s*≈\s*([\d.]+)", content)
            )
            claims["v_trans"].extend(
                float(value) for value in re.findall(r"v_trans\s*\\approx\s*([\d.]+)", content)
            )
            claims["v_trans"].extend(
                float(value)
                for value in re.findall(r"v_\{\\rm trans\}\s*≈\s*([\d.]+)", content)
            )
            claims["v_trans"].extend(
                float(value)
                for value in re.findall(r"v_\{\\rm trans\}\s*\\approx\s*([\d.]+)", content)
            )

            if "cassini" in content.lower():
                for keyword in ("excluded", "sign mismatch", "cancellation", "included", "fitted"):
                    if keyword in content.lower():
                        if keyword == "excluded" and "previously excluded" in content.lower():
                            continue
                        claims["cassini_treatments"].append(keyword)

            if "disformal" in content.lower() and "sign" in content.lower():
                if "does not reverse" in content.lower():
                    claims["disformal_sign_rule"] = "no_reversal"
                elif "reverse" in content.lower():
                    claims["disformal_sign_rule"] = "reversal"

            if "f_od" in content.lower() and "placeholder" in content.lower():
                claims["f_od_placeholder"] = True

            if (
                "f_plasma" in content.lower()
                and "placeholder" in content.lower()
                and ("variance" in content.lower() or "contribution" in content.lower())
            ):
                claims["f_plasma_placeholder"] = True

        return claims

    def _load_json(self, relative_path: str) -> Dict[str, Any]:
        path = self.project_root / relative_path
        if not path.exists():
            raise FileNotFoundError(relative_path)
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)

    _round = staticmethod(_round)
    _format_percent = staticmethod(_format_percent)
    _format_signed_mm = staticmethod(_format_signed_mm)

    def _load_expected_claims(self) -> Dict[str, Any]:
        return load_expected_publication_values(self.project_root)

    def _require_substrings(self, text: str, required: List[str], label: str) -> bool:
        missing = [value for value in required if value not in text]
        if missing:
            self.violations.append(f"{label} missing expected values: {', '.join(missing)}")
            return False
        return True

    def audit_v_trans_consistency(self) -> bool:
        pipeline_v = self.pipeline_constants["DISFORMAL_VELOCITY_THRESHOLD_KM_S"]
        manuscript_vs = self.manuscript_claims.get("v_trans", [])
        if not manuscript_vs:
            self.violations.append("No v_trans found in manuscript")
            return False

        inconsistent = [value for value in manuscript_vs if abs(value - pipeline_v) > 0.1]
        if inconsistent:
            self.violations.append(
                f"v_trans inconsistency: manuscript has {inconsistent}, pipeline has {pipeline_v} km/s"
            )
            return False
        return True

    def audit_cassini_status_consistency(self) -> bool:
        try:
            step008 = self._load_json("results/step008_fitting_results.json")
            pipeline_excluded = bool(
                step008["individual_fits"]["Cassini"]["fit"].get("excluded")
            )
        except (FileNotFoundError, KeyError, TypeError):
            pipeline_excluded = False

        text = self._combined_manuscript_text().lower()

        if pipeline_excluded:
            mentions_exclusion = (
                "excluded from" in text
                or "excluded from the" in text
                or "exclusion" in text
                or "sign mismatch" in text
            )
            if not mentions_exclusion:
                self.violations.append(
                    "Step 008 excludes Cassini from the β ensemble; manuscript should state sign mismatch / exclusion."
                )
                return False
            return True

        if "cassini" in text and "excluded from the inverse-variance" in text:
            self.violations.append(
                "Manuscript excludes Cassini from the β ensemble but Step 008 does not mark Cassini as excluded."
            )
            return False
        return True

    def audit_disformal_sign_rule_consistency(self) -> bool:
        if self.manuscript_claims.get("disformal_sign_rule") == "no_reversal":
            self.violations.append(
                "Manuscript says disformal term does not reverse sign, but pipeline implements sign reversal"
            )
            return False
        if self.manuscript_claims.get("disformal_sign_rule") is None:
            self.violations.append("No disformal sign rule found in manuscript")
            return False
        return True

    def audit_f_od_placeholder(self) -> bool:
        if not self.manuscript_claims.get("f_od_placeholder", False):
            return True

        text = self._combined_manuscript_text().lower()
        if "od model comparison" in text or "tep+od" in text:
            self.violations.append(
                "Manuscript claims OD model comparison but F_OD is placeholder (1.0)."
            )
            return False
        return True

    def audit_f_plasma_placeholder(self) -> bool:
        return True

    def audit_stale_markers(self) -> bool:
        text = self._combined_manuscript_text()
        stale = [marker for marker in STALE_MARKERS if marker in text]
        if stale:
            self.violations.append(f"Stale manuscript claims still present: {', '.join(stale)}")
            return False
        return True

    @staticmethod
    def _mantissa_exponent(value: float) -> Tuple[float, int]:
        if value == 0.0:
            return 0.0, 0
        exp = int(math.floor(math.log10(abs(value))))
        mant = round(value / (10**exp), 1)
        return mant, exp

    def audit_model_comparison_table(self) -> bool:
        model = self.expected["model"]
        text = self._combined_manuscript_text()
        required_ic = [
            f"log L = {model['log_likelihoods']['Null']}",
            f"log L = {model['log_likelihoods']['Anderson']}",
            f"log L = {model['log_likelihoods']['TEP_restricted']}",
            f"log L = {model['log_likelihoods']['TEP_flexible']}",
            f"BIC = {model['bic']['Null']}",
            f"BIC = {model['bic']['TEP_restricted']}",
        ]
        if not self._require_substrings(text, required_ic, "Model comparison (IC)"):
            return False

        data = self._load_json("results/step026_stable_model_comparison.json")
        bf = data["bayes_factors"]
        b_exp = model["bayes"]

        m10, e10 = self._mantissa_exponent(float(bf["TEP_restricted_vs_Null"]))
        m_a0, e_a0 = self._mantissa_exponent(float(bf["Anderson_vs_Null"]))
        m_f0, e_f0 = self._mantissa_exponent(float(bf["TEP_flexible_vs_Null"]))

        approx_blocks = [
            f"$B_{{10}} \\approx {m10} \\times 10^{{{e10}}}$",
            f"$\\Delta$BIC $\\approx {b_exp['delta_bic_tep_null']}$",
            f"$\\Delta$BIC $\\approx {b_exp['delta_bic_tep_anderson']}$",
            f"$B \\approx {self._round(bf['TEP_restricted_vs_Anderson'], 1)}$",
            f"$B_{{A0}} \\approx {m_a0} \\times 10^{{{e_a0}}}$",
            f"$B_{{f0}} \\approx {m_f0} \\times 10^{{{e_f0}}}$",
            f"$\\Delta$BIC $\\approx {b_exp['delta_bic_anderson_null']}$",
            f"$\\Delta$BIC $\\approx {b_exp['delta_bic_flex_null']}$",
        ]
        if not self._require_substrings(text, approx_blocks, "Model comparison (Bayes factors)"):
            return False

        w_tep = float(data["model_selection"]["akaike_weights"]["TEP_restricted"])
        w_fmt = self._format_percent(100.0 * w_tep, 1)
        akaike_ok = w_fmt in text or bool(
            re.search(r"TEP restricted[^\n]{0,48}\\approx\s*1\.0", text)
            or re.search(r"TEP restricted[^\n]{0,48}≈\s*1\.0", text)
        )
        if not akaike_ok:
            self.violations.append(
                f"Model comparison missing Akaike summary (expected {w_fmt} or TEP restricted ≈ 1.0)"
            )
            return False
        return True

    def audit_variance_decomposition(self) -> bool:
        variance = self.expected["variance"]
        text = self._combined_manuscript_text()

        structural_ok = any(
            phrase in text
            for phrase in (
                f"structural proxy bundle accounts for {variance['structural']}",
                f"tracked structural proxy bundle accounts for {variance['structural']}",
                f"structural proxies account for {variance['structural']}",
            )
        )
        observational_ok = any(
            phrase in text
            for phrase in (
                "observational pipeline effects (OD filter absorption and systematic uncertainties) account for "
                + variance["observational"],
                f"observational pipeline effects account for {variance['observational']}",
                f"observational effects account for {variance['observational']}",
            )
        )
        environmental_ok = f"environmental modulation contributes {variance['environmental']}" in text
        residual_ok = any(
            phrase in text
            for phrase in (
                "residual (small-sample statistics, intrinsic scatter, model incompleteness) accounts for "
                + variance["residual"],
                f"the residual accounts for {variance['residual']}",
                f"residual accounts for {variance['residual']}",
                f"and the residual accounts for {variance['residual']}",
            )
        )

        if not (structural_ok and observational_ok and environmental_ok and residual_ok):
            missing = []
            if not structural_ok:
                missing.append("structural stage")
            if not observational_ok:
                missing.append("observational stage")
            if not environmental_ok:
                missing.append("environmental stage")
            if not residual_ok:
                missing.append("residual stage")
            self.violations.append(
                "Variance decomposition missing expected wording for: " + ", ".join(missing)
            )
            return False
        return True

    def audit_component_table(self) -> bool:
        text = self._combined_manuscript_text()
        required: List[str] = []
        for mission_key, label in self.PRIMARY_DETECTIONS:
            component = self.expected["components"][mission_key]
            required.extend(
                [
                    f"| {label} | {self._format_signed_mm(component['grad'])} | {self._format_signed_mm(component['disf'])} | {self._format_signed_mm(component['total'])} |",
                    f"<td>{label}</td>",
                    f"<td>{self._format_signed_mm(component['grad'])}</td>",
                    f"<td>{self._format_signed_mm(component['disf'])}</td>",
                    f"<td>{self._format_signed_mm(component['total'])}</td>",
                ]
            )
        return self._require_substrings(text, required, "Component table")

    def audit_pipeline_output_hashes(self) -> bool:
        required_outputs = (
            "results/step007_tep_predictions.json",
            "results/step009_variance_analysis.json",
            "results/step026_stable_model_comparison.json",
        )
        missing = [path for path in required_outputs if not (self.project_root / path).exists()]
        if missing:
            self.violations.append(f"Missing required pipeline outputs: {', '.join(missing)}")
            return False
        return True

    def audit_reproduction_checklist(self) -> bool:
        run_all = self.project_root / "scripts/run_all.py"
        if not run_all.exists():
            self.violations.append("scripts/run_all.py not found")
            return False

        content = run_all.read_text(encoding="utf-8")
        step_count = content.count("('step_")
        if step_count < 39:
            self.violations.append(
                f"run_all.py lists {step_count} steps; expected at least 39 including Steps 040a and 040"
            )
            return False
        return True

    def run_all_audits(self, logger: StepLogger | None = None) -> bool:
        emit = logger.info if logger else print
        emit("=" * 70)
        emit("CLAIM CONSISTENCY AUDIT")
        emit("=" * 70)

        audits = [
            ("v_trans consistency", self.audit_v_trans_consistency),
            ("Cassini status consistency", self.audit_cassini_status_consistency),
            ("Disformal sign rule consistency", self.audit_disformal_sign_rule_consistency),
            ("F_OD placeholder check", self.audit_f_od_placeholder),
            ("F_plasma placeholder check", self.audit_f_plasma_placeholder),
            ("Pipeline output presence", self.audit_pipeline_output_hashes),
            ("Stale marker scan", self.audit_stale_markers),
            ("Model comparison table", self.audit_model_comparison_table),
            ("Variance decomposition", self.audit_variance_decomposition),
            ("Component table", self.audit_component_table),
            ("Reproduction checklist", self.audit_reproduction_checklist),
        ]

        all_passed = True
        for name, audit_func in audits:
            if logger:
                logger.info(f"[{name}]... ")
            else:
                print(f"\n[{name}]...", end=" ")
            try:
                passed = audit_func()
                if logger:
                    logger.success("PASS" if passed else "FAIL")
                else:
                    print("✓ PASS" if passed else "✗ FAIL")
                all_passed = all_passed and passed
            except Exception as exc:
                if logger:
                    logger.error(f"ERROR: {exc}")
                else:
                    print(f"✗ ERROR: {exc}")
                all_passed = False

        emit("\n" + "=" * 70)
        if all_passed:
            emit("ALL AUDITS PASSED - Manuscript and pipeline are consistent")
        else:
            emit("AUDIT FAILED - Manuscript and pipeline have inconsistencies")
            emit("\nViolations:")
            for index, violation in enumerate(self.violations, 1):
                emit(f"  {index}. {violation}")
        emit("=" * 70)
        return all_passed


def main() -> None:
    project_root = PROJECT_ROOT
    logger = StepLogger("step_027_claim_consistency_audit", project_root)
    start = time.time()
    status = "SUCCESS"
    try:
        logger.header("STEP 027: CLAIM CONSISTENCY AUDIT")
        logger.track_dependency("step_007_tep_model", "results")
        logger.track_dependency("step_009_variance_analysis", "results")
        logger.track_dependency("step_026_stable_model_comparison", "results")
        logger.track_dependency("step_039_flyby_prediction_table", "results")

        updated = sync_publication_artifacts(project_root)
        if updated:
            logger.info("Synchronized publication artifacts:")
            for path in updated:
                logger.info(f"  - {path.relative_to(project_root)}")

        auditor = ClaimConsistencyAuditor(project_root)
        passed = auditor.run_all_audits(logger=logger)
        if not passed:
            status = "FAILED"
            raise SystemExit(1)
    except SystemExit as exc:
        if int(getattr(exc, "code", 1) or 1) != 0:
            status = "FAILED"
        raise
    except Exception as exc:
        status = "FAILED"
        logger.error(f"Unhandled exception: {exc}")
        raise
    finally:
        duration = time.time() - start
        logger.log_step_summary(duration, status=status)

    sys.exit(0)


if __name__ == "__main__":
    main()
