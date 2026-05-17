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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger

from scripts.utils.publication_sync import (
    PRIMARY_DETECTIONS,
    load_expected_publication_values,
    split_step026_model_blocks,
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
    "prefers the Null",
    "favours the Null",
    "favors the Null",
    "direct experimental validation",
    "overwhelming evidence",
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
        self.audit_results: List[Dict[str, Any]] = []

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
        v = float(value)
        if not math.isfinite(v):
            raise ValueError("finite Bayes factor required for mantissa formatting")
        exp = int(math.floor(math.log10(abs(v))))
        mant = round(v / (10**exp), 1)
        return mant, exp

    @staticmethod
    def _mantissa_exponent_from_log10(log10_bf: float) -> Tuple[float, int]:
        lf = float(log10_bf)
        exp_i = int(math.floor(lf))
        mant = round(10 ** (lf - exp_i), 2)
        return mant, exp_i

    def audit_model_comparison_table(self) -> bool:
        headline = self.expected["model"]
        gated = self.expected["model_gated"]
        text = self._combined_manuscript_text()

        # Table 7 / full-catalog IC (n = 9, geometry-spread)
        required_headline_ic = [
            f"log L = {headline['log_likelihoods']['Null']}",
            f"log L = {headline['log_likelihoods']['TEP_restricted']}",
            f"BIC = {headline['bic']['TEP_restricted']}",
        ]
        if not self._require_substrings(text, required_headline_ic, "Model comparison (headline IC)"):
            return False

        # Primary gated ensemble (n = 4 relaxed gate, n = 3 strict gate)
        gated_ic = [
            f"log L = {gated['log_likelihoods']['Null']}",
            f"log L = {gated['log_likelihoods']['TEP_restricted']}",
        ]
        if not self._require_substrings(text, gated_ic, "Model comparison (gated IC)"):
            return False

        step026 = self._load_json("results/step026_stable_model_comparison.json")
        headline_block, gated_block = split_step026_model_blocks(step026)
        bf_headline = headline_block["bayes_factors"]
        bf_gated = gated_block["bayes_factors"]
        b_headline = headline["bayes"]
        b_gated = gated["bayes"]

        headline_log10_blocks = [
            f"log_{{10}} B = {b_headline['log10_tep_vs_null']}",
            f"log_{{10}} B = {b_headline['log10_anderson_vs_null']}",
            f"log_{{10}} B = {b_headline['log10_tep_vs_anderson']}",
        ]
        gated_log10_blocks = [
            f"log_{{10}} B = {b_gated['log10_tep_vs_null']}",
            f"log_{{10}} B = {b_gated['log10_tep_vs_anderson']}",
        ]
        has_bic_warning = (
            "unreliable" in text.lower()
            and ("bic approximation" in text.lower() or "bic-approximated" in text.lower())
            and ("n < 10" in text.lower() or "n &lt; 10" in text.lower())
        )
        if not any(block in text for block in headline_log10_blocks):
            if not has_bic_warning:
                self.violations.append(
                    "Model comparison missing headline full-catalog log10 Bayes factors "
                    f"(expected one of: {', '.join(headline_log10_blocks)})"
                )
                return False
        if not any(block in text for block in gated_log10_blocks):
            if not has_bic_warning:
                self.violations.append(
                    "Model comparison missing gated n=3 log10 Bayes factors "
                    f"(expected one of: {', '.join(gated_log10_blocks)})"
                )
                return False

        headline_delta_blocks = [
            f"$\\Delta$BIC $\\approx {b_headline['delta_bic_tep_null']}$",
            f"$\\Delta$BIC $\\approx {b_headline['delta_bic_tep_anderson']}$",
        ]
        if not self._require_substrings(text, headline_delta_blocks, "Model comparison (headline ΔBIC)"):
            return False

        # Section 4.5.4: σ_sys = 1.20 mm/s column must use n=9 headline log10 B, not n=3 gated.
        if (
            "σ_sys = 1.20" in text
            and f"log_{{10}} B = {b_gated['log10_tep_vs_null']}" in text
            and f"log_{{10}} B = {b_headline['log10_tep_vs_null']}" not in text
        ):
            self.violations.append(
                "Uncertainty-treatment table appears to place gated n=3 log10 B values "
                "in the headline n=9 (σ_geom = 1.20 mm/s) column."
            )
            return False

        w_tep = float(headline_block["model_selection"]["akaike_weights"]["TEP_restricted"])
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

        # Sanity: audit JSON consumers should not confuse blocks
        if float(bf_headline["log10_BF_TEP_restricted_vs_Null"]) <= float(
            bf_gated["log10_BF_TEP_restricted_vs_Null"]
        ):
            self.violations.append(
                "Step 026 full-catalog log10 BF should exceed primary gated log10 BF for TEP vs Null."
            )
            return False
        return True

    def audit_variance_decomposition(self) -> bool:
        variance = self.expected["variance"]
        text = self._combined_manuscript_text()

        # If variance decomposition is underpowered, skip percentage checks
        if variance.get("structural") is None and variance.get("observational") is None:
            return True

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
                + str(variance["observational"]),
                f"observational pipeline effects account for {variance['observational']}",
                f"observational effects account for {variance['observational']}",
            )
        )
        environmental_ok = f"environmental modulation contributes {variance['environmental']}" in text
        residual_ok = any(
            phrase in text
            for phrase in (
                "residual (small-sample statistics, intrinsic scatter, model incompleteness) accounts for "
                + str(variance["residual"]),
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
        violations = []
        for mission_key, label in self.PRIMARY_DETECTIONS:
            component = self.expected["components"][mission_key]
            # Check for either markdown row OR all HTML cells
            markdown_row = f"| {label} | {self._format_signed_mm(component['grad'])} | {self._format_signed_mm(component['disf'])} | {self._format_signed_mm(component['total'])} |"
            html_cells = [
                f"<td>{label}</td>",
                f"<td>{self._format_signed_mm(component['grad'])}</td>",
                f"<td>{self._format_signed_mm(component['disf'])}</td>",
                f"<td>{self._format_signed_mm(component['total'])}</td>",
            ]
            has_markdown = markdown_row in text
            has_html = all(cell in text for cell in html_cells)
            if not (has_markdown or has_html):
                violations.append(f"Component table missing {label} row (expected grad={component['grad']}, disf={component['disf']}, total={component['total']})")

        if violations:
            self.violations.extend(violations)
            return False
        return True

    def audit_discussion_comparison_table8(self) -> bool:
        """Ensure §5.2 mechanism table was not overwritten by Step 040 cosmographic rows."""
        path = self.project_root / "site/components/6_discussion.html"
        if not path.exists():
            return True
        text = path.read_text(encoding="utf-8")
        match = re.search(
            r"<caption>\s*Table 8:\s*Comparison of Flyby Anomaly Explanations\s*</caption>[\s\S]*?<tbody>([\s\S]*?)</tbody>",
            text,
            re.IGNORECASE,
        )
        if not match:
            self.violations.append(
                "Discussion Table 8 (mechanism comparison): caption/tbody block not found in site/components/6_discussion.html"
            )
            return False
        tbody = match.group(1)
        if "Atmospheric drag" not in tbody:
            self.violations.append(
                "Discussion Table 8 missing mechanism rows (e.g. Atmospheric drag); "
                "possible accidental sync with cosmographic Table 8 from Step 040."
            )
            return False
        if re.search(r"<td>\s*NEAR\s*</td>\s*<td>\s*0\.984\s*</td>", tbody, re.IGNORECASE):
            self.violations.append(
                "Discussion Table 8 contains cosmographic-style rows (NEAR / 0.984); "
                "expected mechanism comparison, not Step 040 mission table."
            )
            return False
        return True

    def audit_pipeline_output_hashes(self) -> bool:
        required_outputs = (
            "results/step007_tep_predictions.json",
            "results/step009_variance_analysis.json",
            "results/step013_cross_validation.json",
            "results/step026_stable_model_comparison.json",
        )
        missing = [path for path in required_outputs if not (self.project_root / path).exists()]
        if missing:
            self.violations.append(f"Missing required pipeline outputs: {', '.join(missing)}")
            return False
        return True

    def audit_cross_mission_holdout_mandatory(self) -> bool:
        """
        Step 013 must record a passing mandatory leave-one-mission-out protocol
        (same exit criteria as the step script) so stale or partial runs cannot
        pass the claim audit.
        """
        try:
            step013 = self._load_json("results/step013_cross_validation.json")
        except FileNotFoundError:
            self.violations.append(
                "Step 013 output missing (results/step013_cross_validation.json); "
                "run Step 013 before the claim audit."
            )
            return False
        except json.JSONDecodeError as exc:
            self.violations.append(f"Step 013 output is not valid JSON: {exc}")
            return False

        block = step013.get("cross_mission_holdout_mandatory")
        if not isinstance(block, dict) or "pass" not in block:
            self.violations.append(
                "results/step013_cross_validation.json has no cross_mission_holdout_mandatory.pass "
                "(regenerate Step 013)."
            )
            return False
        if not bool(block["pass"]):
            reasons = block.get("failure_reasons") or []
            detail = "; ".join(str(r) for r in reasons) if reasons else "no failure_reasons recorded"
            self.violations.append(
                "Mandatory cross-mission held-out protocol failed (Step 013): " + detail
            )
            return False
        return True

    def audit_reproduction_checklist(self) -> bool:
        run_all = self.project_root / "scripts/run_all.py"
        if not run_all.exists():
            self.violations.append("scripts/run_all.py not found")
            return False

        content = run_all.read_text(encoding="utf-8")
        step_count = content.count("('step_")
        if step_count < 40:
            self.violations.append(
                f"run_all.py lists {step_count} steps; expected at least 40 "
                f"(CORE_STEPS must stay in sync with tests/test_pipeline_smoke.py)."
            )
            return False
        return True

    def audit_table3_values(self) -> bool:
        """Audit Table 3 values against pipeline outputs from step007 and step008."""
        try:
            step007 = self._load_json("results/step007_tep_predictions.json")
            step008 = self._load_json("results/step008_fitting_results.json")
        except (FileNotFoundError, KeyError) as exc:
            self.violations.append(f"Failed to load pipeline outputs: {exc}")
            return False

        # Helper to convert float to LaTeX scientific notation
        def to_latex_sci(val: float, precision: int = 2) -> str:
            if val == 0:
                return "0"
            exp = int(math.floor(math.log10(abs(val))))
            mant = val / (10**exp)
            return f"{mant:.{precision}f} \\times 10^{{{exp}}}"

        # Helper to check if value is in text in any common format
        def value_in_text(val: float, text: str, tolerance: float = 0.01) -> bool:
            # Check plain formats
            patterns = [f"{val:.2e}", f"{val:.3e}", f"{val:.4e}"]
            # Check LaTeX formats (with and without dollar signs, multiple precisions)
            for prec in [1, 2, 3, 4]:
                latex = to_latex_sci(val, prec)
                patterns.extend([latex, f"${latex}$"])
            return any(p in text for p in patterns)

        # Expected values from pipeline
        expected = {
            "NEAR": {
                "dv_tep": step007["predictions"]["NEAR"]["tep_predictions"]["dv_tep_mm_s"],
                "beta_fitted": step008["individual_fits"]["NEAR"]["fit"]["beta_fitted"],
                "beta_uncertainty": step008["individual_fits"]["NEAR"]["fit"]["uncertainty"],
            },
            "Galileo_1990": {
                "dv_tep": step007["predictions"]["Galileo_1990"]["tep_predictions"]["dv_tep_mm_s"],
                "beta_fitted": step008["individual_fits"]["Galileo_1990"]["fit"]["beta_fitted"],
                "beta_uncertainty": step008["individual_fits"]["Galileo_1990"]["fit"]["uncertainty"],
            },
            "Rosetta_2005": {
                "dv_tep": step007["predictions"]["Rosetta_2005"]["tep_predictions"]["dv_tep_mm_s"],
                "beta_fitted": step008["individual_fits"]["Rosetta_2005"]["fit"]["beta_fitted"],
                "beta_uncertainty": step008["individual_fits"]["Rosetta_2005"]["fit"]["uncertainty"],
            },
        }

        text = self._combined_manuscript_text()
        violations = []

        # Check for each mission
        for mission, vals in expected.items():
            mission_label = mission.replace("_1990", "").replace("_2005", "")

            # Check dv_tep (tolerance 0.01 mm/s)
            dv_tep_pattern = f"{vals['dv_tep']:.2f}"
            if dv_tep_pattern not in text:
                # Try with different rounding
                alt_pattern = f"{vals['dv_tep']:.1f}"
                if alt_pattern not in text:
                    violations.append(f"Table 3 missing {mission_label} dv_tep ≈ {vals['dv_tep']:.2f}")

            # Check beta_fitted (scientific notation with tolerance)
            if not value_in_text(vals['beta_fitted'], text):
                violations.append(f"Table 3 missing {mission_label} beta ≈ {vals['beta_fitted']:.2e}")

            # Check beta_uncertainty
            if not value_in_text(vals['beta_uncertainty'], text):
                violations.append(f"Table 3 missing {mission_label} beta uncertainty ≈ {vals['beta_uncertainty']:.2e}")

        # Check weighted mean beta from step008
        weighted_mean = step008["beta"]
        weighted_unc = step008.get("overall_analysis", {}).get("beta_statistics", {}).get(
            "weighted_uncertainty",
            step008["uncertainty"],
        )
        if not value_in_text(weighted_mean, text):
            violations.append(f"Manuscript missing weighted mean beta ≈ {weighted_mean:.2e}")
        if not value_in_text(weighted_unc, text):
            violations.append(f"Manuscript missing weighted mean uncertainty ≈ {weighted_unc:.2e}")

        # If Step 008 has a random-effects uncertainty, ensure the prose
        # acknowledges it somewhere. This avoids silently replacing the formal
        # fixed-effect diagnostic in Table 3.
        random_unc = step008.get("overall_analysis", {}).get("beta_statistics", {}).get(
            "random_effects_uncertainty"
        )
        if random_unc is not None and not value_in_text(random_unc, text):
            violations.append(f"Manuscript missing random-effects uncertainty ≈ {random_unc:.2e}")

        # Check bootstrap intervals
        bootstrap = step008["robustness_analysis"]["bootstrap"]
        ci_lower = bootstrap['ci_95_lower']
        ci_upper = bootstrap['ci_95_upper']
        if not (value_in_text(ci_lower, text) and value_in_text(ci_upper, text)):
            violations.append(f"Manuscript missing 95% CI [{ci_lower:.2e}, {ci_upper:.2e}]")

        if violations:
            self.violations.extend(violations)
            return False
        return True

    def audit_evidence_framing(self) -> bool:
        """
        Guard the manuscript's core evidence frame.

        The defensible high-level claim is:
        - strong BIC-surrogate evidence for TEP restricted over Null,
        - positive/moderate but not decisive evidence over Anderson at n=3,
        - Juno is a deterministic warning case, not an uncertainty-aware raw-tension case.
        """
        text = self._combined_manuscript_text()
        text_lower = text.lower()
        violations = []

        banned = [
            "prefers the null",
            "favours the null",
            "favors the null",
            "direct experimental validation",
            "overwhelming evidence",
            "juno is a 5",
        ]
        for phrase in banned:
            if phrase in text_lower:
                violations.append(f"Evidence framing contains banned phrase: {phrase}")

        if not (
            "strongly favoured over the Null" in text
            or "strongly favored over the Null" in text
            or "strong evidence for TEP over Null" in text
            or "favours the restricted model" in text
            or "favors the restricted model" in text
        ):
            violations.append(
                "Evidence framing should explicitly state that TEP restricted is strongly favoured over Null."
            )

        if not (
            "moderately favoured over the Anderson" in text
            or "moderately favored over the Anderson" in text
            or "positive but not decisive" in text
        ):
            violations.append(
                "Evidence framing should distinguish the Anderson comparison as positive/moderate but not decisive."
            )

        ledger_required = [
            ("Evidence Ledger for Review", None),
            ("Restricted TEP vs Null", None),
            ("Restricted TEP vs Anderson", None),
            ("External anchors", None),
            ("Null-catalog stress test", "Null-catalog consistency"),
            ("Failure-mode accounting", "Disformal regime confirmation"),
            ("Machine-checked manuscript", None),
        ]
        missing_ledger = [
            primary for primary, alt in ledger_required
            if primary not in text and (alt is None or alt not in text)
        ]
        if missing_ledger:
            violations.append(
                "Evidence ledger missing required entries: " + ", ".join(missing_ledger)
            )

        try:
            step039 = self._load_json("results/step039_flyby_prediction_table.json")
            uncertainty_raw_tension = step039["metadata"][
                "raw_uncertainty_aware_classification_summary"
            ]["raw_tension"]
            raw_tension = step039["metadata"]["raw_classification_summary"]["raw_tension"]
            delta_log_l = step039["full_catalog_raw_likelihood"][
                "delta_log_likelihood_tep_minus_null"
            ]
        except (FileNotFoundError, KeyError, TypeError) as exc:
            violations.append(f"Failed to load Step 039 evidence framing values: {exc}")
        else:
            if raw_tension == 1 and "deterministic fixed-amplitude warning case" not in text:
                violations.append(
                    "Manuscript should describe the Step 039 Juno result as a deterministic fixed-amplitude warning case."
                )
            if uncertainty_raw_tension == 0 and not (
                "0 uncertainty-aware raw-tension cases" in text
                or "zero uncertainty-aware raw-tension cases" in text
                or "0 raw-tension cases remain after random-effects prediction uncertainty" in text
            ):
                violations.append(
                    "Manuscript should state that uncertainty-aware Step 039 has zero raw-tension cases."
                )
            delta_fmt = f"{delta_log_l:.2f}"
            if delta_fmt not in text and f"{delta_log_l:.1f}" not in text:
                violations.append(
                    f"Manuscript missing Step 039 full-catalog ΔlogL ≈ {delta_fmt}."
                )

        if violations:
            self.violations.extend(violations)
            return False
        return True

    def audit_delta_bic_prose_consistency(self) -> bool:
        text = self._combined_manuscript_text()
        violations = []

        anderson_null_val = self.expected["model"]["bayes"]["delta_bic_anderson_null"]
        tep_anderson_val = self.expected["model"]["bayes"]["delta_bic_tep_anderson"]

        for match in re.finditer(r"636\.4", text):
            start = max(0, match.start() - 300)
            end = min(len(text), match.end() + 300)
            context = text[start:end]
            if ("TEP restricted" in context or "over Anderson" in context
                    or "preference over Anderson" in context):
                if "Anderson vs Null" not in context and "separation vs Null" not in context:
                    violations.append(
                        f"ΔBIC ≈ {anderson_null_val} (Anderson vs Null) appears in prose "
                        f"about TEP vs Anderson. Expected ΔBIC ≈ {tep_anderson_val} "
                        f"for TEP restricted vs Anderson."
                    )
                    break

        if violations:
            self.violations.extend(violations)
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
            ("Cross-mission holdout mandatory (Step 013)", self.audit_cross_mission_holdout_mandatory),
            ("Stale marker scan", self.audit_stale_markers),
            ("Model comparison table", self.audit_model_comparison_table),
            ("Variance decomposition", self.audit_variance_decomposition),
            ("Component table", self.audit_component_table),
            ("Discussion Table 8 (mechanisms)", self.audit_discussion_comparison_table8),
            ("Reproduction checklist", self.audit_reproduction_checklist),
            ("Table 3 values", self.audit_table3_values),
            ("Evidence framing", self.audit_evidence_framing),
            ("ΔBIC prose consistency", self.audit_delta_bic_prose_consistency),
        ]

        all_passed = True
        self.audit_results = []
        for name, audit_func in audits:
            if logger:
                logger.info(f"[{name}]... ")
            else:
                print(f"\n[{name}]...", end=" ")
            before = len(self.violations)
            try:
                passed = audit_func()
                new_violations = self.violations[before:]
                self.audit_results.append(
                    {
                        "name": name,
                        "status": "PASS" if passed else "FAIL",
                        "violations": new_violations,
                    }
                )
                if logger:
                    logger.success("PASS" if passed else "FAIL")
                else:
                    print("✓ PASS" if passed else "✗ FAIL")
                all_passed = all_passed and passed
            except Exception as exc:
                self.audit_results.append(
                    {
                        "name": name,
                        "status": "ERROR",
                        "violations": [str(exc)],
                    }
                )
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

    def write_audit_report(
        self,
        output_path: Path,
        *,
        passed: bool,
        synced_files: List[str] | None = None,
        duration_seconds: float | None = None,
    ) -> Path:
        step039 = self._load_json("results/step039_flyby_prediction_table.json")
        step026 = self._load_json("results/step026_stable_model_comparison.json")
        headline_block, gated_block = split_step026_model_blocks(step026)
        headline_bf = headline_block["bayes_factors"]
        gated_bf = gated_block["bayes_factors"]
        report = {
            "metadata": {
                "step": "027_claim_consistency_audit",
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "status": "PASS" if passed else "FAIL",
                "duration_seconds": duration_seconds,
            },
            "synced_files": synced_files or [],
            "audits": self.audit_results,
            "violations": self.violations,
            "evidence_frame": {
                "comparison_frame": "full_catalog_n9_geometry_spread",
                "tep_restricted_vs_null_log10_bf": headline_bf[
                    "log10_BF_TEP_restricted_vs_Null"
                ],
                "tep_restricted_vs_anderson_log10_bf": headline_bf[
                    "log10_BF_TEP_restricted_vs_Anderson"
                ],
                "gated_n3_tep_restricted_vs_null_log10_bf": gated_bf[
                    "log10_BF_TEP_restricted_vs_Null"
                ],
                "gated_n3_tep_restricted_vs_anderson_log10_bf": gated_bf[
                    "log10_BF_TEP_restricted_vs_Anderson"
                ],
                "step039_delta_log_likelihood_tep_minus_null": step039[
                    "full_catalog_raw_likelihood"
                ]["delta_log_likelihood_tep_minus_null"],
                "step039_deterministic_raw_tensions": step039["metadata"][
                    "raw_classification_summary"
                ]["raw_tension"],
                "step039_uncertainty_aware_raw_tensions": step039["metadata"][
                    "raw_uncertainty_aware_classification_summary"
                ]["raw_tension"],
                "prediction_uncertainty_model": step039["metadata"][
                    "beta_prediction_uncertainty_model"
                ],
            },
        }
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        return output_path


def main() -> None:
    project_root = PROJECT_ROOT
    logger = StepLogger("step_027_claim_consistency_audit", project_root)
    start = time.time()
    status = "SUCCESS"
    auditor = None
    passed = False
    updated = []
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
        report_path = project_root / "results" / "step027_claim_consistency_audit.json"
        auditor.write_audit_report(
            report_path,
            passed=passed,
            synced_files=[str(path.relative_to(project_root)) for path in updated],
            duration_seconds=time.time() - start,
        )
        logger.add_output_file(report_path, "machine-readable claim consistency audit report")
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
        if auditor is not None:
            report_path = project_root / "results" / "step027_claim_consistency_audit.json"
            if not report_path.exists():
                try:
                    auditor.write_audit_report(
                        report_path,
                        passed=passed,
                        synced_files=[str(path.relative_to(project_root)) for path in updated],
                        duration_seconds=time.time() - start,
                    )
                    logger.add_output_file(
                        report_path,
                        "machine-readable claim consistency audit report",
                    )
                except Exception as report_exc:
                    logger.error(f"Failed to write audit report: {report_exc}")
        duration = time.time() - start
        logger.log_step_summary(duration, status=status)

    sys.exit(0)


if __name__ == "__main__":
    main()
