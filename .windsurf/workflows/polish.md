---
description: Polish and finalize manuscript before PDF generation
---

# Polish Workflow

Run this before PDF generation to ensure all metadata and content is consistent.

## CITATION.cff Metadata (REQUIRED before PDF generation)

Verify and update CITATION.cff with latest manuscript metadata:

1. **Read CITATION.cff** and check all fields are current:
   ```bash
   cat CITATION.cff
   ```

2. **Verify these fields match the manuscript**:
   - `title`: Must match manuscript title exactly
   - `version`: Format as `vX.Y (Codename)` - e.g., `v0.1 (Sintra)`
   - `date-released`: Update to current date if modified
   - `doi`: Must match Zenodo DOI
   - `abstract`: Must match manuscript abstract (first ~500 chars)
   - `keywords`: Must include all relevant topic keywords

3. **Check consistency across files**:
   - CITATION.cff `title` ↔ manuscript title
   - CITATION.cff `abstract` ↔ zenodo.txt abstract
   - CITATION.cff `version` ↔ PDF filename version
   - CITATION.cff `doi` ↔ zenodo.txt DOI

4. **Verify version in CITATION.cff** - read the current version, do not auto-increment:
   - Format should be: `vX.Y (Codename)` - e.g., `v0.1 (Sintra)`
   - Only manually change version when manuscript content changes significantly
   - Keep codename consistent for the same paper iteration
   - The PDF generator reads this version directly from CITATION.cff

## Manuscript Polish Checklist

1. Do not use 'we'.
2. Ensure appropirate humble tone throughout, no overly self-aggrandizing language or hyperbole.
3. Make sure the paper is storng and presents a good compelling argument.
4. Ensure the writing is logical and flows beautifully
5. The manuscirpt should be a pleasure to read and easy to understand and engaging.
6. Ensure consistency throughout paper
7. Ensure consistency throghout entire codebase
8. The zenodo.txt and readme must match the absract from the manuscript
9. we should not set timelines for future years and future research. we are one independent researcher without insitutional backing so should not set unachievable plans. 
10. Nothing should be fabricated, fake or made up. 
11. Everything must be real. 
12. All data and numbers must be consistent with analysis pipeline outputs.
13. **Integrity**: Verify every citation exists and supports the specific claim made (no hallucinated references).
14. **Integrity**: Ensure all data is correct and any referenced values are consistent with literature.
15. **Integrity**: Ensure no text is copied from other sources without proper attribution.
16. **Writing**: Remove filler words and "weasel words" (e.g., "clearly", "obviously", "very") - let the data speak.
17. Verify that units and notation are used consistently throughout and consistent with other /manuscripts
18. Ensure all references are correct and up to date and all titles and meta data are correct. 
19. Ensure there is no bold or strong midway through a sentence. 

## After Polish - Generate PDF

Once CITATION.cff is verified and manuscript is polished:

```bash
/generate-pdf
```

ALWAYS MAKE CHANGES TO THE MANUSCRIPT SOURCE @COMPONENTS HTML FILES AND NOT THE MARKDOWN DIRECTLY. 