<!--
Thanks for contributing to the K* verification suite!

Please confirm each of the checkboxes below before requesting review.
-->

## Summary

<!-- Two-to-three sentences describing the change and its motivation. -->

## Type of change

- [ ] Bug fix (non-breaking)
- [ ] New verification check or Lean theorem
- [ ] Machine-generator / tooling improvement
- [ ] Documentation only
- [ ] Registry / manifest data change

## Verification checklist

- [ ] `lake build` still exits 0 (Lean4).
- [ ] `python lean4/scripts/check_sorry.py --layer1` passes.
- [ ] `python lean4/scripts/generate_all.py` regenerates without drift.
- [ ] `python verify_registry.py` passes or diff is documented.
- [ ] `pdflatex manuscript.tex` and `pdflatex supplemental_material.tex`
      build clean (if LaTeX changed).
- [ ] If any universal-n Lean theorem added or modified, axiom footprint
      was confirmed via `#print axioms` and the registry updated.

## Backwards-compatibility

<!--
If this PR changes the schema of verification_manifest.json, axiom_report,
or proofs_registry.yaml, describe the migration path.  Otherwise leave as
"No schema change."
-->

No schema change.

## Linked issues

<!-- Closes #N, Relates to #M, etc. -->
