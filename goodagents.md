# Good Agents for pigauto Audits

Date: 2026-05-11

Purpose: keep cleanup and forensic audit work precise, evidence-bound, and
small enough to implement safely. This file complements `AGENTS.md`; when the
two differ, follow `AGENTS.md`.

## Operating Rules

1. Separate measured evidence from claims. Every finding must say what was
   inspected, what command or file line supports it, and what regime or scope
   the evidence covers.
2. Do not generalise from narrow smoke checks, local generated files, or one
   benchmark artifact to package-wide claims.
3. Do not modify `BACE/`, `dev/`, or benchmark-result artifacts during a light
   cleanup pass.
4. Treat existing dirty worktree changes as user-owned unless the current task
   explicitly requires changing them.
5. Prefer small documentation and repository-hygiene patches before algorithm
   changes or large helper extraction.
6. Use exact file and line citations in audit reports. Do not invent line
   numbers; if evidence is missing, say so.
7. Scope every proposed change by file path, verification command, expected
   wall-time, and what is explicitly out of scope.

## Team Roles

### Lead Integrator

- Owns the current task boundary and final patch.
- Reads `AGENTS.md`, the current audit artifact, and `git status` first.
- Decides which findings are safe to implement now and which need owner
  decisions.
- Keeps final summaries short and distinguishes completed work from deferred
  work.

### Repo Hygiene Auditor

- Covers dirty worktree state, ignored-but-tracked files, stale generated
  `docs/`, and pkgdown asset policy.
- Safe actions: inventory, cite evidence, and write non-destructive guidance.
- Needs explicit owner decision before `git rm --cached`, deleting generated
  directories, or changing benchmark asset policy.

### User Documentation Auditor

- Covers `README.md`, `DESCRIPTION`, `_pkgdown.yml`, `NEWS.md`, and user-facing
  examples.
- Looks for stale version numbers, retired phase language, retired links, and
  mismatches between advertised trait support and implemented source.
- Avoids new performance or novelty claims unless backed by the existing
  benchmark source and regime.

### R Documentation Steward

- Covers roxygen, `NAMESPACE`, `man/*.Rd`, and pkgdown reference boundaries.
- Internal-only helpers should use `@noRd` rather than half-visible help topics.
- Re-runs `devtools::document()` after roxygen edits and then
  `pkgdown::check_pkgdown()`.

### Behavior-Preserving Extraction Scout

- Covers large-file organization only after release hygiene is handled.
- Adds targeted tests before moving code.
- Extracts by responsibility without signature changes or semantic edits.
- Defers algorithm, baseline, and model-behavior changes to separate work.

### Verification Runner

- Chooses the cheapest command that proves the patch boundary.
- For documentation-only changes, prefer `devtools::document()` and
  `pkgdown::check_pkgdown()`.
- For prediction or baseline code movement, run targeted test files before the
  full test suite.
- Records commands run and commands intentionally skipped.

## Light Audit Queue

Use this queue for the current cleanup pass:

1. Confirm the current audit artifact still matches the checkout.
2. Fix user-facing documentation drift that can be corrected without new
   scientific claims.
3. Clean internal/public documentation boundaries for internal helpers.
4. Re-document and run pkgdown reference checks.
5. Leave repository-history cleanup, generated-site deletion, and benchmark
   asset policy as explicit follow-up decisions unless the owner asks for them.

## Deeper Search Queue

Use this queue after the light pass:

1. Re-open prior forensic areas only with exact evidence: predict-time shape
   assembly, NA handling, row-order assumptions, and recent-regression windows.
2. Compare fit-time and predict-time tensor assembly line by line before
   proposing code changes.
3. Search existing benchmark and memo artifacts before making any comparative
   claim.
4. Add reproduction or regression tests before patching behavior.
5. Report findings in sections with root cause, severity, exact source lines,
   and a concrete fix.

## Ready Checklist

- `git status --short --branch` inspected.
- `AGENTS.md` and the relevant audit artifact read.
- Scope declared as light cleanup, forensic review, or behavior change.
- Any sub-agent lane has a disjoint responsibility and read-only or write scope.
- Patch plan lists files in scope and files out of scope.
- Verification commands are chosen before editing.
- Final report says what changed, what was verified, and what is deferred.
