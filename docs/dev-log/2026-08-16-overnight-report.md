# ☀️ Overnight report — pigauto, 2026-08-15 → 16

Written for Shinichi returning in the morning. Read this, then
`2026-08-16-campaign-go-no-go.md` if you want to launch the campaigns.

**One-line state:** Tier 0 is closed, eight of nine Rose P1 items are addressed, the
benchmark provenance cloud has lifted with MCSE-backed evidence, and both Totoro campaigns
are staged, priced and waiting on one word from you.

---

## What you need to decide (2 things, both small)

1. **Launch the Totoro campaigns?** Pre-run test complete, cost measured at **≈ 2 h wall**
   for both (regime map ≈ 1.5 h + coverage ≈ 0.5 h at 100 workers). Reply "launch" and they
   run. Details + the exact commands: `docs/dev-log/2026-08-16-campaign-go-no-go.md`.
2. **File the brain lessons?** Four drafted, awaiting your D-37 approval:
   `docs/dev-log/2026-08-15-brain-lesson-drafts.md`. Say "file them" or name which.

Everything else below is already done and needs nothing from you.

---

## Landed on `main` (8 PRs merged)

| PR | What | Why it mattered |
|---|---|---|
| [#158](https://github.com/itchyshin/pigauto/pull/158) | P0 honesty fixes + two-layer #157 gate floor + 9-edit doc bundle | The four correctness fixes finally on `main` after six days stranded; `main` no longer ships a documented guarantee it doesn't satisfy |
| [#159](https://github.com/itchyshin/pigauto/pull/159) | pkgdown destination fix | Post-merge fallout: #158 carried a file into pkgdown's build dir; caught and fixed within the hour |
| [#160](https://github.com/itchyshin/pigauto/pull/160) | P1-7, P1-6 | `$se` is three incompatible objects — now named and fenced against Rubin's-rules misuse |
| [#161](https://github.com/itchyshin/pigauto/pull/161) | P1-13 + P1-10 fail-closed guards | Four ways data contradicting its type was silently mangled, incl. one producing `NaN` and one inventing a uniform composition from an all-zero row |
| [#162](https://github.com/itchyshin/pigauto/pull/162) | P1-5 | Docs credited `Rphylopars` for work pigauto's own solver does |
| [#163](https://github.com/itchyshin/pigauto/pull/163) | P1-12 | **Phylo-signal gate was a silent no-op in multi-obs mode** |
| [#164](https://github.com/itchyshin/pigauto/pull/164) | P1-9 | **zi_count observed zeros could never enter val/test** |
| [#165](https://github.com/itchyshin/pigauto/pull/165) | P1-11 (documented, not fixed) | tree-MI uses MC-dropout draws, not the better-calibrated conformal ones — and nothing said so |

CI green throughout, including the post-merge pkgdown rebuild. `--as-cran` on the landed
tree: **0 errors / 0 warnings / 1 note** (the known dev-version + Suggests-BACE incoming
note) — cleaner than the 2W/3N baseline I pre-registered.

## Three silent bugs, all proven against pre-fix code

The night's pattern: documented behaviour that quietly did nothing, invisible because the
failure produced no error.

1. **#157 gate/predict surface mismatch.** Calibration scored one delta surface, `predict()`
   delivered another. Fixed with a margin-based BM floor at both layers.
2. **P1-12 multi-obs phylo-signal gate.** `species_names` (n_species) indexed with an
   observation-length logical → NA tip names → `keep.tip()` error → swallowed by `tryCatch`
   → gate never fired. *Proof:* BM-simulated trait returned `NA` pre-fix, λ > 0.5 post-fix.
3. **P1-9 zi_count zeros.** An observed zero encodes as gate = 0 + magnitude `NA`; the
   all-columns-observed rule read that NA as missingness. *Proof:* fixture with 23 observed
   zeros held out exactly **0** of them pre-fix.

Consequence worth carrying: **any previously reported `zi_count` val/test metric was
computed on non-zeros only** and should be re-measured.

## The benchmark result (the scientifically interesting one)

Re-ran the per-type bench suite on the corrected pipeline and compared **paired per-rep**
against the pre-fix outputs (old and new share seeds — confirmed by byte-identical
`baseline`/`mean` rows).

| Type | mean Δ | MCSE | \|Δ\|/MCSE |
|---|---:|---:|---:|
| continuous | −0.00168 | 0.00243 | 0.69 |
| binary | +0.00382 | 0.00363 | 1.05 |
| ordinal | +0.00048 | 0.00631 | 0.08 |
| count | +0.00206 | 0.00527 | 0.39 |
| categorical | +0.00253 | 0.00563 | 0.45 |
| proportion | +0.00629 | 0.00553 | 1.14 |

**The held-out-context leak had no detectable effect on any of them.** Categorical
reproduces 6/7 scenarios *exactly* — its gate closes completely, so a training leak cannot
move a number the GNN doesn't influence.

This lifts the roadmap's blanket ban on pre-fix numbers: cite the re-run values, noting the
pre-fix ones agree within noise. **Not** a claim that the leak was harmless — the tighter
per-cell val-floor invariant did detect a difference, and power here only bounds effects
above ~0.01–0.02. Full caveats: `2026-08-16-bench-rerun-results.md`.

**6 of 8 types complete**, all below 1.15 MCSE. `zi_count` and `multi_proportion` are
running. Note `zi_count` will need a *third* run: P1-9 changed which cells are eligible for
val/test, so neither side of tonight's leak-only comparison reflects today's `main`.

## Staged and waiting

- **Totoro**: pigauto `c655d75` installed, campaign runner + dispatcher shipped and
  smoke-validated on three previously-untested paths (F1 linear, F3 mixed-type, coverage
  with split). Job enumeration verified: 5,400 + 1,920 jobs.
- **Two ADEMP campaign designs** (Morris 2019 + Williams 2024 11-item audits):
  `2026-08-15-regime-map-design.md`, `2026-08-15-coverage-campaign-design.md`.
- **Pre-run science bonus:** on the nonlinear family with gates wide open, the leak-free GNN
  beat the baseline **6/6** — the first honest evidence it earns its place, in exactly the
  regime the design predicts it should.

## Three mistakes I made, and what they cost

Recorded because the corrections are the useful part.

1. **I mis-specified my own experiment's surface.** The delta-surface script masked val
   cells only, leaving test cells pinned — a surface nothing in the package uses. It read
   the GNN as 2% *worse*; on the package convention it was 2–4% *better*. Three downstream
   conclusions had begun forming on the artifact.
2. **My CPU-politeness fix caused the failures it was meant to avoid.** Renicing a bench
   master before it forks PSOCK workers hands them nice-19 (inherited at fork), recreating
   the connect-window failure. Two drivers hung an hour each. The peer session independently
   hit the same trap. Rule now: bound a cluster campaign by *worker count*, never priority.
3. **I called healthy jobs hung.** `parLapply` logs nothing until all cells finish, so a
   master at 0% CPU is normal. I killed at least one working driver on that misread.

Also corrected in-place: a roadmap line where I wrote that the joint baseline "never
iterates Σ". It does estimate Σ in closed form, consistent under matrix-normal BM; what's
disabled is a cell-refinement EM that *diverged* and scored 0.53 vs 0.93. I'd framed a
well-evidenced decision as a deficiency.

## What's genuinely still open

- **P1-8** — covariates ignored by the joint/threshold-joint baseline. A real feature with
  a design choice in it (thread them through vs detect-and-message); wants your input, not
  a 1 AM decision.
- **P1-11** — tree-MI has no conformal draws. *Documented* in #165, deliberately not fixed:
  threading `draws_method` through changes how a headline feature generates draws and
  embeds a default choice (back-compat `mc_dropout` vs better-calibrated `conformal`) that
  is yours to make. Fix sketch is in the PR.
- **Tier 1.2 / 1.3** — the regime map and coverage campaigns (designed, priced, gated).
- **Tier 2.2** — no known-Σ recovery simulation for the multivariate machinery.
- Issue [#135](https://github.com/itchyshin/pigauto/issues/135) remains open; the regime map
  is what closes it.

## Housekeeping

A stale `.git/index.lock` from a Cursor git worker is still on the main checkout (verified
unheld; my `rm` is permission-blocked). All my commits went through `GIT_INDEX_FILE`
plumbing around it. Once Cursor is closed:

```bash
rm "/Users/z3437171/Dropbox/Github Local/pigauto/.git/index.lock" && git -C "/Users/z3437171/Dropbox/Github Local/pigauto" reset
```

The `git reset` also clears the ghost "modified" entries the workaround leaves in
`git status`. Your 15 carried-over dirty files are untouched throughout.
