# Strand7Tools

Reads Strand7 (R3.1.1) beam results, builds design members from chains of beams, and checks them
to EN 1993-1-1 (UK National Annex): cross-sections (6.2) and member buckling (6.3.1, 6.3.2.2).
The user is a structural engineer: every result must be checkable by hand, through the worked
calculations.

## Layout

- `run_ec3.py`: the main entry point. Set PROJECT and the worked calculations at the top, or
  run `python run_ec3.py Test\test_model3.toml`.
- `strand7tools/`: the package.
  - `ec3/`: pure calculation (no Strand7 and no printing).
    - `sections.py`: the **Section interface** and `ISection`, plus `Steel`.
    - `cross_section.py`: 6.2.
    - `buckling.py`: 6.3, with Mcr/C1 from the moment diagram.
    - `design.py`: combining cases, per-member summaries, and status ranking.
  - `strand7/`: everything that calls the Strand7 API (read-only).
    - `api.py`: loads the API lazily, so the rest imports without Strand7.
    - `members.py`, `forces.py`, `sections.py`.
  - `app.py`: `run_checks(project)` returns a `RunResult`. It does all the calculation and
    prints nothing, so a future GUI can use it.
  - `report.py`: console tables and CSV files.
  - `worked.py`: worked-calculation text built from result objects.
  - `project.py`, `settings.py`, `design_cases.py`, `member_table.py`: the inputs.
  - `data/settings.toml`, `data/steel_grades.csv`: the default design settings (UK NA values).
- `Test/`: the example project folder. It holds the project files `test_model2.toml` and
  `test_model3.toml` (portal frame), the member tables and design cases.
  - The `.st7` and result files are gitignored.
  - Output goes to `Test/output*/`, which is also gitignored.
- `tools/`: helper scripts. They take an optional project file argument.
  - `list_members`, `list_forces`, `list_sections`, `show_inputs`.
  - `check_portal.py`: compares the force assembly with an independent frame solver on
    `test_model3.st7`.
- `tests/`: pytest. Needs no Strand7.

## Commands

```
python -m pytest                      # unit tests (no Strand7 needed)
python run_ec3.py                     # checks for Test/test_model2.toml (Strand7 must not have the model open)
python tools/check_portal.py          # independent check of the force assembly (needs test_model3.st7 + .LSA)
```

## Rules

- **Never modify or save a Strand7 model.** Open models read-only only.
- **Nothing is passed silently.**
  - Anything not implemented is reported in `not_checked`, which makes the result INCOMPLETE
    or "FAILS (incomplete)". It is never OK.
  - Informational caveats go in `notes`.
- **Calculation and presentation are separate.**
  - `ec3/` returns result objects (`util`, `values`, states). It never prints.
  - `worked.py` turns those results into text.
  - When you add a check, add both, plus a test.
- **Units.**
  - `ec3/` uses SI throughout: N, m, Pa, N·m.
  - Member force diagrams are in kN and kNm.
  - Models must be in metres, with forces in N, kN or MN.
  - E and G are converted from the model's stress unit.
- **Axes.**
  - y-y is the major axis = Strand7 principal axis 1, so My = BM2 and Vz = SF2.
  - z-z is the minor axis, so Mz = BM1 and Vy = SF1.
  - Forces are in member axes, tension positive.
  - The beam-to-member sign factors are in `strand7/members.py`. They are verified by
    `check_portal.py`.
- **E and G come from the Strand7 material by default** (`source = "model"`). This keeps the
  analysis and the design consistent. A warning appears if they differ from 210/81 GPa.
- **New profiles:** subclass `Section` in `ec3/sections.py` and implement every abstract
  method. Then extend `strand7/sections.py` (`read_section` and `make_section`) and the
  shape-specific text in `worked.py`. Leave `cross_section.py` and `buckling.py` alone unless a
  rule is not yet behind the interface.
- **Dependencies.** Pure Python 3.13 with no numpy. Only pytest is used for tests.
- **Before and after a refactor**, run `run_ec3.py` on both projects with every worked
  calculation switched on, and diff the output. The numbers and the text should not change.
- **Git.** Commit as PBRimbault with the GitHub noreply email. Commit or push only when asked.

## Open points

- 6.3.3 (Annex B, Method 2 per the UK NA) is not built yet. Compression with bending is reported
  as INCOMPLETE.
- C1 for "direct" design cases comes from Strand7's own stations, and is about 1% high.
- The UK NA values for 6.3.2.2 have not been confirmed.
- Both test models use E = 200 GPa and G = 76.9 GPa, so each run gives a warning. This is
  intentional: the model values are used.
