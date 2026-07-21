# pydna architecture review

Scope: package *structure*, *maintainability*, and *hygiene* — plus a verified read on
the core *algorithms*. Not correctness of the science, not dev tooling.

Every claim below carries a `file:line` or a number checked directly against the tree.

## TL;DR

The algorithms are sound (core suite passes; hot path delegated to a C suffix-array). The
problems are one level up — organization, not math.

- **Structure** — 37 modules flat in one directory; a 44-line `types.py` leaf drags CRISPR
  into the import root, and 54 in-function imports hide real dependency cycles. Fix: layer
  the tree (`core / io / assembly / methods / pcr / …`) and enforce acyclicity in CI.
- **Maintainability** — god-classes (`Dseq` 51 methods, `Dseqrecord` 41), a 3,285-line
  `assembly2.py` that conflates every cloning method with no registry, and CI logic
  duplicated across four workflows with no single task entrypoint.
- **Hygiene** — `common_sub_strings.py` is 74% commented-out dead code; plus a dead module,
  a stale Plotly block, and inconsistent deprecation warnings.

---

## 1. Structure — where code lives, what imports what

**37 modules, 21,435 lines, flat in a single directory.** No layering, so everything can
import everything — and does.

**Dependency tangle:**

- **`types.py` (44 lines) is not a leaf.** It runs `from pydna.crispr import _cas`
  (`types.py:17`). Every core module imports `types` for annotations, so
  `import pydna.dseq` → `pydna.types` → `pydna.crispr` at runtime: a type-alias leaf
  drags a *feature* module into the import root of the core sequence class.
- **54 `import pydna.*` statements sit inside functions** — deferred imports that hide
  existing cycles (`dseqrecord ↔ assembly2`, `opencloning_models ↔ dseqrecord`).
- **`utils.py` reaches upward** — imports `alphabet, codon, types`, and via `types`
  transitively reaches `crispr`. A bottom-layer helper depending on the top layer.
- **`opencloning_models.py` (1,421 lines) is a bidirectional hub** — imports seven
  modules and is imported back by six.

**Proposed target tree** — re-home existing modules into layers; break cycles in the
process. Not a rewrite.

```
pydna/
  types.py   # pure aliases — zero pydna imports (breaks the crispr cycle)
  core/      seq, dseq, dseqrecord, seqrecord, amplicon, contig, alphabet, codon
  io/        parsers, readers, genbank, genbankfixer, snapgene_history_parser
  assembly/  graph (assembly2 algorithm), common_sub_strings   # algorithm only
  methods/   gibson, goldengate, gateway, crispr, cre_lox, recombinase, fusionpcr  # + registry
  pcr/       amplify, design, primer, tm, oligonucleotide_hybridization, primer_screen
  viz/       gel, ladders, _pretty
  models/    opencloning_models
  utils.py
```

**Invariant to enforce (CI via `import-linter`):** `core` and `types` import nothing
above them; `methods` never import each other; no cycles. This single rule is the
difference between a navigable tree and the current tangle.

---

## 2. Maintainability — the API a contributor must work against

- **God-classes.** `Dseq` carries **51 methods** in 3,169 lines (`dseq.py`); `Dseqrecord`
  carries **41** in 1,639 (`dseqrecord.py`). Neither fits in a contributor's head.
- **No method registry.** `assembly2.py` (3,285 lines) conflates the assembly *graph
  algorithm* with *every cloning technique*, importing `cre_lox, crispr, gateway,
  recombinase`. Adding a method means editing the god-module rather than dropping a file
  into `methods/`.
- **Ambiguous public surface.** Five overlapping sequence types with no namespace:
  `Seq`, `Dseq`, `SeqRecord`, `Dseqrecord`, `FakeSeq`. Nothing signals which is canonical.
- **No `__all__` discipline.** `assembly2.py` exposes **56 top-level functions**; the
  public API is "whatever isn't underscore-prefixed."
- **Versioned module name in a public tree.** `assembly.py` and `assembly2.py` both define
  a class named `Assembly`; `pydna.all:48` still exports the **deprecated** one.
- **No single task entrypoint.** Project verbs are scattered across `run_test.py`, the
  `.github/workflows/` files, and dotfiles — there is no `Makefile` or task runner, so a
  contributor cannot discover how to test/lint/build from one place. Add universal targets
  (`setup / test / test-all / lint / fmt / build / clean`).
- **CI/CD duplicates logic and pins toolchain inline.** Four GitHub Actions workflows
  (test+coverage, PyPI publish, TestPyPI publish, docs) over a 5-Python × 3-OS matrix. Test
  invocation is duplicated between `run_test.py` and the workflow rather than sharing one
  entrypoint, and the toolchain version is pinned inside the YAML
  (`pydna_test_and_coverage_workflow.yml`). Consolidate CI onto the single task entrypoint
  so local and CI run identical commands.

---

## 3. Hygiene — removable cruft and inconsistency

- **`common_sub_strings.py` is 74% dead code** — lines 21–312 are an abandoned pure-Python
  suffix-array implementation, fully commented out. The live algorithm is 30 lines.
- **`threading_timer_decorator_exit.py` (109 lines)** is imported only by the deprecated
  `assembly.py` — dead once that module goes.
- **`contig.py:400–465`** — ~70 lines of commented-out Plotly with `FIXME "I see no reason
  for it"`.
- **Inconsistent deprecation class** — some sites raise stdlib `DeprecationWarning`
  (`assembly.py:69`, `design.py:786`), others `_PydnaDeprecationWarning` (`dseqrecord.py:444`).
- **`Dseq.join()` (dseq.py:3082)** — O(n²) fold (`result = result + self + x`); harmless at
  small N, quadratic at large.
