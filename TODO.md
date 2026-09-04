# TODO

Open work on this pipeline. Status tags: **OPEN** / **DEFERRED** (deliberately not
now) / **BLOCKED** (waiting on something) / **WATCH** (verify, not a code change).

---

## Licensing: user-provided tools — kmernorm, Pheniqs

Two tools this pipeline runs cannot be safely bundled or redistributed in a
container image:

- **kmernorm** (Mingkun Li) — **DONE** (2026-09-04, `main.nf` commit `ea0dac2`).
  No license anywhere (SourceForge page or source archive), i.e.
  all-rights-reserved by default copyright, and `brwnj/kmernorm:v1.0.0` had no
  published Dockerfile, build provenance, or license either. `KMERNORM_v1_0_0`
  now declares `container null` — the public repo no longer bundles or
  references any container for it. Installing kmernorm and wiring it up (e.g. a
  per-process `module`/`beforeScript` in your own gitignored `nextflow.config`)
  is the user's own responsibility now. charlie's own working setup still runs
  it, via the natively-installed `kmernorm/1.0.5` module referenced only from
  the gitignored `nextflow.config` — not exposed on the remote.
- **Pheniqs** (`biosails/pheniqs`) — **OPEN**. The repository's top-level
  `LICENSE` is a restrictive NYU research license (internal/non-commercial
  only; no redistribution, modification, or sublicensing without a signed
  agreement), which conflicts with the AGPL-3.0-or-later headers on individual
  source files. `PHENIQS_SAMPLE_DEMULTIPLEX` / `PHENIQS_DEMULTIPLEX` currently
  reference `quay.io/biocontainers/pheniqs:2.1.0--py39ha79081e_6` (a
  Bioconda/biocontainers-built image, not `biosails/pheniqs` directly, but
  still built from that same restrictively-licensed source).

**Plan for Pheniqs, matching what was just done for kmernorm:** require the
user to install Pheniqs themselves, on `PATH` / in their own env, entirely at
their own discretion — the pipeline does not bundle, vouch for, or warrant it.
This is the approach the sibling public repo
**[BigelowLab/GORG-Dark-SAG-assembly](https://github.com/BigelowLab/GORG-Dark-SAG-assembly)**
already takes for kmernorm.

1. **`README.md`** — add a `## Requirements` bullet ("not bundled in a
   container; you provide your own build/install") and an
   `## Installing Pheniqs` section mirroring GORG-Dark's `## Installing
   kmernorm` (what it is, why it's not bundled, how to get it, then the
   standing disclaimer: *"Installing and running this is entirely at your own
   discretion — confirm license terms and satisfy yourself of the tool's
   suitability independently; this pipeline doesn't warrant or vouch for it.
   See `THIRD_PARTY_LICENSES.md` for more."*). Update the "every process runs
   in its own container except ..." line to include Pheniqs alongside kmernorm.

2. **`THIRD_PARTY_LICENSES.md`** — apply the same "runs it via a user-provided
   binary rather than depending on any third-party container … entirely at the
   discretion of the user" framing already used for kmernorm to the Pheniqs
   rows, and update the "Summary of obligations" accordingly.

3. **`main.nf`** — drop the `quay.io/biocontainers/pheniqs:...` container from
   `PHENIQS_SAMPLE_DEMULTIPLEX` / `PHENIQS_DEMULTIPLEX` (declare `container
   null`, same pattern as `KMERNORM_v1_0_0`); charlie's own `nextflow.config`
   (gitignored) would need a matching `module`/`beforeScript` override to keep
   running locally. Note: changing a process's container re-hashes its tasks
   and busts the `-resume` cache for it and everything downstream — batch this
   when no run is mid-flight, not during one.
