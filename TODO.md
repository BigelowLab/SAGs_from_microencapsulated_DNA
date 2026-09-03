# TODO

Open work on this pipeline. Status tags: **OPEN** / **DEFERRED** (deliberately not
now) / **BLOCKED** (waiting on something) / **WATCH** (verify, not a code change).

---

## Licensing: user-provided tools — kmernorm, Pheniqs, ViralRecall — **OPEN**

Three tools this pipeline runs cannot be safely bundled or redistributed in a
container image:

- **kmernorm** (Mingkun Li) — no license anywhere (SourceForge page or source
  archive), i.e. all-rights-reserved by default copyright. Currently still
  references the third-party image `brwnj/kmernorm:v1.0.0`, which itself has no
  published Dockerfile, build provenance, or license.
- **Pheniqs** (`biosails/pheniqs`) — the repository's top-level `LICENSE` is a
  restrictive NYU research license (internal/non-commercial only; no
  redistribution, modification, or sublicensing without a signed agreement),
  which conflicts with the AGPL-3.0-or-later headers on individual source files.
- **ViralRecall** (`faylward/viralrecall`) — no `LICENSE` file or license
  statement in the repository. (An unrelated MIT-licensed fork exists but is not
  the code this pipeline runs.) Also currently runs via a `beforeScript`-activated
  host conda env, not a container.

**Plan for the final (public, archival) version:** the pipeline requires the user
to install kmernorm, Pheniqs, and ViralRecall themselves, on `PATH` / in the
expected env, entirely at their own discretion — the pipeline does not bundle,
vouch for, or warrant any of them. This is the approach the sibling public repo
**[BigelowLab/GORG-Dark-SAG-assembly](https://github.com/BigelowLab/GORG-Dark-SAG-assembly)**
already takes for kmernorm.

Copy that repo's existing language and extend it to all three tools:

1. **`README.md`**
   - `## Requirements`: add a bullet per tool — "not bundled in a container; you
     provide your own build/install" — pointing at an install section.
   - Add `## Installing kmernorm` / `## Installing Pheniqs` / `## Installing
     ViralRecall` sections (mirror GORG-Dark's `## Installing kmernorm`: what it
     is, why it's not bundled, how to get it, then the standing disclaimer —
     *"Installing and running this is entirely at your own discretion — confirm
     license terms and satisfy yourself of the tool's suitability independently;
     this pipeline doesn't warrant or vouch for it. See `THIRD_PARTY_LICENSES.md`
     for more."*). Keep GORG-Dark's macOS `kmernorm` build warning.
   - Note in the intro / `## Requirements` that every process runs in its own
     container **except** these (like GORG-Dark's "every process except
     `KMERNORM_v1_0_0` executes in its own container").

2. **`THIRD_PARTY_LICENSES.md`** — the kmernorm row already carries GORG-Dark's
   "runs it via a user-provided binary rather than depending on any third-party
   container … entirely at the discretion of the user" wording; apply the same
   framing to the Pheniqs and ViralRecall rows once the pipeline actually stops
   shipping them, and update the "Summary of obligations" accordingly.

3. **`main.nf`** — drop `container 'brwnj/kmernorm:v1.0.0'` from `KMERNORM_v1_0_0`
   (declare `container null`, as GORG-Dark does); give `PHENIQS_*` and
   `VIRALRECALL2` the same treatment (no container / documented user-provided
   env). Note: changing a process's container re-hashes its tasks and busts the
   `-resume` cache for it and everything downstream — batch these together, not
   mid-run.
