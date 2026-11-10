# Third-party notices

The skill itself is MIT-licensed (`LICENSE`). The files below are third-party code under their own licence.

## DeepMind science-skills

Source: https://github.com/google-deepmind/science-skills, commit `370f8f33e5bb27dffe8289ea9bdec436686613ce`.
Copyright 2026 Google LLC. Licensed under the Apache License, Version 2.0; full text in `scripts/LICENSE-APACHE-2.0`.

| File in this repo | Upstream file | Changes |
|---|---|---|
| `scripts/gnomad_variant_frequency.py` | `skills/gnomad_database/scripts/get_variant_frequency.py` | renamed; pinned `polite-http==0.1.3`; `--output` may be a bare filename; requests `filters` on every callset and `populations` on joint |
| `scripts/gnomad_gene_constraint.py` | `skills/gnomad_database/scripts/get_gene_constraint.py` | renamed; pinned `polite-http==0.1.3`; `--output` may be a bare filename; also returns `lof_z` and `mis_z` |
| `scripts/clinvar_api.py` | `skills/clinvar_database/scripts/clinvar_api.py` | pinned `polite-http==0.1.3` and `python-dotenv==1.2.4`; error messages name the endpoint, never the URL with the API key and email; reads only the NCBI settings from `~/.env` |
| `scripts/ucsc_conservation.py` | `skills/ucsc_conservation_and_tfbs/scripts/get_conservation.py` | renamed; pinned `polite-http==0.1.3` |
| `scripts/alphagenome_avi.py` | `skills/alphagenome_variant_impact_score/scripts/alphagenome_atlas_avi.py` | renamed; dependencies moved from the upstream `pyproject.toml` into an inline script block, pinned |

Each file keeps Google's licence header and lists its changes below it.

## Runtime dependencies

`uv` installs these from PyPI on first use; none is bundled. `polite-http` (Apache-2.0) is DeepMind's rate-limited HTTP client published as a package. The AVI script also uses the `alphagenome` client (Apache-2.0) and common data libraries, pinned in its inline block.

## Data terms

The scripts query public services, each with its own terms of use, which you are responsible for following:

- gnomAD: https://gnomad.broadinstitute.org/policies
- NCBI ClinVar and E-utilities: https://www.ncbi.nlm.nih.gov/home/about/policies/
- UCSC Genome Browser API: https://genome.ucsc.edu/conditions.html
- AlphaGenome API (key required): https://deepmind.google.com/science/alphagenome/
