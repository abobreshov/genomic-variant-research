# Genomic Variant Research Skill

This repository turns a real variant-classification workflow into a reusable, auditable skill for genetic variant analysis. A specialist should review every output.

An AI agent skill for systematic computational analysis of human genetic variants. Feed it a gene + variant → get an ACMG evidence table and classification from free public bioinformatics APIs.

Built for [Claude Code](https://docs.anthropic.com/en/docs/claude-code) and compatible AI agents.

## What It Does

```
Input: Gene + Variant (e.g., STRC c.4976A>C)
    │
    ├─ Step 1: Coordinates & annotation (Ensembl VEP)
    ├─ Step 2: Computational evidence (REVEL, with AlphaMissense and CADD as context)
    ├─ Step 3: Population frequency (gnomAD + ClinVar)
    ├─ Step 4: Protein structure (AlphaFold DB + AF3)
    ├─ Step 5: Splicing & regulatory (SpliceAI + AlphaGenome)
    ├─ Step 6: Conservation (UniProt orthologs + Clustal Omega + OrthoDB)
    ├─ Step 7: Disorder (IUPred3)
    ├─ Step 8: Clinical context (ClinGen + OMIM + ClinicalTrials.gov)
    ├─ Fresh-context check of the evidence table
    │
Output: ACMG evidence table with classification
```

## Tools Covered

The workflow calls a subset of these; the rest are listed for follow-up research.

| Category | Tools |
|----------|-------|
| **Structural Biology** | AlphaFold DB, AlphaFold 3, AlphaGenome, SWISS-MODEL, ESMFold, IUPred3, Mol*, PyMOL |
| **Variant Scoring** | AlphaMissense, REVEL, CADD, SpliceAI, dbNSFP, VarSome, InterVar, Franklin |
| **Population Genetics** | gnomAD, ClinVar, seqr, DECIPHER, ClinGen |
| **Conservation** | UniProt, UCSC Genome Browser, Ensembl REST API, OrthoDB, Clustal Omega |
| **Disease-Specific** | OMIM, Deafness Variation Database, Hereditary Hearing Loss Homepage |
| **Gene Therapy** | ClinicalTrials.gov, AAV Capsid Database, ASGCT |
| **Protein Analysis** | SignalP, NetGPI, NetNGlyc |
| **Literature** | PubMed/PMC |

All tools are free. No wet lab required. Total compute cost: ~$50-100 in LLM API calls.

## Requirements

`curl`, `python3` and [`uv`](https://docs.astral.sh/uv/); details under *Prerequisites* in `SKILL.md`. The gnomAD, ClinVar, UCSC conservation and AlphaGenome AVI steps run scripts in `scripts/`, vendored from [DeepMind's science-skills](https://github.com/google-deepmind/science-skills). AVI needs an `ALPHAGENOME_API_KEY`.

## Install

Install the whole folder, since `SKILL.md` calls the scripts beside it:

```bash
# For Claude Code / Hermes
git clone https://github.com/lyfar/genomic-variant-research.git ~/.claude/skills/genomic-variant-research
```

Then ask your agent: *"Score variant BRCA1 c.5266dupC using the genomic variant research workflow"*

## Testing

Run these prompts in a fresh session with the skill on and off, and compare the evidence tables:

1. "Classify STRC c.4976A>C (p.Glu1659Ala). It is in trans with a pathogenic 98 kb deletion." (missense, family data, pseudogene caveat)
2. "Score BRCA1 c.5266dupC." (frameshift: PVS1, well-known pathogenic)
3. "Is GJB2 c.-23+1G>A pathogenic?" (splice-site, SpliceAI path)

## Hard-Won Pitfalls

These cost us weeks to discover. They're all in the skill, but the critical ones:

1. **OrthoDB sequences are UNALIGNED** — never check positional conservation by taking position N from raw FASTA. Sequences have indels. Must pairwise-align a window around the target. We got "3.8% conservation" for a 100%-conserved position before catching this.

2. **PP1 ≠ conservation** — PP1 (co-segregation) requires affected family members, not cross-species conservation. Conservation is context for PP3, not a separate criterion. Mixing these up changes classification.

3. **SIFT/PolyPhen fail for genes with pseudogenes** — pseudogene sequences contaminate the alignment, making conserved positions appear variable. Use REVEL instead (AlphaMissense only when REVEL has no score).

4. **Don't double-count PP3** — AlphaMissense, REVEL and CADD measure overlapping signal, and conservation is already one of their inputs. Pick one calibrated predictor (REVEL) before looking at the scores and use only its level; picking the strongest after the fact inflates the evidence.

5. **CADD direct API returns "Not found"** — use Ensembl VEP with `?CADD=1` instead.

6. **gnomAD constraint: low pLI is EXPECTED for recessive genes** — check missense Z-score instead.

7. **Pseudogenes break gnomAD too** — reads from a near-identical paralog (STRC and STRCP1) mismap, so absence from gnomAD (PM2) and constraint scores are unreliable there. Check coverage and paralogs before applying PM2.

## Background

This skill was developed during a real variant classification project for STRC-related hearing loss (DFNB16). The worked example in the skill (STRC E1659A → VUS) is a real case that led to contact with researchers at Harvard Medical School and Fudan University.

The full research is published at [strc.egor.lol](https://strc.egor.lol).

## License

MIT, except the vendored scripts in `scripts/`, which stay Apache-2.0 (`THIRD_PARTY_NOTICES.md`). Use it, fork it, improve it. If you find a new pitfall, please open an issue or PR.
