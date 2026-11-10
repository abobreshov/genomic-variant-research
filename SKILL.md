---
name: genomic-variant-research
description: Runs a reproducible computational workflow for a human genetic variant, from HGVS or genomic coordinates to an ACMG/AMP evidence table and classification, using free public APIs (Ensembl VEP, REVEL, AlphaMissense, CADD, SpliceAI, gnomAD, ClinVar, AlphaFold, OrthoDB, ClinGen). Use when asked to analyse, score or classify a variant, assess pathogenicity, assign ACMG criteria such as PP3 or PM2, interpret a VUS, or research the gene behind a genetic condition such as hearing loss.
license: MIT
metadata:
  tags: "genomics, variant-classification, acmg, dna, protein, gene-therapy"
---

# Genomic Variant Research

Input: a gene and variant, as HGVS on a named transcript or as GRCh38 `chr-pos-ref-alt`. Output: an ACMG evidence table and a classification, with every value traced to a call made in this run.

This is research support, not a clinical report. Say so in every output: a clinical geneticist must review it before it informs care, because computational evidence alone cannot carry a diagnosis.

## Rules for the whole run

- **Take identifiers from API responses in this run, never from memory**: coordinates, transcript, strand, UniProt ID (a remembered Q7RTU1 for Q7RTU9 once slipped through).
- **Check the numbering before any protein lookup**: the UniProt residue at VEP's `protein_start` must equal VEP's reference amino acid. VEP numbers the MANE protein; AlphaMissense, AlphaFold and IUPred number the UniProt canonical one, and an isoform mismatch silently shifts every lookup.
- **Log each value with the call that produced it and the data release**, as you go. Scores change between releases, and the evidence table is built from this log.
- **Count each evidence code once.** Tools that share inputs measure the same signal, so stacking them inflates the evidence.

## Prerequisites

Steps 1-4, 7 and 8 need only `curl` and `python3` with network access. In addition:

- Step 6 alignment: Biopython (`pip install biopython`), or Clustal Omega (`clustalo`, packaged as `clustal-omega`).
- Step 5 AlphaGenome, optional: `pip install alphagenome python-dotenv` and an `ALPHAGENOME_API_KEY` from Google DeepMind, set in the environment or in `~/.env` (the file DeepMind's own skills read). Never print the key.

## Checklist

Copy this and tick it off:

```text
- [ ] 1 Coordinates: vcf_string, strand, MANE transcript, protein position from VEP; residue matches UniProt
- [ ] 2 PP3/BP4: one predictor's level, tool named
- [ ] 3 Frequency: gnomAD AF, or absence plus the coverage and paralog check; ClinVar entries
- [ ] 4 Structure: pLDDT at the position
- [ ] 5 Splicing: SpliceAI delta scores
- [ ] 6 Conservation: aligned column across species
- [ ] 7 Disorder: IUPred score at the position
- [ ] 8 Clinical: ClinGen gene validity, OMIM inheritance, trials
- [ ] 9 Evidence table filled from the template
- [ ] 10 Fresh-context check passed
```

A step is done when its values and their sources are in the log. If step 10 finds an error, return to Step 1-8 that produced the value, fix it, and run step 10 again.

## Step 1: Coordinates and annotation

One VEP call also returns the scores for steps 2 and 5. Use the MANE Select row.

```bash
curl -s "https://rest.ensembl.org/vep/human/hgvs/TRANSCRIPT:c.VARIANT?mane=1&vcf_string=1&AlphaMissense=1&REVEL=1&CADD=1&SpliceAI=1" \
  -H "Content-Type:application/json" | python3 -c "
import sys, json
for v in json.load(sys.stdin):
    print(f'chr{v[\"seq_region_name\"]}:{v[\"start\"]} strand {v.get(\"strand\")} transcript alleles {v.get(\"allele_string\")} genomic {v.get(\"vcf_string\")}')
    for c in v.get('transcript_consequences', []):
        if not c.get('mane_select'):
            continue
        print(c['transcript_id'], c.get('gene_symbol'), c['consequence_terms'], c.get('amino_acids'), c.get('protein_start'))
        print('  REVEL', c.get('revel'), 'AM', c.get('alphamissense', {}).get('am_pathogenicity'),
              'CADD', c.get('cadd_phred'), 'SpliceAI', c.get('spliceai'))"
```

For GRCh38 `chr-pos-ref-alt` input, use `https://rest.ensembl.org/vep/human/region/CHR:POS-END:1/ALT` with the same parameters, where END = POS + length(REF) - 1: the span must cover the whole reference allele, or a deletion errors and a longer REF is scored as a different one-base change (`17-43124027-ACT-A` is `17:43124027-43124029:1/A`; an insertion such as `17-43124027-A-ACT` keeps `POS-POS`). Take the gnomAD id, the AlphaGenome alleles and the ClinVar window position from the `vcf_string` it returns. CADD's own API often answers "Not found", so take CADD from VEP.

**Use `vcf_string` as the variant's genomic id** for gnomAD and AlphaGenome: it is GRCh38, forward strand, VCF-style with the anchor base for indels (`15-43600551-T-G`, `17-43124027-ACT-A`). For HGVS input, `allele_string` is on the transcript strand: on a minus-strand gene (`strand` -1) it is the complement of the genomic change, and a gnomAD lookup with it finds nothing and reads as absence. HBB c.20A>T is `allele_string` A/T but `vcf_string` 11-5227002-T-A (rs334, common in African-ancestry populations).

PVS1 (a null variant: nonsense, frameshift, canonical ±1,2 splice, initiation codon, whole-exon deletion) applies only where loss of function is the gene's known disease mechanism. Set its strength with the ClinGen PVS1 decision tree (Abou Tayoun 2018): a truncation that escapes nonsense-mediated decay, in the last exon or the last 50 nt of the penultimate one, is weaker.

## Step 2: Computational evidence (PP3 / BP4)

Use one predictor, chosen before you look at the scores: **REVEL**, from the VEP call. Use AlphaMissense instead only when REVEL has no score, and say so. Report the others, and cross-species conservation, as context, not as extra evidence: they overlap with REVEL, and conservation is one of its inputs.

| Tool | PP3 Supporting | PP3 Moderate | PP3 Strong | BP4 Supporting | Source |
|---|---|---|---|---|---|
| REVEL | ≥0.644 | ≥0.773 | ≥0.932 | ≤0.290 | Pejaver 2022 |
| AlphaMissense | ≥0.792 | ≥0.906 | ≥0.990 | ≤0.169 | Bergquist 2025 |
| CADD (PHRED) | ≥25.3 | ≥28.1 | none | ≤22.7 | Pejaver 2022 |

SIFT and PolyPhen-2 reach weaker evidence levels than REVEL in the same calibration, and their alignments can pull in paralog sequences, which makes conserved positions look variable. Leave them out.

The full AlphaMissense table for the protein (`https://alphafold.ebi.ac.uk/files/AF-UNIPROT_ID-F1-aa-substitutions.csv`; if it returns 403, the hg38 file on Zenodo record 8208688) shows regional hotspots in a sliding window. That is context for PM1, which needs a known hotspot or functional domain without benign variation. PP3 and PM1 together count at most as Strong (Pejaver 2022).

## Step 3: Population frequency (BA1 / BS1 / PM2)

```bash
curl -s -X POST 'https://gnomad.broadinstitute.org/api' -H 'Content-Type:application/json' \
  -d '{"query":"{variant(variantId:\"VCF_STRING\",dataset:gnomad_r4){genome{ac an af filters faf95{popmax popmax_population}} exome{ac an af filters faf95{popmax popmax_population}} flags}}"}'
```

`VCF_STRING` is the `vcf_string` from step 1, never the transcript-strand alleles.

- BA1: frequency above 5% in any continental population with at least 2,000 alleles, unless the gene has a ClinGen BA1 exception (Ghosh 2018). Use `faf95.popmax`, not the overall `af`, which can hide the population that crosses 5%, and only from a callset whose `filters` list is empty: a filtered callset's frequency is not trusted. If every callset is filtered, record the filters and apply none of BA1, BS1 or PM2. For BS1, use the gene's ClinGen expert panel threshold when one exists.
- Absent gives PM2_Supporting (ClinGen downgraded PM2 to Supporting), but only once you have checked that gnomAD could have seen the variant:
  - coverage at the position is adequate (gnomAD browser coverage track). Low coverage is not absence.
  - the gene has no near-identical paralog or pseudogene: check `https://rest.ensembl.org/homology/symbol/human/GENE?type=paralogues;format=condensed`, the UCSC segmental-duplication track, and OMIM or the literature, since pseudogenes are often missing from the Ensembl list. Short reads from near-identical copies mismap, so absence there is weak evidence. STRC and its pseudogene STRCP1 are the standing example.

Gene constraint:

```bash
curl -s -X POST 'https://gnomad.broadinstitute.org/api' -H 'Content-Type:application/json' \
  -d '{"query":"{gene(gene_symbol:\"GENE\",reference_genome:GRCh38){gnomad_constraint{pli oe_lof oe_lof_upper oe_mis oe_mis_upper lof_z mis_z}}}"}'
```

For a recessive gene a low pLI is expected, because one broken copy is tolerated; a high missense Z says missense changes are under selection. The paralog caveat above applies to constraint too.

ClinVar: `curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=GENE[gene]&retmode=json"`, then look for the variant's HGVS in the results.

## Step 4: Protein structure

`curl -s "https://alphafold.ebi.ac.uk/api/prediction/UNIPROT_ID"` returns the file URLs; read per-residue confidence from `plddtDocUrl` instead of building the file name, because the model version changes. pLDDT above 90 at the position means a structured region, where a missense change is more likely to matter.

AlphaFold 3 (alphafoldserver.com) can fold the wild type and the mutant, but it is not built to predict the effect of a point change: near-identical pTM is expected and does not show the structure is unaffected.

## Step 5: Splicing

SpliceAI comes from the VEP call in step 1. Outside the canonical ±1,2 splice positions, a delta score of 0.2 or more supports PP3 for a splicing effect and 0.1 or less supports BP4 (ClinGen SVI splicing guidance, Walker 2023). Do not apply PP3 for both a missense and a splicing effect of the same variant. If VEP has no SpliceAI score, use https://spliceailookup.broadinstitute.org/ (its pip package fails to install on current Python).

AlphaGenome, optional; its training data has no inner-ear tissue:

```python
import os
from dotenv import load_dotenv
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

load_dotenv(os.path.expanduser('~/.env'))  # an exported variable still wins
model = dna_client.create(os.environ['ALPHAGENOME_API_KEY'])
variant = genome.Variant(chromosome='chrN', position=POS, reference_bases='REF', alternate_bases='ALT')  # from vcf_string (step 1)
interval = variant.reference_interval.resize(dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_100KB'])
scores = model.score_variant(interval=interval, variant=variant,
                             variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()))
df = variant_scorers.tidy_scores(scores)
```

## Step 6: Conservation

Fetch orthologs from UniProt (`https://rest.uniprot.org/uniprotkb/search?query=gene:GENE+AND+reviewed:true&format=fasta&size=20`) or OrthoDB (`https://data.orthodb.org/current/search?query=GENE&species=9606`, then `https://data.orthodb.org/current/fasta?id=GROUP_ID`).

These sequences are **unaligned**. Indels shift positions, so residue N in one species is not residue N in another; reading position N straight from the raw FASTA once reported 3.8% conservation for a position conserved in all 40 mammals. Align a window around the variant to each sequence and read the aligned column:

```python
# needs Biopython: pip install biopython
import sys
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

# usage: python3 conservation.py human.fasta orthologs.fasta POSITION
human = str(next(SeqIO.parse(sys.argv[1], "fasta")).seq)
pos = int(sys.argv[3])                      # 1-based protein position from step 1
start = max(0, pos - 31)
window, q = human[start:pos + 30], pos - 1 - start   # +-30 aa gives a unique local match

aligner = Align.PairwiseAligner(mode="local", open_gap_score=-10, extend_gap_score=-0.5)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
for rec in SeqIO.parse(sys.argv[2], "fasta"):
    aln = aligner.align(window, str(rec.seq))[0]
    hit = next((rec.seq[ts + q - qs] for (qs, qe), (ts, te) in zip(*aln.aligned) if qs <= q < qe), "-")
    print(rec.id, hit)
```

With Clustal Omega: `clustalo -i orthologs.fasta -o aligned.fasta`, then read the column of the human residue.

## Step 7: Disorder

`curl -s "https://iupred3.elte.hu/iupred3/UNIPROT_ID.json"`: the score for residue N is `iupred2[N-1]`. Above 0.5 is disordered; a variant in an ordered region is more likely to matter.

## Step 8: Clinical context

- ClinGen gene-disease validity: https://search.clinicalgenome.org/kb/genes/HGNC:ID
- OMIM: MIM number, phenotype, inheritance. Its API needs approval, so use the web pages.
- Trials: `curl "https://clinicaltrials.gov/api/v2/studies?query.term=GENE+gene+therapy&format=json"`

## Evidence table

Use this exact structure:

```markdown
| Code | Level | Value | Source (call, release) |
|---|---|---|---|
| PP3 | Moderate | REVEL 0.81 | VEP REST, release N |

Context, not counted: <other scores, conservation, structure>
Classification: <class> (rule: <combining rule from the table below>)
Research output, not clinically reviewed.
```

Use family data (PM3 in trans with a pathogenic variant, PP1 co-segregation) only when the user provides it. PP1 means co-segregation with disease in a family; cross-species conservation never counts as PP1.

## Combining rules

Richards et al. 2015, with ClinGen SVI's PM2 change applied.

| Class | Any one of |
|---|---|
| Pathogenic | 1 Very strong + (≥1 Strong, or ≥2 Moderate, or 1 Moderate + 1 Supporting, or ≥2 Supporting); ≥2 Strong; 1 Strong + (≥3 Moderate, or 2 Moderate + ≥2 Supporting, or 1 Moderate + ≥4 Supporting) |
| Likely pathogenic | 1 Very strong + 1 Moderate; 1 Very strong + 1 Supporting (ClinGen SVI, with PM2_Supporting); 1 Strong + 1-2 Moderate; 1 Strong + ≥2 Supporting; ≥3 Moderate; 2 Moderate + ≥2 Supporting; 1 Moderate + ≥4 Supporting |
| Benign | BA1; ≥2 Strong benign |
| Likely benign | 1 Strong benign + 1 Supporting benign; ≥2 Supporting benign |
| Uncertain (VUS) | none of the above, or pathogenic and benign evidence that conflict |

## Fresh-context check

Give the evidence table and the call log to a subagent that did not gather the evidence. It checks that:

- every value matches its logged source;
- each code is counted once: one PP3 tool, conservation not added to it, missense and splicing PP3 not both applied;
- PM2 is not applied where coverage or a paralog makes absence unreliable;
- the classification follows the combining rules.

## Worked example: STRC c.4976A>C (p.Glu1659Ala)

| Code | Level | Value | Source |
|---|---|---|---|
| PP3 | Supporting | REVEL 0.65 (≥0.644) | VEP REST, Ensembl 116, ENST00000450892 |
| PM3 | Moderate | in trans with a pathogenic 98 kb deletion | family data |

Context, not counted: absent from gnomAD (15-43600551-T-G not found), but PM2 is not applied because reads from the near-identical STRCP1 pseudogene mismap, so absence is unreliable here; AlphaMissense 0.9016 (Supporting range); CADD 25.5 (Supporting range); pLDDT 95.69 at E1659; AlphaFold 3 pTM 0.63 wild type, 0.64 mutant; gnomAD pLI 0.14 and missense Z 7.62 (recessive and constrained, with the pseudogene caveat); E1659 sits in the AlphaMissense region with the most pathogenic-range scores (69%); conserved in 40 of 40 mammals (aligned OrthoDB); IUPred ordered; absent from ClinVar.

Classification: VUS (1 Moderate + 1 Supporting meets no likely-pathogenic rule). The original analysis recorded REVEL 0.789 (PP3_Moderate) and CADD 27.5 from an unrecorded source, and counted PM2_Supporting; with all of those (2 Moderate + 1 Supporting) the class is still VUS. Scores differ between sources and releases, which is why every score is logged with its call and release.
