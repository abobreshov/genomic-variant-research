---
name: genomic-variant-research
description: Systematic computational analysis of genetic variants — from coordinates to ACMG classification. Chains 36 tools into a reproducible workflow for any gene/variant.
tags: [genomics, variant-classification, acmg, dna, protein, gene-therapy]
triggers:
  - variant analysis
  - variant classification
  - ACMG classification
  - score variant
  - pathogenicity
  - genetic variant research
  - gene research
  - hearing loss genetics
---

# Genomic Variant Research

Systematic workflow for computational analysis of any genetic variant. Input: gene + variant → Output: full evidence package with ACMG classification.

## Workflow Overview

```
Gene + Variant
    │
    ├─ Step 1: Coordinates & Annotation (Ensembl VEP)
    ├─ Step 2: Pathogenicity Scores (AlphaMissense + REVEL + CADD)
    ├─ Step 3: Population Frequency (gnomAD + ClinVar)
    ├─ Step 4: Protein Structure (AlphaFold DB + AF3)
    ├─ Step 5: Splicing (SpliceAI + AlphaGenome)
    ├─ Step 6: Conservation (UniProt orthologs + Clustal Omega)
    ├─ Step 7: Disorder & Domains (IUPred3 + ESMFold)
    ├─ Step 8: Clinical Context (ClinGen + OMIM + ClinicalTrials)
    └─ Output: ACMG Evidence Summary
```

## Step 1: Coordinates & Annotation

```bash
# Ensembl VEP — variant effect prediction
curl -s "https://rest.ensembl.org/vep/human/hgvs/NM_TRANSCRIPT:c.VARIANT?" \
  -H "Content-Type:application/json" | python3 -c "
import sys,json
d = json.load(sys.stdin)
for v in d:
    print(f'Location: chr{v[\"seq_region_name\"]}:{v[\"start\"]}')
    for csq in v.get('transcript_consequences', []):
        if 'missense' in str(csq.get('consequence_terms', [])):
            print(f'Gene: {csq[\"gene_symbol\"]} | AA: {csq.get(\"amino_acids\",\"\")} | Pos: {csq.get(\"protein_start\",\"\")}')"
```

Extract: chromosome, position, strand, amino acid change, protein position, transcript ID.

## Step 2: Pathogenicity Scores

### AlphaMissense (PP3 evidence)
```bash
curl -O "https://alphafold.ebi.ac.uk/files/AF-UNIPROT_ID-F1-aa-substitutions.csv"
grep "E1659A" AF-Q7RTU9-F1-aa-substitutions.csv
```
Thresholds (Pejaver 2022): ≥0.840 → PP3_Moderate, ≥0.949 → PP3_Strong

### REVEL (PP3 evidence)
Via dbNSFP or Ensembl VEP.
Thresholds (Pejaver 2022): ≥0.773 → PP3_Moderate, ≥0.932 → PP3_Strong

### CADD
```bash
curl -s "https://rest.ensembl.org/vep/human/region/CHR:POS-POS:STRAND/ALT?CADD=1" \
  -H "Content-Type:application/json"
```
PHRED ≥20 = top 1%, ≥25 = strong, ≥30 = very strong.

**Don't double-count PP3** — use strongest single predictor.

## Step 3: Population Frequency (PM2 / BA1 / BS1)

### gnomAD
```bash
curl -s -X POST 'https://gnomad.broadinstitute.org/api' \
  -H 'Content-Type:application/json' \
  -d '{"query":"{variant(variantId:\"CHR-POS-REF-ALT\",dataset:gnomad_r4){genome{ac af an}}}"}'
```
Absent → PM2_Supporting. AF >5% → BA1.

### gnomAD Gene Constraint (important for recessive genes)
```bash
curl -s -X POST 'https://gnomad.broadinstitute.org/api' \
  -H 'Content-Type:application/json' \
  -d '{"query":"{gene(gene_symbol:\"GENE\",reference_genome:GRCh38){gnomad_constraint{pli oe_lof oe_lof_upper oe_mis oe_mis_upper lof_z mis_z}}}"}'
```
For recessive genes: pLI will be LOW (expected — one broken copy is tolerated). Check **missense Z-score** — high value means missense variants are under selection (functionally important).

### ClinVar
```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=GENE[gene]&retmode=json"
```

## Step 4: Protein Structure

### AlphaFold DB
```bash
curl -s "https://alphafold.ebi.ac.uk/api/prediction/UNIPROT_ID"
```
pLDDT >90 at variant position = structured region = missense more likely damaging.

### AlphaFold 3 Server
Submit WT and mutant at https://alphafoldserver.com — compare pTM/ipTM.

## Step 5: Splicing & Regulatory

### AlphaGenome
```python
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
dna_model = dna_client.create(os.environ['ALPHAGENOME_API_KEY'])
variant = genome.Variant(chromosome='chrN', position=POS,
    reference_bases='REF', alternate_bases='ALT')
interval = variant.reference_interval.resize(
    dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_100KB'])
scores = dna_model.score_variant(interval=interval, variant=variant,
    variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()))
df = variant_scorers.tidy_scores(scores)
```

### SpliceAI
Web: https://spliceailookup.broadinstitute.org/ — enter CHR-POS-REF-ALT

## Step 6: Conservation

```bash
# UniProt orthologs
curl -s "https://rest.uniprot.org/uniprotkb/search?query=gene:GENE+AND+reviewed:true&format=fasta&size=20" > orthologs.fasta
# Align
clustalo -i orthologs.fasta -o aligned.fasta --outfmt=fasta
```
100% conservation across >8 mammals = strong evidence.

### Extended Conservation via OrthoDB
```bash
# Get ortholog group
curl -s 'https://data.orthodb.org/v11/search?query=GENE&species=9606'
# Download FASTA (WARNING: these are UNALIGNED)
curl -s 'https://data.orthodb.org/v11/fasta?id=GROUP_ID' > orthologs_raw.fasta
```
**CRITICAL**: OrthoDB FASTA are raw unaligned sequences. You MUST align before checking positional conservation. Use BioPython pairwise alignment of a window around the variant position against each species to verify the correct aligned column. Do NOT just check position N in each raw sequence — indels shift positions and you'll get garbage results.

## Step 7: Disorder & Domains

### IUPred3
```bash
curl -s "https://iupred3.elte.hu/plot_json/UNIPROT_ID" > disorder.json
```
Variant in ordered region → missense more likely damaging.

## Step 8: Clinical Context

- **ClinGen**: https://search.clinicalgenome.org/kb/genes/HGNC:ID
- **OMIM**: Check MIM number, phenotype, inheritance
- **ClinicalTrials.gov**: `curl "https://clinicaltrials.gov/api/v2/studies?query.term=GENE+gene+therapy&format=json"`

## ACMG Evidence Assembly

| Criteria | Evidence | Source |
|----------|----------|--------|
| PVS1 | Null variant in LoF gene | VEP consequence |
| PM1 | Functional domain / hot spot | AlphaFold pLDDT + domains |
| PM2 | Absent from population DBs | gnomAD + ClinVar |
| PM3 | In trans with known pathogenic | Family segregation |
| PP3 | Computational evidence | AlphaMissense/REVEL/CADD |
| PP4 | Specific phenotype match | Clinical data |

**Classification (Richards 2015):**
- Pathogenic: ≥1 VeryStrong + ≥1 Strong, OR ≥2 Strong
- Likely Pathogenic: 1 Strong + 1-2 Moderate, OR ≥3 Moderate
- VUS: doesn't meet LP or LB

## Critical Pitfalls

1. **OrthoDB sequences are UNALIGNED** — NEVER check positional conservation by taking position N from raw FASTA. Sequences have indels so position 1659 in species A ≠ position 1659 in species B. MUST pairwise-align a window (±30 aa) around the target position for each species, then check the aligned column. Naive position-based check gave "3.8% conservation" for a 100%-conserved position. Use BioPython pairwise2 if clustalo/muscle not installed.
2. **PP1 ≠ conservation** — PP1 = co-segregation in family, NOT cross-species conservation
3. **Don't double-count PP3** — AlphaMissense + REVEL measure overlapping signal
4. **Pseudogene contamination** — SIFT/PolyPhen fail for genes with pseudogenes. Use AlphaMissense/REVEL
5. **Minus strand** — verify strand before cDNA→genomic conversion
6. **REVEL thresholds** — use Pejaver 2022 ClinGen-calibrated only
7. **gnomAD coverage** — low coverage ≠ absent variant
8. **AlphaGenome** — no inner ear tissue in training data
9. **CADD direct API returns "Not found"** — use Ensembl VEP with `?CADD=1` parameter instead
10. **UniProt ID verification** — always verify gene→UniProt mapping via search API, never from memory
11. **AlphaFold confidence JSON** — v4 may 404, try v6
12. **gnomAD constraint for recessive genes** — low pLI is EXPECTED. High missense Z-score is the signal.
13. **AlphaMissense heatmap** — download full CSV to find regional hotspots via sliding window.
14. **SpliceAI API times out** — web interface works, CLI pip package has broken deps (pkg_resources). Use Ensembl VEP with SpliceAI plugin or web lookup.
15. **OMIM/VarSome APIs need approval** — registration doesn't give immediate access. Use web interfaces + get scores via VEP/gnomAD/AlphaFold direct APIs instead.

## Reference: STRC E1659A (Worked Example)

| Tool | Result | ACMG |
|------|--------|------|
| VEP | missense_variant, E/A, pos 1659 | — |
| AlphaMissense | 0.9016 (≥0.840) | PP3_Moderate |
| REVEL | 0.789 (≥0.773) | PP3_Moderate |
| CADD | PHRED 27.5 (top 0.18%) | Supports PP3 |
| gnomAD | Absent (0/251K) | PM2_Supporting |
| ClinVar | Absent | PM2 support |
| AlphaFold | pLDDT 95.69 at E1659 | Structured region |
| AF3 | WT pTM 0.63, E1659A pTM 0.64 | Structure intact |
| gnomAD constraint | pLI=0.14, missZ=7.62 | Recessive + constrained |
| AM heatmap | E1659 region = 69% pathogenic (highest) | Hotspot domain |
| Conservation | 100% across 40 mammals (OrthoDB aligned) | Strong conservation |
| IUPred3 | 0.05 at 1659 (ordered) | Structured |
| PM3 | In trans with 98kb deletion (P) | PM3_Moderate |
| **Final** | **PM3(M) + PP3(M) + PM2(S) = VUS** | **VUS-high** |
