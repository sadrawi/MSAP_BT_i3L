# Sequence-only in silico MSAP primer screening — proof of concept

This Python 3 workflow simulates the **potential** EcoRI–CCGG (HpaII/MspI) MSAP fragments from an unmethylated reference. It does **not** predict methylation or guarantee observed gel bands.

## Reference rice genome
Use *Oryza sativa* Japonica Group, cultivar Nipponbare, IRGSP-1.0, NCBI **GCF_001433935.1**. Record release/accession and checksum. Download an **unmasked nuclear genome FASTA** from [NCBI assembly](https://www.ncbi.nlm.nih.gov/assembly/578461) or [Ensembl Plants](https://plants.ensembl.org/Oryza_sativa/Info/Index). Genome file should contain full chromosomal sequences (plain `.fasta` or `.fasta.gz`), rather than CDS/transcripts, a protein FASTA, or a repeat-masked assembly. Exclude organelle sequences for a nuclear-only analysis.

## Run with included reproducible demo
```bash
python make_demo.py
python test_msap.py
python msap_sim.py --fasta synthetic_demo.fasta --e-bases 2 --hm-bases 2 --min-bp 100 --max-bp 500 --bin-bp 5 --out demo_results
```

## Run a real rice genome
```bash
python msap_sim.py --fasta rice_IRGSP-1.0.fa.gz --e-bases 3 --hm-bases 3 --min-bp 100 --max-bp 500 --bin-bp 5 --out rice_results
```
This runs **4096 pairs** from three selective nucleotides on EcoRI and HM ends. For a quicker initial genome test use 2/2 (256 pairs). For large genomes the input FASTA is read one contig at a time; however output size and temporary memory increase with the number of qualifying fragments.

## Output
- `predicted_fragments.csv`: genomic cut positions, length, selective extensions, orientation.
- `primer_ranking.csv`: fragment count, approximately distinguishable size bins, collisions, contig spread per primer pair. Ranking penalizes crowded bins; **not a calibrated quality score**.

## Assumptions/limitations
- Complete digestion of unmethylated DNA by EcoRI and HpaII/MspI; both are treated as the same CCGG cutter **only in this sequence-only model**.
- Only adjacent heterotypic cuts are possible valid fragments. Other cuts disrupt them.
- No methylation, enzyme kinetics, DNA polymorphism, partial digestion, heterogeneity, amplification efficiency, adapter/ligase biases, assay noise or tissue specificity.
- Fragments are sized from actual enzyme cleavage positions without adapters; electrophoresis may display adapter/primer-dependent size offsets. Validate against laboratory standards.
- Size-bin occupancy is only a crude co-migration approximation; fragment sequences of the same size may appear as one band.
- **Selective bases** are matched immediately internal to the recognition site and are reverse complemented when necessary. This assumes a standard EcoRI/HM adapter-primer structure; check laboratory primer and adapter designs before ordering.
- This is a *teaching/research prototype*, not a validated assay-design pipeline. A high-scoring combination may be unsuitable experimentally.

## Streamlit app
```bash
pip install -r requirements.txt
streamlit run app.py
```
- **Genome**: the synthetic demo, an upload (plain/`.gz` FASTA or an NCBI Datasets `.zip`; up to 1 GB via `.streamlit/config.toml`), or a FASTA path on your own machine, read in place. Organelle-like contigs (Mt/Pt, mitochondrion, chloroplast/plastid) start excluded.
- The FASTA is scanned once (every fragment up to 1,500 bp); selective bases, size window and bin width then update instantly. A synthetic 120 Mb genome scanned in ~4 s, so a 370 Mb rice assembly takes roughly 10–20 s.
- **Tabs**: ranking, E × HM heatmap, virtual gel (ladders, co-migrating bands brighter), per-pair electropherogram with genomic coordinates, and a contig summary.
- **Downloads**: `predicted_fragments.csv` and `primer_ranking.csv` are byte-identical to the CLI's for the same inputs; `run_info.json` records the input SHA-256, parameters and the equivalent CLI command. If contigs are excluded, the CLI (which reads every contig) will differ.
- Works on old local Streamlit installs (tested 1.16 and 1.26) and current releases (1.66).
- `python test_app.py` checks app ↔ CLI parity on the demo and randomized genomes (plain, `.gz`, `.zip`) and drives the UI headlessly (UI part needs Streamlit ≥ 1.28).

### Deploying on Streamlit Community Cloud
Push these files to the **root** of a GitHub repo (keep `.streamlit/config.toml`), then deploy with main file `app.py`; `requirements.txt` is picked up automatically. Upload the genome **gzipped** (`*.fa.gz`): apps get roughly 0.7–2.7 GB of memory, and a plain rice FASTA alone is ~380 MB. The "Local path" source is hidden when Streamlit runs headless, as on Cloud, because it would read the server's disk; set `MSAP_LOCAL_PATH=1` (or `0`) to override.
