# phyloseq-taxid-resolver
Batch-resolve NCBI Taxonomy IDs (TaxIDs) from a phyloseq-style taxonomy CSV using Biopython Entrez, with a live tqdm progress bar, caching, resume, safe partial writes, and optional species-rank verification.

> Script entrypoint: **`add_taxids.py`**.  
> Input: CSV with columns `Taxon, Domain, Phylum, Class, Order, Family, Genus, Species`.  
> Output: same CSV **+ `NCBI_TaxID`** appended.

---

## Why this exists
When working with **microbiome** data (e.g. `read_phyloseq()` in R), you often need canonical **NCBI TaxIDs** for downstream joins (KEGG/RefSeq/GTDB), deduping, and QC. This tool converts “Genus species” at scale—politely and reproducibly.

---

## Features
-  **Fast & robust**: JSON `esearch` + XML `efetch` (optional species-rank check).
-  **Resume safe**: JSON cache persists lookups; `--resume` picks up where you left off.
-  **Partial writes**: Streams to `*.part`, then atomically renames to final CSV.
-  **Progress**: `tqdm` live progress bar and `--verbose` per-name logs.
-  **Synonym smart**: Built-in fixes (e.g., `Propionibacterium acnes` → `Cutibacterium acnes`; `Clostridium difficile` → `Clostridioides difficile`).
-  **Reproducible**: Rate-limit aware and deterministic per input row.

---

## Install

```bash
# Python 3.8+
python -m pip install -r requirements.txt
```

## Quickstart

# macOS/Linux
python add_taxids.py taxonomy_phyloseq.csv taxonomy_with_taxid.csv \
  --email "you@example.com" \
  --api-key "YOUR_NCBI_API_KEY" \
  --resume --verbose

# Windows PowerShell
python .\add_taxids.py .\taxonomy_phyloseq.csv .\taxonomy_with_taxid.csv `
  --email 'you@example.com' `
  --api-key 'YOUR_NCBI_API_KEY' `
  --resume --verbose

Tip: first run a tiny batch to validate everything:

python add_taxids.py taxonomy_phyloseq.csv taxonomy_with_taxid.csv --email "you@example.com" --api-key "KEY" --limit 20 --resume --verbose

## Input format

A CSV with these headers (order doesn’t matter, but Species must exist):
Taxon, Domain, Phylum, Class, Order, Family, Genus, Species

Species should be Genus + specific epithet (e.g., Abiotrophia defectiva), or Genus sp. if that’s truly all you have.

## Output

taxonomy_with_taxid.csv — same columns as input plus NCBI_TaxID.

A lookup cache taxid_cache.json (created automatically when --resume is used).

## Options (CLI)

--email           (required) your contact email for NCBI

--api-key         (recommended) boosts rate limit (NCBI)

--sleep           override polite delay (auto: ~0.34s no key; ~0.12s with key)

--resume          use/update a JSON cache (default: taxid_cache.json)

--limit N         resolve only first N unique species (debug)

--no-rank-check   skip species-rank verification (faster, but less strict)

--no-progress     disable tqdm progress bar (CI logs, etc.)

--verbose         log each resolution to stdout


# Performance & API etiquette

NCBI recommends ~3 requests/sec without a key and ~10 requests/sec with an API key. The script auto-throttles accordingly; you can tune --sleep. Exceeding limits may cause errors or temporary blocks.

## Notes & troubleshooting

Ambiguous or old names: we normalize a few common synonyms. You can extend the SYNONYMS dict in add_taxids.py. 

Placeholders (Genus sp.): NCBI often has placeholder taxa (e.g., Acidiphilium sp.) with their own TaxIDs; this is expected behavior.

Proxies/firewalls: If responses look like HTML (login pages), configure your proxy or run on an open network.

Security: Don’t commit your API key. Prefer passing it at runtime or via an env var (then read it in Python and pass to Entrez).


