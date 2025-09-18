#!/usr/bin/env python3
"""
add_taxids.py — Resolve NCBI TaxIDs for a taxonomy CSV (phyloseq-style)
- Uses JSON for ESearch (parsed via json.loads)
- Species-rank verification via EFetch XML (decoded safely)
- Real-time progress bar (tqdm), verbose logs
- Checkpoint cache (resume-safe)
- Safe partial writes (.part file flushed continuously)

Input CSV must contain columns:
  Taxon, Domain, Phylum, Class, Order, Family, Genus, Species

Output CSV is identical + NCBI_TaxID
"""
import argparse, csv, json, os, sys, time, re
from typing import Dict, List
from urllib.error import HTTPError, URLError
from Bio import Entrez

# ---------- Utilities ----------
def clean_name(name: str) -> str:
    name = (name or "").strip()
    name = re.sub(r"\s+", " ", name.replace("_", " "))
    return name

def load_cache(path: str) -> Dict[str,str]:
    if path and os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}

def save_cache(path: str, cache: Dict[str,str]):
    if not path:
        return
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)

def open_out_writer(out_path: str, fieldnames):
    out_tmp = out_path + ".part"
    g = open(out_tmp, "w", newline="", encoding="utf-8")
    w = csv.DictWriter(g, fieldnames=fieldnames)
    w.writeheader()
    g.flush()
    return g, w, out_tmp

# Preferred synonyms (old -> new) to improve hit rates
SYNONYMS = {
    "Acholeplasma modicum": "Haploplasma modicum",
    "Propionibacterium acnes": "Cutibacterium acnes",
    "Clostridium difficile": "Clostridioides difficile",
}

def esearch_ids_json(term: str) -> List[str]:
    # Use JSON and parse manually to avoid Entrez.read XML requirements
    h = Entrez.esearch(db="taxonomy", term=term, retmode="json", retmax=20)
    raw = h.read()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", "ignore")
    data = json.loads(raw)
    return data.get("esearchresult", {}).get("idlist", [])

def efetch_xml_raw(taxid: str) -> str:
    h = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    raw = h.read()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", "ignore")
    return raw

def rank_is_species(xml_text: str) -> bool:
    return "<Rank>species</Rank>" in xml_text

def best_taxid_for_name(name: str, sleep_s: float, verbose: bool=False, rank_check: bool=True) -> str:
    name0 = clean_name(name)
    name1 = SYNONYMS.get(name0, name0)
    terms = (f'"{name1}"[SCIN]', f"{name1}[SCIN]", name1)
    for term in terms:
        try:
            ids = esearch_ids_json(term)
            if ids:
                if rank_check:
                    for tid in ids:
                        xml = efetch_xml_raw(tid)
                        if rank_is_species(xml):
                            if verbose:
                                print(f"[species] {name0} -> {tid}", flush=True)
                            time.sleep(sleep_s)
                            return tid
                # if no species-ranked match, accept first id
                if verbose:
                    print(f"[first]   {name0} -> {ids[0]}", flush=True)
                time.sleep(sleep_s)
                return ids[0]
        except (HTTPError, URLError) as e:
            if verbose:
                print(f"[warn] HTTP error for {name0} term={term!r}: {e}", flush=True)
            time.sleep(max(0.8, sleep_s*2))
        except Exception as e:
            if verbose:
                print(f"[warn] parse/error for {name0} term={term!r}: {e}", flush=True)
            time.sleep(max(0.5, sleep_s*2))
        finally:
            time.sleep(sleep_s)
    if verbose:
        print(f"[miss]    {name0} -> ''", flush=True)
    return ""

def main():
    pa = argparse.ArgumentParser()
    pa.add_argument("in_csv", help="Input taxonomy CSV (phyloseq-style)")
    pa.add_argument("out_csv", help="Output CSV with NCBI_TaxID")
    pa.add_argument("--email", required=True, help="Entrez.email")
    pa.add_argument("--api-key", dest="api_key", default=None, help="NCBI API key (recommended)")
    pa.add_argument("--cache", default="taxid_cache.json", help="Cache file for Species->TaxID")
    pa.add_argument("--sleep", type=float, default=None, help="Seconds to sleep between calls (auto if omitted)")
    pa.add_argument("--resume", action="store_true", help="Resume using existing cache/out file if present")
    pa.add_argument("--limit", type=int, default=None, help="Resolve only the first N species (debug)")
    pa.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bar")
    pa.add_argument("--verbose", action="store_true", help="Print each resolution")
    pa.add_argument("--no-rank-check", action="store_true", help="Skip species-rank verification (faster)")
    args = pa.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    # Auto sleep based on NCBI guidance: ~3 rps w/o key, ~10 rps w/ key
    sleep_s = args.sleep if args.sleep is not None else (0.34 if not args.api_key else 0.12)

    # Read input
    with open(args.in_csv, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        rows = list(r)
    if not rows:
        print("Input CSV has no rows.", file=sys.stderr); sys.exit(1)

    fieldnames = list(rows[0].keys())
    if "Species" not in fieldnames:
        print("Input CSV must have a 'Species' column.", file=sys.stderr); sys.exit(2)
    if "NCBI_TaxID" not in fieldnames:
        fieldnames.append("NCBI_TaxID")

    # Build unique species
    uniq = []
    for row in rows:
        nm = clean_name(row.get("Species",""))
        if nm:
            uniq.append(nm)
    uniq = sorted(set(uniq))
    if args.limit:
        uniq = uniq[:args.limit]

    # Cache
    cache = load_cache(args.cache) if args.resume else {}
    if args.verbose:
        print(f"Unique species: {len(uniq)}; cached: {len(cache)}; sleep={sleep_s}s", flush=True)

    # Output writer
    g, w, out_tmp = open_out_writer(args.out_csv, fieldnames)

    # Progress bar
    try:
        from tqdm import tqdm
        pb = None if args.no_progress else tqdm(total=len(uniq), unit="sp.", desc="Resolving NCBI TaxIDs")
    except Exception:
        pb = None
        if not args.no_progress:
            print("[note] tqdm not installed; pip install tqdm", flush=True)

    # Resolve + stream write
    try:
        for row in rows:
            nm = clean_name(row.get("Species",""))
            if nm and nm not in cache:
                cache[nm] = best_taxid_for_name(nm, sleep_s=sleep_s, verbose=args.verbose, rank_check=(not args.no_rank_check))
                if pb: pb.update(1)
                if len(cache) % 25 == 0:
                    save_cache(args.cache, cache)
            row["NCBI_TaxID"] = cache.get(nm, "")
            w.writerow(row)
            g.flush()
    except KeyboardInterrupt:
        print("\n[Interrupted] Writing partial output and saving cache...", flush=True)
        save_cache(args.cache, cache)
        g.flush(); g.close()
        os.replace(out_tmp, args.out_csv)
        print(f"Partial file saved to: {args.out_csv}", flush=True)
        sys.exit(130)
    except Exception as e:
        print(f"\n[Error] {e}", file=sys.stderr)
        save_cache(args.cache, cache)
        g.flush(); g.close()
        os.replace(out_tmp, args.out_csv)
        print(f"[Recovered] Partial file saved to: {args.out_csv}", flush=True)
        sys.exit(1)
    finally:
        if pb: pb.close()

    # Finalize
    save_cache(args.cache, cache)
    g.flush(); g.close()
    os.replace(out_tmp, args.out_csv)
    resolved = sum(1 for s in uniq if cache.get(s))
    print(f"Done. Resolved {resolved}/{len(uniq)} unique species. Output: {args.out_csv}", flush=True)

if __name__ == "__main__":
    main()
