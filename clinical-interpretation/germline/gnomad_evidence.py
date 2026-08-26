"""
Population allele frequency evidence via gnomAD's public GraphQL API.
Uses per-variant queries rather than downloading a local gnomAD file
(the full dataset runs into hundreds of GB). Rate-limiting is real
(reports of blocking after ~10 rapid queries), so this module uses
retry-with-backoff and disk caching.
"""
import json
import time
from pathlib import Path

import requests

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
DEFAULT_DATASET = "gnomad_r4"
REQUEST_DELAY_SECONDS = 3.0
MAX_RETRIES = 4
INITIAL_BACKOFF_SECONDS = 5.0


def _cache_path(cache_dir):
    path = Path(cache_dir) / "gnomad_af_cache.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _load_cache(cache_dir):
    path = _cache_path(cache_dir)
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def _save_cache(cache_dir, cache):
    with open(_cache_path(cache_dir), "w") as f:
        json.dump(cache, f, indent=2)


def _to_gnomad_variant_id(chrom, pos, ref, alt):
    bare_chrom = chrom.replace("chr", "")
    return f"{bare_chrom}-{pos}-{ref}-{alt}"


def query_gnomad_af(chrom, pos, ref, alt, cache_dir="resources", dataset=DEFAULT_DATASET):
    cache = _load_cache(cache_dir)
    variant_id = _to_gnomad_variant_id(chrom, pos, ref, alt)

    if variant_id in cache:
        return cache[variant_id]

    query = """
    query VariantAF($variantId: String!, $dataset: DatasetId!) {
        variant(variantId: $variantId, dataset: $dataset) {
            genome { af }
            exome { af }
        }
    }
    """

    result = None
    backoff = INITIAL_BACKOFF_SECONDS
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = requests.post(
                GNOMAD_API_URL,
                json={"query": query, "variables": {"variantId": variant_id, "dataset": dataset}},
                headers={"Content-Type": "application/json"},
                timeout=15,
            )
            if response.status_code == 429:
                retry_after = response.headers.get("Retry-After")
                wait = float(retry_after) if retry_after else backoff
                print(f"  gnomAD rate-limited on {variant_id} (attempt {attempt}/{MAX_RETRIES}), waiting {wait:.0f}s...")
                time.sleep(wait)
                backoff *= 2
                continue
            response.raise_for_status()
            result = response.json()
            break
        except requests.RequestException as e:
            print(f"  gnomAD query failed for {variant_id} (attempt {attempt}): {e}")
            time.sleep(backoff)
            backoff *= 2

    time.sleep(REQUEST_DELAY_SECONDS)

    if result is None:
        print(f"  gnomAD lookup for {variant_id} could not be resolved after {MAX_RETRIES} attempts")
        return None

    data = result.get("data", {}).get("variant")
    if data is None:
        cache[variant_id] = None
        _save_cache(cache_dir, cache)
        return None

    af = None
    if data.get("genome") is not None:
        af = data["genome"].get("af")
    if af is None and data.get("exome") is not None:
        af = data["exome"].get("af")

    cache[variant_id] = af
    _save_cache(cache_dir, cache)
    return af


def get_frequency_evidence(af):
    if af is None:
        return ["PM2"]
    if af > 0.05:
        return ["BA1"]
    if af > 0.01:
        return ["BS1"]
    return []
