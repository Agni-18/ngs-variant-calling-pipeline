"""
PVS1 (null variant) evidence from VEP consequence annotations, gated
by a curated gene-mechanism list. PVS1 only applies where LOF is an
established disease mechanism -- silence is the safe default.
"""
NULL_VARIANT_CONSEQUENCES = {
    "stop_gained", "frameshift_variant", "splice_acceptor_variant",
    "splice_donor_variant", "start_lost", "transcript_ablation",
}

# JAG1: haploinsufficiency established for Alagille syndrome.
# ASXL1: truncating variants established in myeloid malignancies / Bohring-Opitz.
# GNAS: deliberately excluded -- imprinting/parent-of-origin-dependent mechanism.
LOF_MECHANISM_GENES = {"JAG1", "ASXL1"}


def get_pvs1_evidence(gene_symbol: str, consequence: str) -> list[str]:
    if gene_symbol not in LOF_MECHANISM_GENES:
        return []
    consequence_terms = set(consequence.split("&"))
    if consequence_terms & NULL_VARIANT_CONSEQUENCES:
        return ["PVS1"]
    return []
