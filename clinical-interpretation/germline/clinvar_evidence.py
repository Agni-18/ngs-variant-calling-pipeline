"""
Evidence extraction from ClinVar. Cross-references variants against
a local ClinVar VCF subset to derive PP5/BP6 evidence codes.
"""
import gzip
from dataclasses import dataclass


@dataclass
class ClinVarRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    clnsig: str
    clnrevstat: str
    clndn: str


_PATHOGENIC_CLNSIG = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
_BENIGN_CLNSIG = {"Benign", "Likely_benign", "Benign/Likely_benign"}


def load_clinvar_index(clinvar_vcf_path: str) -> dict:
    index = {}
    opener = gzip.open if clinvar_vcf_path.endswith(".gz") else open
    with opener(clinvar_vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
            info = fields[7]

            info_dict = {}
            for entry in info.split(";"):
                if "=" in entry:
                    key, _, value = entry.partition("=")
                    info_dict[key] = value

            record = ClinVarRecord(
                chrom=chrom, pos=pos, ref=ref, alt=alt,
                clnsig=info_dict.get("CLNSIG", ""),
                clnrevstat=info_dict.get("CLNREVSTAT", ""),
                clndn=info_dict.get("CLNDN", ""),
            )
            key = (f"chr{chrom}" if not chrom.startswith("chr") else chrom, pos, ref, alt)
            index[key] = record
    return index


def get_clinvar_evidence(chrom, pos, ref, alt, clinvar_index):
    key = (chrom, pos, ref, alt)
    record = clinvar_index.get(key)
    if record is None:
        return [], None

    if record.clnsig in _PATHOGENIC_CLNSIG:
        return ["PP5"], record
    if record.clnsig in _BENIGN_CLNSIG:
        return ["BP6"], record

    return [], record
