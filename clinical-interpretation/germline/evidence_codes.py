"""
ACMG/AMP 2015 evidence codes and strength levels.

Reference: Richards S, Aziz N, Bale S, et al. Standards and guidelines
for the interpretation of sequence variants: a joint consensus
recommendation of the American College of Medical Genetics and
Genomics and the Association for Molecular Pathology. Genet Med.
2015;17(5):405-424.
"""
from dataclasses import dataclass
from enum import Enum


class Strength(Enum):
    STAND_ALONE = "stand_alone"
    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"


class Direction(Enum):
    PATHOGENIC = "pathogenic"
    BENIGN = "benign"


@dataclass(frozen=True)
class EvidenceCode:
    code: str
    strength: Strength
    direction: Direction
    description: str


PVS1 = EvidenceCode("PVS1", Strength.VERY_STRONG, Direction.PATHOGENIC,
    "Null variant in a gene where loss of function is a known mechanism of disease.")
PS1 = EvidenceCode("PS1", Strength.STRONG, Direction.PATHOGENIC,
    "Same amino acid change as a previously established pathogenic variant.")
PS2 = EvidenceCode("PS2", Strength.STRONG, Direction.PATHOGENIC,
    "De novo (confirmed maternity/paternity) in a patient with the disease.")
PS3 = EvidenceCode("PS3", Strength.STRONG, Direction.PATHOGENIC,
    "Well-established functional studies show a damaging effect.")
PS4 = EvidenceCode("PS4", Strength.STRONG, Direction.PATHOGENIC,
    "Prevalence in affected individuals significantly increased vs controls.")
PM1 = EvidenceCode("PM1", Strength.MODERATE, Direction.PATHOGENIC,
    "Located in a mutational hot spot / critical functional domain.")
PM2 = EvidenceCode("PM2", Strength.MODERATE, Direction.PATHOGENIC,
    "Absent (or extremely low frequency) in population databases.")
PM3 = EvidenceCode("PM3", Strength.MODERATE, Direction.PATHOGENIC,
    "For recessive disorders, detected in trans with a pathogenic variant.")
PM4 = EvidenceCode("PM4", Strength.MODERATE, Direction.PATHOGENIC,
    "Protein length change in a non-repeat region or stop-loss.")
PM5 = EvidenceCode("PM5", Strength.MODERATE, Direction.PATHOGENIC,
    "Novel missense at a residue with a different established pathogenic missense.")
PM6 = EvidenceCode("PM6", Strength.MODERATE, Direction.PATHOGENIC,
    "Assumed de novo without confirmed parentage.")
PP1 = EvidenceCode("PP1", Strength.SUPPORTING, Direction.PATHOGENIC,
    "Cosegregation with disease in multiple affected family members.")
PP2 = EvidenceCode("PP2", Strength.SUPPORTING, Direction.PATHOGENIC,
    "Missense in a gene with low benign missense rate.")
PP3 = EvidenceCode("PP3", Strength.SUPPORTING, Direction.PATHOGENIC,
    "Multiple computational evidence lines support a deleterious effect.")
PP4 = EvidenceCode("PP4", Strength.SUPPORTING, Direction.PATHOGENIC,
    "Patient phenotype highly specific for a gene with single known etiology.")
PP5 = EvidenceCode("PP5", Strength.SUPPORTING, Direction.PATHOGENIC,
    "Reputable source reports pathogenic, evidence unavailable for review.")

BA1 = EvidenceCode("BA1", Strength.STAND_ALONE, Direction.BENIGN,
    "Allele frequency >5% in a population database.")
BS1 = EvidenceCode("BS1", Strength.STRONG, Direction.BENIGN,
    "Allele frequency greater than expected for the disorder.")
BS2 = EvidenceCode("BS2", Strength.STRONG, Direction.BENIGN,
    "Observed in healthy adults for a fully penetrant early-onset disease.")
BS3 = EvidenceCode("BS3", Strength.STRONG, Direction.BENIGN,
    "Well-established functional studies show no damaging effect.")
BS4 = EvidenceCode("BS4", Strength.STRONG, Direction.BENIGN,
    "Lack of segregation in affected family members.")
BP1 = EvidenceCode("BP1", Strength.SUPPORTING, Direction.BENIGN,
    "Missense in a gene where primarily truncating variants cause disease.")
BP2 = EvidenceCode("BP2", Strength.SUPPORTING, Direction.BENIGN,
    "Observed in trans/cis with a pathogenic variant.")
BP3 = EvidenceCode("BP3", Strength.SUPPORTING, Direction.BENIGN,
    "In-frame indel in a repetitive region without known function.")
BP4 = EvidenceCode("BP4", Strength.SUPPORTING, Direction.BENIGN,
    "Multiple computational evidence lines suggest no deleterious effect.")
BP5 = EvidenceCode("BP5", Strength.SUPPORTING, Direction.BENIGN,
    "Variant found in a case with an alternate molecular basis for disease.")
BP6 = EvidenceCode("BP6", Strength.SUPPORTING, Direction.BENIGN,
    "Reputable source reports benign, evidence unavailable for review.")
BP7 = EvidenceCode("BP7", Strength.SUPPORTING, Direction.BENIGN,
    "Synonymous variant with no predicted splicing impact, not conserved.")

ALL_CODES = {c.code: c for c in [
    PVS1, PS1, PS2, PS3, PS4, PM1, PM2, PM3, PM4, PM5, PM6,
    PP1, PP2, PP3, PP4, PP5, BA1, BS1, BS2, BS3, BS4,
    BP1, BP2, BP3, BP4, BP5, BP6, BP7,
]}
