"""
ACMG/AMP 2015 classification engine.
Implements the evidence-combining rules from Table 5 of Richards et al. 2015.
"""
from collections import Counter
from dataclasses import dataclass
from enum import Enum

from .evidence_codes import ALL_CODES, Direction, Strength


class Classification(Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    UNCERTAIN_SIGNIFICANCE = "Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


@dataclass
class ClassificationResult:
    classification: Classification
    applied_codes: list[str]
    rule_matched: str


def _counts_by_strength(codes: list[str], direction: Direction) -> Counter:
    counts = Counter()
    for code in codes:
        evidence = ALL_CODES[code]
        if evidence.direction == direction:
            counts[evidence.strength] += 1
    return counts


def classify(codes: list[str]) -> ClassificationResult:
    for code in codes:
        if code not in ALL_CODES:
            raise ValueError(f"Unknown evidence code '{code}'.")

    path_result = _apply_pathogenic_rules(codes)
    benign_result = _apply_benign_rules(codes)

    if path_result is not None and benign_result is not None:
        return ClassificationResult(
            classification=Classification.UNCERTAIN_SIGNIFICANCE,
            applied_codes=codes,
            rule_matched=(
                f"Conflicting evidence: pathogenic rule "
                f"'{path_result.rule_matched}' and benign rule "
                f"'{benign_result.rule_matched}' both satisfied"
            ),
        )
    if path_result is not None:
        return path_result
    if benign_result is not None:
        return benign_result

    return ClassificationResult(
        classification=Classification.UNCERTAIN_SIGNIFICANCE,
        applied_codes=codes,
        rule_matched="No combining rule satisfied (insufficient evidence)",
    )


def _apply_pathogenic_rules(codes: list[str]):
    c = _counts_by_strength(codes, Direction.PATHOGENIC)
    very_strong = c[Strength.VERY_STRONG]
    strong = c[Strength.STRONG]
    moderate = c[Strength.MODERATE]
    supporting = c[Strength.SUPPORTING]

    if very_strong >= 1 and (
        strong >= 1 or moderate >= 2
        or (moderate >= 1 and supporting >= 1) or supporting >= 2
    ):
        return ClassificationResult(Classification.PATHOGENIC, codes,
            "1 Very Strong + (>=1 Strong OR >=2 Moderate OR 1 Moderate+1 Supporting OR >=2 Supporting)")
    if strong >= 2:
        return ClassificationResult(Classification.PATHOGENIC, codes, ">=2 Strong")
    if strong >= 1 and (
        moderate >= 3 or (moderate >= 2 and supporting >= 2) or (moderate >= 1 and supporting >= 4)
    ):
        return ClassificationResult(Classification.PATHOGENIC, codes,
            "1 Strong + (>=3 Moderate OR 2 Moderate+2 Supporting OR 1 Moderate+4 Supporting)")

    if very_strong >= 1 and moderate >= 1:
        return ClassificationResult(Classification.LIKELY_PATHOGENIC, codes, "1 Very Strong + 1 Moderate")
    if strong >= 1 and 1 <= moderate <= 2:
        return ClassificationResult(Classification.LIKELY_PATHOGENIC, codes, "1 Strong + 1-2 Moderate")
    if strong >= 1 and supporting >= 2:
        return ClassificationResult(Classification.LIKELY_PATHOGENIC, codes, "1 Strong + >=2 Supporting")
    if moderate >= 3:
        return ClassificationResult(Classification.LIKELY_PATHOGENIC, codes, ">=3 Moderate")
    if moderate >= 2 and supporting >= 2:
        return ClassificationResult(Classification.LIKELY_PATHOGENIC, codes, "2 Moderate + >=2 Supporting")
    if moderate >= 1 and supporting >= 4:
        return ClassificationResult(Classification.LIKELY_PATHOGENIC, codes, "1 Moderate + >=4 Supporting")

    return None


def _apply_benign_rules(codes: list[str]):
    c = _counts_by_strength(codes, Direction.BENIGN)
    stand_alone = c[Strength.STAND_ALONE]
    strong = c[Strength.STRONG]
    supporting = c[Strength.SUPPORTING]

    if stand_alone >= 1:
        return ClassificationResult(Classification.BENIGN, codes, "1 Stand-Alone (BA1)")
    if strong >= 2:
        return ClassificationResult(Classification.BENIGN, codes, ">=2 Strong")

    if strong >= 1 and supporting >= 1:
        return ClassificationResult(Classification.LIKELY_BENIGN, codes, "1 Strong + 1 Supporting")
    if supporting >= 2:
        return ClassificationResult(Classification.LIKELY_BENIGN, codes, ">=2 Supporting")

    return None
