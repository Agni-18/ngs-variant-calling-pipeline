"""Tests for the ACMG/AMP 2015 combining-rules classifier."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from germline.classifier import Classification, classify


class TestPathogenic:
    def test_very_strong_plus_strong(self):
        assert classify(["PVS1", "PS1"]).classification == Classification.PATHOGENIC

    def test_very_strong_plus_two_moderate(self):
        assert classify(["PVS1", "PM1", "PM2"]).classification == Classification.PATHOGENIC

    def test_very_strong_plus_moderate_plus_supporting(self):
        assert classify(["PVS1", "PM2", "PP3"]).classification == Classification.PATHOGENIC

    def test_very_strong_plus_two_supporting(self):
        assert classify(["PVS1", "PP1", "PP3"]).classification == Classification.PATHOGENIC

    def test_two_strong(self):
        assert classify(["PS1", "PS3"]).classification == Classification.PATHOGENIC

    def test_strong_plus_three_moderate(self):
        assert classify(["PS1", "PM1", "PM2", "PM5"]).classification == Classification.PATHOGENIC

    def test_strong_plus_two_moderate_plus_two_supporting(self):
        assert classify(["PS1", "PM1", "PM2", "PP1", "PP3"]).classification == Classification.PATHOGENIC

    def test_strong_plus_moderate_plus_four_supporting(self):
        assert classify(["PS1", "PM1", "PP1", "PP2", "PP3", "PP4"]).classification == Classification.PATHOGENIC


class TestLikelyPathogenic:
    def test_very_strong_plus_one_moderate(self):
        assert classify(["PVS1", "PM2"]).classification == Classification.LIKELY_PATHOGENIC

    def test_strong_plus_one_moderate(self):
        assert classify(["PS1", "PM2"]).classification == Classification.LIKELY_PATHOGENIC

    def test_strong_plus_two_moderate(self):
        assert classify(["PS1", "PM1", "PM2"]).classification == Classification.LIKELY_PATHOGENIC

    def test_strong_plus_two_supporting(self):
        assert classify(["PS1", "PP1", "PP3"]).classification == Classification.LIKELY_PATHOGENIC

    def test_three_moderate(self):
        assert classify(["PM1", "PM2", "PM5"]).classification == Classification.LIKELY_PATHOGENIC

    def test_two_moderate_plus_two_supporting(self):
        assert classify(["PM1", "PM2", "PP1", "PP3"]).classification == Classification.LIKELY_PATHOGENIC

    def test_one_moderate_plus_four_supporting(self):
        assert classify(["PM2", "PP1", "PP2", "PP3", "PP4"]).classification == Classification.LIKELY_PATHOGENIC


class TestBenign:
    def test_stand_alone(self):
        assert classify(["BA1"]).classification == Classification.BENIGN

    def test_two_strong_benign(self):
        assert classify(["BS1", "BS2"]).classification == Classification.BENIGN


class TestLikelyBenign:
    def test_strong_plus_supporting(self):
        assert classify(["BS1", "BP1"]).classification == Classification.LIKELY_BENIGN

    def test_two_supporting_benign(self):
        assert classify(["BP1", "BP4"]).classification == Classification.LIKELY_BENIGN


class TestUncertainSignificance:
    def test_insufficient_evidence(self):
        assert classify(["PP3"]).classification == Classification.UNCERTAIN_SIGNIFICANCE

    def test_empty_evidence(self):
        assert classify([]).classification == Classification.UNCERTAIN_SIGNIFICANCE

    def test_conflicting_evidence_reported_as_vus(self):
        result = classify(["PS1", "PS3", "BS1", "BS2"])
        assert result.classification == Classification.UNCERTAIN_SIGNIFICANCE
        assert "Conflicting" in result.rule_matched


class TestInputValidation:
    def test_unknown_code_raises(self):
        with pytest.raises(ValueError, match="Unknown evidence code"):
            classify(["PZ99"])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
