import pytest
import mhctools
from mhctools.base_predictor import BasePredictor

from varlens import mhc_binding


@pytest.fixture
def binding_predictor(monkeypatch):
    """Exercise real Topiary variant processing with deterministic MHC scores."""
    class Predictor(BasePredictor):
        calls = []

        def predict_peptides(self, peptides):
            self.calls.append(tuple(peptides))
            return mhctools.BindingPredictionCollection([
                mhctools.BindingPrediction(
                    allele=allele, peptide=peptide,
                    affinity=(20 if allele.endswith("02:02") else 100) + len(peptide),
                    prediction_method_name="test",
                )
                for peptide in peptides for allele in self.alleles
            ])

    mhc_binding.CACHED_BINDING_AFFINITIES.clear()
    mhc_binding.BINDING_PREDICTORS.clear()
    monkeypatch.setattr(mhctools, "NetMHCpan", Predictor)
    yield Predictor
    mhc_binding.CACHED_BINDING_AFFINITIES.clear()
    mhc_binding.BINDING_PREDICTORS.clear()
