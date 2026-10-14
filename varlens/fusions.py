"""Supplied fusion RNA -> Isovar evidence -> Topiary protein fragments.

This adapter does not discover fusions from BAMs or infer a coding frame from
gene names. Isovar validates the supplied sequence, mappings, references and
read observations, retaining every supported frame and unresolved outcome.
"""

from hashlib import sha256
import json

from .errors import InputError
from .version import __version__


def analyze_fusion(payload, *, alleles=(), epitope_lengths=(8, 9, 10, 11),
                   min_fragments=2, predictor=None):
    """Return a JSON-serializable reconstruction, fragments and predictions.

    ``payload`` follows Isovar's supplied-fusion input contract (``fusion``,
    ``references``, ``reads``, optional ``reference_names``). The Sid fixture
    envelope fields ``note`` and ``original_records`` are preserved but not
    passed to Isovar as evidence. Unknown fields are rejected by Isovar.

    All coding hypotheses, including ambiguous/partial-CDS hypotheses, become
    separately identified fragments. Their status is carried onto prediction
    rows; no hypothesis is selected as the true protein. No alleles/predictor
    means reconstruction only. ``predictor`` may be a configured TopiaryPredictor
    for a Python caller; otherwise supplied alleles select mhctools.NetMHCpan.
    """
    try:
        import isovar
        from isovar.fusion import fusion_from_dict, reconstruct_fusion
    except ImportError as error:
        raise InputError("Fusion reconstruction requires pip install 'varlens[fusion]'") from error
    import topiary

    if not isinstance(payload, dict):
        raise InputError("Expected one supplied-fusion JSON object")
    if predictor is not None and alleles:
        raise InputError("Configure alleles on the predictor or pass alleles, not both")
    # Copy so caller mutations cannot change the recorded provenance.
    supplied = json.loads(json.dumps(payload, allow_nan=False))
    native = {k: v for k, v in supplied.items() if k not in ("note", "original_records")}
    try:
        reconstruction = reconstruct_fusion(
            *fusion_from_dict(native), peptide_lengths=epitope_lengths,
            min_fragments=min_fragments)
    except (ValueError, TypeError) as error:
        raise InputError("Invalid supplied fusion: %s" % error) from error
    if reconstruction.get("schema") != "isovar.fusion_rna.v3":
        raise InputError("Unsupported Isovar fusion result schema")

    fragments = _fragments(reconstruction)
    result = {
        "schema": "varlens.fusion.v1",
        "software": {"varlens": __version__, "isovar": isovar.__version__,
                     "topiary": topiary.__version__},
        "input": supplied,
        "reconstruction": reconstruction,
        "fragments": [fragment.to_dict() for fragment in fragments],
        "prediction_status": "not_requested",
        "predictions": [],
    }
    if alleles or predictor is not None:
        if not fragments:
            result["prediction_status"] = "no_coding_products"
        else:
            from .mhc_binding import fragment_predictions
            predictions = fragment_predictions(
                fragments, alleles=alleles, epitope_lengths=epitope_lengths,
                predictor=predictor)
            # pandas encodes missing predictions as JSON null, never NaN or 0.
            result["predictions"] = json.loads(predictions.to_json(orient="records"))
            result["prediction_status"] = "completed"
    return result


def _fragments(result):
    """Adapt only the validated result produced above, without retranslation."""
    from topiary import ProteinFragment

    fragments = []
    for path in result["paths"]:
        junction, = path["junctions"]
        support = junction["direct_support"]
        for index, translation in enumerate(path["translations"]):
            # A boundary inside a codon targets that residue. At a codon
            # boundary the zero-width interval requires a peptide to cross it.
            # Keep both sides of an inserted sequence as separate boundaries.
            boundaries = translation["junction_in_translated_cds"]
            targets = sorted({(b // 3, (b + 2) // 3) for b in boundaries})
            identity = [result["sample_id"], result["event_id"], path["path_id"],
                        index, path["sequence_sha256"], translation["amino_acids"]]
            identifier = "fusion:" + sha256(json.dumps(identity).encode()).hexdigest()
            fragments.append(ProteinFragment(
                fragment_id=identifier, source_type="sv:fusion",
                sequence=translation["amino_acids"], target_intervals=targets,
                variant=result["event_id"], sample_name=result["sample_id"],
                n_rna_alt_reads=support["reads"],
                n_rna_alt_fragments=support["fragments"],
                annotations={
                    "fusion_reconstruction_status": result["status"],
                    "fusion_reconstruction_reasons": result["reasons"],
                    "fusion_path_id": path["path_id"],
                    "fusion_translation_index": index,
                    "fusion_complete_5prime": translation["complete_5prime"],
                    "fusion_translation_observed": translation["translation_observed"],
                    "fusion_transcript_ids": translation["transcript_ids"],
                    "fusion_frame_evidence": translation["frame_evidence"],
                    "fusion_provenance": result["provenance"],
                    "fusion_rna_evidence_scope": "supplied_direct_junction_observations",
                    "fusion_junction_in_translated_cds": boundaries,
                },
            ))
    return fragments
