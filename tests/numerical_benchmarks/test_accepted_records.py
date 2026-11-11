"""CLI and bundle path for records that already have a frozen reference."""

import json
import subprocess
from pathlib import Path

import pytest
from examples.iv_bolus_case import EXAMPLE_DIR

from pharmacoweave.accepted_records import evaluate_accepted_record
from pharmacoweave.cli import main
from pharmacoweave.errors import AnalysisError
from pharmacoweave.qsp.turnover import simulate_linear_turnover


def test_accepted_record_matches_the_direct_call_and_rejects_an_engine_response(tmp_path: Path) -> None:
    request = {
        "kind": "linear_turnover",
        "kin": 2,
        "kout": 1,
        "r0": 0,
        "times_h": [0, 1],
        "context_of_use": "Research check through the accepted-record command.",
    }
    direct = simulate_linear_turnover(
        kin=2,
        kout=1,
        r0=0,
        times_h=[0, 1],
        context_of_use="Research check through the accepted-record command.",
    )
    evaluated = evaluate_accepted_record(request)
    assert evaluated.curve[1].response == pytest.approx(direct.curve[1].response, rel=1e-12, abs=1e-12)
    with pytest.raises(AnalysisError, match="engine response"):
        evaluate_accepted_record({"kind": "engine_response", "text": "not interpreted"})

    request_path = tmp_path / "request.json"
    output_path = tmp_path / "result.json"
    bundle = tmp_path / "bundle"
    request_path.write_text(json.dumps(request), encoding="utf-8")
    assert (
        main(
            [
                "accepted-record",
                str(request_path),
                "--output",
                str(output_path),
                "--study",
                str(EXAMPLE_DIR),
                "--bundle",
                str(bundle),
                "--run-id",
                "accepted-record",
            ]
        )
        == 0
    )
    written = json.loads(output_path.read_text(encoding="utf-8"))
    stored = json.loads((bundle / "analysis" / "accepted_record.json").read_text(encoding="utf-8"))
    assert written["curve"][1]["response"] == pytest.approx(direct.curve[1].response, rel=1e-12, abs=1e-12)
    assert stored["curve"][1]["response"] == pytest.approx(direct.curve[1].response, rel=1e-12, abs=1e-12)


def test_population_table_request_does_not_start_a_subprocess(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_if_called(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("an accepted population table started a subprocess")

    monkeypatch.setattr(subprocess, "run", fail_if_called)
    loaded = evaluate_accepted_record(
        {
            "kind": "population_result_table",
            "rows": [
                {
                    "subject_id": "SUBJ-01",
                    "parameter_name": "CL",
                    "estimate": 0.2,
                    "unit": "L/h",
                    "estimation_method": "caller",
                    "evidence_grade": "estimated",
                }
            ],
        }
    )
    assert loaded.rows[0].parameter_name == "CL"
