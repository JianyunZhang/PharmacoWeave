"""A caller executable starts only with a file, a version, and license text."""

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import pytest

from pharmacoweave.cli import main
from pharmacoweave.errors import AnalysisError, PlannedCapabilityError
from pharmacoweave.integrations.external_backends import run_backend, run_caller_engine, supported_backends

_LICENSE = "Caller-supplied license text for this research run.\n"


def test_missing_credentials_do_not_start_a_subprocess(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fail_if_called(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("subprocess was started without a complete caller engine slot")

    monkeypatch.setattr(subprocess, "run", fail_if_called)
    executable = _echo_executable(tmp_path)
    with pytest.raises(AnalysisError, match="all present"):
        run_caller_engine(tmp_path / "missing.cmd", "request", engine_version="1.0", license_text=_LICENSE)
    with pytest.raises(AnalysisError, match="all present"):
        run_caller_engine(executable, "request", engine_version="  ", license_text=_LICENSE)
    with pytest.raises(AnalysisError, match="all present"):
        run_caller_engine(executable, "request", engine_version="1.0", license_text=" \n")
    empty = tmp_path / "empty.cmd"
    empty.write_bytes(b"")
    with pytest.raises(AnalysisError, match="all present"):
        run_caller_engine(empty, "request", engine_version="1.0", license_text=_LICENSE)
    assert supported_backends() == ["pharmacoweave-analytical"]


def test_complete_slot_records_file_hash_and_leaves_the_response_uninterpreted(tmp_path: Path) -> None:
    executable = _echo_executable(tmp_path)
    request = "slot-request-v1\n"
    recorded = run_caller_engine(
        executable,
        request,
        engine_version=" 1.2.3 ",
        license_text=_LICENSE,
    )
    assert recorded.engine_version == "1.2.3"
    assert recorded.license_text == _LICENSE
    assert recorded.executable_sha256 == hashlib.sha256(executable.read_bytes()).hexdigest()
    assert recorded.request_sha256 == hashlib.sha256(request.encode("utf-8")).hexdigest()
    assert recorded.response_text == request
    assert recorded.response_sha256 == hashlib.sha256(request.encode("utf-8")).hexdigest()
    assert recorded.exit_code == 0
    assert recorded.constitutes_dosing_recommendation is False
    assert any("not a PharmacoWeave estimate" in item for item in recorded.limitations)
    assert "nlmixr2" in " ".join(recorded.limitations)
    assert supported_backends() == ["pharmacoweave-analytical"]


def test_failed_child_is_not_stored_as_a_result(tmp_path: Path) -> None:
    executable = _failing_executable(tmp_path)
    with pytest.raises(AnalysisError, match="exited 3"):
        run_caller_engine(executable, "request", engine_version="1.0", license_text=_LICENSE)


def test_catalog_engines_stay_separate_from_the_slot(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_if_called(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("catalog engine started a subprocess")

    monkeypatch.setattr(subprocess, "run", fail_if_called)
    with pytest.raises(PlannedCapabilityError):
        run_backend("nlmixr2")
    with pytest.raises(PlannedCapabilityError):
        run_backend("pk-sim")


def test_cli_slot_writes_the_recorded_run(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    executable = _echo_executable(tmp_path)
    request = tmp_path / "request.txt"
    license_file = tmp_path / "LICENSE.txt"
    output = tmp_path / "run.json"
    request.write_text("cli-request\n", encoding="utf-8")
    license_file.write_text(_LICENSE, encoding="utf-8")
    assert (
        main(
            [
                "engine-slot",
                str(request),
                "--executable",
                str(executable),
                "--engine-version",
                "0.9",
                "--license-file",
                str(license_file),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    printed = capsys.readouterr().out
    assert f"executable_sha256={hashlib.sha256(executable.read_bytes()).hexdigest()}" in printed
    assert "engine_version=0.9" in printed
    assert "dosing_recommendation=false" in printed
    assert "not a PharmacoWeave estimate" in output.read_text(encoding="utf-8")


def _echo_executable(directory: Path) -> Path:
    if os.name == "nt":
        path = directory / "echo_engine.cmd"
        python = str(Path(sys.executable))
        path.write_text(
            "@echo off\r\n"
            f'"{python}" -c "import pathlib,sys; pathlib.Path(sys.argv[2]).write_bytes('
            'pathlib.Path(sys.argv[1]).read_bytes())" "%~1" "%~2"\r\n'
            "exit /b 0\r\n",
            encoding="utf-8",
        )
        return path
    path = directory / "echo_engine.sh"
    path.write_text('#!/bin/sh\ncp "$1" "$2"\n', encoding="utf-8")
    path.chmod(0o755)
    return path


def _failing_executable(directory: Path) -> Path:
    if os.name == "nt":
        path = directory / "fail_engine.cmd"
        path.write_text("@echo off\r\nexit /b 3\r\n", encoding="utf-8")
        return path
    path = directory / "fail_engine.sh"
    path.write_text("#!/bin/sh\nexit 3\n", encoding="utf-8")
    path.chmod(0o755)
    return path
