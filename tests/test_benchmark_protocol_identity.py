"""Frozen protocol byte identities must survive Git checkout and source archives."""
import hashlib
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("relative, expected", [
    ("benchmarks/zymo-even/protocol.json", "30c859006c71a974a99baefba3644161afa1fc8fd695fd439851548f27b07ea0"),
    ("benchmarks/comparison/protocol.json", "44b236dd6e3bbd97c0cab31d1f59e0e61c3351cfc4ccfc939881e2adb6253b5a"),
    ("benchmarks/comparison/tied-null-protocol.json", "5905dfa525aa219dfd383cc6b7802180fc6ca5c857032b71213a7995e6a2cef5"),
])
def test_frozen_protocol_preserves_recorded_byte_identity(relative, expected):
    assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected
