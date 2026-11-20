"""Configuration errors must not silently select different pipeline settings."""

from dataclasses import asdict

import pytest
import yaml

from metaquest.exceptions import ConfigError
from metaquest.settings import get_config, load_config, reset_config


@pytest.fixture(autouse=True)
def clear_config_cache():
    reset_config()
    yield
    reset_config()


@pytest.mark.parametrize("section", [
    "databases", "assembly", "classification", "annotation", "abundance",
    "preprocessing", "screening", "reporting",
])
def test_unknown_section_setting_is_rejected(tmp_path, section):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({section: {"typo": 1}}), encoding="utf-8")

    with pytest.raises(ConfigError, match=rf"{section}\.typo"):
        load_config(path)


@pytest.mark.parametrize("section", ["classificaton", "ml", 123])
def test_unknown_top_level_section_is_rejected(tmp_path, section):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({section: {}}), encoding="utf-8")

    with pytest.raises(ConfigError, match=str(section)):
        load_config(path)


@pytest.mark.parametrize("content", ["[]", "false", "0", "settings", "[broken"])
def test_invalid_yaml_document_is_a_configuration_error(tmp_path, content):
    path = tmp_path / "config.yaml"
    path.write_text(content, encoding="utf-8")

    with pytest.raises(ConfigError, match="config.yaml"):
        load_config(path)


@pytest.mark.parametrize("section", ["databases", "classification", "reporting"])
@pytest.mark.parametrize("value", [None, False, 0, [], "threads=2"])
def test_section_must_be_a_mapping(tmp_path, section, value):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({section: value}), encoding="utf-8")

    with pytest.raises(ConfigError, match=section):
        load_config(path)


def test_empty_config_preserves_defaults(tmp_path):
    expected = asdict(load_config())
    path = tmp_path / "config.yaml"
    path.write_text("# No overrides\n", encoding="utf-8")

    assert asdict(load_config(path)) == expected


def test_valid_partial_overrides_and_database_precedence(tmp_path, monkeypatch):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({
        "classification": {"threads": 2},
        "reporting": {"plot_formats": ["pdf"]},
        "databases": {"base_dir": str(tmp_path / "yaml-db")},
    }), encoding="utf-8")
    monkeypatch.setenv("METAQUEST_DB_DIR", str(tmp_path / "env-db"))

    config = load_config(path)
    assert config.classification.threads == 2
    assert config.classification.bracken_threshold == 10
    assert config.reporting.plot_formats == ("pdf",)
    assert config.databases.base_dir == tmp_path / "env-db"
    assert load_config(path, db_dir=tmp_path / "cli-db").databases.base_dir == tmp_path / "cli-db"


def test_failed_load_does_not_replace_cached_configuration(tmp_path):
    previous = load_config()
    path = tmp_path / "config.yaml"
    path.write_text("classification:\n  typo: 2\n", encoding="utf-8")

    with pytest.raises(ConfigError):
        load_config(path)
    assert get_config() is previous


def test_directory_config_path_is_a_configuration_error(tmp_path):
    with pytest.raises(ConfigError, match="Cannot read configuration"):
        load_config(tmp_path)
