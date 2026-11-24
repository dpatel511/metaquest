"""Host-filter accounting and native controls; SAM text POS is 1-based."""
from dataclasses import asdict, replace
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import signal
import shutil
import subprocess
import sys
import time

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pytest

from metaquest.exceptions import ConfigError, MetaQuestError
from metaquest.settings import load_config


def host_module():
    from metaquest.core import host_removal
    return host_removal


def write_reads(path, records):
    with path.open('w') as handle:
        for name, sequence in records:
            handle.write(f'@{name}\n{sequence}\n+\n{"I" * len(sequence)}\n')
    return path


def read_output(path):
    with gzip.open(path, 'rt') as handle:
        return list(FastqGeneralIterator(handle))


@pytest.mark.parametrize('large', [False, True])
def test_host_index_hashes_every_shard_and_detects_changes(tmp_path, large):
    host = host_module()
    prefix = tmp_path / 'host with spaces'
    suffix = 'bt2l' if large else 'bt2'
    for number in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
        Path(f'{prefix}.{number}.{suffix}').write_bytes(number.encode())
    before = host.index_identity(prefix)
    assert len(before['sha256']) == 6
    Path(f'{prefix}.rev.2.{suffix}').write_bytes(b'changed')
    assert host.index_identity(prefix) != before


@pytest.mark.parametrize('kind', ['missing', 'partial', 'mixed', 'empty'])
def test_host_index_preflight_rejects_incomplete_or_ambiguous_prefix(tmp_path, kind):
    host = host_module()
    prefix = tmp_path / 'host'
    if kind != 'missing':
        for number in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
            Path(f'{prefix}.{number}.bt2').write_bytes(b'index')
        if kind == 'partial':
            Path(f'{prefix}.rev.2.bt2').unlink()
        elif kind == 'mixed':
            Path(f'{prefix}.1.bt2l').write_bytes(b'other')
        else:
            Path(f'{prefix}.3.bt2').write_bytes(b'')
    with pytest.raises(ConfigError, match='index'):
        host.index_files(prefix)


def test_host_removal_is_disabled_by_default_and_cli_overrides_yaml(tmp_path):
    config = load_config()
    assert config.preprocessing.host_index is None
    path = tmp_path / 'config.yaml'
    path.write_text('preprocessing:\n  host_index: configured-prefix\n')
    assert load_config(path).preprocessing.host_index == Path('configured-prefix')
    assert load_config(path, host_index=tmp_path / 'cli-prefix').preprocessing.host_index == tmp_path / 'cli-prefix'


@pytest.mark.parametrize('value', ['false', '42', '[]'])
def test_invalid_host_index_type_is_configuration_error(tmp_path, value):
    config = tmp_path / 'config.yaml'
    config.write_text(f'preprocessing:\n  host_index: {value}\n')
    with pytest.raises(ConfigError, match='host_index'):
        load_config(config)


NATIVE = bool(shutil.which('bowtie2') and shutil.which('bowtie2-build') and shutil.which('samtools'))


@pytest.fixture
def native_reference(tmp_path):
    rng = random.Random(20261123)
    host = ''.join(rng.choices('ACGT', k=12000))
    microbe = ''.join(rng.choices('ACGT', k=12000))
    reference = tmp_path / 'reference.fa'
    reference.write_text(f'>host\n{host}\n')
    prefix = tmp_path / 'index with spaces'
    subprocess.run(['bowtie2-build', str(reference), str(prefix)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return prefix, host, microbe


def control_reads(tmp_path, host, microbe, mode, *, duplicates=False):
    def reverse(sequence):
        return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    left, right = [], []
    for number, (first, second) in enumerate(((host, host), (host, microbe), (microbe, host), (microbe, microbe))):
        name = 'same' if duplicates else f'fragment_{number}'
        left.append((name + '/1', first[100:250]))
        right.append((name + '/2', reverse(second[450:600])))
    r1 = write_reads(tmp_path / 'R1.fastq', left)
    if mode == 'single':
        return [r1]
    return [r1, write_reads(tmp_path / 'R2.fastq', right)]


@pytest.mark.skipif(not NATIVE, reason='native Bowtie2/Samtools controls')
@pytest.mark.parametrize('mode,expected', [('paired', 1), ('single', 2)])
def test_native_host_filter_removes_any_matching_mate_and_accounts(tmp_path, native_reference, mode, expected):
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, host, microbe, mode)
    output, summary = host_module().run_host_removal(reads, tmp_path / 'result', prefix,
                                                   paired=mode != 'single', threads=2)
    assert summary['input'] == 4
    assert summary['retained'] == expected
    assert summary['removed'] == 4 - expected
    assert summary['count_unit'] == ('reads' if mode == 'single' else 'fragments')
    assert all(len(read_output(path)) == expected for path in output)
    expected_ids = {'fragment_2', 'fragment_3'} if mode == 'single' else {'fragment_3'}
    assert {row[0].removesuffix('/1') for row in read_output(output[0])} == expected_ids
    assert len(summary['state']['index']['sha256']) == 6
    assert summary['state']['pair_policy'] == 'remove_pair_if_either_maps'
    assert summary['state']['tool_versions']['bowtie2']
    assert not summary['reused']
    previous = [(p, p.read_bytes()) for p in output]
    _, reused = host_module().run_host_removal(reads, tmp_path / 'result', prefix,
                                              paired=mode != 'single', threads=2, reuse=True)
    assert reused['reused']
    assert all(path.read_bytes() == raw for path, raw in previous)
    output[0].write_bytes(b'corrupted')
    _, rebuilt = host_module().run_host_removal(reads, tmp_path / 'result', prefix,
                                               paired=mode != 'single', threads=2, reuse=True)
    assert not rebuilt['reused']
    assert len(read_output(output[0])) == expected


@pytest.mark.skipif(not NATIVE, reason='native Bowtie2/Samtools controls')
def test_duplicate_name_conversion_loss_fails_before_completion(tmp_path, native_reference):
    prefix, host, microbe = native_reference
    # Two consecutive non-host pairs share the same SAM query name.
    reads = control_reads(tmp_path, microbe, microbe, 'paired', duplicates=True)
    with pytest.raises(MetaQuestError, match='count|identifier'):
        host_module().run_host_removal(reads, tmp_path / 'result', prefix, paired=True, threads=2)
    assert not (tmp_path / 'result/preprocessing/host_removal/summary.json').exists()


@pytest.mark.skipif(not NATIVE, reason='native Bowtie2/Samtools controls')
def test_preprocessing_disabled_can_still_remove_host_and_reports_zero_retained(tmp_path, native_reference):
    from metaquest.pipeline.context import PipelineContext
    from metaquest.pipeline.stages.preprocessing import run_preprocessing_stage
    prefix, host, _ = native_reference
    reads = control_reads(tmp_path, host, host, 'paired')
    config = load_config()
    config = replace(config, preprocessing=replace(config.preprocessing, enabled=False, host_index=prefix))
    ctx = PipelineContext(config, reads, tmp_path / 'result')
    with pytest.raises(MetaQuestError, match='no reads|No reads'):
        run_preprocessing_stage(ctx)
    summary = json.loads((ctx.output_dir / 'preprocessing/host_removal/summary.json').read_text())
    assert summary['retained'] == 0 and summary['removed'] == summary['input'] == 4


def resume_fixture(tmp_path, config):
    from metaquest.core.analysis import _jsonable
    from metaquest.pipeline.runner import _database_identity
    reads = write_reads(tmp_path / 'input.fastq', [('read', 'ACGT' * 25)])
    flags = {'taxonomy_only': True, 'skip_functional': True}
    metadata = {'status': 'completed', 'metaquest_version': '2.0.0b4',
                'input_files': [{'path': str(reads.resolve()), 'size_bytes': reads.stat().st_size,
                                 'sha256': hashlib.sha256(reads.read_bytes()).hexdigest()}],
                'effective_config': _jsonable(asdict(config)), 'workflow': flags,
                'database_identity': _database_identity(config, **flags), 'tool_provenance': {}}
    path = tmp_path / 'analysis_metadata.json'
    path.write_text(json.dumps(metadata))
    return reads, flags, metadata, path


def test_legacy_resume_without_host_setting_is_unchanged(tmp_path):
    from metaquest.core.analysis import _validate_resume
    config = load_config(db_dir=tmp_path / 'db')
    reads, flags, metadata, path = resume_fixture(tmp_path, config)
    del metadata['effective_config']['preprocessing']['host_index']
    path.write_text(json.dumps(metadata))
    original = path.read_bytes()
    _validate_resume(tmp_path, [reads], config, flags)
    assert path.read_bytes() == original
    enabled = replace(config, preprocessing=replace(config.preprocessing, host_index=tmp_path / 'host'))
    with pytest.raises(ConfigError, match='configuration'):
        _validate_resume(tmp_path, [reads], enabled, flags)


def test_resume_detects_index_changes_even_when_size_and_mtime_match(tmp_path):
    from metaquest.core.analysis import _validate_resume
    prefix = tmp_path / 'host'
    for number in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
        Path(f'{prefix}.{number}.bt2').write_bytes(b'AAA')
    config = load_config(db_dir=tmp_path / 'db', host_index=prefix)
    reads, flags, _, _ = resume_fixture(tmp_path, config)
    path = Path(f'{prefix}.1.bt2')
    old = path.stat()
    path.write_bytes(b'TTT')
    os.utime(path, ns=(old.st_atime_ns, old.st_mtime_ns))
    with pytest.raises(ConfigError, match='database identity'):
        _validate_resume(tmp_path, [reads], config, flags)


def host_study(tmp_path):
    from test_comparison import study
    sheet = study(tmp_path)
    hashes = {part + '.bt2': hashlib.sha256(part.encode()).hexdigest()
              for part in ('1', '2', '3', '4', 'rev.1', 'rev.2')}
    for directory in sorted(tmp_path.glob('[0-9]*')):
        path = directory / 'analysis_metadata.json'
        metadata = json.loads(path.read_text())
        metadata['effective_config']['preprocessing']['host_index'] = f'/local/{directory.name}/host'
        metadata['effective_config']['preprocessing']['threads'] = int(directory.name) + 1
        metadata['tool_provenance']['host_removal'] = {
            'status': 'completed', 'schema_version': 1, 'input': 150, 'retained': 120, 'removed': 30,
            'count_unit': 'fragments', 'state': {'index': {'sha256': hashes}, 'preset': 'sensitive', 'seed': 0,
                'tool_versions': {'bowtie2': '2.5.4', 'samtools': '1.18'}, 'pair_policy': 'remove_pair_if_either_maps'}}
        path.write_text(json.dumps(metadata))
    return sheet


def test_cohort_accepts_identical_host_indexes_at_different_paths(tmp_path):
    from metaquest.core.comparison import run_comparison
    sheet = host_study(tmp_path)
    assert run_comparison(sheet, tmp_path / 'comparison').is_file()


@pytest.mark.parametrize('change', ['hash', 'policy', 'preset', 'version', 'disabled', 'missing', 'counts'])
def test_cohort_rejects_host_reference_methods_and_count_mismatches(tmp_path, change):
    from metaquest.core.comparison import ComparisonError, run_comparison
    sheet = host_study(tmp_path)
    path = tmp_path / '000/analysis_metadata.json'
    metadata = json.loads(path.read_text())
    host = metadata['tool_provenance']['host_removal']
    if change == 'hash':
        host['state']['index']['sha256']['1.bt2'] = hashlib.sha256(b'changed').hexdigest()
    elif change == 'policy':
        host['state']['pair_policy'] = 'both_mates_only'
    elif change == 'preset':
        host['state']['preset'] = 'very-sensitive'
    elif change == 'version':
        host['state']['tool_versions']['bowtie2'] = '2.5.3'
    elif change == 'disabled':
        metadata['effective_config']['preprocessing']['host_index'] = None
    elif change == 'missing':
        del metadata['tool_provenance']['host_removal']
    else:
        host['retained'] = 119
    path.write_text(json.dumps(metadata))
    with pytest.raises(ComparisonError):
        run_comparison(sheet, tmp_path / 'comparison')


def test_legacy_and_explicitly_disabled_host_profiles_remain_compatible(tmp_path):
    from test_comparison import study
    from metaquest.core.comparison import run_comparison
    sheet = study(tmp_path)
    path = tmp_path / '000/analysis_metadata.json'
    metadata = json.loads(path.read_text())
    metadata['effective_config']['preprocessing']['host_index'] = None
    path.write_text(json.dumps(metadata))
    assert run_comparison(sheet, tmp_path / 'comparison').is_file()


@pytest.mark.skipif(not NATIVE, reason='native Bowtie2/Samtools controls')
def test_host_native_evidence_tampering_invalidates_cache(tmp_path, native_reference):
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, host, microbe, 'paired')
    _, evidence = host_module().run_host_removal(reads, tmp_path / 'result', prefix, paired=True, threads=2)
    log = tmp_path / 'result' / evidence['attempt_dir'] / 'flagstat.json'
    log.write_text('{}')
    _, rerun = host_module().run_host_removal(reads, tmp_path / 'result', prefix, paired=True, threads=2, reuse=True)
    assert not rerun['reused'] and rerun['retained'] == 1
    assert log.read_text() == '{}'


@pytest.mark.skipif(not NATIVE, reason='native Bowtie2/Samtools controls')
@pytest.mark.parametrize('enabled', [False, True])
def test_interleaved_host_removal_and_eight_stage_reports(tmp_path, native_reference, enabled, monkeypatch):
    from metaquest.pipeline.context import PipelineContext
    from metaquest.pipeline.stages.preprocessing import run_preprocessing_stage
    from metaquest.reporting.stable_reporter import generate_stable_reports
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, host, microbe, 'paired')
    interleaved = tmp_path / 'interleaved.fastq'
    with reads[0].open() as first, reads[1].open() as second, interleaved.open('w') as out:
        for left, right in zip(FastqGeneralIterator(first), FastqGeneralIterator(second), strict=True):
            for name, sequence, quality in (left, right):
                out.write(f'@{name}\n{sequence}\n+\n{quality}\n')
    config = load_config()
    config = replace(config, preprocessing=replace(config.preprocessing, enabled=enabled, host_index=prefix, threads=2))
    output = tmp_path / 'result'
    output.mkdir()
    ctx = PipelineContext(config, [interleaved], output, read_mode='interleaved')
    run_preprocessing_stage(ctx)
    assert ctx.preprocessing['host_removal']['retained'] == 1
    assert len(ctx.analysis_input_files) == 2
    generate_stable_reports(ctx)
    summary = json.loads((output / 'analysis_summary.json').read_text())
    assert len(summary['stage_reports']) == 8
    for path in (output / 'report.html', output / 'reports/01_preprocessing.html'):
        text = path.read_text()
        assert 'Host removal' in text and 'Remove the whole pair' in text
        assert 'host_removal/summary.json' in text
    # Detailed provenance remains downloadable without burying the host counts.
    assert 'host removal / state / inputs' not in (output / 'reports/01_preprocessing.html').read_text()
    assert 'host removal / state / inputs' in (output / 'reports/01_preprocessing.tsv').read_text()
    if not enabled:
        from metaquest.pipeline.runner import _stage_output_hashes
        assert {'split_R1.fastq.gz', 'split_R2.fastq.gz'} <= set(_stage_output_hashes(ctx, 'Read Preprocessing'))
    ctx.resume = True
    monkeypatch.setattr('metaquest.io.utils.split_interleaved', lambda *a, **k: pytest.fail('unchanged split was regenerated'))
    run_preprocessing_stage(ctx)
    assert ctx.preprocessing['host_removal']['reused']


def test_host_html_escapes_reference_and_keeps_zero_fraction():
    from metaquest.reporting.host_removal import host_html
    evidence = {'input': 4, 'retained': 0, 'removed': 4, 'count_unit': 'fragments', 'retained_fraction': 0.0,
                'state': {'preset': 'sensitive', 'index': {'prefix': '<script>unsafe</script>'},
                          'tool_versions': {'bowtie2': '2.5.4', 'samtools': '1.18'}}, 'interpretation': []}
    text = host_html(evidence)
    assert '<script>' not in text and '&lt;script&gt;' in text
    assert '0.00%' in text


@pytest.mark.parametrize('tool,version', [('bowtie2', '2.4.5'), ('bowtie2', '3.0.0'), ('samtools', '1.17'), ('samtools', '2.0')])
def test_host_tool_versions_reject_unsupported_native_versions(monkeypatch, tool, version):
    versions = {'bowtie2': '2.5.4', 'samtools': '1.18', tool: version}
    def probe(command, **kwargs):
        name = command[0]
        text = f'bowtie2-align-s version {versions[name]}' if name == 'bowtie2' else f'samtools {versions[name]}'
        return subprocess.CompletedProcess(command, 0, text, '')
    monkeypatch.setattr(host_module().subprocess, 'run', probe)
    with pytest.raises(ConfigError, match='version'):
        host_module().tool_versions()


@pytest.mark.skipif(not NATIVE, reason='native Bowtie2/Samtools controls')
def test_empty_host_input_has_null_fraction_and_all_unmapped_reads_are_preserved(tmp_path, native_reference):
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, microbe, microbe, 'paired')
    output, summary = host_module().run_host_removal(reads, tmp_path / 'retained', prefix, paired=True, threads=2)
    assert summary['retained'] == summary['input'] == 4 and summary['removed'] == 0
    assert all(len(read_output(path)) == 4 for path in output)
    for path in reads:
        path.write_text('')
    output, summary = host_module().run_host_removal(reads, tmp_path / 'empty', prefix, paired=True, threads=2)
    assert summary['retained_fraction'] is None
    assert summary['input'] == summary['retained'] == summary['removed'] == 0
    assert all(read_output(path) == [] for path in output)


@pytest.mark.skipif(not NATIVE, reason='native Bowtie2/Samtools controls')
def test_shared_reference_sequence_is_removed_regardless_of_biological_label(tmp_path, native_reference):
    prefix, host, _ = native_reference
    reads = [write_reads(tmp_path / 'shared.fastq', [('microbial_label_shared_sequence', host[100:250])])]
    _, summary = host_module().run_host_removal(reads, tmp_path / 'result', prefix, paired=False, threads=2)
    assert summary['removed'] == 1 and summary['retained'] == 0
    assert any('Shared microbial/reference sequence' in line for line in summary['interpretation'])


@pytest.mark.skipif(not NATIVE, reason='native Bowtie2/Samtools controls')
def test_failed_mapper_keeps_previous_reads_and_failed_attempt_logs(tmp_path, native_reference, monkeypatch):
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, host, microbe, 'paired')
    output, _ = host_module().run_host_removal(reads, tmp_path / 'result', prefix, paired=True, threads=2)
    previous = [p.read_bytes() for p in output]
    tools = tmp_path / 'tools'
    tools.mkdir()
    mapper = tools / 'bowtie2'
    mapper.write_text(f'#!{sys.executable}\nimport sys\nif "--version" in sys.argv: print("bowtie2 version 2.5.4")\nelse: print("controlled mapper failure", file=sys.stderr); sys.exit(7)\n')
    mapper.chmod(0o755)
    monkeypatch.setenv('PATH', str(tools) + os.pathsep + os.environ['PATH'])
    with pytest.raises(MetaQuestError, match='Host-filter'):
        host_module().run_host_removal(reads, tmp_path / 'result', prefix, paired=True, threads=2)
    assert [p.read_bytes() for p in output] == previous
    logs = list((tmp_path / 'result/preprocessing/host_removal').glob('attempt-*/bowtie2.log'))
    assert any('controlled mapper failure' in p.read_text() for p in logs)
    assert not (tmp_path / 'result/preprocessing/host_removal/summary.json').exists()


@pytest.mark.skipif(not NATIVE or os.name != 'posix', reason='native POSIX output safety')
def test_host_output_fifo_fails_without_opening_it(tmp_path, native_reference):
    import metaquest
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, host, microbe, 'paired')
    output = tmp_path / 'result/preprocessing/host_removal'
    output.mkdir(parents=True)
    os.mkfifo(output / 'retained_R1.fastq.gz')
    env = os.environ.copy()
    env['PYTHONPATH'] = str(Path(metaquest.__file__).resolve().parents[1])
    code = ('from pathlib import Path;from metaquest.core.host_removal import run_host_removal;'
            'from metaquest.exceptions import ConfigError;import sys\n'
            'try: run_host_removal([Path(sys.argv[1]),Path(sys.argv[2])],Path(sys.argv[3]),Path(sys.argv[4]),paired=True,threads=2)\n'
            'except ConfigError as e: assert "regular file" in str(e)\n'
            'else: raise AssertionError("FIFO accepted")')
    subprocess.run([sys.executable, '-c', code, *map(str, reads), str(tmp_path / 'result'), str(prefix)],
                   env=env, check=True, timeout=10)


@pytest.mark.skipif(not NATIVE or os.name != 'posix', reason='native POSIX output safety')
def test_host_output_symlink_never_replaces_external_directory(tmp_path, native_reference):
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, host, microbe, 'paired')
    outside = tmp_path / 'outside'
    outside.mkdir()
    marker = outside / 'keep.txt'
    marker.write_text('keep')
    result = tmp_path / 'result'
    (result / 'preprocessing').mkdir(parents=True)
    (result / 'preprocessing/host_removal').symlink_to(outside, target_is_directory=True)
    with pytest.raises(ConfigError, match='output path'):
        host_module().run_host_removal(reads, result, prefix, paired=True, threads=2)
    assert marker.read_text() == 'keep' and list(outside.iterdir()) == [marker]


@pytest.mark.skipif(not NATIVE or os.name != 'posix', reason='native POSIX output safety')
def test_host_output_file_symlink_promotion_preserves_external_file(tmp_path, native_reference):
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, host, microbe, 'paired')
    result = tmp_path / 'result'
    output = result / 'preprocessing/host_removal'
    output.mkdir(parents=True)
    external = tmp_path / 'external.fastq.gz'
    external.write_bytes(b'preserve external bytes')
    target = output / 'retained_R1.fastq.gz'
    target.symlink_to(external)
    _, summary = host_module().run_host_removal(reads, result, prefix, paired=True, threads=2)
    assert external.read_bytes() == b'preserve external bytes'
    assert not target.is_symlink() and summary['retained'] == 1


@pytest.mark.skipif(not NATIVE or sys.platform != 'linux', reason='native Linux process groups')
@pytest.mark.parametrize('interrupt', [signal.SIGINT, signal.SIGTERM])
def test_host_streaming_interrupt_reaps_mapper_converter_and_counter(tmp_path, native_reference, interrupt):
    import metaquest
    from test_process_recovery import alive
    prefix, host, microbe = native_reference
    reads = control_reads(tmp_path, host, microbe, 'paired')
    tools = tmp_path / 'tools'
    tools.mkdir()
    pid_file = tmp_path / 'mapper.pid'
    mapper = tools / 'bowtie2'
    mapper.write_text(
        f'#!{sys.executable}\nimport os,subprocess,sys,time\nfrom pathlib import Path\n'
        'if "--version" in sys.argv: print("bowtie2 version 2.5.4"); sys.exit(0)\n'
        'child=subprocess.Popen([sys.executable,"-c","import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)"])\n'
        f'Path({str(pid_file)!r}).write_text(f"{{os.getpid()}} {{child.pid}}")\n'
        'print("@HD\\tVN:1.6\\n@SQ\\tSN:host\\tLN:12000",flush=True)\ntime.sleep(60)\n')
    mapper.chmod(0o755)
    driver = tmp_path / 'driver.py'
    driver.write_text('from metaquest.cli import main\nfrom metaquest.io import utils\n'
                      'utils.run_system_check=lambda *a,**k:None\nmain()\n')
    env = os.environ.copy()
    env['PYTHONPATH'] = str(Path(metaquest.__file__).resolve().parents[1])
    env['PATH'] = str(tools) + os.pathsep + env['PATH']
    stderr_path = tmp_path / 'driver.stderr'
    stderr_log = stderr_path.open('wb')
    process = subprocess.Popen([sys.executable, str(driver), '--no-color', 'run', '--host-index', str(prefix), '--min-sequences', '1',
                                '--paired', *map(str, reads), '--taxonomy-only', '--output', str(tmp_path / 'result')],
                               env=env, stdout=subprocess.DEVNULL, stderr=stderr_log, start_new_session=True)
    pids = []
    try:
        deadline = time.monotonic() + 10
        while not pid_file.exists() and process.poll() is None and time.monotonic() < deadline:
            time.sleep(.05)
        assert pid_file.exists(), 'mapper did not start: ' + stderr_path.read_text()
        pids = [int(x) for x in pid_file.read_text().split()]
        while time.monotonic() < deadline:
            children = Path(f'/proc/{process.pid}/task/{process.pid}/children').read_text().split()
            pids = [int(x) for x in pid_file.read_text().split()] + [int(x) for x in children]
            if len(set(pids)) >= 4:
                break
            time.sleep(.02)
        assert len(set(pids)) >= 4  # mapper, its child, Samtools converter and counter
        os.kill(process.pid, interrupt)
        assert process.wait(timeout=12) == 128 + interrupt
        deadline = time.monotonic() + 3
        while any(alive(pid) for pid in set(pids)) and time.monotonic() < deadline:
            time.sleep(.05)
        assert not any(alive(pid) for pid in set(pids))
    finally:
        for pid in pids:
            if alive(pid):
                os.kill(pid, signal.SIGKILL)
        if process.poll() is None:
            process.kill()
        process.wait()
        stderr_log.close()
