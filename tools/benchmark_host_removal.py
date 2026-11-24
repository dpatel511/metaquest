#!/usr/bin/env python3
"""Reproduce small host-filter integration/resource controls, not accuracy claims."""
import argparse
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import random
import resource
import statistics
import subprocess
import sys
import time

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import metaquest
from metaquest.config import __version__
from metaquest.core.host_removal import run_host_removal, tool_versions, index_identity
from metaquest.io.processes import managed_process


def reverse(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def prepare(inputs, protocol):
    inputs.mkdir()
    rng = random.Random(protocol['seed'])
    host = ''.join(rng.choices('ACGT', k=protocol['host_reference_bases']))
    other = ''.join(rng.choices('ACGT', k=protocol['unrelated_reference_bases']))
    (inputs / 'host.fa').write_text(f'>synthetic_reference\n{host}\n')
    truth = []
    with (inputs / 'R1.fastq').open('w') as r1, (inputs / 'R2.fastq').open('w') as r2:
        for group, count in protocol['groups'].items():
            first = host if group in ('both_reference_matching', 'first_mate_matching') else other
            second = host if group in ('both_reference_matching', 'second_mate_matching') else other
            for number in range(count):
                # Synthetic reference slices use 0-based Python offsets.
                position = number * 17 % (min(len(host), len(other)) - 450)
                name = f'{group}_{number}'
                length = protocol['read_length_bp']
                for mate, handle, sequence in ((1, r1, first[position:position+length]),
                                                (2, r2, reverse(second[position+300:position+300+length]))):
                    handle.write(f'@{name}/{mate}\n{sequence}\n+\n{"I" * length}\n')
                truth.append({'fragment_id': name, 'group': group, 'retained': group == 'neither_mate_matching'})
    assert len(truth) == protocol['paired_fragments']
    (inputs / 'truth.json').write_text(json.dumps(truth, indent=2))


def records(path):
    with (gzip.open if str(path).endswith('.gz') else open)(path, 'rt') as handle:
        yield from FastqGeneralIterator(handle)


def check_output(inputs, output):
    truth = {row['fragment_id'] for row in json.loads((inputs / 'truth.json').read_text()) if row['retained']}
    observed = set()
    expected = {}
    for left, right in zip(records(inputs / 'R1.fastq'), records(inputs / 'R2.fastq'), strict=True):
        key = left[0].removesuffix('/1')
        if key in truth:
            expected[key] = (left[1:], right[1:])
    for left, right in zip_longest(records(output[0]), records(output[1])):
        assert left is not None and right is not None, 'paired counts differ'
        key = left[0].removesuffix('/1')
        assert key == right[0].removesuffix('/2') and key not in observed
        assert key in expected and (left[1:], right[1:]) == expected[key]
        observed.add(key)
    assert observed == truth
    return len(observed)


def worker(kind, output, inputs, protocol):
    output.mkdir()
    reads = [inputs / f'R{i}.fastq' for i in (1, 2)]
    if kind == 'candidate':
        paths, summary = run_host_removal(reads, output, inputs / 'index', paired=True, threads=2)
        assert summary['input'] == protocol['paired_fragments']
        assert summary['state']['preset'] == 'sensitive' and summary['state']['seed'] == 0
    else:
        paths = [output / f'R{i}.fastq.gz' for i in (1, 2)]
        command = ['bowtie2', '-x', str(inputs / 'index'), *protocol['native_settings'], '-1', str(reads[0]), '-2', str(reads[1])]
        conversion = ['samtools', 'fastq', '-@', '0', '-c', '6', '-F', '2304', '-f', '12', '-N',
                      '-1', str(paths[0]), '-2', str(paths[1]), '-0', str(output / 'other.fastq'), '-s', str(output / 'singletons.fastq'), '-']
        with (output / 'bowtie2.log').open('w') as log, (output / 'fastq.log').open('w') as converted:
            with managed_process(command, stdout=subprocess.PIPE, stderr=log) as align:
                with managed_process(conversion, stdin=align.stdout, stdout=subprocess.DEVNULL, stderr=converted) as fastq:
                    align.stdout.close()
                    assert fastq.wait() == 0
                    assert align.wait() == 0
    retained = check_output(inputs, paths)
    (output / 'verified.json').write_text(json.dumps({'retained_fragments': retained, 'decoded_records_match_truth': True}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, default=Path(__file__).resolve().parents[1] / 'benchmarks/host-removal/protocol.json')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--worker', choices=('native', 'candidate'), help=argparse.SUPPRESS)
    parser.add_argument('--inputs', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    if args.worker:
        worker(args.worker, args.output, args.inputs, protocol)
        return
    assert not args.output.exists(), 'Use a fresh evidence directory'
    assert Path(metaquest.__file__).is_relative_to(sys.prefix), metaquest.__file__
    assert __version__ == protocol['candidate_version']
    output = args.output.resolve()
    output.mkdir()
    inputs = output / 'inputs'
    prepare(inputs, protocol)
    command = ['bowtie2-build', str(inputs / 'host.fa'), str(inputs / 'index')]
    with (output / 'index-build.log').open('w') as log:
        subprocess.run(['/usr/bin/time', '-v', '-o', str(output / 'index-build-resources.txt'), *command], stdout=log, stderr=subprocess.STDOUT, check=True)
    data = {'protocol_sha256': hashlib.sha256(args.protocol.read_bytes()).hexdigest(), 'protocol': protocol,
            'metaquest_version': __version__, 'module': metaquest.__file__, 'tool_versions': tool_versions(),
            'index': index_identity(inputs / 'index'), 'inputs_sha256': {name: hashlib.sha256((inputs / name).read_bytes()).hexdigest()
                                                                       for name in ('host.fa', 'R1.fastq', 'R2.fastq', 'truth.json')}, 'runs': []}
    (output / 'results.json').write_text(json.dumps(data, indent=2))
    for i, kind in enumerate(protocol['resource_order']):
        directory = output / f'{i+1}-{kind}'
        resource_log = output / f'{i+1}-{kind}-resources.txt'
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        started = time.perf_counter_ns()
        command = [sys.executable, str(Path(__file__).resolve()), '--protocol', str(args.protocol.resolve()), '--output', str(directory), '--worker', kind, '--inputs', str(inputs)]
        with (output / f'{i+1}-{kind}.log').open('w') as log:
            result = subprocess.run(['/usr/bin/time', '-v', '-o', str(resource_log), *command], stdout=log, stderr=subprocess.STDOUT)
        elapsed = (time.perf_counter_ns()-started)/1e9
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        row = {'kind': kind, 'exit': result.returncode, 'elapsed_seconds': elapsed,
               'process_tree_cpu_seconds': after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
               'gnu_time_maximum_rss_kib': int(next(line.split(':', 1)[1] for line in resource_log.read_text().splitlines() if 'Maximum resident set size' in line))}
        data['runs'].append(row)
        (output / 'results.json').write_text(json.dumps(data, indent=2))
        result.check_returncode()
        print(kind, f'{elapsed:.3f}s', flush=True)
    medians = {kind: statistics.median(row['elapsed_seconds'] for row in data['runs'] if row['kind'] == kind) for kind in ('native', 'candidate')}
    limit = max(2*medians['native'], medians['native']+1)
    data['budget'] = {'medians_seconds': medians, 'limit_seconds': limit, 'passed': medians['candidate'] <= limit}
    (output / 'results.json').write_text(json.dumps(data, indent=2))
    assert data['budget']['passed'], data['budget']
    print(json.dumps(data['budget'], indent=2))


if __name__ == '__main__':
    main()
