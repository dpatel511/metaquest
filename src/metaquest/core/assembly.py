"""
Metagenomic Assembly Module — MEGAHIT orchestration.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List

from Bio import SeqIO

from ..exceptions import AssemblyError
from ..io.output_formatter import get_formatter


def assemble_reads_to_fasta(
    reads: List[str],
    output_dir: Path,
    formatter=None,
    *,
    threads: int = 8,
) -> Path:
    """
    Perform metagenomic assembly using MEGAHIT.

    Args:
        reads: List of input read file paths.
        output_dir: Output directory.
        formatter: Output formatter (uses global if None).
        threads: Number of threads.

    Returns:
        Path to assembled contigs FASTA file.

    Raises:
        AssemblyError: If assembly fails or produces invalid output.
    """
    if formatter is None:
        formatter = get_formatter()

    megahit_outdir = output_dir / "megahit_assembly"
    final_contigs = megahit_outdir / "final.contigs.fa"

    cmd = ["megahit", "-t", str(threads)]
    if len(reads) == 1:
        cmd.extend(["-r", str(reads[0])])
    elif len(reads) == 2:
        cmd.extend(["-1", str(reads[0]), "-2", str(reads[1])])

    cmd.extend(["-o", str(megahit_outdir), "--force"])

    # MEGAHIT creates POSIX FIFOs for compressed reads. Windows-mounted output
    # directories under WSL cannot provide them; keep scratch on the system temp filesystem.
    with TemporaryDirectory(prefix="metaquest-megahit-") as temporary:
        cmd.extend(["--tmp-dir", temporary])
        returncode, _, stderr = formatter.run_subprocess(
            cmd, "Running MEGAHIT assembler", show_command=True
        )

    if returncode != 0:
        raise AssemblyError(f"MEGAHIT failed (exit {returncode}). Check log: {megahit_outdir}/log")

    if not final_contigs.exists():
        raise AssemblyError(f"MEGAHIT did not produce output: {final_contigs}")

    if final_contigs.stat().st_size == 0:
        raise AssemblyError("MEGAHIT produced empty contigs file")

    contig_count = sum(1 for _ in SeqIO.parse(final_contigs, "fasta"))
    if contig_count == 0:
        raise AssemblyError("Assembly contains no valid contigs")

    formatter.debug(f"Assembly: {contig_count} contigs produced")
    return final_contigs
