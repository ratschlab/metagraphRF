"""readfish aligner plugin backed by an in-process MetaGraph index (pymetagraph).

A read counts as mapped if MetaGraph reports anything for it in the selected
mode, and every reported label becomes one alignment with that label as the
contig, so readfish's targets and multi_*/single_* conditions apply to labels:

- ``query``: label query (``metagraph query``), the labels sharing at least
  ``discovery_fraction`` of the read's k-mers.
- ``map``: seed-only label mapping (``metagraph align --align-map-only``),
  the labels whose seeds cover the most of the read.
- ``align``: sequence-to-graph alignment (``metagraph align``). With
  ``labeled = false`` alignments carry no label and are reported with the
  contig ``*``, as the March 2026 runs did through ``metagraph server_query``.

Reference coordinates are not known (label annotations have none), so every
alignment spans the read's aligned or mapped query range on the target.
"""
import os
import sys
from dataclasses import dataclass
from typing import Any, Iterable, Optional

import attrs

# settings each mode accepts, besides graph, annotation, threads, mode, debug_log
MODE_SETTINGS = {
    "query": {"num_top_labels": int, "discovery_fraction": float, "presence_fraction": float},
    "map": {"min_exact_match": float, "min_seed_length": int},
    "align": {"labeled": bool, "min_exact_match": float, "max_alternative_alignments": int,
              "max_num_nodes_per_seq_char": float, "min_seed_length": int},
}
# defaults where they differ from pymetagraph's (which follow the metagraph CLI)
MODE_DEFAULTS = {
    "query": {"num_top_labels": 1, "discovery_fraction": 0.0},
}
UNLABELED = "*"


@dataclass
class Alignment:
    ctg: str
    r_st: int
    r_en: int
    strand: int


@attrs.define
class Result:
    """Result holder

    This should be progressively filled with data from the basecaller,
    barcoder, and then the aligner.

    :param channel: The channel that this read is being sequenced on
    :param read_id: The read ID assigned to this read by MinKNOW
    :param seq: The basecalled sequence for this read
    :param barcode: The barcode that has been assigned to this read
    :param basecall_data: Any extra data that the basecaller may want to send to the aligner
    :param alignment_data: Any extra alignment data
    """

    channel: int
    read_id: str
    seq: str
    barcode: Optional[str] = attrs.field(default=None)
    basecall_data: Optional[Any] = attrs.field(default=None)
    alignment_data: Optional[list[Alignment]] = attrs.field(default=None)


def _to_bool(value) -> bool:
    if isinstance(value, str):
        if value.lower() in ("true", "1", "yes"):
            return True
        if value.lower() in ("false", "0", "no"):
            return False
        raise ValueError(f"not a boolean: {value!r}")
    return bool(value)


class Aligner:
    """
    :param graph: Path to the graph (``.dbg``).
    :param annotation: Path to the annotation (e.g. ``.column.annodbg``).
    :param threads: Threads used by pymetagraph.
    :param mode: ``query``, ``map`` or ``align``; the remaining settings
        depend on it (see ``MODE_SETTINGS``) and are passed to pymetagraph.
    :param debug_log: File (or ``stdout``/``stderr``) to log MAP/NOMAP/SKIP per read to.
    """

    def __init__(self, graph: str, annotation: str, threads: int = 1, mode: str = "query",
                 debug_log: Optional[str] = None, **kwargs):
        if mode not in MODE_SETTINGS:
            raise ValueError(f"Unknown mode '{mode}', expected one of {', '.join(MODE_SETTINGS)}")
        unknown = set(kwargs) - set(MODE_SETTINGS[mode])
        if unknown:
            raise ValueError(f"Settings not supported in mode '{mode}': {', '.join(sorted(unknown))}"
                             f" (supported: {', '.join(MODE_SETTINGS[mode])})")
        self.mode = mode
        self.params = dict(MODE_DEFAULTS.get(mode, {}))
        for key, value in kwargs.items():
            convert = _to_bool if MODE_SETTINGS[mode][key] is bool else MODE_SETTINGS[mode][key]
            self.params[key] = convert(value)

        if debug_log == "stdout":
            self.logfile = sys.stdout
        elif debug_log == "stderr":
            self.logfile = sys.stderr
        elif debug_log:
            self.logfile = open(debug_log, "w")
        else:
            self.logfile = open(os.devnull, "w")

        import pymetagraph
        self.graph_path = graph
        self.graph = pymetagraph.AnnotatedGraph(graph, annotation, threads=int(threads))

    def validate(self) -> None:
        pass

    @property
    def initialised(self) -> bool:
        return True

    def describe(self, regions: list, barcodes: dict) -> str:
        return (f"metagraphRF: {self.graph_path} (k={self.graph.k}, {len(self.graph.labels)} labels), "
                f"mode {self.mode} {self.params}")

    def _alignments(self, sequences: list) -> list:
        """Per sequence, a list of Alignment, empty if nothing was found."""
        if self.mode == "query":
            return [[Alignment(label, 0, len(seq), 1) for label, _ in hits]
                    for seq, hits in zip(sequences, self.graph.query(sequences, **self.params))]
        if self.mode == "map":
            return [[] if m is None else
                    [Alignment(label, m.query_begin, m.query_end, -1 if m.orientation else 1)
                     for label in m.labels]
                    for m in self.graph.map(sequences, **self.params)]
        alignments = []
        for paths in self.graph.align(sequences, **self.params):
            alignments.append([Alignment(label, a.query_begin, a.query_end, -1 if a.orientation else 1)
                               for a in paths for label in (a.labels or [UNLABELED])])
        return alignments

    def map_reads(self, calls: Iterable[Result]) -> Iterable[Result]:
        skipped = []
        metadata = []
        sequences = []
        for result in calls:
            if result.seq:
                metadata.append(result)
                sequences.append(result.seq)
            else:
                skipped.append(result)

        for result, alignments in zip(metadata, self._alignments(sequences)):
            result.alignment_data = alignments
            print(f"{'MAP' if alignments else 'NOMAP'}\t{result.read_id}\t{result.seq}", file=self.logfile)
            yield result

        for result in skipped:
            result.alignment_data = []
            print(f"SKIP\t{result.read_id}\tN", file=self.logfile)
            yield result

    def disconnect(self):
        if self.logfile and self.logfile not in [sys.stdout, sys.stderr]:
            self.logfile.close()
