import pandas as pd
from metagraph.client import GraphClient
import os
import sys
from pathlib import Path
from typing import Optional, Any, Iterable, List, Dict
from dataclasses import dataclass, fields
import attrs
import argparse

from pandas.core.sample import sample


@dataclass
class Params:
    port: int = 5555
    discovery_threshold: float = 0.2


def get_params_from_kwargs(kwargs):
    # Filter kwargs to only include fields that are in the Param dataclass
    param_fields = {f.name for f in fields(Params)}
    filtered_kwargs = {k: v for k, v in kwargs.items() if k in param_fields}
    return Params(**filtered_kwargs)


def get_params_from_args():
    """
    Get parameters from sys.argv
    :return: A Params dataclass instance
    """
    parser = argparse.ArgumentParser(description="Minknow server simulator")

    # Automatically add arguments based on the Config dataclass fields
    for field in fields(Params):
        field_type = field.type

        # For fields that are required, explicitly set 'required=True'
        if field.default is field.default_factory:  # Check if default is absent
            is_required = True
        else:
            is_required = False

        if field_type == List[str]:
            # For List fields, handle them as an optional list of values
            parser.add_argument(
                f'--{field.name}',
                type=str,
                nargs='+',  # Accepts multiple values
                required=is_required,  # Only set required if field has no default
                help=f"List of {field.name} values"
            )
        elif field_type == bool:
            # Boolean flags, no value expected
            parser.add_argument(
                f'--{field.name}',
                action='store_true',  # Presence means True
                required=is_required,  # Only set required if field has no default
                help=f"Set {field.name} to True"
            )
        else:
            # For other types, handle with default or required flag
            parser.add_argument(
                f'--{field.name}',
                type=field_type,
                default=getattr(Params, field.name, None),  # Use None if no default
                required=is_required,  # Only set required if field has no default
                help=f"Set the {field.name} (default: {getattr(Params, field.name, None)})"
            )

    parsed_args = parser.parse_args()
    params = Params(**vars(parsed_args))
    return params


@dataclass
class Alignment:
    ctg: str
    r_st: int
    r_en: int
    strand: int


def mgresults2alignmentdict(results: pd.DataFrame) -> Dict[int, Alignment]:
    alignments = {}
    for idx, row in results.iterrows():
        start = row['kmer_coords'][0].split('-')[1]
        end = row['kmer_coords'][-1].split('-')[1]
        ctg = row['sample'][:-1]
        strand = 1 if row['sample'][-1] == '+' else -1
        alignments[int(row['seq_description'])] = Alignment(ctg=ctg, r_st=start, r_en=end, strand=strand)
    return alignments


@attrs.define
class Result:
    """Result holder

    This should be progressively filled with data from the basecaller,
    barcoder, and then the aligner.

    :param channel: The channel that this read is being sequenced on
    :param read_id: The read ID assigned to this read by MinKNOW
    :param seq: The basecalled sequence for this read
    :param decision: The ``Decision`` that has been made, this will by used to determine the ``Action``
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

class Aligner:

    def __init__(self, debug_log: str | None = None, **kwargs):
        if debug_log:
            if debug_log == 'stdout':
                self.logfile = sys.stdout
            elif debug_log == 'stderr':
                self.logfile = sys.stderr
            else:
                self.logfile = open(debug_log, 'w')
        else:
            self.logfile = None
        self.kwargs = kwargs
        self.params = get_params_from_kwargs(kwargs)
        self.client = GraphClient('localhost', port=self.params.port, api_path='')

    def validate(self) -> None:
        pass

    @property
    def initialised(self) -> bool:
        return True

    def describe(self, regions: list, barcodes: dict) -> str:
        return "Running metagraph client. Hopefully, I will return a more meaningful description in the future"

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

        results = self.client.search(sequences, discovery_threshold=self.params.discovery_threshold, top_labels=1, query_coords=True)
        alignments = mgresults2alignmentdict(results)
        for i in range(len(metadata)):
            result = metadata[i]
            if i in alignments:
                result.alignment_data = [alignments[i]]
            else:
                result.alignment_data = []
            yield result

        for result in skipped:
            result.alignment_data = []
            yield result


    def disconnect(self):
        if self.logfile and self.logfile not in [sys.stdout, sys.stderr]:
            self.logfile.close()

