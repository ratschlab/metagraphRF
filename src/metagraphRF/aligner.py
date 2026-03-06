from metagraph.client import GraphClient
import sys
from typing import Optional, Any, Iterable, List, Dict
from dataclasses import dataclass, fields
import attrs
import argparse

@dataclass
class Params:
    port: int = 5555
    min_exact_match: float = 0.5


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
    parser = argparse.ArgumentParser(description="Readfish-compatible Metagraph API")

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
        try:
            if self.client.ready():
                print("Metagraph client is ready")
                return
        except ConnectionRefusedError:
            print("Metagraph server is not running.")
            print("Please run `metagraph server_query -i <graph.dbg> -a <grapg.column.annodbg>`")
        finally:
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

        results = self.client.align(sequences,min_exact_match=self.params.min_exact_match)
        aligned_sequences = set(results['seq_descritpion'])

        for i in range(len(metadata)):
            result = metadata[i]
            if i in aligned_sequences:
                # we send a dummy alignment
                result.alignment_data = [Alignment('Found', 0, len(result.seq), 1)]
            else:
                result.alignment_data = []
            yield result

        for result in skipped:
            result.alignment_data = []
            yield result


    def disconnect(self):
        if self.logfile and self.logfile not in [sys.stdout, sys.stderr]:
            self.logfile.close()


def test_aligner():
    sequences = [x.strip() for x in """
            TCAGCGCTTCAGGTGCAATGGCTACCGGCTGGGTTCGCGATGGATCCACCTGGTACTACCTCACGCCCTCGGT.
            CTAATAACTGAAATAGAGGCTGTTATAAATGCAGATAAATTATGCTGAATGAAGAACGATATGAATATAAAGG.
            AAATGAGCCAATTGTCCAAGCAAGCAAATTGAAAATATTACAAAGGCCCAAGAGGCACTTGCTAAAGAAGTAG.
            GCCACGCCGTCCAATGAAGTTTTCACGATACAAACGGCGTTCACAAAATCGAGCTGATCGTCATAGCCGACCG.
            ACAATTTGTTGAAATGATGGTGCAGGATTGAAAGAAAGCCATCCTCATGATGTGCAATCACAAACGAGGCTCT.
            AAACGCCCAACTGCAAGCCCAAAGCGCATACGAGACGGATAATGATCGTCTCGTAACCGACGGCTCGAGGAAA.
            ACCCGATTTCTTCTCACTCTGTAGCCATAATAGCCAAAAAACTGATGAATCTTTCGATTCATCAGTTTTAAAG.
            ACCATAAGGTACCTCCCTATTAAATATCAGTCATCATCTTCTTCAGCATGATGGCTGAAACTTCAAGGCTGTA
        """.split('.\n')]

    aligner = Aligner()
    aligner.validate()
    df = aligner.client.align(sequences)
    print(df)


if __name__ == "__main__":
    test_aligner()