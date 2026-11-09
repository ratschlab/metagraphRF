# metagraphRF

A [readfish](https://github.com/LooseLab/readfish) aligner plugin that makes adaptive-sampling
decisions with a MetaGraph index loaded into the readfish process through
[pymetagraph](https://github.com/ratschlab/metagraph/tree/master/metagraph/api/python/pymetagraph).
No `metagraph server_query` is needed.

## Install

```
pip install .
```

This builds pymetagraph, and so MetaGraph, from source (see its README for the build
requirements). To use a pymetagraph that is already installed, run `pip install --no-deps .`

## Configuration

```toml
[mapper_settings.metagraphRF]
graph = "/path/to/graph.dbg"
annotation = "/path/to/graph.column.annodbg"
threads = 8
mode = "query"              # query | map | align
discovery_fraction = 0.03   # mode-specific settings, see below
```

Every label MetaGraph reports for a read becomes one alignment with that label as the contig,
so readfish's `targets` and `single_*`/`multi_*` conditions apply to labels. A read with nothing
reported is `no_map`.

| mode | MetaGraph operation | settings (defaults) |
|---|---|---|
| `query` | label query (`metagraph query`): labels sharing at least `discovery_fraction` of the read's k-mers | `num_top_labels` (1), `discovery_fraction` (0.0), `presence_fraction` (0.0) |
| `map` | seed-only label mapping (`metagraph align -a --align-map-only`): labels whose seeds cover the most of the read | `min_exact_match` (0.7), `min_seed_length` (min(k, 19)), `max_seed_length` |
| `align` | alignment (`metagraph align`) | `labeled` (true), `min_exact_match` (0.7), `max_alternative_alignments` (1), `min_seed_length` (min(k, 19)), `max_seed_length`, `max_nodes_per_seq_char` (5.0) |

With `labeled = false`, `align` behaves like the `/align` endpoint of `metagraph server_query`
that earlier versions of this plugin used: alignments have no labels and are reported with the
contig `*`, so they only count as mapped. Settings a mode does not support raise an error.

`debug_log = "<file>"` (or `stdout`/`stderr`) logs `MAP`/`NOMAP`/`SKIP` per read.
