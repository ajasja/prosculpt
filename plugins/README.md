# Backbone filter plugins

Plugins are the backbone filters prosculpt applies to RFdiffusion's output,
right after `run_rfdiff` and before anything is rechained and passed on to
ProteinMPNN. Rejecting a bad backbone here saves the whole MPNN + structure
prediction + scoring pass it would otherwise trigger.

Each plugin is a plain Python file exposing a single function. Prosculpt loads
it by path, so a plugin can live in this folder, next to a job's own config, or
anywhere else on disk.

## The contract

```python
def filter_backbone(pdb: str | Path, kwargs: dict) -> bool:
    ...
```

- `pdb` - the backbone PDB to judge; the function is called once per file.
- `kwargs` - this filter's own settings, straight from the job config.
- Return `True` to keep the backbone, `False` to reject it.

That's all that's required - no base class, no registration, no imports from
prosculpt. A file that doesn't define `filter_backbone` is rejected when the
plugin is loaded (an `AttributeError`), before any design work starts, as is a
`filter_script` path that doesn't exist (a `FileNotFoundError`).

PyRosetta is already initialised by the time a filter runs (`prosculpt_run.py`
calls `init("-mute all")` on import), so a plugin can use it directly without
initialising it itself.

## Configuring filters

In a job's YAML config:

```yaml
rfdiff_backbone_filters:
  - filter_name: num_neighbors              # label, used in the log
    filter_script: plugins/num_neighbors.py
    delete_failed: false                    # false: move rejects aside. true: delete them
    min_num_neighbors: 7                    # everything below is passed to the plugin as kwargs
    max_num_neighbors: 50
    neighbor_distance: 8.0
```

`filter_name`, `filter_script` and `delete_failed` are consumed by prosculpt
itself; **every other key becomes an entry in the plugin's `kwargs`**. There are
no other reserved names and no default arguments - each plugin decides what it
needs and what to do when a key is missing (`kwargs.get("neighbor_distance", 8)`
for an optional one, `kwargs["base_pdb"]` for a mandatory one).

Filters run in the order listed, and a backbone is rejected as soon as one of
them returns `False` - the remaining filters aren't run on it.

`filter_script` is resolved as: an absolute path is used as-is; a relative path
is checked against the current working directory first, then against
prosculpt's own installation directory (which is why `plugins/SS_filter.py`
works no matter where the job runs from).

## What happens to rejected backbones

All of it inside RFdiffusion's own output directory:

- **Passed** - left where it is.
- **Rejected**, `delete_failed: false` - the `.pdb` and its `.trb` are moved
  into a `failed_filters/` subdirectory, so you can still see what was thrown
  out and why (the log names the filter that rejected it).
- **Rejected**, `delete_failed: true` - both files are deleted.
- **Plugin raised an exception** - logged as an error and the backbone is
  **kept**. A broken filter can't silently discard a run.

If *every* backbone is rejected, the pipeline aborts with an error rather than
carrying on with nothing to design - usually a sign the thresholds are too
strict for what RFdiffusion is producing.

## The included plugins

All three use PyRosetta.

| Plugin | Passes when | Config keys |
| --- | --- | --- |
| `SS_filter.py` | Chain A contains more than 2 DSSP secondary structure elements (helices + strands) | none |
| `num_neighbors.py` | The number of non-chain-B residues within `neighbor_distance` of chain A's first helix falls in `[min_num_neighbors, max_num_neighbors]` | `min_num_neighbors` (default 0), `max_num_neighbors` (default unbounded), `neighbor_distance` (default 8) |
| `clash_check.py` | The design has no inter-chain clashes once assembled against a reference structure | `base_pdb` (required), `residues_to_align` (default 85), `skip_last_residues` (default 200) |

`clash_check.py` is written for a specific assembly geometry (it superimposes
two copies of the design onto a base structure by chain A's C-terminal
residues) - treat it as a worked example of a non-trivial filter rather than
something to point at an arbitrary target.

Usage examples live in [`Examples/filters.yaml`](../Examples/filters.yaml), and
`config/run.yaml` carries a commented-out block for each.

## Writing your own

```python
from pathlib import Path
import logging

logger = logging.getLogger(__name__)   # goes to the job's own log


def filter_backbone(pdb: str | Path, kwargs: dict) -> bool:
    """One-line description of what passes."""
    threshold = kwargs.get("my_threshold", 10)

    value = ...  # measure something about `pdb`
    logger.info(f"{Path(pdb).name}: value={value}, threshold={threshold}")

    return value <= threshold
```

Then point a config at it:

```yaml
rfdiff_backbone_filters:
  - filter_name: my_filter
    filter_script: plugins/my_filter.py
    delete_failed: false
    my_threshold: 12
```

Worth keeping in mind:

- The filter runs once per backbone, in the same process as the rest of the
  pipeline - an expensive filter is paid for on every structure RFdiffusion
  produced.
- Log what the measured value actually was, not just pass/fail. It's the only
  way to tell "the filter works and these backbones are bad" from "the
  threshold is wrong" when a run comes back empty.
- The dashboard's **Run job** tab can upload a filter script and set its
  arguments for you (RFdiffusion backbone filters module), which is a quick way
  to try one out without hand-editing YAML.
