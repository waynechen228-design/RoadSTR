# Example data

Three injections from one 3730xl run (CLA Identifiler Plus, POP-7, GS500 LIZ, dye set G5):

| File | Content |
|---|---|
| `allelic_ladder.fsa` | allelic ladder |
| `positive_control_9947A.fsa` | positive control DNA 9947A (commercial reference DNA) |
| `negative_control.fsa` | amplification negative control |

No donor or patient material is included. Try:

```bash
str-typer validate-kit --kit cla_identifiler_plus --ladder examples/allelic_ladder.fsa
str-typer run --kit cla_identifiler_plus --ladder examples/allelic_ladder.fsa \
    --positive examples/positive_control_9947A.fsa --negative examples/negative_control.fsa \
    --samples examples/positive_control_9947A.fsa --plots --out example_results
```
