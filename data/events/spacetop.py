from pathlib import Path

import polars as pl


evs = []
for ev in Path("sourcedata/design/spacetop").rglob("*events.tsv"):
    try:
        tmp = pl.read_csv(ev, separator="\t", null_values="n/a")
        if "event_type" in tmp.columns:
            tmp = tmp.rename({"event_type": "trial_type"})
        evs.append(
            tmp.select("onset", "duration", "trial_type").with_columns(
                pl.col("onset").cast(pl.Float64),
                pl.col("duration").cast(pl.Float64),
                pl.col("trial_type").cast(pl.Utf8),
                src=pl.lit(str(ev)),
            )
        )
    except pl.exceptions.ColumnNotFoundError as e:
        print(ev)
        raise e

out: pl.DataFrame = pl.concat(evs)

out2 = out.with_columns(
    sub=pl.col("src").str.extract(r"sub-(\w+)_ses"),
    ses=pl.col("src").str.extract(r"ses-(\d+)_task"),
    task=pl.col("src").str.extract(r"task-(\w+)_acq"),
    run=pl.col("src").str.extract(r"run-(\d+)_"),
).drop("src")

out2.write_csv("derivatives/spacetop.tsv", separator="\t")
out2.write_parquet("derivatives/spacetop.parquet")
