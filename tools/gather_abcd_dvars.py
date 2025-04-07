import polars as pl
import pyarrow.dataset as ds

dset = ds.dataset(
    "/Users/psadil/git/manuscripts/motion/data/dvars/dataset=abcd",
    format="ipc",
    partitioning=["sub", "ses", "src"],
)

pl.scan_pyarrow_dataset(dset).select("t", "DPD", "sub", "ses", "src").with_columns(
    task=pl.col("src").str.extract(r"task-(\w+)_run"),
    run=pl.col("src").str.extract(r"run-(\w+)_bold"),
    sub=pl.col("sub").str.extract(r"sub=(\w+)"),
    ses=pl.col("ses").str.extract(r"ses=(\w+)"),
).select("t", "DPD", "sub", "ses", "task", "run").collect().write_parquet(
    "/Users/psadil/git/manuscripts/motion/data/abcd_dvars.parquet"
)


dset = ds.dataset(
    "/Users/psadil/git/manuscripts/motion/data/dvars/dataset=ukb",
    format="ipc",
    partitioning=["sub", "ses", "src"],
)

pl.scan_pyarrow_dataset(dset).select("t", "DPD", "sub", "ses", "src").with_columns(
    task=pl.when(pl.col("src").str.contains("20227"))
    .then(pl.lit("rest"))
    .otherwise(pl.lit("faces/shape")),
    sub=pl.col("sub").str.extract(r"sub=(\w+)"),
    ses=pl.col("ses").str.extract(r"ses=(\w+)").cast(pl.Int8),
).filter(pl.col("src").str.contains("20227-filtered_func_data_clean$")).select(
    "t", "DPD", "sub", "ses", "task"
).collect().write_parquet("/Users/psadil/git/manuscripts/motion/data/ukb_dvars.parquet")
