from pathlib import Path
import re
import polars as pl
from scipy import signal
import numpy as np

from pymrimisc import motion
import nitime.algorithms as tsa


def find_peak_freq(
    series: pl.Series, tr: float, lower: float = 0.2, upper: float = 0.6
) -> float:
    try:
        x, y, _ = tsa.multi_taper_psd(
            signal.detrend(series.to_numpy()),
            Fs=1 / tr,
            adaptive=True,
            jackknife=False,
            NW=8,
            NFFT=512,
        )
        xi = np.where(np.logical_and(x >= lower, x <= upper))
        yy = y[xi]
        xx = x[xi]
        out = xx[np.where(yy == np.max(yy))][0]
    except:
        return np.nan

    return out


def add_filtered(
    d: pl.DataFrame,
    tr: float,
    cols: list[str],
    col_for_peak: list[str],
    ped_trans: str = "trans_y",
) -> pl.DataFrame:
    peaks = []
    for _, group in d.group_by(cols):
        peaks.append(
            group.sort("t")
            .group_by(cols)
            .agg(
                pl.col(ped_trans).map_batches(
                    lambda x: find_peak_freq(x, tr=tr, lower=0.1, upper=0.6),
                    return_dtype=pl.Float64,
                    returns_scalar=True,
                )
            )
        )

    tmp: pl.DataFrame = pl.concat(peaks)
    peaks_out = (
        tmp.filter(pl.col(ped_trans).is_not_nan())
        .group_by("sub", *col_for_peak)
        .agg(pl.col(ped_trans).median())
        .group_by(col_for_peak)
        .agg(
            med=pl.col(ped_trans).median(),
            lower=pl.col(ped_trans).quantile(0.5),
            upper=pl.col(ped_trans).quantile(0.75),
        )
    )

    d_long = d.drop("rmsd", strict=False).unpivot(index=[*cols, "t"])

    filtered = []
    for _, group in d_long.group_by(*cols, "variable"):
        peaks_to_use = peaks_out.join(
            group, on=col_for_peak, how="semi"
        ).to_dicts()[0]
        _filtered = motion.filter_band_stop(
            group.sort("t").select("value").to_series(),
            tr=tr,
            w0=(peaks_to_use["upper"] + peaks_to_use["lower"]) / 2,
            bw=peaks_to_use["upper"] - peaks_to_use["lower"],
        )
        filtered.append(
            pl.DataFrame(
                group.select(*cols, "variable", "t").with_columns(
                    value=_filtered
                )
            )
        )

    tmp2: pl.DataFrame = pl.concat(filtered)
    return (
        tmp2.pivot(on="variable", index=[*cols, "t"])
        .unpivot(index=[*cols, "t"])
        .with_columns(pl.col("variable") + pl.lit("_filtered"))
        .pivot(on="variable", index=[*cols, "t"])
        .join(d, on=[*cols, "t"], how="right")
    )


def _add_fd(d: pl.DataFrame, filtered: bool, cols: list[str]) -> pl.DataFrame:
    if filtered:
        d = d.drop(pl.selectors.ends_with("_x", "_y", "_z")).rename(
            lambda x: x.removesuffix("_filtered")
        )
        out_col = "framewise_displacement_filtered"
    else:
        out_col = "framewise_displacement"

    with_fd = []
    for _, group in d.group_by(cols):
        with_fd.append(motion.add_fd(group.sort("t")))

    return (
        pl.concat(with_fd)
        .select([*cols, "t", "framewise_displacement"])
        .rename({"framewise_displacement": out_col})
    )


def add_fd(d: pl.DataFrame, cols: list[str]) -> pl.DataFrame:
    d_filtered_fd = _add_fd(d, filtered=True, cols=cols)
    d_fd = _add_fd(d, filtered=False, cols=cols).join(
        d_filtered_fd, on=[*cols, "t"]
    )

    return d.join(d_fd, on=[*cols, "t"])


root = Path("rawdata")
derivatives = Path("derivatives")


hcpaging = (
    pl.scan_ipc(root / "dataset=hcpa")
    .with_columns(
        ped=pl.col("src").str.extract(r"(AP|PA)"),
        task=pl.col("src").str.extract(
            r"(REST1|REST2|CARIT|FACENAME|VISMOTOR)"
        ),
        sub=pl.col("sub").cast(pl.Utf8),
        ses=pl.col("ses").cast(pl.Utf8),
    )
    .drop("src", pl.selectors.contains("deriv"))
    .collect()
    .pipe(
        add_filtered,
        tr=0.72,
        cols=["sub", "ses", "task", "ped"],
        col_for_peak=["task"],
    )
    .pipe(add_fd, cols=["sub", "ses", "task", "ped"])
)

hcpaging.write_parquet(derivatives / "HCPAgingRec.parquet")
hcpaging.write_csv(derivatives / "HCPAgingRec.tsv", separator="\t")
del hcpaging

hcpd = (
    pl.scan_ipc(root / "dataset=hcpd")
    .with_columns(
        ped=pl.col("src").str.extract(r"(AP|PA)"),
        task=pl.col("src").str.extract(
            r"(REST1a|REST1b|REST1|REST2a|REST2b|REST2|CARIT|EMOTION|GUESSING)"
        ),
        sub=pl.col("sub").cast(pl.Utf8),
        ses=pl.col("ses").cast(pl.Utf8),
    )
    .drop("src", pl.selectors.contains("deriv"))
    .collect()
    .fill_null(0)
    .pipe(
        add_filtered,
        tr=0.72,
        cols=["sub", "ses", "task", "ped"],
        col_for_peak=["task"],
    )
    .pipe(add_fd, cols=["sub", "ses", "task", "ped"])
)

hcpd.write_parquet(derivatives / "HCPDevelopmentRec.parquet")
hcpd.write_csv(derivatives / "HCPDevelopmentRec.tsv", separator="\t")
del hcpd

hcpya = (
    pl.scan_ipc(root / "dataset=hcpya")
    .with_columns(
        ped=pl.col("src").str.extract(r"(LR|RL)"),
        task=pl.col("src").str.extract(
            r"(REST1|REST2|RELATIONAL|EMOTION|MOTOR|GAMBLING|LANGUAGE|SOCIAL|WM)"
        ),
        sub=pl.col("sub").cast(pl.Utf8),
        ses=pl.col("ses").cast(pl.Utf8),
    )
    .drop("src", pl.selectors.contains("deriv"))
    .fill_null(0)
    .collect()
    .pipe(
        add_filtered,
        tr=0.72,
        cols=["sub", "ses", "task", "ped"],
        col_for_peak=["task"],
        ped_trans="trans_x",
    )
    .pipe(add_fd, cols=["sub", "ses", "task", "ped"])
)
hcpya.write_csv(
    derivatives / "human-connectome-project-openaccess.tsv", separator="\t"
)

hcpya.write_parquet(
    derivatives / "human-connectome-project-openaccess.parquet",
)
del hcpya

# some abcd scans are just way too short
abcd_too_short = (
    pl.scan_ipc(root / "dataset=abcd")
    .with_columns(
        task=pl.col("src").str.extract(r"(rest|mid|sst|nback)"),
        run=pl.col("src").str.extract(r"run-(\d+)"),
    )
    .drop("src")
    .group_by(["sub", "ses", "task", "run"])
    .len()
    .filter(pl.col("len") < 10)
)

abcd = (
    pl.scan_ipc(root / "dataset=abcd")
    .with_columns(
        task=pl.col("src").str.extract(r"(rest|mid|sst|nback)"),
        run=pl.col("src").str.extract(r"run-(\d+)"),
    )
    .join(abcd_too_short, on=["sub", "ses", "task", "run"], how="anti")
    .drop("src")
    .collect()
    .pipe(
        add_filtered,
        tr=0.8,
        cols=["sub", "ses", "task", "run"],
        col_for_peak=["task", "ses"],
    )
    .pipe(add_fd, cols=["sub", "ses", "task", "run"])
)
abcd.write_parquet(derivatives / "abcd.parquet")
abcd.write_csv(derivatives / "abcd.tsv", separator="\t")
del abcd

ukb = (
    pl.scan_ipc(root / "dataset=ukb")
    .with_columns(
        task=pl.when(pl.col("src") == 20227)
        .then(pl.lit("rest"))
        .otherwise(pl.lit("faces/shapes")),
        sub=pl.col("sub").cast(pl.Utf8),
        ses=pl.col("ses").cast(pl.Utf8),
    )
    .drop("src")
    .collect()
    .pipe(
        add_filtered,
        tr=0.735,
        cols=["sub", "ses", "task"],
        col_for_peak=["task"],
    )
    .pipe(add_fd, cols=["sub", "ses", "task"])
)
ukb.write_parquet(derivatives / "ukb.parquet")
ukb.write_csv(derivatives / "ukb.tsv", separator="\t")
del ukb

# spacetop
tsvs = []
for subdir in (
    Path("sourcedata")
    / "dcs04/smart/data/spatialtopology/fmriprep/results/fmriprep"
).glob("sub*"):
    print(subdir)
    sub = re.findall(r"(?<=sub-)\d{4}", subdir.name)
    for sesdir in subdir.glob("ses*"):
        ses = re.findall(r"(?<=ses-)\d{2}", sesdir.name)
        for tsv in sesdir.rglob("*tsv"):
            task = re.findall(r"(?<=task-)[a-zA-Z]+", tsv.name)
            run = re.findall(r"(?<=run-)\d+", tsv.name)
            tsvs.append(
                pl.read_csv(tsv, separator="\t", null_values="n/a")
                .select(pl.selectors.starts_with("rot", "trans", "rmsd"))
                .drop(
                    pl.selectors.ends_with("power2"),
                    pl.selectors.contains("deriv"),
                )
                .with_row_index("t")
                .with_columns(
                    sub=pl.lit(sub[0]),
                    ses=pl.lit(ses[0]),
                    task=pl.lit(task[0]),
                    run=pl.lit(run[0]),
                )
            )

tmp: pl.DataFrame = pl.concat(tsvs)
spacetop = (
    tmp.pipe(
        add_filtered,
        tr=0.46,
        cols=["sub", "ses", "task", "run"],
        col_for_peak=["task"],
    )
    .pipe(add_fd, cols=["sub", "ses", "task", "run"])
    .fill_null(0)
)

spacetop.write_parquet(derivatives / "spacetop.parquet")
spacetop.write_csv(derivatives / "spacetop.tsv", separator="\t")
