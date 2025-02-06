from pathlib import Path
import re
import polars as pl

from scipy import signal
import nitime.algorithms as tsa


def get_peaks(d: pl.DataFrame, tr: float, cols: list[str]) -> pl.DataFrame:
    peaks = []
    for param in ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]:
        for names, group in d.group_by(cols):
            if group.shape[0] < 100:
                continue
            print(names)
            try:
                x, y, _ = tsa.multi_taper_psd(
                    signal.detrend(group.select(param).to_series().to_numpy()),
                    Fs=1 / tr,
                    adaptive=True,
                    jackknife=False,
                    NW=8,
                    NFFT=512,
                )
                peaks.append(
                    group.select(cols)
                    .unique()
                    .with_columns(key=1, param=pl.lit(param))
                    .join(
                        pl.DataFrame({"freq": x, "pxx": y, "key": 1}), on="key"
                    )
                    .drop("key"),
                )
            except Exception as e:
                print(e)
    return pl.concat(peaks)


root = Path("rawdata")
derivatives = Path("derivatives")

(
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
    .pipe(get_peaks, tr=0.72, cols=["sub", "ses", "task", "ped"])
    .write_parquet(derivatives / "hcpa_spectrum.parquet")
)


(
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
    .pipe(get_peaks, tr=0.72, cols=["sub", "ses", "task", "ped"])
    .write_parquet(derivatives / "hcpd_spectrum.parquet")
)


(
    pl.scan_ipc(root / "dataset=hcpya")
    .collect()
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
    .pipe(get_peaks, tr=0.72, cols=["sub", "ses", "task", "ped"])
    .write_parquet(
        derivatives / "hcpya_spectrum.parquet",
    )
)

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

(
    pl.scan_ipc(root / "dataset=abcd")
    .with_columns(
        task=pl.col("src").str.extract(r"(rest|mid|sst|nback)"),
        run=pl.col("src").str.extract(r"run-(\d+)"),
    )
    .join(abcd_too_short, on=["sub", "ses", "task", "run"], how="anti")
    .drop("src")
    .collect()
    .pipe(get_peaks, tr=0.8, cols=["sub", "ses", "task", "run"])
    .write_parquet(derivatives / "abcd_spectrum.parquet")
)

(
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
    .pipe(get_peaks, tr=0.735, cols=["sub", "ses", "task"])
    .write_parquet(derivatives / "ukb_spectrum.parquet")
)


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
(
    tmp.pipe(get_peaks, tr=0.46, cols=["sub", "ses", "task", "run"])
    .fill_null(0)
    .write_parquet(derivatives / "spacetop_spectrum.parquet")
)
