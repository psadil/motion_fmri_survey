from pathlib import Path
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
                    group
                    .select(cols)
                    .unique()
                    .with_columns(key=1, param=pl.lit(param))
                    .join(pl.DataFrame({"freq": x, "pxx": y, "key": 1}), on="key")
                    .drop("key"),
                )
            except Exception as e:
                print(e)
    return pl.concat(peaks)


(
    pl
    .scan_csv(
        "../abcc-4-0-0/abcc-xcp_d_v0.13.0/**/*run*motion.tsv",
        separator="\t",
        include_file_paths="src",
    )
    .with_columns(
        sub=pl.col("src").str.extract(r"sub-(\w+)"),
        ses=pl.col("src").str.extract(r"ses-(\w+)"),
        task=pl.col("src").str.extract(r"(rest|mid|sst|nback)"),
        run=pl.col("src").str.extract(r"run-(\d+)").cast(pl.Int8),
    )
    .with_columns(
        ses=pl.col("ses").replace_strict({
            "00A": "Baseline",
            "02A": "Year2",
            "04A": "Year4",
            "06A": "Year6",
        })
    )
    .drop("src", pl.selectors.contains("filter"), pl.selectors.contains("frame"))
    .collect()
    .pipe(get_peaks, tr=0.8, cols=["sub", "ses", "task", "run"])
    .write_parquet(Path("derivatives") / "abcc_spectrum.parquet")
)
