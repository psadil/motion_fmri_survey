from pathlib import Path
import re

import polars as pl
import pandas as pd


def main(src: Path, dst_stem: Path, dst_root: Path) -> None:
    evs = []
    for sub in src.glob("*"):
        sub_id = re.findall(r"\w*\d{6,7}(?=_*)", str(sub.name))[0]
        for task in (sub / "MNINonLinear" / "Results").glob("tfMRI*"):
            task_id = re.findall(r"(?<=tfMRI_)\w+(?=_)", str(task.name))[0]
            dir_id = re.findall("(?<=_)(AP|PA|LR|RL)$", task.name)[0]
            for ev in (task / "EVs").rglob("*.txt"):
                try:
                    evs.append(
                        pl.from_pandas(
                            pd.read_csv(
                                ev,
                                sep=r"\s+",
                                names=["onset", "duration", "amplitude"],
                                dtype={
                                    "onset": float,
                                    "duration": float,
                                    "amplitude": float,
                                },
                            )
                        )
                        .drop("amplitude")
                        .with_columns(
                            sub=pl.lit(sub_id),
                            task=pl.lit(task_id),
                            dir=pl.lit(dir_id),
                            trial_type=pl.lit(ev.stem),
                        )
                    )
                except pl.exceptions.NoDataError:
                    continue

    out: pl.DataFrame = pl.concat(evs)
    out.write_parquet(dst_root / f"{dst_stem}.parquet")
    out.write_csv(dst_root / f"{dst_stem}.tsv", separator="\t")


dst_root = Path("derivatives")
for ds in ["HCPAgingRec", "HCPDevelopmentRec"]:
    main(
        src=Path(
            f"sourcedata/design/{ds}/dcs07/smart/data/{ds}/sourcedata/fmriresults01/"
        ),
        dst_stem=ds,
        dst_root=dst_root,
    )

ds = "human-connectome-project-openaccess"
main(
    src=Path(f"sourcedata/design/{ds}/dcs07/smart/data/{ds}/HCP"),
    dst_stem=ds,
    dst_root=dst_root,
)
