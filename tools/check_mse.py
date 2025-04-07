import argparse
import logging
from pathlib import Path

import nibabel as nb
import polars as pl

logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    level=logging.INFO,
)


def get_measures(src: Path, dst: Path) -> None:
    if dst.exists():
        logging.info(f"Skipping {dst}")
    elif not (parent := dst.parent).exists():
        parent.mkdir(parents=True)

    logging.info(f"Creating {dst}")
    nii = nb.loadsave.load(src)

    ts = (
        pl.from_numpy(nii.get_fdata())
        .rename(lambda col: col.removeprefix("column_"))
        .with_row_index(name="t", offset=0)
        .unpivot(index="t")
    )

    avg = ts.group_by("variable").agg(
        med=pl.col("value").median(), avg=pl.col("value").mean()
    )

    ts.join(avg, on=["variable"]).group_by("t").agg(
        mse=(pl.col("value") - pl.col("avg")).pow(2).mean(),
        mad=(pl.col("value") - pl.col("med")).abs().median(),
    ).sort("t").write_ipc(dst, compression="zstd")


def main(i: int, src: Path, dst_root: Path):
    subs = [sub for sub in src.read_text().splitlines()]
    to_process = Path(subs[i])
    logging.info(f"Processing {to_process}")
    for run in ["1", "2"]:
        for ped in ["LR", "RL"]:
            parent = to_process / "MNINonLinear" / "Results" / f"rfMRI_REST{run}_{ped}"
            img = parent / f"rfMRI_REST{run}_{ped}_Atlas_MSMAll.dtseries.nii"
            if not img.exists():
                continue
            dst = (
                dst_root
                / f"sub={to_process.name}"
                / f"src={img.stem.removesuffix('.dtseries')}"
                / "outliers.arrow"
            )
            get_measures(img, dst)
    logging.info("finished")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("i", type=int)
    parser.add_argument("src", type=Path)
    parser.add_argument("dst_root", type=Path)

    args = parser.parse_args()
    main(i=args.i, src=args.src, dst_root=args.dst_root)
