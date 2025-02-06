import argparse
from pathlib import Path
import re
import logging

import pydantic
import polars as pl
import numpy as np


logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s", level=logging.INFO
)


class MotionProcesser(pydantic.BaseModel):
    src: Path
    dst_root: Path

    @property
    def sub_id(self) -> str:
        sub = re.findall(r"(?<=sub-)[a-zA-Z0-9]+", self.src.name)
        if not (len(sub) == 1):
            msg = f"Unexpected length of sub_id for {self.src}"
            raise RuntimeError(msg)

        return sub[0]

    @property
    def ses_id(self) -> str:
        ses = re.findall(r"(?<=ses-)[a-zA-Z0-9]+", self.src.name)
        if not (len(ses) == 1):
            msg = f"Unexpected length of sub_id for {self.src}"
            raise RuntimeError(msg)

        return ses[0]

    @property
    def src_id(self) -> str:
        return self.src.stem

    @property
    def dst(self) -> Path:
        return (
            self.dst_root
            / f"sub={self.sub_id}"
            / f"ses={self.ses_id}"
            / f"src={self.src_id}"
            / "motion.arrow"
        )

    def process_run(self) -> None:
        if self.dst.exists():
            logging.info(f"{self.dst} already exists--skipping")
            return
        if not (parent := self.dst.parent).exists():
            parent.mkdir(parents=True)

        (
            pl.read_csv(
                self.src,
                separator="\t",
                schema_overrides={"t_indx": pl.Int64},
            )
            .with_columns(
                pl.selectors.contains("_x", "_y", "_z").cast(pl.Float64)
            )
            .rename({"t_indx": "t"})
            .with_columns(pl.selectors.starts_with("rot") * np.pi / 180)
            .write_ipc(self.dst, compression="zstd")
        )


def main(src: Path, dst_root: Path):
    for to_process in src.glob("sub*"):
        logging.info(f"Processing {to_process}")
        for txt in to_process.glob("ses*/func/*motion.tsv"):
            motion_processor = MotionProcesser(src=txt, dst_root=dst_root)
            motion_processor.process_run()

    logging.info("completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("src", type=Path)
    parser.add_argument("dst_root", type=Path)

    args = parser.parse_args()
    main(src=args.src, dst_root=args.dst_root)
