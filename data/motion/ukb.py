import argparse
from pathlib import Path
import re
import logging

import pydantic
import polars as pl
import pandas as pd

from pymrimisc import motion

logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s", level=logging.INFO
)


class MotionProcesser(pydantic.BaseModel):
    src: Path
    dst_root: Path

    @property
    def sub_id(self) -> str:
        sub = re.findall(r"\d{7}", str(self.src))
        if not (len(sub) == 1):
            msg = f"Unexpected length of sub_id for {self.src}"
            raise RuntimeError(msg)

        return sub[0]

    @property
    def ses_id(self) -> str:
        ses = re.findall(r"(?<=ses-)\d+", str(self.src))
        if not (len(ses) == 1):
            msg = f"Unexpected length of ses_id for {self.src}"
            raise RuntimeError(msg)

        return ses[0]

    @property
    def src_id(self) -> str:
        if "20227" in str(self.src):
            datatype = "20227"
        elif "20249" in str(self.src):
            datatype = "20249"
        else:
            msg = "Unknown file"
            raise RuntimeError(msg)
        return datatype

    @property
    def dst(self) -> Path:
        return (
            self.dst_root
            / f"sub={self.sub_id}"
            / f"ses={self.ses_id}"
            / f"src={self.src_id}"
            / "motion.arrow"
        )

    @property
    def rmsd(self) -> Path:
        return self.src.with_name("prefiltered_func_data_mcf_rel.rms")

    def process_run(self) -> None:
        if self.dst.exists():
            logging.info(f"{self.dst} already exists--skipping")
            return
        if not (parent := self.dst.parent).exists():
            parent.mkdir(parents=True)

        rmsd_tail = pl.read_csv(
            self.rmsd,
            separator="\n",
            has_header=False,
            new_columns=["rmsd"],
        )
        rmsd = pl.concat(
            [pl.DataFrame({"rmsd": 0.0}), rmsd_tail]
        ).with_row_index(name="t")

        # polars unable to read handle variable whitespace delim
        tmp = pd.read_csv(
            self.src,
            sep=r"\s+",
            names=["rot_x", "rot_y", "rot_z", "trans_x", "trans_y", "trans_z"],
        )

        pl.from_pandas(tmp).with_row_index("t").join(rmsd, on="t").write_ipc(
            self.dst, compression="zstd"
        )


def main(src: Path, dst_root: Path):
    for to_process in src.glob("sub*"):
        logging.info(f"Processing {to_process}")
        for ses in [2, 3]:
            for datatype in ["20249", "20227"]:
                for data in (
                    to_process / f"ses-{ses}" / "non-bids" / datatype / "fMRI"
                ).glob("*MRI*/mc/prefiltered_func_data_mcf.par"):
                    motion_processor = MotionProcesser(
                        src=data, dst_root=dst_root
                    )
                    logging.info(f"Working on {motion_processor.dst}")
                    motion_processor.process_run()

    logging.info("completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("src", type=Path)
    parser.add_argument("dst_root", type=Path)

    args = parser.parse_args()
    main(src=args.src, dst_root=args.dst_root)
