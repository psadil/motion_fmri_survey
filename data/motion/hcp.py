import argparse
from pathlib import Path
import re
import logging

import pydantic
import numpy as np
import polars as pl
import pandas as pd

from pymrimisc import motion


logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s", level=logging.INFO
)


class MotionProcesser(pydantic.BaseModel):
    src: Path
    dst_root: Path
    sessions: Path | None = None

    @property
    def sub_id(self) -> str:
        sub = re.findall(r"\d{6,7}", str(self.src))
        if not (len(sub) == 1):
            msg = f"Unexpected length of sub_id for {self.src}"
            raise RuntimeError(msg)

        return sub[0]

    @property
    def ses_id(self) -> str:
        if self.sessions:
            ses = (
                pl.read_csv(self.sessions / f"{self.sub_id}.csv")
                .with_columns(ses=pl.col("Session Day").str.extract(r"(\d)"))
                .filter(pl.col("Scan Description") == self.src_id)
                .select("ses")
                .to_series()
                .to_list()
            )
            if not len(ses) == 1:
                # missing sessions (likely due to scans being unusable)
                # https://wiki.humanconnectome.org/docs/HCP%20Data%20Release%20Updates%20Known%20Issues%20and%20Planned%20fixes.html
                if (
                    self.sub_id == "809252"
                    and self.src_id == "tfMRI_SOCIAL_RL"
                ):
                    ses = ["2"]
                elif self.sub_id == "196952" and self.src_id == "tfMRI_WM_LR":
                    ses = ["1"]
                else:
                    msg = f"Unexpected length of ses_id for {self.src}"
                    raise RuntimeError(msg)
            ses_id = ses[0]
        else:
            ses_id = "1"
        return ses_id

    @property
    def src_id(self) -> str:
        return self.src.parent.name

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
        return self.src.with_name("Movement_RelativeRMS.txt")

    def process_run(self) -> None:
        if self.dst.exists():
            logging.info(f"{self.dst} already exists--skipping")
            return
        if not (parent := self.dst.parent).exists():
            parent.mkdir(parents=True)

        rmsd = pl.read_csv(
            self.rmsd,
            separator="\n",
            has_header=False,
            new_columns=["rmsd"],
            schema={"rmsd": pl.Float64},
        ).with_row_index("t")

        # polars unable to read handle variable whitespace delim
        tmp = pd.read_csv(
            self.src,
            sep=r"\s+",
            names=[
                "trans_x",
                "trans_y",
                "trans_z",
                "rot_x",
                "rot_y",
                "rot_z",
                "trans_x_derivative1",
                "trans_y_derivative1",
                "trans_z_derivative1",
                "rot_x_derivative1",
                "rot_y_derivative1",
                "rot_z_derivative1",
            ],
        )

        pl.from_pandas(tmp).with_row_index("t").with_columns(
            pl.selectors.starts_with("rot") * np.pi / 180
        ).join(rmsd, on="t").write_ipc(self.dst, compression="zstd")


def main(src: Path, dst_root: Path, sessions: Path | None = None):
    for to_process in src.glob("*"):
        logging.info(f"Processing {to_process}")
        for txt in (to_process / "MNINonLinear" / "Results").glob(
            "*MRI_*_*/Movement_Regressors.txt"
        ):
            if "_7T_" in str(txt):
                continue
            motion_processor = MotionProcesser(
                src=txt, dst_root=dst_root, sessions=sessions
            )
            motion_processor.process_run()

    logging.info("completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("src", type=Path)
    parser.add_argument("dst_root", type=Path)
    parser.add_argument("--sessions", type=Path)

    args = parser.parse_args()
    main(src=args.src, dst_root=args.dst_root, sessions=args.sessions)
