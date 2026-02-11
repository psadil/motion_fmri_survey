from pathlib import Path
import typing
import logging

import numpy as np
import polars as pl
import pydantic
from nilearn.glm import first_level


logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    level=logging.INFO,
)


class Header(pydantic.BaseModel):
    data_type: str
    dim1: int
    dim2: int
    dim3: int
    dim4: int
    datatype: int
    pixdim1: float
    pixdim2: float
    pixdim3: float
    pixdim4: float
    cal_max: float
    cal_min: float
    file_type: str

    @classmethod
    def from_fslinfo(cls, src: Path) -> typing.Self:
        raw = src.read_text().splitlines()
        items = [item.split() for item in raw]
        header = {x[0]: x[1] for x in items}
        return cls(**header)  # type: ignore


def write_events(events: pl.DataFrame, dst: Path) -> None:
    events.with_columns(
        sub=pl.col("src").str.extract(r"sub-(\w+)_ses"),
        ses=pl.col("src").str.extract(r"ses-(\w+)_task"),
        task=pl.col("src").str.extract(r"task-(\w+)_run"),
        run=pl.col("src").str.extract(r"run-(\d{2})"),
    ).drop("src").write_parquet(dst)


root = Path("sourcedata/design/abcd")
for task in ["nback", "sst", "mid"]:
    logging.info(f"Finding {task}")
    designs = []
    for s, sub in enumerate(
        (root / "dcs07" / "smart" / "data" / "abcd" / "rawdata").glob("sub*")
    ):
        logging.info(f"Processing {s}: {sub.name}")
        for events in sub.rglob(f"*{task}*events.tsv"):
            header_file = root / events.name.replace("events.tsv", "bold")
            if not header_file.exists():
                logging.warning(f"Unable to find {header_file}")
                continue

            header = Header.from_fslinfo(header_file)

            tsv = pl.read_csv(events, separator="\t")

            design = first_level.make_first_level_design_matrix(
                np.arange(
                    0,
                    header.dim4 * header.pixdim4,
                    header.pixdim4,
                ),
                events=tsv.to_pandas(),
                hrf_model=None,
                drift_model=None,
            ).rename_axis("time")

            designs.append(
                pl.from_pandas(design, include_index=True)
                .with_columns(src=pl.lit(events.name))
                .with_row_index("t")
            )

    write_events(
        pl.concat(designs, how="diagonal"),
        Path("derivatives") / "abcd" / f"{task}.parquet",
    )

pl.scan_csv(
    root
    / "dcs07"
    / "smart"
    / "data"
    / "abcd"
    / "rawdata"
    / "*"
    / "*"
    / "func"
    / "*events.tsv",
    separator="\t",
    include_file_paths="src",
    null_values="n/a",
).with_columns(
    sub=pl.col("src").str.extract(r"sub-(\w+)_ses"),
    ses=pl.col("src").str.extract(r"ses-(\w+)_task"),
    task=pl.col("src").str.extract(r"task-(\w+)_run"),
    run=pl.col("src").str.extract(r"run-(\w+)_event"),
).sink_parquet(
    "derivatives/abcd.parquet"
)
