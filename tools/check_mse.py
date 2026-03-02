import re
from pathlib import Path

import numpy as np
import polars as pl

# $  mapfile -t subs < <(find /dcs07/smart/data/ukb/rawdata -maxdepth 1 -type d -name "sub*")
# $  for s in "${subs[@]}"; do src=$s/ses-2/non-bids/25752; if [[ -d $src ]]; then ss=$(echo $s | grep -oP "[0-9]{7}"); echo ${src}/${ss}_25752_2_0.txt >> 25752_2_0; fi; done

files = Path("25752_2_0").read_text().splitlines()

for f in files:
    sub = re.findall(r"\d{7}", f)[0]
    if not (dst := Path(f"data/{sub=}/data.parquet")).exists():
        print(f"{sub=}")
        pl.DataFrame(
            {"connectivity": np.loadtxt(f)}
        ).with_row_index().lazy().sink_parquet(dst, mkdir=True)

cols = pl.scan_parquet("data").select("index").unique().collect().sort("index")
# now consolidate into a single file
pl.scan_parquet("data").pivot(on="index", on_columns=cols, index="sub").sort(
    pl.col("sub")
).sink_parquet("25752_2_0.parquet")
