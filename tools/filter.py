import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy import signal

d = pl.scan_parquet(
    "/Users/psadil/Library/CloudStorage/OneDrive-JohnsHopkins/data/motion/derivatives/human-connectome-project-openaccess.parquet"
).collect()

tr = 0.72
bw = 0.41
w0 = 0.31
b, a = signal.iirnotch(w0=0.31, Q=w0 / bw, fs=1 / tr)

s = d.select("trans_x").to_series().to_numpy()
n = range(len(s))

fig = plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot(n, s, color="r", linewidth=2)
plt.xlabel("Time", fontsize=20)
plt.ylabel("Magnitude", fontsize=18)
plt.title("Noisy Signal", fontsize=20)

# Apply notch filter to the noisy signal using signal.filtfilt
out = signal.filtfilt(b, a, s)

# Plot output signal of notch filter
plt.plot(n, out)
plt.xlabel("Time", fontsize=20)
plt.ylabel("Magnitude", fontsize=18)
plt.title("Filtered Signal", fontsize=20)
plt.subplots_adjust(hspace=0.5)
fig.tight_layout()
plt.show()


fig = plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot(n[100:200], s[100:200], color="r", linewidth=2)
plt.xlabel("Time", fontsize=20)
plt.ylabel("Magnitude", fontsize=18)
plt.title("Noisy Signal", fontsize=20)
# Plot output signal of notch filter
plt.plot(n[100:200], out[100:200])
plt.xlabel("Time", fontsize=20)
plt.ylabel("Magnitude", fontsize=18)
plt.title("Filtered Signal", fontsize=20)
plt.subplots_adjust(hspace=0.5)
fig.tight_layout()
plt.show()

# now, try finding peaks
x, y = signal.periodogram(s, fs=1 / tr)
fig = plt.figure()
plt.plot(x, np.log(y))
plt.xlabel("Frequency")
plt.ylabel("Log Power Spectral Density")
plt.show()

xi = np.where(np.logical_and(x >= 0.1, x <= 0.4))
xx = x[xi]
yy = y[xi]

peak = np.where(yy == np.max(yy))

fig = plt.figure()
plt.plot(xx, np.log(yy))
plt.plot(xx[peak], np.log(yy[peak]), "ro")
plt.xlabel("Frequency")
plt.ylabel("Log Power Spectral Density")
plt.show()


def find_peak(
    series: pl.Series, lower: float = 0.1, upper: float = 0.4, tr=0.72
) -> float:
    x, y = signal.periodogram(series.to_numpy(), fs=1 / tr)
    xi = np.where(np.logical_and(x >= lower, x <= upper))
    yy = y[xi]
    xx = x[xi]
    return xx[np.where(yy == np.max(yy))][0]


peaks = (
    d.drop(pl.selectors.contains("deriv"), "ses")
    .unpivot(index=["t", "sub", "task", "ped"])
    .sort("sub", "task", "ped", "variable", "t")
    .group_by("sub", "task", "ped", "variable")
    .agg(
        pl.col("value").map_batches(
            find_peak, returns_scalar=True, return_dtype=pl.Float64
        )
    )
)

peaks.group_by("task", "variable").agg(
    med=pl.col("value").median(),
    lower=pl.col("value").quantile(0.25),
    upper=pl.col("value").quantile(0.75),
)

peaks.group_by("task", "variable").agg(med=pl.col("value").median()).pivot(
    on=["variable"], values=["med"]
)
