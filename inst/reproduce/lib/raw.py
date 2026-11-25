import hashlib
import os
import sys
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd


SKAB = "https://raw.githubusercontent.com/waico/SKAB/b2c0d46c2971dcbfe71e26087b6d231998bb91c2/data/"
SOURCES = {
    "normal": (SKAB + "anomaly-free/anomaly-free.csv",
               "7641a8742957100ef424dd97adfb3eeab27d4d180645fe2b3c49b7f7146ae9a6"),
    "rotor": (SKAB + "other/7.csv",
              "1bee14f49bfc817ba2c335fb427baa71153d8f8ab52e72652f15f74e4b77daec"),
    "skin": ("https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/33_skin.npz",
             "d02aa01ffc1b6bf2b8b1ac250c90e8a8eacd0cefc962940ed5cc4bf5deac9c95"),
}
SENSORS = ["Accelerometer1RMS", "Accelerometer2RMS", "Current", "Pressure",
           "Temperature", "Thermocouple", "Voltage", "RateRMS"]


def get_source(record, work):
    url, expected = SOURCES[record]
    local = os.environ.get("FDBAH_DATA_DIR")
    path = (Path(local) if local else work) / url.rsplit("/", 1)[1]
    if not local:
        with urlopen(url, timeout=90) as response:
            path.write_bytes(response.read())
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Source checksum differs for {path.name}")
    return path


def prepare_skab(path, record):
    data = pd.read_csv(path, sep=None, engine="python").rename(
        columns={"Volume Flow RateRMS": "RateRMS"})
    data[SENSORS] = data[SENSORS].apply(pd.to_numeric, errors="coerce")
    labels = (pd.to_numeric(data["anomaly"], errors="coerce") if "anomaly" in data
              else pd.Series(0, index=data.index))
    valid = (np.isfinite(data[SENSORS]).all(axis=1)
             & data["datetime"].notna() & labels.isin([0, 1]))
    sample = data.loc[valid, SENSORS].copy()
    sample["label"] = labels.loc[valid].astype(int)
    if record == "normal":
        if len(sample) != 9405:
            raise ValueError("Unexpected SKAB normal record length")
        sample = sample.iloc[np.arange(1500) * (len(sample) - 1) // 1499]
    return sample


def prepare_skin(path):
    with np.load(path, allow_pickle=False) as archive:
        x = np.asarray(archive["X"], dtype=float)
        y = np.asarray(archive["y"]).reshape(-1).astype(int)
    def row_hash(i):
        return hashlib.sha256(f"Skin:{i + 1}".encode("ascii")).digest()
    rows = np.sort(sorted(range(len(y)), key=row_hash)[:2500])
    sample = pd.DataFrame(x[rows], columns=["B", "G", "R"])
    sample["label"] = y[rows]
    return sample


if __name__ == "__main__":
    record, out = sys.argv[1], Path(sys.argv[2])
    source = get_source(record, out.parent)
    if record == "skin":
        prepare_skin(source).to_csv(out, index=False)
    else:
        prepare_skab(source, record).to_csv(out, index=False, float_format="%.17g")
