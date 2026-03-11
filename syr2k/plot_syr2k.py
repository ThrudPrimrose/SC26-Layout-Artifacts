import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ---- Load CSV ----
df = pd.read_csv("data.csv")

# ---- Build variant name ----
def make_variant(row):
    parts = []

    if row["bi"] != 0 or row["bj"] != 0 or row["bk"] != 0:
        parts.append(f"blocked_{int(row['bi'])}_{int(row['bj'])}_{int(row['bk'])}")

    if row["ti"] != 0 or row["tj"] != 0 or row["tk"] != 0:
        parts.append(f"tiled_{int(row['ti'])}_{int(row['tj'])}_{int(row['tk'])}")

    if not parts:
        return "baseline"

    return "+".join(parts)

df["variant"] = df.apply(make_variant, axis=1)


# ---- Create color grouping (blocked and tiled share color if dims equal) ----
def color_group(v):
    if v == "baseline":
        return "baseline"

    # extract dimensions
    import re
    m = re.search(r'_(\d+)_(\d+)_(\d+)', v)
    if m:
        return f"{m.group(1)}_{m.group(2)}_{m.group(3)}"
    return v

df["color_group"] = df["variant"].apply(color_group)


# ---- Runtime column (adjust if needed) ----
runtime_col = "runtime"


# ---- Plot ----
plt.figure(figsize=(12,6))

sns.violinplot(
    data=df,
    x="variant",
    y=runtime_col,
    hue="color_group",
    dodge=False,
    cut=0
)

plt.xticks(rotation=60)
plt.tight_layout()
plt.show()