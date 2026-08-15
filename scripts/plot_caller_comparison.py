"""
Combine hap.py summary.csv outputs from GATK and DeepVariant into
a single precision / recall / F1 comparison bar chart, split by
variant type (SNP vs INDEL). Called via Snakemake `script:`.
"""
import pandas as pd
import matplotlib.pyplot as plt

input_paths = list(snakemake.input)
output_path = snakemake.output[0]

frames = []
for path in input_paths:
    caller = "gatk" if "gatk" in path else "deepvariant"
    df = pd.read_csv(path)
    # hap.py summary.csv has a 'Type' column (SNP/INDEL) and a
    # 'Filter' column (PASS/ALL) plus METRIC.Precision/Recall/F1_Score
    df = df[df["Filter"] == "PASS"]
    df["caller"] = caller
    frames.append(df[["Type", "caller", "METRIC.Precision", "METRIC.Recall", "METRIC.F1_Score"]])

combined = pd.concat(frames, ignore_index=True)
combined.to_csv(output_path.replace(".png", "_data.csv"), index=False)

metrics = ["METRIC.Precision", "METRIC.Recall", "METRIC.F1_Score"]
labels = ["Precision", "Recall", "F1"]
variant_types = combined["Type"].unique()

fig, axes = plt.subplots(1, len(variant_types), figsize=(5 * len(variant_types), 4), sharey=True)
if len(variant_types) == 1:
    axes = [axes]

for ax, vtype in zip(axes, variant_types):
    subset = combined[combined["Type"] == vtype]
    pivot = subset.set_index("caller")[metrics]
    pivot.columns = labels
    pivot.plot(kind="bar", ax=ax, legend=(ax is axes[0]))
    ax.set_title(vtype)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("")
    ax.set_ylabel("Score")
    ax.tick_params(axis="x", rotation=0)

fig.suptitle("GATK vs DeepVariant — hap.py benchmark against GIAB HG002 truth set")
fig.tight_layout()
fig.savefig(output_path, dpi=150)
