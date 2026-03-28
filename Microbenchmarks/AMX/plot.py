import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
for fname in ["benchmark.csv", "benchmark_2.csv"]:
    df = pd.read_csv(fname)

    # Compute median statistics
    median_stats = df.groupby('kernel').agg(
        median_time_ns=('time_ns', 'median'),
        median_gflops=('gflops', 'median')
    ).reset_index()

    # Sort kernels by median GFLOPS for better visualization
    median_stats = median_stats.sort_values('median_gflops', ascending=False)
    kernel_order = median_stats['kernel'].tolist()

    # Violin plot of GFLOPS
    plt.figure(figsize=(12, 6))
    sns.violinplot(x='kernel', y='gflops', data=df, order=kernel_order, inner='quartile')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('GFLOPS')
    plt.title('Distribution of GFLOPS across 100 runs per kernel')
    plt.tight_layout()
    plt.savefig(f'gflops_violin_{fname}.png', dpi=300)
    plt.show()

    # Print median summary
    print(f"Median runtime and GFLOPS per kernel for {fname}:")
    print(median_stats.to_string(index=False))