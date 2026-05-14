# Partition Survival Analysis V3

## Summary

`partition_survivalV3.py` is a Python script that performs **optimal threshold-based survival analysis** on clinical or biomedical data. Given a dataset with survival times, event indicators, and a continuous variable (e.g., a biomarker or gene expression level), the script exhaustively searches for the threshold value that best partitions patients into two groups with maximally different survival outcomes.

Two statistical tests are used to evaluate candidate thresholds:

- **Log-Rank test** — the standard non-parametric test for comparing survival distributions; equally weights all time points.
- **Gehan-Breslow-Wilcoxon test** — a weighted variant that gives more weight to earlier events, making it more sensitive to differences in early survival.

For each test, the script identifies the threshold that produces the highest test statistic, then generates **Kaplan-Meier survival curves** for the resulting groups.

## Dependencies

| Package    | Version | Purpose                                    |
|------------|---------|--------------------------------------------|
| pandas     | ≥ 1.3   | Data loading and manipulation              |
| numpy      | ≥ 1.20  | Numerical operations                       |
| lifelines  | ≥ 0.27  | Survival analysis (KM fitting, log-rank)   |
| matplotlib | ≥ 3.4   | Plotting Kaplan-Meier curves               |

### Install Dependencies

```bash
pip install pandas numpy lifelines matplotlib
```

## Input File Format

The script expects a **tab-delimited** text file named `Input.txt` in the same directory as the script. The file must contain the following three columns (header names are case-sensitive):

| Column     | Type    | Description                                          |
|------------|---------|------------------------------------------------------|
| `Months`   | Numeric | Survival or follow-up time in months                 |
| `Event`    | Integer | Event indicator: `1` = death/event, `0` = censored   |
| `Variable` | Numeric | Continuous variable to partition (e.g., biomarker)   |

### Example `Input.txt`

```
Months	Event	Variable
12.5	1	3.2
24.0	0	7.8
6.3	1	1.1
36.1	0	5.5
18.7	1	4.0
```

## Usage

1. Place `Input.txt` (tab-delimited, with the required columns) in the same directory as the script.
2. Run the script:

```bash
python partition_survivalV3.py
```

3. Review the console output and generated files.

### Configuration

The minimum group size for each partition is controlled by the `min_group_size` variable on **line 123** of the script. The default value is `10`, meaning any candidate threshold that would produce a group smaller than 10 samples is skipped. Adjust this value as needed for your dataset size.

## Output

The script produces three outputs:

### 1. Console Output

Prints the best threshold, test statistic, and p-value for both the Log-Rank and Gehan-Breslow-Wilcoxon tests.

### 2. `survival_analysis_results.txt`

A text file written to the working directory containing the same results:

```
Log-Rank Test Results:
Best Threshold: <value>
Test Statistic: <value>
P-Value: <value>

Gehan-Breslow-Wilcoxon Test Results:
Best Threshold: <value>
Test Statistic: <value>
P-Value: <value>
```

### 3. Kaplan-Meier Plot Images

- **`KM_LogRank.png`** — Survival curves for the two groups defined by the Log-Rank optimal threshold.
- **`KM_Gehan.png`** — Survival curves for the two groups defined by the Gehan-Breslow-Wilcoxon optimal threshold.

Each plot shows the estimated survival probability over time for both groups (`Variable ≤ threshold` vs. `Variable > threshold`).

## How It Works

1. **Load data** from `Input.txt` and validate required columns and event encoding.
2. **Threshold search** — For each unique value of `Variable`, split patients into two groups (≤ threshold and > threshold). Skip splits where either group has fewer than `min_group_size` samples.
3. **Statistical testing** — For each valid split, compute the Log-Rank (or Gehan-Breslow-Wilcoxon) test statistic. Track the threshold that maximizes the test statistic.
4. **Report results** — Write the best threshold, statistic, and p-value to both the console and `survival_analysis_results.txt`.
5. **Plot** — Fit Kaplan-Meier curves for the two groups at each optimal threshold and save the figures as PNG files.

## Notes

- The p-values reported are **unadjusted** for multiple testing. Because the script tests many candidate thresholds, the reported p-values may be optimistically low. Consider applying a correction (e.g., Bonferroni) or using permutation-based approaches if formal significance testing is required.
- The script uses `plt.show()`, which will open interactive plot windows when run outside of a headless environment. In headless/server environments, the plots are still saved to disk.
- If the input file uses a different delimiter or filename, modify **line 9** of the script accordingly.
