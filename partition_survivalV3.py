import pandas as pd
import numpy as np
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

# Load the data
input_file = "Input.txt"
data = pd.read_csv(input_file, sep="\t")  # Assuming tab-delimited file

# Ensure data contains the required columns
if not {"Months", "Event", "Variable"}.issubset(data.columns):
    raise ValueError("Input file must contain 'Months', 'Event', and 'Variable' columns.")

# Extract relevant columns
time = data["Months"]
event = data["Event"]
variable = data["Variable"]

# Ensure proper conventions for Event column (1 = death, 0 = censored)
if not set(event.unique()).issubset({0, 1}):
    raise ValueError("Event column must contain only 0 (censored) and 1 (death).")

# Function to find best threshold using a specified test
def find_best_threshold(time, event, variable, min_group_size=5, test_type="logrank"):
    """
    Finds the best threshold for splitting the continuous variable into two groups,
    using either the Log-Rank test or Gehan-Breslow-Wilcoxon test.
    
    Parameters:
    - time: Series of survival times.
    - event: Series of event indicators (1=death, 0=censored).
    - variable: Series of the continuous variable.
    - min_group_size: Minimum number of samples required in each group.
    - test_type: Type of test to use ("logrank" or "gehan").
    
    Returns:
    - best_threshold: The optimal threshold value.
    - best_statistic: The highest test statistic.
    - best_p_value: The p-value for the best threshold.
    """
    best_threshold = None
    best_statistic = -np.inf
    best_p_value = 1.0

    unique_values = sorted(variable.unique())
    
    # Loop through potential thresholds
    for threshold in unique_values:
        # Split data into two groups based on the threshold
        group1 = (variable <= threshold)
        group2 = (variable > threshold)

        # Check if both groups meet the minimum size requirement
        if group1.sum() < min_group_size or group2.sum() < min_group_size:
            continue

        if test_type == "logrank":
            # Perform Log-Rank test
            results = logrank_test(
                time[group1], time[group2],
                event[group1], event[group2]
            )
        elif test_type == "gehan":
            # Perform Gehan-Breslow-Wilcoxon test
            results = multivariate_logrank_test(
                time, 
                np.where(group1, "Group1", "Group2"), 
                event, 
                weightings="wilcoxon"
            )
        else:
            raise ValueError("Invalid test_type. Use 'logrank' or 'gehan'.")

        # Update best threshold if the statistic improves
        if results.test_statistic > best_statistic:
            best_statistic = results.test_statistic
            best_p_value = results.p_value
            best_threshold = threshold

    return best_threshold, best_statistic, best_p_value

# Function to plot Kaplan-Meier curves based on a threshold and optionally save the plot as a PNG file
def plot_km_curve(threshold, time, event, variable, title_suffix="", save_filename=None):
    """
    Plots Kaplan-Meier survival curves for two groups defined by a threshold.
    
    Parameters:
    - threshold: The value used to split the variable.
    - time: Series of survival times.
    - event: Series of event indicators.
    - variable: Series of the continuous variable.
    - title_suffix: Suffix to add to the plot title for clarification.
    - save_filename: Optional filename to save the plot as a .png image.
    """
    kmf_low = KaplanMeierFitter()
    kmf_high = KaplanMeierFitter()

    # Split the data into two groups based on the threshold
    group_low = (variable <= threshold)
    group_high = (variable > threshold)

    # Fit Kaplan-Meier curves
    kmf_low.fit(time[group_low], event[group_low], label=f"Variable <= {threshold}")
    kmf_high.fit(time[group_high], event[group_high], label=f"Variable > {threshold}")

    # Plot the survival curves
    plt.figure(figsize=(10, 6))
    kmf_low.plot_survival_function()
    kmf_high.plot_survival_function()
    plt.title(f"Kaplan-Meier Survival Curves {title_suffix}")
    plt.xlabel("Months")
    plt.ylabel("Survival Probability")
    plt.legend()
    plt.grid()

    if save_filename:
        plt.savefig(save_filename)  # Save the figure as a PNG file

    plt.show()

# Set minimum group size (adjustable)
min_group_size = 10

# Find best threshold using Log-Rank test
logrank_results = find_best_threshold(time, event, variable, min_group_size, test_type="logrank")

# Find best threshold using Gehan-Breslow-Wilcoxon test
gehan_results = find_best_threshold(time, event, variable, min_group_size, test_type="gehan")

# Save results to a text file
with open("survival_analysis_results.txt", "w") as f:
    f.write("Log-Rank Test Results:\n")
    f.write(f"Best Threshold: {logrank_results[0]}\n")
    f.write(f"Test Statistic: {logrank_results[1]}\n")
    f.write(f"P-Value: {logrank_results[2]}\n\n")
    
    f.write("Gehan-Breslow-Wilcoxon Test Results:\n")
    f.write(f"Best Threshold: {gehan_results[0]}\n")
    f.write(f"Test Statistic: {gehan_results[1]}\n")
    f.write(f"P-Value: {gehan_results[2]}\n")

# Display results in console
print("Log-Rank Test Results:")
print(f"Best Threshold: {logrank_results[0]}")
print(f"Test Statistic: {logrank_results[1]}")
print(f"P-Value: {logrank_results[2]}")

print("\nGehan-Breslow-Wilcoxon Test Results:")
print(f"Best Threshold: {gehan_results[0]}")
print(f"Test Statistic: {gehan_results[1]}")
print(f"P-Value: {gehan_results[2]}")

# Plot and save Kaplan-Meier curves using the best Log-Rank threshold
plot_km_curve(logrank_results[0], time, event, variable, 
              title_suffix="(Log-Rank Threshold)", 
              save_filename="KM_LogRank.png")

# Plot and save Kaplan-Meier curves using the best Gehan-Breslow-Wilcoxon (Wilcoxon) threshold
plot_km_curve(gehan_results[0], time, event, variable, 
              title_suffix="(Gehan/Wilcoxon Threshold)", 
              save_filename="KM_Gehan.png")
