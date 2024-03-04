import math
import os
import sys
import argparse
import json
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Process input and output directories')
parser.add_argument('--file', help='Specify the input file', required=False)

args = parser.parse_args()
fname = args.file

with open(fname, 'r') as file:
    data = json.load(file)



def get_all_radii(data):
    radii = set()

    for prot_id in data:
        for chain_results in data[prot_id]:
            results = chain_results['results']

            for record in results:
                radius = record['radius']
                radii.add(radius)

    return sorted(list(radii))

def add(dict, key, value):
    if key not in dict:
        dict[key] = 0
    dict[key] += value

def expand_histogram(histogram, new_values):
    for val in new_values:
        add(histogram, val, 1)


def get_histogram_for(radius, data):
    histogram = {}

    for prot_id in data:
        for chain_results in data[prot_id]:
            results = chain_results['results']
            for record in results:
                curr_radius = record['radius']

                if radius != curr_radius:
                    continue

                expand_histogram(histogram, record['protrusion'])

    return histogram

def calculate_mean_and_coefficient_of_variation(histogram):
    total_occurrences = sum(histogram.values())
    mean_value = sum(key * freq for key, freq in histogram.items()) / total_occurrences

    variance = sum(((key - mean_value) ** 2) * freq for key, freq in histogram.items()) / total_occurrences
    standard_deviation = math.sqrt(variance)

    coefficient_of_variation = (standard_deviation / mean_value)

    return mean_value, coefficient_of_variation

def calculate_mean_and_variance(histogram):
    total_occurrences = sum(histogram.values())
    mean_value = sum(key * freq for key, freq in histogram.items()) / total_occurrences
    variance = sum(((key - mean_value) ** 2) * freq for key, freq in histogram.items()) / total_occurrences

    return mean_value, variance

histograms = []    
radii = get_all_radii(data)

for radius in radii:
    histogram = get_histogram_for(radius, data)

    # mean, var = calculate_mean_and_variance(histogram)
    # print(f'r={radius}\tvar={var}\tmean={mean}')

    mean, coef_var = calculate_mean_and_coefficient_of_variation(histogram)
    print(f'r={radius}\tcoef var={coef_var}\tmean={mean}')

    histograms.append(histogram)


values_list = []
occurrences_list = []

for my_dict in histograms:
    values = list(my_dict.keys())
    occurrences = list(my_dict.values())

    values_list.append(values)
    occurrences_list.append(occurrences)

# Plotting

fig, ax = plt.subplots()

for i in range(len(histograms)):
    ax.bar(values_list[i], occurrences_list[i], label=f'Radius {radii[i]}', alpha=0.4)

ax.set_xlabel('Protrusion')
ax.set_ylabel('Number of Occurrences')
ax.legend()
plt.title('Comparison of protrusion using different radii')

ax.set_xlim(0, 60)
ax.set_ylim(0, 700000)

plt.savefig(f'comparison-for-{get_file_name_without_extension(fname)}.png')
