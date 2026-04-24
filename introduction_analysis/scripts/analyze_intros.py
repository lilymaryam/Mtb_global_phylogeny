import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from pathlib import Path

assert Path.cwd().name == "introduction_analysis", "Not in introduction_analysis directory"

nodes = {}
country = {}
count = 0
with open('data/introduce_output.tsv') as f:
    line = f.readline().strip().split()
    print([line[1],line[6],line[12], 'lineage'])
    #for country 
    for line in f:
        line = line.strip().split()
        #print(len(line))
        if len(line) == 17:
            #print(line)
            lineage= line[15].split('.')[0].replace('lineage','L')
            if lineage != 'none':
                line = [line[1],line[6],line[12], lineage]
            else:
                line = [line[1],line[6],line[12], 'None']
            #print(line)
            #break
        else:
            line = [line[1],line[6],line[12], 'None']
        assert len(line) == 4
        if line[3] == 'None' or ';' in line[3]:
            #print(line)'][]
            #print(line[0])
            #count+=1
            continue
        if line[0] not in nodes and line[3]!='None':
            nodes[line[0]] = [line[3]]
        else:
            if line[3] != 'None' and line[3] not in nodes[line[0]]:
                nodes[line[0]].append(line[3])
        if line[2] not in country:
            country[line[2]] = {line[0]: line[1]}
        else:
            country[line[2]][line[0]] = line[1]

rank = {}
for c in country:
    rank[c] = 0
    for n in country[c]: 
        if len(nodes[n]) == 1:
            rank[c] += 1
            print(c, n, nodes[n][0], country[c][n])

# Get the top 10 key-value pairs as a list of tuples
top_10 = sorted(rank.items(), key=lambda item: item[1], reverse=True)[1:11]

print(top_10)

data = {}
for t in top_10:
    data[t[0]] = {'singles': {}, 'multiples': {}}
    for n in country[t[0]]:
        lin = nodes[n][0] 
        #print(country[t[0]][n])
        key = 'multiples'

        if country[t[0]][n] == '1':
            print('key', key)
            key = 'singles'
        #else = 
        data[t[0]][key][lin] = data[t[0]][key].get(lin, 0) + 1
        '''
        if lin in data[t[0]][key]:
            data[t[0]][key][lin] += 1
        if lin not in data[t[0]][key]:
            data[t[0]][key][lin] = 1
        '''

'''
data = {
    'GBR': {'singles': {'L4': 1399, 'L1': 589, 'La1': 34, 'L2': 359, 'L3': 986, 'L6': 11, 'La3': 27, 'L5': 22, 'L9': 1, 'L7': 1}, 
            'multiples': {'L4': 869, 'L3': 526, 'L1': 252, 'La1': 30, 'L2': 181, 'L5': 8, 'L6': 7, 'La3': 10, 'L7': 1}},
    'AUS': {'singles': {'L2': 830, 'L1': 937, 'L4': 657, 'L3': 581, 'La1': 2, 'L6': 3, 'La3': 15, 'L5': 3}, 
            'multiples': {'L4': 105, 'L1': 158, 'L2': 159, 'L3': 84, 'La1': 2, 'L9': 1}},
    'IND': {'singles': {'L4': 348, 'L2': 208, 'L3': 851, 'L1': 576, 'La3': 10}, 
            'multiples': {'L2': 101, 'L4': 118, 'L1': 223, 'L3': 373, 'La3': 9, 'La1': 2}},
    'USA': {'singles': {'L5': 8, 'L6': 18, 'L1': 440, 'L4': 884, 'L3': 261, 'L2': 543, 'La1': 48, 'La3': 7, 'La2': 2, 'L7': 1}, 
            'multiples': {'L4': 295, 'L2': 117, 'La1': 30, 'La3': 5, 'L1': 77, 'L3': 36, 'L6': 4, 'L5': 2, 'L9': 1}},
    'CHN': {'singles': {'L2': 777, 'L4': 264, 'L3': 7, 'L1': 24, 'La1': 2}, 
            'multiples': {'L2': 590, 'L4': 125, 'L3': 4, 'La2': 1, 'La1': 4, 'L1': 1}},
    'ZAF': {'singles': {'L2': 206, 'L4': 557, 'L1': 36, 'L3': 31, 'L5': 1, 'La3': 1}, 
            'multiples': {'L4': 372, 'L2': 153, 'L1': 28, 'L3': 9, 'La1': 3}},
    'CAN': {'singles': {'L4': 207, 'L3': 88, 'L1': 345, 'L2': 264, 'La1': 11, 'L5': 4, 'La3': 29, 'L6': 2}, 
            'multiples': {'L4': 71, 'L2': 55, 'L1': 60, 'La1': 2, 'La3': 12, 'L3': 6}},
    'GEO': {'singles': {'L2': 267, 'L4': 341, 'L3': 4, 'La1': 4}, 
            'multiples': {'L2': 240, 'L4': 226, 'L3': 1}},
    'DEU': {'singles': {'L2': 174, 'La3': 3, 'L4': 598, 'L3': 111, 'L1': 52, 'La1': 9, 'L9': 1, 'L6': 7, 'L5': 2}, 
            'multiples': {'L2': 21, 'L4': 91, 'L3': 4, 'L1': 2, 'La2': 1, 'L5': 1}},
    'ITA': {'singles': {'L2': 76, 'L4': 555, 'L3': 51, 'L5': 6, 'L1': 66, 'La1': 8, 'L6': 14, 'La2': 1, 'La3': 2}, 
            'multiples': {'L2': 37, 'L3': 10, 'L4': 194, 'L1': 13, 'L6': 4, 'La1': 3}}
}
'''
# Extract categories
categories = sorted(set(k for c in data.values() for m in c.values() for k in m))


print(f"Number of lineages: {len(categories)}")
print(f"Lineages: {categories}")
print("Script currently expects 11 lineages, if more are detected, more colors will need to be added.")


tol_vibrant = ["#EE6677", "#4477AA", "#228833", "#CCBB44", "#66CCEE", 
               "#AA3377", "#BBBBBB", "#E69F00", "#56B4E9", "#009E73", "#F0E442"]

color_map = {cat: tol_vibrant[i] for i, cat in enumerate(categories)}
# Assign colors
#colors = plt.cm.get_cmap("Set3", len(categories))
#color_map = {cat: colors(i) for i, cat in enumerate(categories)}

#color_map['La2'] = '#D62728'  # Red
#color_map['La3'] = '#9467BD' 
# Setup figure
fig, ax = plt.subplots(figsize=(10, 6))

bar_width = 0.35  # Adjusted width for spacing
gap = 0.05        # Space between singles/multiples
x_labels = list(data.keys())
x = np.arange(len(x_labels))

bottom_singles = np.zeros(len(x_labels))
bottom_multiples = np.zeros(len(x_labels))

bars = {}  # Store bars for legend

for category in categories:
    singles_values = [data[country]['singles'].get(category, 0) for country in x_labels]
    multiples_values = [data[country]['multiples'].get(category, 0) for country in x_labels]

    bars_singles = ax.bar(x - (bar_width + gap) / 2, singles_values, bar_width, 
                           bottom=bottom_singles, color=color_map[category], label=category if category not in bars else "")
    bars_multiples = ax.bar(x + (bar_width + gap) / 2, multiples_values, bar_width, 
                            bottom=bottom_multiples, color=color_map[category])

    bottom_singles += singles_values
    bottom_multiples += multiples_values

    bars[category] = bars_singles[0]  # Store one bar per category for legend

# Add * for singles and ** for multiples above the highest bar segment
'''
for i in range(len(x_labels)):
    ax.annotate("1", xy=(x[i] - (bar_width + gap) / 2, bottom_singles[i] + 5), 
                ha='center', fontsize=12, color='black', fontweight='bold')
    ax.annotate(">1", xy=(x[i] + (bar_width + gap) / 2, bottom_multiples[i] + 5), 
                ha='center', fontsize=12, color='black', fontweight='bold')
'''
y_offset = -5  # Adjust as needed to position below bars

for i in range(len(x_labels)):
    ax.annotate("1", xy=(x[i] - (bar_width + gap) / 2, 0), 
                ha='center', va='top', fontsize=10, color='black',
                xytext=(0, y_offset), textcoords="offset points")
    
    ax.annotate(">1", xy=(x[i] + (bar_width + gap) / 2, 0), 
                ha='center', va='top', fontsize=10, color='black',
                xytext=(0, y_offset), textcoords="offset points")

# Adjust x-axis label positions to be lower
ax.set_xlabel("Country", fontsize=12, labelpad=15)  # Increase labelpad to move x-axis label lower
ax.set_ylabel("Count", fontsize=12)

# Lower x-tick labels
ax.tick_params(axis="x", which="both", pad=10)  # Increase padding to move them down


ax.set_xticks(x)
ax.set_xticklabels(x_labels)
ax.set_ylabel("Count")
ax.set_title("MTB Introductions by Country and Lineage")

# Restore legend for mutation types inside the plot area
legend1 = ax.legend(bars.values(), bars.keys(), title="Lineages",bbox_to_anchor=(1, 1),  loc='upper right')

# Create second legend for * and ** markers
star_patch = mpatches.Patch(color='white', label="1 Single Introductions")
double_star_patch = mpatches.Patch(color='white', label=">1 Multiple Transmissions")
fig.legend(handles=[star_patch, double_star_patch], loc="upper right", bbox_to_anchor=(.86, .941), title="Bar Labels")

# Tighten layout and ensure both legends are visible
plt.tight_layout()
#plt.subplots_adjust(bottom=0.2)

plt.savefig('Full_Tree_Introduction_Analysis.png', format='png', dpi=600)