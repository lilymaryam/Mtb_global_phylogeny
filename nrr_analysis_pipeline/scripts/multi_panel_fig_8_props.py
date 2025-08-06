
import matplotlib.pyplot as plt
import sys
import numpy as np
import statistics
import os

'''
This is the final version of the multipanel simulated data figure. This one tracks proportion of scores rather than histogram. 
'''

def figs(datafile):
    """Extract Jaccard similarity scores from a file."""
    jaccard_scores = []
    total = 0
    ones = 0
    try:
        with open(datafile) as df:
            for line in df:
                line = line.strip().split()
                if len(line) < 3:
                    print(f"Skipping malformed line: {line}")
                    continue
                if line[2] != 'None':
                    #print(line[2])
                    total += 1
                    if line[2] == '1.0':
                        ones += 1

                    try:
                        jaccard_scores.append(float(line[2]))
                    except ValueError:
                        print(f"Skipping invalid score: {line[2]}")
    except FileNotFoundError:
        print(f"File not found: {datafile}")
    except Exception as e:
         print(f"Error processing file {datafile}: {e}")
    #print(datafile, ones, total, ones/total)
    #print(datafile, statistics.median(jaccard_scores))
    #print(jaccard_scores)
    
    return jaccard_scores

def proportions(jaccard_scores):
    total = (len(jaccard_scores))
    hundos = 0
    for i in jaccard_scores:
        if i == 1.0:
            hundos += 1
    print('hundos', hundos, total, hundos/total)

    bin_edges = np.linspace(0, 1, 21)
    # Compute the histogram
    hist, bin_edges = np.histogram(jaccard_scores, bins=bin_edges)
    props = [i/total for i in hist]
    #print(type(props))
    #print(props)
    #print(bin_edges.tolist())
    return props, bin_edges.tolist()[:-1]

def main():
    # Define input files and corresponding subplot titles
    
    if not os.getcwd().endswith('nrr_analysis_pipeline'):
        print(f"Error: This script must be run from the 'nrr_analysis_pipeline' directory.")
        print(f"Current directory: {os.getcwd()}")
        print(f"Please navigate to the nrr_analysis_pipeline directory and run the script again.")
        sys.exit(1)
    inputs = [
        ('results/true1.sim1', 'True 1 Sim 1'),
        ('results/true3.sim1', 'True 3 Sim 1'),
        ('results/true5.sim1', 'True 5 Sim 1'),
        ('results/true8.sim1', 'True 8 Sim 1'),
        ('results/true10.sim1', 'True 10 Sim 1'),
        ('results/true1.sim5', 'True 1 Sim 5'),
        ('results/true3.sim5', 'True 3 Sim 5'),
        #('/home/lily/scripts/Mtb_global_phylogeny/202500805-simulatedresults-moved/true3.sim5', 'True 3 Sim 5'),
        ('results/true5.sim5', 'True 5 Sim 5'),
        ('results/true8.sim5', 'True 8 Sim 5'),
        ('results/true10.sim5', 'True 10 Sim 5'),
        ('results/true1.sim10', 'True 1 Sim 10'),
        ('results/true3.sim10', 'True 3 Sim 10'),
        ('results/true5.sim10', 'True 5 Sim 10'),
        ('results/true8.sim10', 'True 8 Sim 10'),
        ('results/true10.sim10', 'True 10 Sim 10'),
    ]
    
    nrows, ncols = 3, 5  # Adjust rows/cols as needed
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16, 8))
    axs = axs.flat
    # Define row and column titles
    row_titles = ['1 SNP Neighborhood', '5 SNP Neighborhood', '10 SNP Neighborhood']
    col_titles = ['Divergence Age 1 Yrs', 'Divergence Age 3 Yrs', 'Divergence Age 5 Yrs', 'Divergence Age 8 Yrs', 'Divergence Age 10 Yrs']
    N = 20
    ind = np.arange(N) 
    width = 0.05
    for i, (filepath, title) in enumerate(inputs):
        ax = axs[i]
        data = figs(filepath)
        print(title)
        if not data:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center')
            ax.set_title(title, fontsize=8)
            ax.axis('off')
            continue
        binned = proportions(data)
        for k in range(len(binned[1])):
        #ax.hist(data, bins=np.linspace(0, 1, 21), alpha=0.5, color='blue', edgecolor='black')
            ax.bar(binned[1][k], binned[0][k], width = .05,  color='blue', edgecolor='black', alpha=0.5, align='edge')
            #print(binned[1][k], binned[0][k])
        if i < 5:
            ax.set_title(col_titles[i], fontsize=12, pad =10)
        if i > 9:
            ax.set_xlabel('Neighborhood Recovery Rate', fontsize=10)
        if i == 0 or i == 5 or i == 10:
            ax.set_ylabel('Proportion of Total Scores', fontsize=10)
        ax.set_ylim(0, 1.0)
        #ax.grid(True)
    
    # Hide unused subplots
    for j in range(len(inputs), len(axs)):
        axs[j].axis('off')

    # Adjust layout to fit titles and subplots
    #fig.subplots_adjust(left=0.05, wspace=0.3, hspace=0.5)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.08, wspace=0.2, hspace=0.2)

    for row_idx in range(nrows):
        # Get the position of the first subplot in the row
        ax = axs[row_idx * ncols]  # First subplot in each row
        pos = ax.get_position()  # Get the Bbox (Bounding Box) of the subplot
        
        # Calculate the vertical position based on the subplot's y0 value (bottom of the subplot)
        y_pos = pos.y0 + (pos.y1 - pos.y0) / 2  # Center of the row
        
        # Place the row title at the calculated position
        fig.text(
            0.015,  # X position (adjust this for alignment)
            y_pos,  # Y position (center of the row)
            row_titles[row_idx],
            va='center', ha='left', fontsize=12, rotation='vertical'
        )
    
    plt.savefig('results/nrr_visual.svg', format='svg', dpi=600)
    

if __name__ == "__main__":
    main()



