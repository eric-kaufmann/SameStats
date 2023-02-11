import csv
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str) 
args = parser.parse_args()

filename = str(args.filename)

# read the CSV file and create a list of (x, y) coordinates
best = []
with open(args.filename+'_generated_data.csv') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        best.append((float(row[0]), float(row[1])))

target = []
with open(args.filename+'_target_data.csv') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        target.append((float(row[0]), float(row[1])))

# create the scatter plot using matplotlib
x, y = zip(*target)
plt.scatter(x, y, c='r')

x, y = zip(*best)
plt.scatter(x, y, c='b')

# save the plot to a file
print('saving scatter plot...')
plt.savefig(args.filename+'_scatter_plot.png')