import csv
import matplotlib.pyplot as plt

# read the CSV file and create a list of (x, y) coordinates
best = []
with open('best_matrix.csv') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        best.append((float(row[0]), float(row[1])))

target = []
with open('target_matrix.csv') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        target.append((float(row[0]), float(row[1])))

# create the scatter plot using matplotlib
x, y = zip(*best)
plt.scatter(x, y, c='b')

x, y = zip(*target)
plt.scatter(x, y, c='r')

# save the plot to a file
print('saving scatter plot...')
plt.savefig('scatter_plot.png')