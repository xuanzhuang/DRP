#!/usr/bin/python

import numpy as np
import sys
import matplotlib.pyplot as plt
import csv
import os.path
from collections import defaultdict
from scipy.optimize import curve_fit

# outputs
# 1. one figure (two lines)
# 2. when y = 0.5, the difference of HSD-LSD in X
# 3. when y = 0.5, the x of HSD and LSD for all samples

# configurable inputs
input_filename_prefix = "line"
input_file_number = 64
output_image_filename_prefix = "plot_"
output_filename = "sigmod_difference.csv"
output_filename_2 = "sigmod_summary.csv"
polyfunc_fit_degree = 3
y_value = 0.5
fig_x_label = "hours since egg deposition"
fig_y_label = "% pariated"
LSD_label = "LSD"
HSD_label = "HSD"
x_lim_min = 0
x_lim_max = 350
x_lim_min = 0
x_lim_max = 1

def sigmoid(x, x0, k):
   y = 1 / (1 + np.exp(-k*(x-x0)))
   return y

out_file = open(output_filename, "a")
out_file_2 = open(output_filename_2, "a")

# input files
LSD_input_file = sys.argv[1]
HSD_input_file = sys.argv[2]

fname_out_png = output_image_filename_prefix + LSD_input_file.split('_')[0] + ".png"

# open file
if (not os.path.isfile(LSD_input_file)):
    print(LSD_input_file + " does not exist.")
    sys.exit()

if (not os.path.isfile(LSD_input_file)):
    print(LSD_input_file + " does not exist.")
    sys.exit()

print("Process " + LSD_input_file + " and " + HSD_input_file)

# read file
columns = defaultdict(list) # each value in each column is appended to a list
with open(LSD_input_file) as f:
	reader = csv.DictReader(f) 		# read rows into a dictionary format
	for row in reader: 			    # read a row as {column1: value1, column2: value2,...}
		for (k,v) in row.items(): 	# go over each column name and value
			columns[k].append(v) 	# append the value into the appropriate list
						            # based on column name k

with open(HSD_input_file) as f:
	reader = csv.DictReader(f)
	for row in reader:
		for (k,v) in row.items():
			columns[k].append(v)

x1=[]
y1=[]
err1=[]
x2=[]
y2=[]
err2=[]

for j in range(0, len(columns['LSD_hrs'])):
	x1.append(float(columns['LSD_hrs'][j]))
	y1.append(float(columns['LSD'][j]))
	if columns['LSD'][j] == "1.0":
		break

for j in range(0, len(columns['HSD_hrs'])):
        x2.append(float(columns['HSD_hrs'][j]))
        y2.append(float(columns['HSD'][j]))
        if columns['HSD'][j] == "1.0":
                break

p1_x0 = (x1[0] + x1[-1])/2
p1_k = 12/(x1[-1] - x1[0])
popt1, pcov1 = curve_fit(sigmoid, x1, y1, p0=[p1_x0, p1_k])

p2_x0 = (x2[0] + x2[-1])/2
p2_k = 12/(x2[-1] - x2[0])
popt2, pcov2 = curve_fit(sigmoid, x2, y2, p0=[p2_x0, p2_k])

# calculate new x's and y's, for plotting
xs1 = np.arange(min(x1), max(x1), 0.5)
ys1 = sigmoid(xs1, *popt1)
xs2 = np.arange(min(x2), max(x2), 0.5)
ys2 = sigmoid(xs2, *popt2)

# solve x for y = 0.5
new_x1 = popt1[0] - np.log(1/y_value - 1)/popt1[1]
new_x2 = popt2[0] - np.log(1/y_value - 1)/popt2[1]

print("new_x1 = " + str(new_x1))
print("new_x2 = " + str(new_x2))

diff = new_x2 - new_x1
print("diff = " + str(diff))

temp = LSD_input_file.split('_')[0].split('L')
text_diff = temp[0]+"_"+temp[1]+","+str(diff)+"\n"
out_file.write(text_diff)

text_value = temp[0]+"_"+temp[1]+","+str(new_x1)+","+str(new_x2)+"\n"
out_file_2.write(text_value)

# plot and save as png file
plt.plot(x1, y1, 'o', color='k')
plt.plot(xs1, ys1, color='b', label=LSD_label)
plt.plot(x2, y2, 'x', color='k')
plt.plot(xs2, ys2, color='r', label=HSD_label)
plt.legend(loc='best')
plt.ylabel(fig_y_label)
plt.xlabel(fig_x_label)
plt.grid()
print("output to " + fname_out_png)
plt.savefig(fname_out_png)
plt.clf()

out_file.close()
out_file_2.close()
