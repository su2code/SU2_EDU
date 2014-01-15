#!/usr/bin/env python

## \file plot_pressure.py
#  \brief python package for gradients
#  \author Thomas D. Economon, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 3.0.0 "eagle"
#
# SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, time
import pandas as pd
from optparse import OptionParser
from numpy import *
from matplotlib import pyplot as plt
from matplotlib import mlab

parser = OptionParser()
parser.add_option("-f", "--file", dest="file",
                  help="surface flow csv file", metavar="FILE")
(options, args)=parser.parse_args()

# Store the file name

filename = options.file

# Load the csv file with the airfoil coordinate and pressure data (sorted)

data = pd.read_csv(filename)
data = data.sort('Global_Index')
data.to_csv(filename, index=False)
data = mlab.csv2rec(filename, comments='#', skiprows=0, checkrows=0)

# Plot the airfoil shape and pressure distribution

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(data.x_coord[:],data.pressure_coefficient[:],'-b',linewidth = 2.0)
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_xlabel(r'$x/c$', fontsize=20)
ax1.set_ylabel(r'$C_p$', fontsize=20)

ax2 = ax1.twinx()
ax2.plot(data.x_coord[:],data.y_coord[:],'-k',linewidth = 1.5)
ax2.axis('equal')
ax2.axis('off')
ax2.set_ylim([-0.1,0.7])
ax2.set_ylabel(r'$y / c$', fontsize=20)
plt.savefig('pressure_distribution.png',format='png')
