# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2025 Swiss Tropical and Public Health Institute
# Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
# Copyright (C) 2020-2025 University of Basel
# Copyright (C) 2025 The Kids Research Institute Australia
#
# OpenMalaria is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

output_txt_path = "/home/bs/sync/code/fork/openMalaria-48.0/output.txt"

from matplotlib import pyplot as plt
import pandas as pd

mm = {
	0: 'nHost',
	1: 'nInfect',
}

mmi = {v: k for k, v in mm.items()}

df = pd.read_csv(output_txt_path, sep="\t", header=1)
df.columns = ['survey', 'age-group', 'measure', 'value']

nHosts_ageGp1 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 1)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp2 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 2)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp3 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 3)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp4 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 4)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp5 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 5)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp6 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 6)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp1001 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 1001)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp1002 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 1002)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp1003 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 1003)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp1004 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 1004)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp1005 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 1005)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts_ageGp1006 = df[(df["measure"] == mmi['nHost']) & (df['age-group'] == 1006)].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nHosts = df[(df["measure"] == mmi['nHost'])].groupby(['survey', 'measure']).sum().value[1:].reset_index()
nInfect = df[(df["measure"] == mmi['nInfect'])].groupby(['survey', 'measure']).sum().value[1:].reset_index()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(nHosts.index, nHosts.value, marker="v", markersize=1, label=f"nHosts")
ax.plot(nInfect.index, nInfect.value, marker="v", markersize=1, label=f"nInfect")
ax.plot(nHosts_ageGp1.index, nHosts_ageGp1.value, marker="v", markersize=1, label=f"nHosts_ageGp1")
ax.plot(nHosts_ageGp2.index, nHosts_ageGp2.value, marker="v", markersize=1, label=f"nHosts_ageGp2")
ax.plot(nHosts_ageGp3.index, nHosts_ageGp3.value, marker="v", markersize=1, label=f"nHosts_ageGp3")
ax.plot(nHosts_ageGp4.index, nHosts_ageGp4.value, marker="v", markersize=1, label=f"nHosts_ageGp4")
ax.plot(nHosts_ageGp5.index, nHosts_ageGp5.value, marker="v", markersize=1, label=f"nHosts_ageGp5")
ax.plot(nHosts_ageGp6.index, nHosts_ageGp6.value, marker="v", markersize=1, label=f"nHosts_ageGp6")
ax.plot(nHosts_ageGp1001.index, nHosts_ageGp1001.value, marker="v", markersize=1, label=f"nHosts_ageGp1001")
ax.plot(nHosts_ageGp1002.index, nHosts_ageGp1002.value, marker="v", markersize=1, label=f"nHosts_ageGp1002")
ax.plot(nHosts_ageGp1003.index, nHosts_ageGp1003.value, marker="v", markersize=1, label=f"nHosts_ageGp1003")
ax.plot(nHosts_ageGp1004.index, nHosts_ageGp1004.value, marker="v", markersize=1, label=f"nHosts_ageGp1004")
ax.plot(nHosts_ageGp1005.index, nHosts_ageGp1005.value, marker="v", markersize=1, label=f"nHosts_ageGp1005")
ax.plot(nHosts_ageGp1006.index, nHosts_ageGp1006.value, marker="v", markersize=1, label=f"nHosts_ageGp1006")
ax.legend(loc="upper left")

plt.show()
