import netCDF4
import matplotlib.pyplot as plt
import datetime
import numpy as np
import sys

inSource = sys.argv[1]
outBaseName = sys.argv[2]

din = netCDF4.Dataset(inSource, "r")
timeSeries = din["temporalEOFs"]
dates = map("".join, din["coldates"])

for pcnum in np.arange(timeSeries.shape[0]):
    datetimeLabels = map(lambda datestr: datetime.datetime.strptime(datestr, "%Y%m%d%H"), dates)
    plt.plot(datetimeLabels, timeSeries[pcnum, :])
    plt.axes().set_aspect(.3*1./plt.axes().get_data_ratio())
    plt.xlabel('Date')
    plt.title("time series for eof {0}".format(pcnum + 1))
    plt.savefig("{0}TimeSeries{1}.pdf".format(outBaseName, pcnum + 1))
    plt.clf()
    print("Wrote time series {0}/{1}".format(pcnum + 1, timeSeries.shape[0]))

