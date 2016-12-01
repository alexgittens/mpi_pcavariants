# reshapes the thermocline EOFs returned from the PCA on the matrix to the original 3D grid
# note that the thermocline data is a submatrix of the original CSFRO matrix, so need to 
# interpret the observedLocations as offsets into the original full CSFRO 3D grid

import numpy as np
from numba import jit
import csv, sys, h5py, netCDF4

FILLVALUE = -999
OFFSET = 9*360*720 # thermocline extraction skipped the first 9 levels

def main(inSource, preprocessMethod, metadataFname, outDest):

    fin = h5py.File(inSource, "r")
    U = fin["U"][:]
    V = fin["V"][:]
    S = fin["S"][:]
    rowMeans = fin["rowMeans"][:]

    md = np.load(metadataFname)
    latGrid = np.array(md["latList"])
    lonGrid = np.array(md["lonList"])
    depths = np.array(md["depthList"])
    dates = map("".join, md["timeStamps"])
    mapToLocations = np.array(md["observedLocations"])

    # reorder time series so columns are increasing in time, left-to-right
    increasingDateIndices = np.argsort(dates)
    V = V[increasingDateIndices, :]
    dates = np.array(dates)[increasingDateIndices]

    writeEOFs(outDest, latGrid, lonGrid, depths, dates, mapToLocations,
            preprocessMethod, rowMeans, U, S, V)

def writeEOFs(outDest, latGrid, lonGrid, depths, dates, mapToLocations, preprocessMethod, rowMeans, U, S, V):

    @jit
    def toEOF(vector):
        eof = -FILL_VALUE*np.ones((len(depths), len(latGrid), len(lonGrid)), vector.dtype)
        curIndex = OFFSET
        curWriteIndexOffset = 0

        for depth in range(len(depths)):
            for lat in range(len(latGrid)):
                for lon in range(len(lonGrid)):
                    if (curWriteIndexOffset == (len(vector) + OFFSET)):
                        return eof
                    elif (curIndex == mapToLocations[curWriteIndexOffset]):
                        eof[depth, lat, lon] = vector[curWriteIndexOffset]
                        curWriteIndexOffset += 1
                    curIndex += 1
        return eof

    rootgrp = netCDF4.Dataset(outDest, "w", format="NETCDF4")

    latDim = rootgrp.createDimension("lat", len(latGrid))
    lonDim = rootgrp.createDimension("lon", len(lonGrid))
    depthDim = rootgrp.createDimension("depth", len(depths))
    eofsDim = rootgrp.createDimension("eofs", S.shape[0])
    datesDim = rootgrp.createDimension("dates", dates.shape[0])
    datelengthDim = rootgrp.createDimension("lengthofdates", len(dates[0]))

    rootgrp.createVariable("lat", latGrid.dtype, ("lat",))[:] = latGrid
    rootgrp.createVariable("lon", lonGrid.dtype, ("lon",))[:] = lonGrid
    rootgrp.createVariable("depth", depths.dtype, ("depth",))[:] = depths
    rootgrp.createVariable("coldates", 'S1', ("dates", "lengthofdates"))[:] = map(list, dates)
    rootgrp.createVariable("temporalEOFs", V.dtype, ("eofs", "dates"))[:] = V.T
    rootgrp.createVariable("singvals", S.dtype, ("eofs",))[:] = S
    rootgrp.createVariable("meanTemps", rowMeans.dtype, ("depth", "lat", "lon"), fill_value = FILLVALUE)[:] = toEOF(rowMeans)
    
    for eofNum in range(U.shape[1]):
        print "Converting EOF {0}".format(eofNum)
        curEOF = rootgrp.createVariable("EOF" + str(eofNum), U.dtype, ("depth", "lat", "lon"), fill_value = FILLVALUE)
        curEOF[:] = toEOF(U[:, eofNum])
        curEOF.type_of_statistical_processing = preprocessMethod
        curEOF.level_type = "Depth below sea level (m)"
        curEOF.grid_type = "Latitude/longitude"
    
    rootgrp.close()

inSource = sys.argv[1]
preprocessMethod = sys.argv[2]
metadataFname = sys.argv[3]
outDest = sys.argv[4]
main(inSource, preprocessMethod, metadataFname, outDest)
