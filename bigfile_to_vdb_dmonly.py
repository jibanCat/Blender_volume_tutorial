"""
Assume on the cluster-

Create a stand-alone env:
```
conda create -n blender_env python=3.7
conda activate blender_env
```

Packages:
```
pip install numpy pyopenvdb cython bigfile
conda install -c bccp nbodykit
```

Need to include the path to pyopenvdb library to the LD_LIBRARY_PATH
```
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.conda/envs/blender_env/lib/:$HOME/.conda/envs/blender_env/lib/python3.7/site-packages/" >> ~/.bashrc
```
Or replace $HOME/.conda to the place you installed your environment.

Note: partially adopted from the script yueying gave me and mahdi's script.
"""
from typing import Optional

import argparse

import numpy as np
import pyopenvdb as vdb
from bigfile import File
from nbodykit.lab import BigFileCatalog


def Scale_data(data):
	return data/(data.max()-data.min())


def get_nonlin_fields(inpath: str, field: int = 1, Nmesh: Optional[int] = None,
    resampler: str = "tsc") -> np.ndarray:
    """
    Parameters:
    ----
    inpath : path for input file.
    field : The numbers 1 - 5 correspond to the different types of particle in
        the simulation: 0 is gas particles (baryons) 1 is CDM, 2 is neutrinos,
        3 is unused, 4 is stars and 5 is black holes.
        Ref: MP-Gadget user manual.
    """
    field_names = ["gas", "cdm", "neutrinos", "unused", "stars", "black holes"]

    bigf = File(inpath)

    header = bigf.open('Header')
    boxsize = header.attrs['BoxSize'][0]
    redshift = 1./header.attrs['Time'][0] - 1

    print("Boxsize {}; Redshift: {:.3g}".format(boxsize, redshift))

    Ng = header.attrs['TotNumPart'][1] ** (1/3)
    Ng = int(np.rint(Ng))
    print("Particle per side {:.3g}".format(Ng))

    if Nmesh is None:
        Nmesh = 2*Ng

    cat = BigFileCatalog(inpath, dataset='{}/'.format(field), header='Header')
    print("Reading", field_names[field], "field ...")

    mesh = cat.to_mesh(Nmesh=Nmesh, resampler=resampler)

    painted_field = mesh.paint(mode="real")

    return painted_field.value # I want ndarray
    

def bigfile2vdb(inpath: str, output_name: str = "MPGadget_tutorial_fullphysics.vdb",
    Nmesh: int = 128, resampler: str = "tsc") -> None:

    data = {
        "gas" : get_nonlin_fields(inpath, field=0, Nmesh=Nmesh, resampler=resampler),
        "cdm" : get_nonlin_fields(inpath, field=1, Nmesh=Nmesh, resampler=resampler),
        "stars" : get_nonlin_fields(inpath, field=4, Nmesh=Nmesh, resampler=resampler),
        "black holes" : get_nonlin_fields(inpath, field=5, Nmesh=Nmesh, resampler=resampler),
    }

    # loop over normalization for vdb
    for key, values in data.items():
        data[key] = Scale_data(values)

    dataCube = []

    dataCube.append(vdb.FloatGrid())
    dataCube[-1].copyFromArray(data["gas"])
    dataCube[-1].name = "gas"
    dataCube[-1].transform = vdb.createLinearTransform(voxelSize=1/(Nmesh))

    dataCube.append(vdb.FloatGrid())
    dataCube[-1].copyFromArray(data["cdm"])
    dataCube[-1].name = "cdm"
    dataCube[-1].transform = vdb.createLinearTransform(voxelSize=1/(Nmesh))

    dataCube.append(vdb.FloatGrid())
    dataCube[-1].copyFromArray(data["stars"])
    dataCube[-1].name = "stars"
    dataCube[-1].transform = vdb.createLinearTransform(voxelSize=1/(Nmesh))

    dataCube.append(vdb.FloatGrid())
    dataCube[-1].copyFromArray(data["black holes"])
    dataCube[-1].name = "black holes"
    dataCube[-1].transform = vdb.createLinearTransform(voxelSize=1/(Nmesh))

    vdb.write(output_name, grids=dataCube)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--inpath", type=str, help="Path to the bigfile snapshot")
    parser.add_argument("--output_name", type=int, default="output.vdb", help="Output .vdb name")

    # input file and out file paths
    parser.add_argument("--Nmesh", type=int, help="The size of 3D painted field you want.")
    parser.add_argument("--resampler", type=str, default="tsc", help="Window mathods in Nbodykit for painting the field.")

    args = parser.parse_args()

    bigfile2vdb(args.inpath, args.output_name, args.Nmesh, args.resampler)
