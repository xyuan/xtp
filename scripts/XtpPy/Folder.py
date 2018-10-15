import h5py
import os
import sys
import numpy as np

class Folder(object):
    def __init__(self, name, h5Handle):
        self._name = name
        self._h5Handle = h5Handle
        self._params = []
        self._data = []
        self._folders = {}

        # populated parameter names
        # attrs

        # Populate data set names
        # datasets


        # populate folders
        for folderName in self._h5Handle.keys():
            try:
                self._folders[folderName] = Folder(folderName,
                                                   self._h5Handle.get(folderName))
            except AttributeError:
                pass


    def __getitem__(self, path):
        path = path.split("/")[0]
        name = path[0]
        s = ""
        for p in path[1:]:
            s += p

        _folders[name][s]

        if name in data:
            dset = h5Handle.get(name)
            temp = np.zeros_like(dset.shape)
            dset.read_direct(temp, temp)
            return temp

    def Print(self, offset=""):
        for folderName in self._folders.keys():
            print(offset + folderName)
            self._folders[folderName].Print(offset+"-")
