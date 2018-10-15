import h5py
import os
import sys

from .Folder import Folder

class OrbReader(object):
    def __init__(self, filePath):
        self._filePath = filePath
        self._h5Handle = h5py.File(self._filePath, 'r')
        self._folders = {}

        # attrs

        # datasets

        #populate tree with folder names
        for folderName in self._h5Handle.keys():
            try:
                self._folders[folderName] = Folder(folderName,
                                                   self._h5Handle.get(folderName))
            except KeyError:
                pass



    def Print(self):
        for folder in self._folders.keys():
            print(folder)
            self._folders[folder].Print("-")


class OrbReader2(object):
    def __init__(self, fileName):
        self._h5Handle = h5py.File(fileName, "r")
        self._rootFolder = Folder("/", self._h5Handle)


    def Print(self):
        self._rootFolder.Print()
