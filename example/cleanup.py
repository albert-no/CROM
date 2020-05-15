# cleanup.py

import os
import shutil

from glob import glob
from subprocess import call


def main():
    file_formats = ['bin', 'o', 'txt', 'subqscores', 'log']
    for file_format in file_formats:
        fnames = glob('./*.'+file_format)
        for fname in fnames:
            os.remove(fname)

    folders = ['bin', 'svd_params', 'logs']
    for folder in folders:
        try:
            shutil.rmtree(folder)
        except:
            print(folder + " does not exsits")


if __name__ == "__main__":
    main()
