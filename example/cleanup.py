# cleanup.py

import os

from glob import glob
from subprocess import call


def main():
    file_formats = ['bin', 'o', 'txt', 'subqscores', 'log']
    for file_format in file_formats:
        fnames = glob('./*.'+file_format)
        for fname in fnames:
            os.remove(fname)


if __name__ == "__main__":
    main()
