# CROMq
CROMq is a lossy compression algorithm for quality values.

## Installing
CROMq can only be used on Linux.
CROMq requires fftw3. To download, please check [here](http://www.fftw.org/download.html).
The following is an example command of installation.

```
./configure --enable-threads --enable-shared --enable-sse2
sudo make
sudo make install
```

## Example data
The example data that we used is [SRR494099](https://www.ebi.ac.uk/ena/data/view/SRR494099).
It is a fastq format file, so we need to extract qscores from it.
The following example is the simple command to extract qscores from fastq file.

```
awk NR%4==0 SRR494099.fastq >> SRR494099.qscore
```

Note that the file must end with `.qscore`.

## Usage
Please check [here](example/README.md) to see how to execute the example code.
