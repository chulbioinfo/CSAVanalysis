# Overview
This script to run 'CSAV' analysis compares the position of CSAV in the multiple sequence alignment (MSA) file with the known SNPs in the chromosomal position. The script first converts the CSAV position in the MSA file to chromosomal position based on the BLAST output file and the chromosomal positions of CSAVs are then compared with population variant dataset (vcf file downloaded from Ensembl) to check if CSAV is conserved within the population data. The script here is unable to run on different dataset.

# System Requirments
## Hardware requirements
This script requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This script is supported for *Windows OS* and *Linux*. The script has been tested on the following systems:
* Windows 10 Education
* Linux: Ubuntu 18.04.3 LTS

## Python Dependencies
This script mainly depends on the Python scientific stack.

    numpy
    scipy


# Running the script
* (option1) If you want to clone the github link:
<pre>
<code>
git clone https://github.com/chulbioinfo/CSAVanalysis/
cd CSAVanalysis/
</code>
</pre>

* (option2) If you use the zip file from the github link:
<pre>
<code>
wget https://github.com/chulbioinfo/CSAVanalysis/archive/master.zip
unzip master.zip
cd CSAVanalysis-master/
</code>
</pre>

* Decompress the demo data of input files (as [fasta](https://en.wikipedia.org/wiki/FASTA_format) format)
<pre>
<code>
cd 0.rawdata/MSA/
unzip MSA.zip
cd ../../
</code>
</pre>

* Run the script
<pre>
<code>
cd 1.CSAV/bin/
python CSAV.py
</code>
</pre>

# License
This project is covered under the Apache 2.0 License.

