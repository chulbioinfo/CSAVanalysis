# Overview
This script to run 'CSAV' analysis to idenify convergent single amino acid variants of target species which share a convergent trait you are interested in.
- - -

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

- - -

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

(Optional step to check inputs before running the script)
  - This sciprt requires following variables:
  1. input path
  2. input file format
  3. whole species list 
  4. target species list
  5. outgroup species list
  6. output path

  - Example variables in the script
<pre>
<code>
    seqPATH = "../../0.rawdata/MSA/pep/"
    seqFORM = ".sate.default.pep.removed.shortname.fasta"
    sID_LIST = ['TAEGU','GEOFO','CORBR','MELUN','NESNO','CALAN','MANVI','FALPE','CARCR','MERNU','PICPU','BUCRH','APAVI','LEPDI','COLST','TYTAL','HALLE','HALAL','CATAU','PELCR','EGRGA','NIPNI','PHACA','FULGL','PYGAD','APTFO','GAVST','PHALE','EURHE','CHAVO','BALRE','OPHHO','CHAPE','CAPCA','CHLUN','TAUER','CUCCA','MESUN','PTEGU','COLLI','PHORU','PODCR','GALGA','MELGA','ANAPL','TINMA','STRCA','HUMAN','ACACH']
    nTargets = 'TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN'
    outgroup_LIST = ['HUMAN','ACACH']
    oPATH = "../output/"
</code>
</pre>


- - -

# Running time
- Based on the demo data (Mutiple peptide sequence alignments of 141 genes with convergent single amino acid variants of vocal learning birds)
- 4 min 52.236422 sec (single thread in a standard computer with AMD Ryzen 5 2400G 3.6GHz and 16GB RAM)

# License
This project is covered under the Apache 2.0 License.

