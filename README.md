# Tigops 

Manipulation of unitigs and contigs. Relatively old code from 2014 but might still be useful for some tasks today:

 - compute coverage from a dataset of reads

 - some other hidden features in the source code

# Usage

    tigops coverage -kmer-size [k value] -tigs [unitigs file]  -reads [reads file]  -out [output file] -name [some_id]

What it does it return the average read k-mer abundance for each sequence of the unitig file. The output will look like:

    >[original unitig header] _cov_[coverage value]_ID_[some_id]
    [sequence]

you can launch it as many times as the number of reads files you can, each time giving the previous output as input, it will append a "_cov_xxx_ID_yyy" at the end of the fasta header.

 
# Project build

For building your project, you should do the following
    
    cd tools/tigops
    mkdir build
    cd build
    cmake ..
    make -j 8

    
Then, you should get a binary holding the name of the project.

Note: the first compilation should take some time since the GATB-CORE library is generated.


Computing coverage for shorter kmers than unitig length
----

Tigops considers only kmers that are unique to the unitig.


Increase supported kmer lengths
-----

Increase the following value in tigops.hpp: 

	#define span KMER_SPAN(2)
