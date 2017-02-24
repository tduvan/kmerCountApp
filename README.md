# kmerCountApp
A command line application that reads given FASTQ files to count the occurance and frequency of desired kmers

There are two different versions which satisfy the same task.
- kmer_count_naive.py is the earlier version of the applicaiton which is a straightforward imlemenation.
- kmer_count_bloomfilt.py is the latter version which employs a bloom filter and some other mods that makes the code more pythonic.

This work is only created for self development purposes. 

Kmer Counting Tool

Supports Python 3.x (has not been tested with older versions)

This Module accepts .TXT or .FASTQ file extensions. The tool creates kmers (word
segments within DNA with the size of k) using the user input. The created kmers
are then grouped and listed according to their frequencies.

The tool is run from the command line with the following positional arguments:

    - input file
    - kmer size
    - number of top kmers to be listed

Also following optional arguments can be used:

    --database_path
    --bloom_itemsize
    --false_pos_prob

Ex/
    Python3 kmer_naive.py input_filepath kmer_size top_count --databasepath --bloomitemsize --false_pos_prob
