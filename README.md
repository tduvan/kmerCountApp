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

========================================================================================================================
kmer_count_bloomfilt.py

The algorithm includes two parameters that are not open to direct alteration by the user:
(HENCE THE PARAMETERS CAN BE MODIFIED WITHIN THE CODE FOR FURTHER OPTIMIZING THE APPLICATION)

    - Dictionary Output Size (600000) : This value sets the number of items to be retrieved from the database when
    merging with the dictionary in each iteration the dictionary fills up.
    
    - Dictionary Size on Memory (25000000) : This value defines the max length of the dictionary that can be used before
    dumping to the database. The length method is used instead of the getsizeof method for performance purposes.
    (The correlation between dictionary length and its  memory usage is roughly 100 MB on RAM ≃ 3000000 items in dictionary)
    
The values can be altered for getting better performance on different system configurations. If the RAM capacity is low,
the second parameter should be kept as small as possible. This will allow the code work on low RAM capacity while compromising process speed.


