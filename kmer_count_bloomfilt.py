"""Kmer Counting Tool

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

The algorithm includes two parameters that are not open to direct alteration by the user:

    - Dictionary Output Size (10 Million) : This value defines the number of items to be retrieved from the database for
    merging with the dictionary in current iteration

    - Dictionary Size on Memory (50000000) : This value defines the max length of the dictionary that can be used before
    dumping to the database. The length method is used instead of the getsizeof method for performance purposes.
    (The correlation between dictionary length and its  memory usage is roughly 100 MB on RAM ≃ 3000000 items in dictionary)

The values can be altered for getting better performance on different system configurations. If the RAM capacity is low,
the second parameter should be kept as small as possible. This will allow the code work on low RAM capacity while comp_
romising process speed.

========================================================================================================================

Disclaimer:

    The bloom filter class is directly used from a version that was found on the internet with minor modifications.
    Detailed information is included in the class explanation.
"""

import os
import sqlite3
import mmh3
import argparse
from collections import defaultdict
from bitstring import BitArray
from math import log, pow, ceil


class BloomFilter(object):
    """
        A generic Bloom Filter implementation.
        A Bloom filter based on Max Burstein's implementation (http://maxburstein.com/blog/creating-a-simple-bloom-fil_ter/)
        Code can be found on https: // gist.github.com / DoggettCK / 4755641
    """

    def __init__(self, size, hash_count):                       # filter class initialization
        self.size = size
        self.hash_count = hash_count
        self.bit_array = BitArray(size)
        self.bit_array.set(0)

    def add(self, string):                                      # add method for the filter
        for seed in range(self.hash_count):
            result = mmh3.hash(string, seed) % self.size
            self.bit_array[result] = 1

    def lookup(self, string):                                   # lookup method for the filter
        for seed in range(self.hash_count):
            result = mmh3.hash(string, seed) % self.size
            if self.bit_array[result] == 0:
                return False
        return True

    @staticmethod
    def suggest_sizes(n, p):
        """Given an expected number of items and probability of false positives, suggests an appropriate size and hash count for a Bloom filter.
    n - expected number of items in filter
    p - probability of false positives (0.0 - 1.0)
    http://en.wikipedia.org/wiki/Bloom_filter#Probability_of_false_positives"""

        if not (0.0 <= p <= 1.0):
            raise ValueError("False probability percentage must be between 0.0 and 1.0")
        if n <= 0:
            raise ValueError("Number of items must be greater than 0")

        l2 = log(2)
        m = -n * log(p) / pow(l2, 2)
        k = (m / n) * l2

        return int(ceil(m)), int(ceil(k))

    @staticmethod
    def create_suggested(n, p):
        """Given an expected number of items and probability of false positives, returns a BloomFilter of appropriate size and hash count.
    n - expected number of items in filter
    p - probability of false positives (0.0 - 1.0)"""

        suggested_sizes = BloomFilter.suggest_sizes(n, p)
        return BloomFilter(suggested_sizes[0], suggested_sizes[1])


class DataBaseConnection(object):
    """
        Database class for handling the SQLite connections
    """

    def __init__(self, file_path):                              # Database Connection manager class initialization
        self.path = file_path
        self.connection = sqlite3.connect(file_path)
        self.cursor = self.connection.cursor()
        self.connect()
        self.create_tbl()

    def connect(self):                                          # method for establishing connection to the database
        return self.connection

    def create_tbl(self):                                       # method for creating a tables in the database
        conn = self.connect()
        conn.execute("CREATE TABLE IF NOT EXISTS Kmers (Kmer TEXT PRIMARY KEY NOT NULL, Kcount INTEGER NOT NULL)")

    def insert(self, dictionary):                               # method for inserting the kmer dictionaries to database
        try:
            for segment in self.get(600000):                  # Gets the first 10M entries in the database (remaining
                dictionary[segment[0]] += segment[1]                #are disregarded)
            self.clear()

        except sqlite3.OperationalError:
            self.create_tbl()

        finally:
            conn = self.connect()
            conn.executemany("INSERT OR IGNORE INTO Kmers (Kmer, Kcount) VALUES (?,?)",
                             ((key, value) for key, value in dictionary.items() if value > 1))  # ignoring unique values
            conn.commit()

    def get(self, rows):                                        # get method for pulling data from database
        return self.cursor.execute("SELECT * FROM Kmers ORDER BY Kcount DESC").fetchmany(rows)

    def clear(self):
        self.connection.execute("DELETE FROM Kmers")            # clears the table items


def save_kmers(segment, dictionary, database_obj):              # buffer function for kmer dictionary in memory
    if len(dictionary) < 25000000:                                  # size limit for the kmer dictionary, 100 MB memory
        dictionary[segment] += 1                                    # usage ≃ len(kmers):3000000 lines
                                                                    # len() method is used instead of getsizeof()method
    else:                                                           # since it works faster
        dictionary[segment] += 1
        database_obj.insert(dictionary)
        dictionary.clear()


def main():

    parser = argparse.ArgumentParser()                              # parsing inputs from the command line
    parser.add_argument("filename", help="Destination of the input FASTQ file")
    parser.add_argument("kmer_size", help="The size of the kmers to search")
    parser.add_argument("top_count", help="The number of most frequent kmers to show")
    parser.add_argument("--database_path", help="Save path for database (default = KmerCount.db)", default="Kmer_count.db")
    parser.add_argument("--bloom_itemsize", help="Estimated item size for the filter (default = 1Billion)",
                        default=10000000)
    parser.add_argument("--false_pos_prob", help="desired probability of false positives (default = 1)", default=1)

    args = parser.parse_args()
    filename = args.filename
    k = int(args.kmer_size)
    top_c = int(args.top_count)
    db_path = args.database_path
    bloomsize = int(args.bloom_itemsize)
    fp_prob = float(args.false_pos_prob)

    kmers = defaultdict(int)

    bf = BloomFilter.create_suggested(bloomsize, fp_prob)           # Bloom filter class instantiation
    db = DataBaseConnection(db_path)                                # Database Connection Manager class instantiation

    with open(filename, 'r') as f:
        for index, line in enumerate(f):

            while line.startswith("@"):
                if index == 0 or len(line) < len(sequence):
                    sequence = f.__next__().rstrip()
                    f.__next__().rstrip()
                    f.__next__().rstrip()
                    break

                else:
                    line = f.__next__().rstrip()

            for segment in range(len(sequence) - k + 1):
                motif = sequence[segment:segment + k]
                if not bf.lookup(motif):
                    bf.add(motif)
                else:
                    save_kmers(motif, kmers, db)

    db.insert(kmers)
    print(db.get(top_c), "\n")
    print("Results are saved to : %s" % db_path, "\n")

    db.connect().commit()
    db.connect().close()

    # os.remove(db.path)                                            # leave commented out for saving the file on disk

if __name__ == '__main__':
    main()

    # TODO: Check input file extensions
    # TODO: Add parallelism to the code
