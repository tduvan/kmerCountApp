# This code contains the answer for the coding challenge question presented by
# SEVEN BRIDGES GENOMICS Biyoteknoloji A.Ş for the second phase of the recruitment
# process.

"""Kmer Counting Tool

Supports Python 3.x (has not been tested with older versions)

This Module accepts .TXT or .FASTQ file extensions. The tool creates kmers (word
segments within DNA with the size of k) using the user input. The created kmers
are then grouped and listed according to their frequencies.

The tool is run from the command line with the following positional arguments:

    - input file
    - kmer size
    - number of top kmers to be listed

Ex/
    Python3 kmer_naive.py input_filepath kmer_size top_count
"""


from collections import defaultdict
from operator import itemgetter
import argparse
import sqlite3
import os


def database_init(directory):
    """ The function initializes the connection to the database file
    """
    return sqlite3.connect(directory)


def database_get(directory, rows):
    """ Returns the requested top rows from the database file
    """
    cursor = database_init(directory).cursor()
    cursor.execute("SELECT * FROM Kmers ORDER BY Kcount DESC")
    read = (cursor.fetchmany(rows))
    database_init(directory).close()
    return read


def database_insert(directory, dictionary):
    """ Inserts the dictionary items into the database file. An existence check is
    made for each insert operation, if the key exists the corresponding value is updated.
    """
    session = database_init(directory)
    session.execute("CREATE TABLE IF NOT EXISTS Kmers (Kmer TEXT PRIMARY KEY NOT NULL, Kcount INTEGER NOT NULL)")

    for name, value in dictionary.items():
        if value > 60:  # semi-arbitrary value for filtering out the kmers with low frequency deduced by trial and error
            try:
                session.execute("INSERT OR ABORT INTO Kmers (Kmer, Kcount) VALUES ('%s', '%d')" % (name, value))

            except sqlite3.IntegrityError:
                cursor = session.cursor()
                cursor.execute("SELECT * FROM Kmers WHERE Kmer = '%s'" % name)
                value += int(cursor.fetchone()[1])
                session.execute("UPDATE Kmers SET Kcount = '%d' WHERE Kmer = '%s'" % (value, name))

    session.commit()
    session.close()


def main():

    parser = argparse.ArgumentParser()  # parsing inputs from the command line
    parser.add_argument("filename", help="Destination of the input FASTQ file")
    parser.add_argument("kmer_size", help="The size of the kmers to search")
    parser.add_argument("top_count", help="The number of most frequent kmers to show")

    args = parser.parse_args()
    filename = args.filename
    k_size = int(args.kmer_size)    # casting the input string to integer
    top_c = int(args.top_count)     # not implemented yet

    db_dir = "temp.db"              # filename for the database file
    kmers = defaultdict(int)

    partition_size = 25000000  # size limit for the kmer dictionary, 100 MB memory usage ≃ len(kmers):3000000 lines
                               # len() method is used instead of getsizeof() method -> works faster
    with open(filename, 'r') as f:
        for index, line in enumerate(f):

            while line.startswith("@"):  # If string starts with "@" and shorter than the sequence length,
                                         # definitely new sequence
                if index == 0 or len(line) < len(sequence):
                    sequence = f.__next__()
                    f.__next__()
                    f.__next__()
                    break

                else:
                    line = f.__next__()

            for segment in range(len(sequence) - k_size + 1):
                kmer = sequence[segment:segment + k_size]

                if len(kmers) < partition_size:  # checking the RAM size in use, if a certain threshold is reached,
                    kmers[kmer] += 1             # the dictionary is saved to database, else continue.

                else:
                    kmers[kmer] += 1             # The dictionary is dumped into the database file.
                    database_insert(db_dir, kmers)
                    kmers.clear()

        if os.path.exists(db_dir):               # if no database file created (where len(kmers) < partition_size )
                                                 # print kmer dictionary else get values from DB and print
            database_insert(db_dir, kmers)
            print(database_get(db_dir, top_c))
            os.remove(db_dir)                    # this part should be commented out for keeping the created database
        else:
            # database_insert(db_dir,kmers)         # this part should be kept commented out to save dictionary to database
            kmers = sorted({key: value for key, value in kmers.items() if value > 1}.items(), reverse=True, key=itemgetter(1))
            print(kmers[:top_c])


if __name__ == '__main__':
    main()
