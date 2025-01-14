#!/usr/bin/python3

import sys
import re
import concurrent.futures

CHUNK_SIZE = 10000  # Number of lines to process at once. Adjust as necessary.

def process_chunk(chunk1, chunk2, flag1, flag2):
    output1, output2 = [], []

    for i in range(0, len(chunk1), 4):
        header1, seq1, strand1, quality1 = chunk1[i:i+4]
        header2, seq2, strand2, quality2 = chunk2[i:i+4]

        match1 = re.search(flag1 + r'[ATC]', seq1)
        match2 = re.search(flag2 + r'[AGC]', seq2)

        if match1 and match2:
            cut1 = match1.end() - 1
            cut2 = match2.end() - 1

            new_seq1 = seq1[cut1:]
            new_seq2 = seq2[cut2:]
            new_quality1 = quality1[cut1:]
            new_quality2 = quality2[cut2:]

            output1.extend([header1, new_seq1, strand1, new_quality1])
            output2.extend([header2, new_seq2, strand2, new_quality2])

    return output1, output2

def main(file1, file2, flag1, flag2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2, \
         open(f"{file1}.trim", "w") as out1, open(f"{file2}.trim", "w") as out2:

        chunk1, chunk2 = [], []
        futures = []

        with concurrent.futures.ThreadPoolExecutor() as executor:
            while True:
                line1, line2 = f1.readline(), f2.readline()
                if not line1 or not line2:
                    if chunk1 and chunk2:
                        futures.append(executor.submit(process_chunk, chunk1, chunk2, flag1, flag2))
                    break
                chunk1.append(line1.strip())
                chunk2.append(line2.strip())

                if len(chunk1) >= CHUNK_SIZE:
                    futures.append(executor.submit(process_chunk, chunk1, chunk2, flag1, flag2))
                    chunk1, chunk2 = [], []

            for future in concurrent.futures.as_completed(futures):
                output1, output2 = future.result()
                out1.write('\n'.join(output1) + '\n')
                out2.write('\n'.join(output2) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <file1> <file2> <Flag1> <Flag2>")
        sys.exit(1)

    file1, file2, flag1, flag2 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    main(file1, file2, flag1, flag2)
