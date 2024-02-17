import sys
import argparse
from msilib import sequence


def count(sequence):
    a_total = 0
    c_total = 0
    g_total = 0
    t_total = 0

    for letter in sequence:
        if letter == 'a':
            a_total += 1
        elif letter == 'c':
            c_total += 1
        elif letter == 't':
            t_total += 1
        elif letter == 'g':
            g_total += 1

    return a_total, c_total, t_total, g_total


def complement(sequence, reverse=False):
    complementary_sequence = ""
    for s in sequence:
        if s.lower() == 'a':
            complementary_sequence += 't'
        elif s.lower() == 't':
            complementary_sequence += 'a'
        elif s.lower() == 'g':
            complementary_sequence += 'c'
        elif s.lower() == 'c':
            complementary_sequence += 'g'

    if reverse:
        new_sequence = reverse_sequence(complementary_sequence)
        return new_sequence
    else:
        return complementary_sequence
def reverse_sequence(sequence):
    new_sequence = sequence[::-1]
    return new_sequence

def calculate_gc_content(sequence):
    a_total, c_total, t_total, g_total = count(sequence)
    x = len(sequence)
    y = c_total + g_total
    gc_content = (y/x)*100
    gc_threshold = gc_mean(gc_content)
    gc_content = round(gc_content, 2)
    return gc_content, gc_threshold

def gc_mean(gc_content):
    gc_threshold = gc_content*1.25
    return gc_threshold


def detect_gc_islands(sequence, gc_threshold, window_size):
    counter = 0
    cpg_islands = []

    while counter + window_size <= len(sequence):
        cpg_window = sequence[counter: counter + window_size]
        current_gc_content = calculate_gc_content(cpg_window)

        if current_gc_content[0] > gc_threshold:
            cpg_islands.append((counter, counter + window_size, current_gc_content))

        counter += 1

    return cpg_islands

def main():
    args = parse_arguments()
    files = args.input
    results = []
    for file in files:

        if len(files) > 3:
            print("error")
        else:
            with open(file, 'r') as f:
                sequence = f.read()
                sequence = sequence.lower()


            if args.basecount:
                a_total, c_total, t_total, g_total = count(sequence)
                results.append(f"A: {a_total}, C: {c_total}, T: {t_total}, G: {g_total}\n")
            if args.reversecomplement:
                complementary_sequence = complement(sequence, reverse=True)
                results.append(f"Reverse complement of {file}:\n{complementary_sequence}\n")
            if args.GCcontent:
                gc_content = calculate_gc_content(sequence)
                results.append(f"The percentage of GC in {file} is {gc_content[0]}\n")
            if args.numberofislands:
                gc_content, gc_threshold = calculate_gc_content(sequence)
                regions = detect_gc_islands(sequence, gc_threshold, window_size=200)
                for start, end, _ in regions:
                    results.append(f"Start: {start}, End: {end}\n")

                results.append(f"No other regions in {file}\n")

            if args.output:
                output_file = 'results.txt'
                with open(output_file, 'w') as output:
                    output.writelines(results)
            else:
                print(results)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate some statistics about DNA sequences.")
    parser.add_argument('--input', nargs='+', required=True, help='Input file names (space-separated)')
    parser.add_argument('--output', action='store_true')
    parser.add_argument('--basecount', action='store_true')
    parser.add_argument('--reversecomplement', action='store_true')
    parser.add_argument('--GCcontent', action='store_true')
    parser.add_argument('--numberofislands', action='store_true')
    return parser.parse_args()

if __name__ == '__main__':
    main()

def test_rev_comp():
    assert complement('ATCG')=='tagc', "complement test"
    assert complement('ATCG', True)=='cgat', "reverse complement test"
    print("Tests passed")


def test_gc_content():
    assert calculate_gc_content('ggggaaaaaaaatttatatatcgcc')[0]==32, "gc_content test"
    print("Tests passed")
