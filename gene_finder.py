# -*- coding: utf-8 -*-
"""
This file contains the functions from soft des project one genefinder. It reads a string of dna and outputs proteins.

@author: MJ McMillen

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'ERROR'
    #end get_complement


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    strand_length = len(dna)   # how many letters are being decoded
    current_dna = 0  # wat letter is being reversed
    curlet= dna[current_dna]
    reversed_letter = get_complement(curlet)  # returns the reversed
    # letter of the current string characters
    reversed_dna = reversed_letter  # adds the letter to the final string
    current_dna = 1
    while current_dna < strand_length:
            curlet= dna[current_dna]
            reversed_letter = get_complement(curlet)
            reversed_dna = reversed_dna + reversed_letter
            current_dna = current_dna + 1
    return reversed_dna[::-1]
    # end get_reverse_complement function


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    SC1 = 'TAG'  # these are the stopcodons
    SC2 = 'TAA'
    SC3 = 'TGA'
    length_dna = len(dna)
    strand_check = 0 #first instance in countdown clock
    if dna[0:3] != 'ATG':
        print('Warning: no start codon') #warns that the beginning is not a start codon
    while strand_check < length_dna:
        cur_cod = dna[strand_check:(strand_check+3)] # takes a slice of 3 nucleotides in frame.
        if cur_cod == SC1 or cur_cod == SC2 or cur_cod == SC3:
            #print (endcodon)
            break
        else:
            strand_check = strand_check+3
    endcodon = strand_check
    return dna[0:endcodon]  # this returns the beginning of string:laststopcod
    # end rest_of_ORF function


def remove_duplicate(values):
    """ This is a function that I added that removed duplicate strings from a
    list. I noticed that my functions were failing doc tests because they gave
    the same piece of dna twice so I decided to write this to remove the
    duplicates while keeping the string in its origional order. The function
    moves from the end to the beginning removing any instance of repetition.

    >>> remove_duplicate(['AAA', 'BBB', 'AAA', 'CCC'])
    ['BBB', 'AAA', 'CCC']
    """
    output = []
    seen = set()
    for value in values:
        if value not in seen:
            output.append(value)
            seen.add(value)
    return (output)



def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # find the atg's
    length_dna = len(dna)
    strand_check = 0 # initializing count
    frames= ['initial'] # list filler so I dont get an error
    while strand_check < length_dna:
        current_codon = dna[strand_check:(strand_check+3)]
        if current_codon == 'ATG':
            ORF = rest_of_ORF(dna[strand_check:])
            frames.append(ORF)
            end_ORF = len(ORF)
            dna = dna[end_ORF:] #this makes sure it does not stop after one ORF or repeat an ORF
            strand_check=0
            #end of if loop
        else:
            strand_check = strand_check+3
        #end of while loop
    return remove_duplicate(frames[1:])
    #end of find_all_ORFs_oneframe function


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    Endy = find_all_ORFs_oneframe(dna)
    dna2 = dna[1:] #dna second reading frame
    Endy= Endy + find_all_ORFs_oneframe(dna2)
    dna3 = dna[2:] #dna third reading frame
    Endy = Endy + find_all_ORFs_oneframe(dna3)
    #print(type(Endy), 'type of find all ORFs')
    return remove_duplicate(Endy)
    # end of find_all_ORFs function


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse_dna = get_reverse_complement(dna)
    first_dna_ORFs = find_all_ORFs(dna)
    reverse_dna_ORFs = find_all_ORFs(reverse_dna)
    return remove_duplicate(first_dna_ORFs + reverse_dna_ORFs)
    # End of find_all_ORFs_both_strands

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_strands = find_all_ORFs_both_strands(dna)
    current_biggest = ['empty']# biggest current ORF
    biggest = 0 #length of current biggest ORF
    for s in all_strands:
        length = len(s)
        if length > biggest:
            biggest = length
            current_biggest = s
            #end of if
        elif length == biggest:
            biggest = length
            current_biggest = current_biggest, biggest
            #end elif
        #end of for
    return current_biggest
    #end of longest_ORF


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    current_biggest = 0
    n = 0
    while n <= num_trials:
        currentdna = shuffle_string(dna)
        longest_in_shuffle = longest_ORF(currentdna)
        biggest_slice = len(longest_in_shuffle)
        if  biggest_slice >= current_biggest:
            current_biggest = biggest_slice
            #end if loop
        #end for loop
        n = n + 1
    return current_biggest


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    lengthy = len(dna)
    counting = 0
    aminosequence = 'Q'
    while counting < lengthy:
        curslice = dna[counting:counting+3]
        if len(curslice)<3:
            break
            #end if statement
        amino_coded = aa_table[curslice]
        aminosequence += amino_coded
        counting = counting +3
        #end of while loop
    return aminosequence[1:]
    #end of function coding_strands_to_AA

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    gene_threshold = longest_ORF_noncoding(dna, 1500)    #how long a sequence must be to be considered a gene
    all_orfs = find_all_ORFs_both_strands(dna)
    #this finds all the ORFs in both strands of dna
    #now I need to filter the strands and throw out the ones less than gene_threshold
    final_amino_acids = ['initial']
    number_of_orfs = len(all_orfs)
    s = 0
    while s < number_of_orfs:
        current_orf = all_orfs[s]
        length_of_slice = len(current_orf)
        if length_of_slice > gene_threshold:
            #find the amino acids.
            aminos_coded = coding_strand_to_AA(current_orf)
            final_amino_acids.append(aminos_coded)
            #end of elif
        s = s + 1
    del final_amino_acids[0]
    #print(len(final_amino_acids)
    return final_amino_acids


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    from load import load_seq
    test_dna = load_seq("./data/X73525.fa")
    print(gene_finder(test_dna))
