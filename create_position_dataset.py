#!/bin/env python

import numpy as np
from numpy import save
import sys
from random import seed
from random import randint
from Bio import SeqIO
import argparse


def create_dataset(negative_file, positive_file,output,total_win_len):
    print("loading files")
    pos_seqs = [x for x in SeqIO.parse(positive_file, 'fasta')]
    neg_seqs = [x for x in SeqIO.parse(negative_file, 'fasta')]
    print("done")
    # total_win_len = 50000

    dataset = np.zeros((len(pos_seqs)*2, 5, total_win_len), dtype=np.int16)
    labels = np.zeros((len(pos_seqs)*2, 2), dtype=np.int32)
    seed(1)
    j = 0
    for i in range(len(pos_seqs)):
        initial_pos = randint(0, total_win_len - len(str(pos_seqs[i].seq)) - 1)

        # to fill the right size of the 50 kb window
        filling_seq_i = randint(0, len(neg_seqs)-1)
        filling_seq = str(neg_seqs[filling_seq_i].seq)
        while len(str(filling_seq)) < initial_pos:
            filling_seq_i = randint(0, len(neg_seqs)-1)
            filling_seq += str(neg_seqs[filling_seq_i].seq)
        final_seq = filling_seq[0:initial_pos]+str(pos_seqs[i].seq)

        # to fill the rigth size of the 50 kb window
        filling_seq_i = randint(0, len(neg_seqs)-1)
        filling_seq = str(neg_seqs[filling_seq_i].seq)
        while len(str(filling_seq)) < total_win_len - len(final_seq):
            filling_seq_i = randint(0, len(neg_seqs)-1)
            filling_seq += str(neg_seqs[filling_seq_i].seq)
        final_seq = final_seq + filling_seq[0:total_win_len - len(final_seq)]
        dataset[j] = one_hot(final_seq)
        labels[j, 0] = initial_pos
        labels[j, 1] = initial_pos + len(pos_seqs[i].seq)
        if labels[j, 0] < 0 or labels[j, 1] < 0:
            print(initial_pos)
            print(initial_pos + len(pos_seqs[i].seq))
            print("----------------------------")
            print(labels[j, 0])
            print(labels[j, 1])
            print(len(pos_seqs[i]))
        j += 1

        # to add negative instances
        filling_seq_i = randint(0, len(neg_seqs)-1)
        filling_seq = str(neg_seqs[filling_seq_i].seq)
        while len(str(filling_seq)) < total_win_len:
            filling_seq_i = randint(0, len(neg_seqs)-1)
            filling_seq += str(neg_seqs[filling_seq_i].seq)
        final_seq = filling_seq[0:total_win_len]
        dataset[j] = one_hot(final_seq)
        labels[j, 0] = 0
        labels[j, 1] = 0
        j += 1

        if i % total_win_len == 0:
            print("doing "+str(i))

    # print(dataset)
    # print(labels)
    print(dataset.shape)
    print(labels.shape)
    if output[-1] == '/':
        np.save(output + "dataset.npy", dataset)
        np.save(output + "dataset_labels.npy", labels)
    else:
        np.save(output + "/" + "dataset.npy", dataset)
        np.save(output + "/" + "dataset_labels.npy", labels)


def one_hot(sequence):
    langu = ['A', 'C', 'G', 'T', 'N']
    rep2d = np.zeros((1, 5, len(sequence)), dtype=np.int16)
    posNucl = 0
    for nucl in sequence:
        posLang = langu.index(nucl.upper())
        rep2d[0][posLang][posNucl] = 1
        posNucl += 1
    return rep2d

def filter(file):
    newFile = open(file+".filtered", "w")
    for te in SeqIO.parse(file, "fasta"):
        seq = str(te.seq)
        filterDna = [x for x in seq if x.upper() in ['A', 'C', 'G', 'T', 'N']]
        newSeq = "".join(filterDna)
        newFile.write(">"+str(te.id)+"\n"+newSeq+"\n")

if __name__ == '__main__':

    # Define parameter
    parser = argparse.ArgumentParser()

    parser.add_argument('--neg',type=str,required=True,help="negative file")
    parser.add_argument('--pos',type=str,required=True,help='positive file')
    parser.add_argument('--output',type=str,required=True,help="output path")
    parser.add_argument('--total_win_len',type=int,default=50000,help="sequence of length")

    arg = parser.parse_args()


    # Read parameter
    total_win_len = arg.total_win_len
    output = arg.output
    pos_file = arg.pos
    neg_file = arg.neg

    # Uniform character case
    filter(pos_file)
    filter(neg_file)

    create_dataset(neg_file+'.filtered', pos_file+'.filtered',output,total_win_len)