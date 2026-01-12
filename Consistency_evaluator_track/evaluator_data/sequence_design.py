###Designing sequences for the consistency evaluator
#May 13th
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pyfaidx
import pysam
import random
import numpy as np
random.seed(10)
from Bio import SeqIO
from Bio.Seq import Seq
def create_fasta():
    accessbile_peaks_ipsc = pd.read_csv("/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Consistency_evaluator/evaluator_data/ENCFF121CAA.bed", sep = "\t", header = None)
    print(accessbile_peaks_ipsc.shape)
    #print(accessbile_peaks_ipsc)


    accessbile_peaks_ipsc.columns = ['chr', 'start', 'end','atac_peak', '4', '5', '6', '7' ,'8', '9']
    accessbile_peaks_ipsc['peak_length'] = accessbile_peaks_ipsc['end'] - accessbile_peaks_ipsc['start']

    accessbile_peaks_ipsc_unique = accessbile_peaks_ipsc.drop_duplicates(subset=['chr', 'start', 'end'], keep='first')
    print(accessbile_peaks_ipsc_unique)

    #pick 150 random rows
    accessbile_peaks_ipsc_unique_150 = accessbile_peaks_ipsc_unique.sample(n=150, random_state=1)
    print(accessbile_peaks_ipsc_unique_150)

    SEQUENCE_LENGTH = 251
    half_len = SEQUENCE_LENGTH // 2

    # Calculate center of each peak
    accessbile_peaks_ipsc_unique_150['midpoint'] = (
        (accessbile_peaks_ipsc_unique_150['start'] + accessbile_peaks_ipsc_unique_150['end']) // 2
    )

    # Compute start and end of the 251 bp region
    accessbile_peaks_ipsc_unique_150['seq_start'] = accessbile_peaks_ipsc_unique_150['midpoint'] - half_len
    accessbile_peaks_ipsc_unique_150['seq_end'] = accessbile_peaks_ipsc_unique_150['midpoint'] + half_len
    accessbile_peaks_ipsc_unique_150['peak_length_250'] = accessbile_peaks_ipsc_unique_150['seq_end'] - accessbile_peaks_ipsc_unique_150['seq_start']
    print(accessbile_peaks_ipsc_unique_150)


    fasta_path = "/arc/project/st-cdeboer-1/iluthra/hg38.fa"
    genome = pysam.FastaFile(fasta_path)


    sequences = pd.DataFrame(columns = ['chr','start','end','sequence'])
    for i in range(0,accessbile_peaks_ipsc_unique_150.shape[0]):
        chrom_current = accessbile_peaks_ipsc_unique_150['chr'].iloc[i]
        start_current = accessbile_peaks_ipsc_unique_150['seq_start'].iloc[i]
        end_current = accessbile_peaks_ipsc_unique_150['seq_end'].iloc[i]
        seq = genome.fetch(chrom_current, start_current,  end_current)
        seq = seq.upper()
        data_2_add = pd.DataFrame({'chr': [chrom_current],'start': [start_current], 'end': [end_current], 'sequence': [seq]})
        sequences = pd.concat([sequences, data_2_add])



    #prep files to run through biasaway
    open_peaks = open("/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Consistency_evaluator/evaluator_data/sequences.fasta", "w")
    for k in range(0, sequences.shape[0]):
        open_peaks.write(">" + sequences['chr'].iloc[k] + ":" + str(sequences['start'].iloc[k]) + '-' + str(sequences['end'].iloc[k]) + "\n" +sequences['sequence'].iloc[k] + "\n")

    # #

### read in all the sequences + shuffled versions
#add mutation to the cm enhancers


open_peaks = SeqIO.parse(open('/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Consistency_evaluator/evaluator_data/sequences.fasta'),'fasta')
open_peaks_seqs = []
open_peaks_names = []

open_peaks_seqs_RC = []
open_peaks_names_RC = []

for seq in open_peaks:
    name, sequence = seq.id, str(seq.seq)
    open_peaks_seqs.append(sequence)
    open_peaks_names.append(name)

    open_peaks_seqs_RC.append(str(Seq(sequence).reverse_complement()))
    open_peaks_names_RC.append(name + '_RC')


open_peaks_mononuc = SeqIO.parse(open('/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Consistency_evaluator/evaluator_data/mononucleotide_shuffled.fa'),'fasta')
open_peaks_mononuc_seqs = []
open_peaks_mononuc_names = []
open_peaks_mononuc_seqs_RC = []
open_peaks_mononuc_names_RC = []

for seq in open_peaks_mononuc:
    name, sequence = seq.description, str(seq.seq)
    open_peaks_mononuc_seqs.append(sequence)
    open_peaks_mononuc_names.append(name.split()[-1]+ '_mononuc_shuffled')

    open_peaks_mononuc_seqs_RC.append(str(Seq(sequence).reverse_complement()))
    open_peaks_mononuc_names_RC.append(name.split()[-1]+ '_mononuc_shuffled' + '_RC')

open_peaks_dinuc = SeqIO.parse(open('/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Consistency_evaluator/evaluator_data/dinucleotide_shuffled.fa'),'fasta')
open_peaks_dinuc_seqs = []
open_peaks_dinuc_names = []

open_peaks_dinuc_seqs_RC = []
open_peaks_dinuc_names_RC = []

for seq in open_peaks_dinuc:
    name, sequence = seq.description, str(seq.seq)
    open_peaks_dinuc_seqs.append(sequence)
    open_peaks_dinuc_names.append(name.split()[-1]+ '_dinuc_shuffled')

    open_peaks_dinuc_seqs_RC.append(str(Seq(sequence).reverse_complement()))
    open_peaks_dinuc_names_RC.append(name.split()[-1]+ '_dinuc_shuffled' + '_RC')

all_data = pd.DataFrame({
    'seq_name': open_peaks_names + open_peaks_names_RC + open_peaks_mononuc_names + open_peaks_mononuc_names_RC+ open_peaks_dinuc_names + open_peaks_dinuc_names_RC,
    'sequence': open_peaks_seqs + open_peaks_seqs_RC + open_peaks_mononuc_seqs + open_peaks_mononuc_seqs_RC + open_peaks_dinuc_seqs + open_peaks_dinuc_seqs_RC
})

print(all_data)
all_data.to_csv('/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Consistency_evaluator/evaluator_data/all_consistency_data.csv', sep = '\t')