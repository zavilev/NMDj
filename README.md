## Description

NMDj is a tool for searching, classification and quantification of local splicing events, which lead to formation of transcripts targeted by nonsense-mediated decay (NMD) system. NMDj takes transcript annotation (in GTF format) and outputs local splicing events in the form of characteristic junctions (CJs). CJs are the set of alternative junctions from the same transcript region, derived either from coding or NMD transcripts. A switch from coding CJs to NMD CJs induces frameshift or creates introns more than 50nt downstream of a stop-codon, which triggers NMD.

NMDj classifies found events in a number of types, including but not limited to poison and essential cassette exons, alternative splice sites and intron retention. Additionally, NMDj can quantify local splicing events in the form of PSI values, based on RNA-Seq junction split-read counts. 

If input GTF originates from transcript assemblers and lacks ORF annotations, NMDj can assign transcripts to genes and annotate ORFs if existing gene and ORF annotations for the same genome were given.

## Installation

To install NMDj [clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system.

    git clone https://github.com/zavilev/NMDj.git

NMDj directory contains `./scripts` and `./example` directories. The `./scripts` directory contains NMDj code (NMDj.py executable script dependent on nmdj_functions.py) and scripts to run example analysis. The `./example` directory contains example input data to run example analysis.

### Dependencies

To run NMDj, you need the following python dependencies:

- python (tested with 3.8)
- pandas (tested with 1.5.2)
- gtfparse (tested with 2.5.0)
- numpy (tested with 1.24.4)
- biopython (tested with 1.81)

To run example script, you additionally need [RSEM](https://github.com/deweylab/RSEM) (tested with 1.3.3), [STAR](https://github.com/alexdobin/STAR) (tested with 2.7.3a) and [pyIPSA](https://github.com/pervouchinelab/pyIPSA)

## Running

Add NMDj.py to your PATH. To get usage options run

    NMDj.py -h

The output is 

    usage: NMDj.py [-h] -g infile.gtf [-o] [-a annotation.gtf] [-G genome.fa] [-n] [-N] [-r transcripts.txt or attr:value] [-p path/to/output/] [-q file.txt]
                   [--threads THREADS] [--no_clustering]

    optional arguments:
      -h, --help            show this help message and exit
      -g infile.gtf, --gtf infile.gtf
                            input GTF file
      -o, --orf             whether to annotate ORFs. Requires --ann and --genome argument to be set
      -a annotation.gtf, --ann annotation.gtf
                            GTF file with annotated genes, transcripts and ORFs
      -G genome.fa, --genome genome.fa
                            single multifasta genome file. Chromosome names must match with GTFs
      -n, --nmd             whether to predict targets of NMD in input.gtf
      -N, --nmdann          whether to predict targets of NMD in annotation.gtf
      -r transcripts.txt or attr:value, --ref transcripts.txt or attr:value
                            Either a file with reference transcript ids, each id on a new line, or a string <attribute>:<value>. In the second case transcript
                            is chosen as reference if the value of its <attribute> contains <value>
      -p path/to/output/, --prefix path/to/output/
                            Prefix for NMDj output files
      -q file.txt, --ipsa_files file.txt
                            File with paths to ipsa files containing counts of RNA-Seq split-reads aligned to junctions
      --threads THREADS     number of threads for parallel processing
      --no_clustering       skip clustering events

## Output

NMDj outputs several files, depending on the options specified. The main output files are:

- junctions.tsv

event_id | transcript_biotype | chromosome | start | end | strand | junction_type | transcripts_with_junction
 --- | --- | --- | --- | --- | --- | --- | --- 
ENSG00000000460_0 | nonsense_mediated_decay | chr1 | 169798958 | 169802621 | + | junction | ENST00000459772, ENST00000466580, ENST00000481744
ENSG00000000460_0 | protein_coding | chr1 | 169798958 | 169800883 | + | junction | ENST00000286031, ENST00000359326,  ENST00000413811, ENST00000472795, ENST00000496973

This file contains CJs. Event_id is a unique ID of event, derived from gene id and an ordinal number (0-based) of event in this gene. Junction_type is either junction (spliced) or intron_retention. Other columns names are self-explanatory.

- classification.tsv

gene_id | event_id | interval_start | interval_end | nmd_id_list | coding_id_list | VIP | event_type | all_CJs
 --- | --- | --- | --- | --- | --- | --- | --- | ---  
ENSG00000000460 | ENSG00000000460_0 | 169798958 | 169802643 | ENST00000459772, ENST00000466580, ENST00000481744 | ENST00000286031, ENST00000359326, ENST00000413811, ENST00000472795, ENST00000496973 | DA:DADA | EE | Yes
ENSG00000000971 | ENSG00000000971_0 | 196673962 | 196675529 | ENST00000695986 | ENST00000359637, ENST00000367429, ENST00000630130, ENST00000695968, ... | DADA:DA | PE | Yes

This file contains classification of events (event_type), each event on a new line. Interval_start and interval_end are borders of event region.

- PSI.tsv - if split-read counts were given

event_id | sample_id | nmd | denom | psi
--- | --- | --- | --- | ---
ENSG00000000460_0 | sample1 | 7.0 | 18.5 | 0.3783783783783784
ENSG00000000971_0 | sample1 | 0.0 | 13.0 | 0.0

- novel_ORF.gtf - if novel transcripts were given

chromosome | source  | type | start | end | score | strand | frame | attributes    
--- | --- | --- | --- | --- | --- | --- | --- | ---
chr1 | NMDj | start_codon | 2412500 | 2412502 | . | - | . | gene_id "ENSG00000157911"; transcript_id "STRG.100.1";
chr1 | NMDj | start_codon | 2412500 | 2412502 | . | - | . | gene_id "ENSG00000157911"; transcript_id "STRG.100.2";
chr1 | NMDj | stop_codon | 62687626 | 62687628 | . | - | . | gene_id "ENSG00000116641"; transcript_id "STRG.1385.8";

## Example

Example script `./scripts/run_example.sh`:
1. runs RSEM simulation of RNA-Seq reads based on GTEx sample expression profile;
2. aligns simulated reads to genome;
3. runs pyIPSA to obtain split-read counts;
4. runs NMDj to find and quantify local events;
5. compares ground truth PSI values with those obtained by NMDj;
6. runs NMDj based on GTF containing novel transcripts obtained by StringTie, to demonstrate ability to annotate ORFs

Modify config section in `./scripts/run_example.sh`, then

    cd ./scripts
    ./run_example.sh

The following content will be created in `./example`:

- `rsem_ref/`. Contains results of rsem-prepare-reference
- `simulation/`. Contains results of rsem-simulate-reads
- `bam/`. Contains result of STAR
- `ipsa_output/`. Contains pyIPSA output
- `nmdj_output/`. Contains NMDj results
- `nmdj_output_stringite/`. Contains NMDj results based on GTF containing novel transcripts obtained by StringTie
- `NMDj_PSI_vs_ground_truth.png` --- a scatter plot of NMDj PSI vs ground truth PSI

Approximate execution time on 8 threads is 10 min
