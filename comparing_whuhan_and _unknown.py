import csv
import json
import sys
import requests
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import Align
from Bio import GenBank
from Bio import SeqIO
from Bio import Entrez
import Bio
pip install biopython

Entrez.email = 'A.N.Other@example.com'
with Entrez.efetch(
    db="nucleotide", rettype="gb", retmode="text", id="NC_045512.2"
) as handle:
    Wuhan_record = SeqIO.read(handle, "genbank")
print("%s with %i features" % (Wuhan_record.id, len(Wuhan_record.features)))

Entrez.email = 'A.N.Other@example.com'
with Entrez.efetch(
    db="nucleotide", rettype="gb", retmode="text", id="MZ077342.1"
) as handle:
    Delta_record = SeqIO.read(handle, "genbank")
print("%s with %i features" % (Delta_record.id, len(Delta_record.features)))

Entrez.email = 'A.N.Other@example.com'
with Entrez.efetch(
    db="nucleotide", rettype="gb", retmode="text", id="MZ007342.1"
) as handle:
    Unknown_record = SeqIO.read(handle, "genbank")
print("%s with %i features" % (Unknown_record.id, len(Unknown_record.features)))

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.mismatch_score = -10
aligner.target_internal_open_gap_score = -20.000000
aligner.target_internal_extend_gap_score = -5.00000
aligner.target_left_open_gap_score = -10.000000
aligner.target_left_extend_gap_score = -0.500000
aligner.target_right_open_gap_score = -10.000000
aligner.target_right_extend_gap_score = -0.500000
aligner.query_internal_open_gap_score = -10.000000
aligner.query_internal_extend_gap_score = -0.500000
aligner.query_left_open_gap_score = -10.000000
aligner.query_left_extend_gap_score = -0.500000
aligner.query_right_open_gap_score = -10.000000
aligner.query_right_extend_gap_score = -0.500000

# two or more sequences are arranged to line up the corresponding nucleotides or amino acids.

spikeW = Wuhan_record.seq[21562:25384].translate()
spikeD = Delta_record.seq[21525:25427].translate()
spikeU = Unknown_record.seq[21523:25427].translate()
alignmentWvsD_aa = aligner.align(spikeW, spikeD)

print(alignmentWvsD_aa[0])
# Now we run the aligner. The output is the amino acid sequence in the single letter code, Wuhan isolate on top.

spikeW_r = Seq(spikeW)
spikeW_record = SeqRecord(spikeW_r)
spikeW_record.id = "Wuhan"
spikeD_r = Seq(spikeD)
spikeD_record = SeqRecord(spikeD_r)
spikeD_record.id = "Delta"
spikeU_r = Seq(spikeU)
spikeU_record = SeqRecord(spikeU_r)
spikeU_record.id = "Unknown"
altogether = spikeW_record.format(
    "fasta") + spikeD_record.format("fasta") + spikeU_record.format("fasta")
print(altogether)
# The following code block prints out the file with three spike sequences. They are in a format called fasta, which is widely used for this sort of job.

# https://rest.uniprot.org/beta/docs/
WEBSITE_API = "https://rest.uniprot.org/beta"

# Helper function to download data


def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response

# Execute the next two code blocks to submit the job. Depending on how busy the EBI computer is, it will take a few seconds to several hours.


r = requests.post("https://www.ebi.ac.uk/Tools/services/rest/clustalo/run", data={
    "email": "lestimpe@gmail.com",
    "iterations": 0,
    "outfmt": "clustal_num",
    "order": "aligned",
    "sequence": altogether
})

job_id = r.text
print(job_id)

# get job status
r = get_url(f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}")
print(r.text)
# If the status is "RUNNING", click the following code block to check. When the job is finished run the last code block to display the aligned sequences.

# get job status
r = get_url(f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}")

print(r.text)
r = get_url(f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal_num")
print(r.text)
