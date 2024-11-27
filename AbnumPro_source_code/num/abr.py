import os
import subprocess
import sys
import tempfile
import gzip
import math

from functools import partial
from textwrap import wrap
from subprocess import Popen, PIPE
from itertools import groupby, islice
from multiprocessing import Pool

from Bio.SearchIO.HmmerIO import Hmmer3TextParser as HMMERParser

def abr(sequences):
    alignments = {}
    alignments["D1"] = run_hmmer(sequences,"D1.hmm")
    alignments["D2"] = run_hmmer(sequences, "D2.hmm")
    alignments["D3"] = run_hmmer(sequences, "D3.hmm")
    alignments["D4"] = run_hmmer(sequences, "D4.hmm")
    results = abr_sequences_from_alignment(sequences,alignments)
    return results

def abr_sequences_from_alignment(sequences, alignments, scheme="imgt", allow=set(["H", "K", "L", "A", "B", "G", "D"]),
                                    assign_germline=False, allowed_species=None):
    
    hit_tables = []
    descriptions = []
    abr_lists = []
    d1_alignments = alignments["D1"]
    d2_alignments = alignments["D2"]
    d3_alignments = alignments["D3"]
    d4_alignments = alignments["D4"]
    for i in range(len(sequences)):
        region = {}
        if d1_alignments[i]==None or d2_alignments[i]==None or d3_alignments[i]==None or d4_alignments[i]==None:
            hit_tables.append({sequences[i][0]:None})
            descriptions.append({"chain_type":None})
            abr_lists.append([])
            continue
        if int(d1_alignments[i]['query_end']) >= int(d2_alignments[i]['query_start']) or int(d2_alignments[i]['query_end']) >= int(d3_alignments[i]['query_start']) or int(d3_alignments[i]['query_end']) >= int(d4_alignments[i]['query_start']):
            hit_tables.append({sequences[i][0]: None})
            descriptions.append({"chain_type":None})
            abr_lists.append([])
            continue
        if d1_alignments[i]['chain_type']!=d2_alignments[i]['chain_type'] or d2_alignments[i]['chain_type']!=d3_alignments[i]['chain_type'] or d3_alignments[i]['chain_type']!=d4_alignments[i]['chain_type']:
            hit_tables.append({sequences[i][0]: None})
            descriptions.append({"chain_type":None})
            abr_lists.append([])
            continue
        
      
        fr1 = sequences[i][1][d1_alignments[i]['query_start']:d1_alignments[i]['query_end'] + 1]
        fr2 = sequences[i][1][d2_alignments[i]['query_start']:d2_alignments[i]['query_end'] + 1]
        fr3 = sequences[i][1][d3_alignments[i]['query_start']:d3_alignments[i]['query_end'] + 1]
        fr4 = sequences[i][1][d4_alignments[i]['query_start']:d4_alignments[i]['query_end'] + 1]
        abr1 = sequences[i][1][d1_alignments[i]['query_end']+1:d2_alignments[i]['query_start']]
        abr2 = sequences[i][1][d2_alignments[i]['query_end'] + 1:d3_alignments[i]['query_start']]
        abr3 = sequences[i][1][d3_alignments[i]['query_end'] + 1:d4_alignments[i]['query_start']]
        region["fr1"] = fr1
        region["abr1"] = abr1
        region["fr2"] = fr2
        region["abr2"] = abr2
        region["fr3"] = fr3
        region["abr3"] = abr3
        region["fr4"] = fr4
        chain_type = ""
        if d1_alignments[i]['chain_type']=="light":
            chain_type="L"
        else:
            chain_type="H"
        descriptions.append({"chain_type":chain_type})
        
        hit_tables.append({sequences[i][0]:region})
        abr_list = []
        abr_list.append(fr1)
        abr_list.append(abr1)
        abr_list.append(fr2)
        abr_list.append(abr2)
        abr_list.append(fr3)
        abr_list.append(abr3)
        abr_list.append(fr4)
        abr_lists.append(abr_list)




    return descriptions,abr_lists,hit_tables
def run_hmmer(sequence_list, hmm_database="D1.hmm", hmmerpath="", ncpu=None, bit_score_threshold=0, hmmer_species=None):
    
    fasta_filehandle, fasta_filename = tempfile.mkstemp(".fasta", text=True)
    with os.fdopen(fasta_filehandle, 'w') as outfile:
        write_fasta(sequence_list, outfile)

    output_filehandle, output_filename = tempfile.mkstemp(".txt", text=True)
    "/cygdrive/c/mypath/myfile"
    # HMM_P1 = HMM.replace("\\", "/").replace("C:", "/cygdrive/c")

    output_filename_P = output_filename.replace(
        "\\", "/").replace("C:", "/cygdrive/c")
    fasta_filename_P = fasta_filename.replace(
        "\\", "/").replace("C:", "/cygdrive/c")
    HMM_P = hmm_database       #"D1.hmm"

    script_path = os.path.abspath(__file__)


    script_directory = os.path.dirname(script_path)
    hmmerpath = script_directory
    
    HMM = HMM_P
    if hmmerpath:
        hmmscan = os.path.join(hmmerpath, "hmmscan")
    else:
        hmmscan = "hmmscan"
    try:
        if ncpu is None:
            command = [hmmscan, "-o", output_filename_P,
                       HMM_P, fasta_filename_P]
        else:
            command = [hmmscan, "-o", output_filename_P,
                       "--cpu", str(ncpu), HMM, fasta_filename]

        working_directory = script_directory
       
        process = subprocess.Popen(
            command, cwd=working_directory, stdout=PIPE, stderr=PIPE)
       
        _, pr_stderr = process.communicate()
        # 关闭进程
        process.terminate()

        if pr_stderr:
           
            _f = os.fdopen(output_filehandle)
            _f.close()

            print(pr_stderr)
        results = parse_hmmer_output(
            output_filehandle, bit_score_threshold=bit_score_threshold, hmmer_species=hmmer_species)

    finally:
        # clear up
        os.remove(fasta_filename)
        os.remove(output_filename)
        print("1")
    print(results)
    best_results = []
    for result in results:
        hit = result[2]
        if hit:
            best_results.append(sorted(hit, key=lambda x: (x['bit_len'],x['bitscore']))[-1])
        else:
            best_results.append(None)

    return best_results

def parse_hmmer_output(filedescriptor="", bit_score_threshold=80, hmmer_species=None):

    results = []
    if type(filedescriptor) is str:
        openfile = open
    elif type(filedescriptor) is int:
        openfile = os.fdopen

    with openfile(filedescriptor) as inputfile:
        p = HMMERParser(inputfile)
        print(p)
        for query in p:
            print(query)
            results.append(_parse_hmmer_query(
                query, bit_score_threshold=bit_score_threshold, hmmer_species=hmmer_species))
    print(results)

    return results

def _parse_hmmer_query(query, bit_score_threshold=80, hmmer_species=None):
    bit_score_threshold = -20
   
    hit_table = [['id', 'description', 'evalue', 'bitscore', 'bias',
                  'query_start', 'query_end','bit_len']]



    top_descriptions, domains, state_vectors = [], [], []

    if query.hsps:  # We have some hits
        print(query.hsps)
       
        hsp_list = query.hsps
        # Iterate over the matches of the domains in order of their e-value (most significant first)
        for hsp in sorted(hsp_list, key=lambda x: x.evalue):
            print(hsp)
            new = True
            
            if hsp.bitscore >= bit_score_threshold:
                # Check to see if we already have seen the domain
                for i in range(len(domains)):
                    if _domains_are_same(domains[i], hsp):
                        new = False
                        break
               
                hit_table.append([hsp.hit_id, hsp.hit_description, hsp.evalue,
                                  hsp.bitscore, hsp.bias, hsp.query_start, hsp.query_end,hsp.hit_end-hsp.hit_start])
                if new:  # It is a new domain and this is the best hit. Add it for further processing.
                    domains.append(hsp)
                    
                    top_descriptions.append(
                        dict(list(zip(hit_table[0], hit_table[-1]))))

       
        ordering = sorted(list(range(len(domains))),
                          key=lambda x: domains[x].query_start)
        domains = [domains[_] for _ in ordering]
        top_descriptions = [top_descriptions[_] for _ in ordering]

    ndomains = len(domains)
  
    for i in range(ndomains):
        domains[i].order = i
        chain_type, region = top_descriptions[i]["id"].split("_")
        domains[i]
        
        top_descriptions[i]["chain_type"] = chain_type # Reparse
        top_descriptions[i]["region"] = region
       
    return hit_table, state_vectors, top_descriptions

def _domains_are_same(dom1, dom2):
   
    dom1, dom2 = sorted([dom1, dom2], key=lambda x: x.query_start)
    if dom2.query_start >= dom1.query_end:
        return False
    return True

def write_fasta(sequences, f):
  
    for name, sequence in sequences:
        print(">%s" % name, file=f)
        print('\n'.join(['\n'.join(wrap(block, width=80))
                         for block in sequence.splitlines()]), file=f)


if __name__ == "__main__":
    # Test and example useage of the abn function.
    sequences = [("12e8:H",
                  "QHLEQSGGGAGGGLVKPGGSLELCCKASGFTFSSYYMCWVRQAPGKGLEWIGCIYAGSSGSRGNTYYGSWVNGRFTLSRDIDQSTGCLQLNSLTVADTAMHYCARLSWXXYWYSGWGTSGAQAP"),
                 ("12e8:L",
                  "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASV"),
                 ("scfv:A",
                  "DIQMTQSPSSLSASVGDRVTITCRTSGNIHNYLTWYQQKPGKAPQLLIYNAKTLADGVPSRFSGSGSGTQFTLTISSLQPEDFANYYCQHFWSLPFTFGQGTKVEIKRTGGGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFDFSRYDMSWVRQAPGKRLEWVAYISSGGGSTYFPDTVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCARQNKKLTWFDYWGQGTLVTVSSHHHHHH"),
                 ("lysozyme:A",
                  "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL")]
    sequences1 = [("1A3R_2|Chain B[auth H]",
                  "QHLEVQLQQSGAELVRPGASVKLSCTTSGFNIKDIYIHWVKQRPEQGLEWIGRLDPANGYTKYDPKFQGKATITVDTSSNTAYLHLSSLTSEDTAVYYCDGYYSYYDMDYWGPGTSVTVSSAKTTAPSVYPLAPVCGDTTGSSVTLGCLVKGYFPEPVTLTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVTSSTWPSQSITCNVAHPASSTKVDKKIEPR"),
                 ("1A3R_1|Chain A[auth L]","DIVMTQSPSSLTVTTGEKVTMTCKSSQSLLNSRTQKNYLTWYQQKPGQSPKLLIYWASTRESGVPDRFTGSGSGTDFTLSISGVQAEDLAVYYCQNNYNYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC")
                 ]
    results = abr(sequences1)
    print(results)