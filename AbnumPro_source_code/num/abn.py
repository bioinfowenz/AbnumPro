

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


from .schemes import *
from .germlines import all_germlines

all_species = list(all_germlines['V']['H'].keys())

amino_acids = sorted(list("QWERTYIPASDFGHKLCVNM"))
set_amino_acids = set(amino_acids)
abn_path = os.path.split(__file__)[0]

scheme_short_to_long = {"m": "martin", "c": "chothia", "k": "kabat", "imgt": "imgt", "kabat": "kabat",
                        "chothia": "chothia", "martin": "martin", "i": "imgt", "a": "aho", "aho": "aho", "wolfguy": "wolfguy", "w": "wolfguy"}

scheme_names = list(scheme_short_to_long.keys())
chain_type_to_class = {"H": "H", "K": "L", "L": "L",
                       "A": "A", "B": "B", "G": "G", "D": "D"}

HMM_path = os.path.join(abn_path, "dat", "HMMs")


all_reference_states = list(range(1, 129))


class HMMscanError(Exception):
    def __init__(self, message):
        
        super(HMMscanError, self).__init__(message)




def read_fasta(filename):
    """
    Read a sequence file and parse as description, string 
    """
    return [r for r in fasta_iter(filename)]


def fasta_iter(fasta_name):

    if fasta_name.endswith('.gz'):  # IOError raised upon iteration if not a real gzip file.
        fh = gzip.open(fasta_name)
    else:
        fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        #header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def write_fasta(sequences, f):

    for name, sequence in sequences:
        print(">%s" % name, file=f)
        print('\n'.join(['\n'.join(wrap(block, width=80))
              for block in sequence.splitlines()]), file=f)


def validate_sequence(sequence):

    assert len(sequence) < 10000, "Sequence too long."
    assert not (set(sequence.upper()) - set_amino_acids), "Unknown amino acid letter found in sequence: %s" % ", ".join(
        list((set(sequence.upper()) - set_amino_acids)))
    return True


def validate_numbering(xxx_todo_changeme, name_seq=[]):

    (numbering, start, end) = xxx_todo_changeme
    name, seq = name_seq
    last = -1
    nseq = ""

    for (index, _), a in numbering:
        assert index >= last, "Numbering was found to decrease along the sequence %s. Please report." % name
        last = index
        nseq += a.replace("-", "")

    assert nseq in seq.replace(
        "-", ""), "The algorithm did not number a contiguous segment for sequence %s. Please report" % name

    return numbering, start, end


def grouper(n, iterable):

    it = iter(iterable)

    def take():
        while 1:
            yield list(islice(it, n))
    return iter(take().__next__, [])


def abn_output(numbered, sequences, alignment_details, outfile, sequence_id=None, domain_id=None):

    assert (sequence_id is not None) or (
        sequence_id is None and domain_id is None), "If domain_id is specified, sequence_id must also be specified."
    for i in range(len(numbered)):
        if sequence_id is None:
            print("# %s" % sequences[i][0], file=outfile)  # print the name
        if numbered[i] is not None:
            if sequence_id is not None:
                if i != sequence_id:
                    continue
            print("# abn numbered", file=outfile)
            for j in range(len(numbered[i])):  # Iterate over domains
                if domain_id is not None:
                    if j != domain_id:
                        continue
                print("# Domain %d of %d" %
                      (j+1, len(numbered[i])), file=outfile)
                print("# Most significant HMM hit", file=outfile)
                print(
                    "#|species|chain_type|e-value|score|seqstart_index|seqend_index|", file=outfile)
                alignment_details[i][j]["evalue"] = str(
                    alignment_details[i][j]["evalue"])
                print("#|%s|%s|%s|%.1f|%d|%d|" % tuple([alignment_details[i][j][field] for field in
                                                        ["species", "chain_type", "evalue", "bitscore"]]
                                                       + [numbered[i][j][1], numbered[i][j][2]]), file=outfile)

                if 'germlines' in alignment_details[i][j]:
                    print('# Most sequence-identical germlines', file=outfile)
                    print(
                        '#|species|v_gene|v_identity|j_gene|j_identity|', file=outfile)
                    (species, vgene), vid = alignment_details[i][j]['germlines'].get(
                        'v_gene', [['', 'unknown'], 0])
                    if vgene is None:
                        vgene, vid = 'unknown', 0
                    (_, jgene), jid = alignment_details[i][j]['germlines'].get(
                        'j_gene', [['', 'unknown'], 0])
                    if jgene is None:
                        jgene, jid = 'unknown', 0
                    print('#|%s|%s|%.2f|%s|%.2f|' %
                          (species, vgene, vid, jgene, jid), file=outfile)
                chain_type = chain_type_to_class[alignment_details[i]
                                                 [j]["chain_type"]]
                print("# Scheme = %s" %
                      alignment_details[i][j]["scheme"], file=outfile)
                if len(numbered[i][j][0]) == 0:
                    print("# Warning: %s scheme could not be applied to this sequence." %
                          alignment_details[i][j]["scheme"], file=outfile)
                for (index, insertion), aa in numbered[i][j][0]:
                    print(chain_type, ("%d" % index).ljust(
                        5), insertion, aa, file=outfile)
        print("//", file=outfile)


def csv_output(sequences, numbered, details, outfileroot):


    chain_types = {}
    pos_ranks = {}
    all_pos = {}
    _lc = {'K': 'KL', 'L': 'KL'}

    # Divide the set into chain types and find how to order the numbering for each type.
    for i in range(len(sequences)):  # Iterate over entries
        if numbered[i] is None:
            continue

        for j in range(len(numbered[i])):  # Iterate over domains.
            # Record the chain type index
            c = details[i][j]['chain_type']
            c = _lc.get(c, c)  # Consider lambda and kappa together.
            chain_types.setdefault(c, []).append((i, j))
            if c not in pos_ranks:
                pos_ranks[c] = {}
                all_pos[c] = set()

            # Update the insertion order for the scheme. i.e. is it A B C or C B A (e.g. imgt 111 and 112 repectively)
            l = -1
            r = 0
            for p, _ in numbered[i][j][0]:
                if p[0] != l:
                    l = p[0]
                    r = 0
                else:
                    r += 1
                pos_ranks[c][p] = max(r, pos_ranks[c].get(p, r))
                all_pos[c].add(p)

    # Write a new file for each chain type. Kappa and lambda are written together as light chains.
    for cts in ['H', 'KL', 'A', 'B', 'G', 'D']:
        if cts in chain_types:
            with open(outfileroot + '_%s.csv' % cts, 'w') as out:

                # Sort the positions by index and insertion order
                positions = sorted(
                    all_pos[cts], key=lambda p: (p[0], pos_ranks[cts][p]))

                # Header line
                fields = ['Id', 'domain_no', 'hmm_species', 'chain_type', 'e-value', 'score', 'seqstart_index', 'seqend_index',
                          'identity_species', 'v_gene', 'v_identity', 'j_gene', 'j_identity']
                fields += [('%d%s' % (p)).strip() for p in positions]
                print(','.join(fields), file=out)

                # Iterate over the domains identified
                for i, j in chain_types[cts]:
                    line = [sequences[i][0].replace(',', ' '),
                            str(j),
                            details[i][j].get('species', ''),
                            details[i][j].get('chain_type', ''),
                            str(details[i][j].get('evalue', '')),
                            str(details[i][j].get('bitscore', '')),
                            str(numbered[i][j][1]),
                            str(numbered[i][j][2]),
                            details[i][j].get('germlines', {}).get(
                                'v_gene', [['', ''], 0])[0][0],
                            details[i][j].get('germlines', {}).get(
                                'v_gene', [['', ''], 0])[0][1],
                            '%.2f' % details[i][j].get('germlines', {}).get(
                                'v_gene', [['', ''], 0])[1],
                            details[i][j].get('germlines', {}).get(
                                'j_gene', [['', ''], 0])[0][1],
                            '%.2f' % details[i][j].get('germlines', {}).get('j_gene', [['', ''], 0])[1]]

                    # Hash the numbering. Insertion order has been preserved in the positions sort.
                    d = dict(numbered[i][j][0])
                    line += [d.get(p, '-') for p in positions]

                    assert len(line) == len(fields)
                    print(','.join(line), file=out)


## Parsing and recognising domain hits from hmmscan ##
def _domains_are_same(dom1, dom2):
    """
    Check to see if the domains are overlapping.
    @param dom1: 
    @param dom2: 

    @return: True or False  
    """
    dom1, dom2 = sorted([dom1, dom2], key=lambda x: x.query_start)
    if dom2.query_start >= dom1.query_end:
        return False
    return True


def _parse_hmmer_query(query, bit_score_threshold=80, hmmer_species=None):

    hit_table = [['id', 'description', 'evalue', 'bitscore', 'bias',
                  'query_start', 'query_end']]

    # Find the best hit for each domain in the sequence.

    top_descriptions, domains, state_vectors = [], [], []

    if query.hsps:  # We have some hits
        # If we have specified a species, check to see we have hits for that species
        # Otherwise revert back to using any species
        if hmmer_species:
            #hit_correct_species = [hsp for hsp in query.hsps if hsp.hit_id.startswith(hmmer_species) and hsp.bitscore >= bit_score_threshold]
            hit_correct_species = []
            for hsp in query.hsps:
                if hsp.bitscore >= bit_score_threshold:
                    for species in hmmer_species:
                        if hsp.hit_id.startswith(species):
                            hit_correct_species.append(hsp)

            if hit_correct_species:
                hsp_list = hit_correct_species
            else:
                print("Limiting hmmer search to species %s was requested but hits did not achieve a high enough bitscore. Reverting to using any species" % (
                    hmmer_species))
                hsp_list = query.hsps
        else:
            hsp_list = query.hsps

        # Iterate over the matches of the domains in order of their e-value (most significant first)
        for hsp in sorted(hsp_list, key=lambda x: x.evalue):
            new = True
            # Only look at those with hits that are over the threshold bit-score.
            if hsp.bitscore >= bit_score_threshold:
                # Check to see if we already have seen the domain
                for i in range(len(domains)):
                    if _domains_are_same(domains[i], hsp):
                        new = False
                        break
                hit_table.append([hsp.hit_id, hsp.hit_description, hsp.evalue,
                                 hsp.bitscore, hsp.bias, hsp.query_start, hsp.query_end])
                if new:  # It is a new domain and this is the best hit. Add it for further processing.
                    domains.append(hsp)
                    # Add the last added to the descriptions list.
                    top_descriptions.append(
                        dict(list(zip(hit_table[0], hit_table[-1]))))

        # Reorder the domains according to the order they appear in the sequence.
        ordering = sorted(list(range(len(domains))),
                          key=lambda x: domains[x].query_start)
        domains = [domains[_] for _ in ordering]
        top_descriptions = [top_descriptions[_] for _ in ordering]

    ndomains = len(domains)
    # If any significant hits were identified parse and align them to the reference state.
    for i in range(ndomains):
        domains[i].order = i
        species, chain = top_descriptions[i]["id"].split("_")
        # Alignment to the reference states.
        state_vectors.append(_hmm_alignment_to_states(
            domains[i], ndomains, query.seq_len))
        top_descriptions[i]["species"] = species  # Reparse
        top_descriptions[i]["chain_type"] = chain
        # Make sure the query_start agree if it was changed
        top_descriptions[i]["query_start"] = state_vectors[-1][0][-1]

    return hit_table, state_vectors, top_descriptions


def _hmm_alignment_to_states(hsp, n, seq_length):


    # Extract the strings for the reference states and the posterior probability strings
    reference_string = hsp.aln_annotation["RF"]
    state_string = hsp.aln_annotation["PP"]

    assert len(reference_string) == len(
        state_string), "Aligned reference and state strings had different lengths. Don't know how to handle"

    _hmm_start = hsp.hit_start
    _hmm_end = hsp.hit_end

    _seq_start = hsp.query_start
    _seq_end = hsp.query_end

    # Extact the full length of the HMM hit
    species, ctype = hsp.hit_id.split('_')
    _hmm_length = get_hmm_length(species, ctype)

 
    if hsp.order == 0 and _hmm_start and _hmm_start < 5:
        n_extend = _hmm_start
        if _hmm_start > _seq_start:
            n_extend = min(_seq_start, _hmm_start - _seq_start)
        state_string = '8'*n_extend + state_string
        reference_string = 'x'*n_extend + reference_string
        _seq_start = _seq_start - n_extend
        _hmm_start = _hmm_start - n_extend

    
    if n == 1 and _seq_end < seq_length and (123 < _hmm_end < _hmm_length):
        n_extend = min(_hmm_length - _hmm_end, seq_length - _seq_end)
        state_string = state_string + '8'*n_extend
        reference_string = reference_string + 'x'*n_extend
        _seq_end = _seq_end + n_extend
        _hmm_end = _hmm_end + n_extend

    
    hmm_states = all_reference_states[_hmm_start: _hmm_end]
    sequence_indices = list(range(_seq_start,  _seq_end))
    h, s = 0, 0  

    state_vector = []
    
    for i in range(len(state_string)):
        if reference_string[i] == "x":  # match state
            state_type = "m"
        else:  # insert state
            state_type = "i"

        if state_string[i] == ".":  # overloading if deleted relative to reference. delete_state
            state_type = "d"
            sequence_index = None
        else:
            sequence_index = sequence_indices[s]
        

        state_vector.append(((hmm_states[h], state_type),  sequence_index))

        # Updates to the indices
        if state_type == "m":
            h += 1
            s += 1
        elif state_type == "i":
            s += 1
        else:  # delete state
            h += 1

    return state_vector


def parse_hmmer_output(filedescriptor="", bit_score_threshold=80, hmmer_species=None):

    results = []
    if type(filedescriptor) is str:
        openfile = open
    elif type(filedescriptor) is int:
        openfile = os.fdopen

    with openfile(filedescriptor) as inputfile:
        p = HMMERParser(inputfile)
        for query in p:
            results.append(_parse_hmmer_query(
                query, bit_score_threshold=bit_score_threshold, hmmer_species=hmmer_species))

    return results


def run_hmmer(sequence_list, hmm_database="ALL", hmmerpath="", ncpu=None, bit_score_threshold=80, hmmer_species=None):


    # Check that hmm_database is available

    assert hmm_database in ["ALL"], "Unknown HMM database %s" % hmm_database
    HMM = os.path.join(HMM_path, "%s.hmm" % hmm_database)


    fasta_filehandle, fasta_filename = tempfile.mkstemp(".fasta", text=True)
    with os.fdopen(fasta_filehandle, 'w') as outfile:
        write_fasta(sequence_list, outfile)

    output_filehandle, output_filename = tempfile.mkstemp(".txt", text=True)
    "/cygdrive/c/mypath/myfile"
    HMM_P1 = HMM.replace("\\", "/").replace("C:", "/cygdrive/c")

    output_filename_P = output_filename.replace(
        "\\", "/").replace("C:", "/cygdrive/c")
    fasta_filename_P = fasta_filename.replace(
        "\\", "/").replace("C:", "/cygdrive/c")
    HMM_P = "ALL.hmm"


    script_path = os.path.abspath(__file__)

    script_directory = os.path.dirname(script_path)
    hmmerpath = script_directory
    if hmmerpath:
        hmmscan = os.path.join(hmmerpath, "hmmscan")
    else:
        hmmscan = "hmmscan"
    try:
        if ncpu is None:
            command = [hmmscan, "-o", output_filename_P,
                       HMM_P,  fasta_filename_P]
        else:
            command = [hmmscan, "-o", output_filename_P,
                       "--cpu", str(ncpu), HMM,  fasta_filename]

        working_directory = script_directory

        process = subprocess.Popen(
            command, cwd=working_directory, stdout=PIPE, stderr=PIPE)
        # process = Popen(command, stdout=PIPE, stderr=PIPE)
        _, pr_stderr = process.communicate()
        # 关闭进程
        process.terminate()

        if pr_stderr:
            
            _f = os.fdopen(output_filehandle)
            _f.close()

            raise HMMscanError(pr_stderr)
        results = parse_hmmer_output(
            output_filehandle, bit_score_threshold=bit_score_threshold, hmmer_species=hmmer_species)

    finally:
        # clear up
        os.remove(fasta_filename)
        os.remove(output_filename)
        print("1")

    return results


def get_hmm_length(species, ctype):

    try:
        return len(list(all_germlines['J'][ctype][species].values())[0].rstrip('-'))
    except KeyError:
        return 128


def number_sequence_from_alignment(state_vector, sequence, scheme="imgt", chain_type=None):

    scheme = scheme.lower()
    if scheme == "imgt":
        return number_imgt(state_vector, sequence)
    elif scheme == "chothia":
        if chain_type == "H":
            return number_chothia_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_chothia_light(state_vector, sequence)
        else:
            raise AssertionError(
                "Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
    elif scheme == "kabat":
        if chain_type == "H":
            return number_kabat_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_kabat_light(state_vector, sequence)
        else:
            raise AssertionError(
                "Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
    elif scheme == "martin":
        if chain_type == "H":
            return number_martin_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_martin_light(state_vector, sequence)
        else:
            raise AssertionError(
                "Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
    elif scheme == "aho":
        # requires the chain type to heuristically put the CDR1 gap in position.
        return number_aho(state_vector, sequence, chain_type)
    elif scheme == "wolfguy":
        if chain_type == "H":
            return number_wolfguy_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_wolfguy_light(state_vector, sequence)
        else:
            raise AssertionError(
                "Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
    else:
        raise AssertionError(
            "Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))


def number_sequences_from_alignment(sequences, alignments, scheme="imgt", allow=set(["H", "K", "L", "A", "B", "G", "D"]),
                                    assign_germline=False, allowed_species=None):
    
    numbered = []
    alignment_details = []
    hit_tables = []
    for i in range(len(sequences)):

       
        hit_table, state_vectors, detailss = alignments[i]

       
        hit_numbered, hit_details = [], []
        for di in range(len(state_vectors)):
            state_vector = state_vectors[di]
            details = detailss[di]
            details["scheme"] = scheme
            details["query_name"] = sequences[i][0]

            
            if state_vector and details["chain_type"] in allow:
                try:
                   
                    hit_numbered.append(validate_numbering(number_sequence_from_alignment(state_vector, sequences[i][1],
                                                                                          scheme=scheme, chain_type=details["chain_type"]), sequences[i]))
                    if assign_germline:
                        details["germlines"] = run_germline_assignment(state_vector, sequences[i][1],
                                                                       details["chain_type"], allowed_species=allowed_species)
                    hit_details.append(details)
                
                except AssertionError as e:
                    print(str(e), file=sys.stderr)
                   
                    raise e
                except Exception as e:
                    print(
                        "Error: Something really went wrong that has not been handled", file=sys.stderr)
                    print(str(e), file=sys.stderr)
                    raise e

        if hit_numbered:
            numbered.append(hit_numbered)
            alignment_details.append(hit_details)
        else:
            numbered.append(None)
            alignment_details.append(None)
        hit_tables.append(hit_table)

    return numbered, alignment_details, hit_tables


def get_identity(state_sequence, germline_sequence):

    assert len(state_sequence) == len(germline_sequence) == 128
    n, m = 0, 0
    for i in range(128):
        if germline_sequence[i] == "-":
            continue
        if state_sequence[i].upper() == germline_sequence[i]:
            m += 1
        n += 1

    if not n:
        return 0
    return float(m)/n


def run_germline_assignment(state_vector, sequence, chain_type, allowed_species=None):
    
    genes = {'v_gene': [None, None],
             'j_gene': [None, None],
             }

   
    state_dict = dict(((i, 'm'), None) for i in range(1, 129))
    state_dict.update(dict(state_vector))
    state_sequence = "".join([sequence[state_dict[(i, 'm')]] if state_dict[(
        i, 'm')] is not None else "-" for i in range(1, 129)])

   
    if chain_type in all_germlines["V"]:
        if allowed_species is not None:
            if not all([sp in all_germlines['V'][chain_type] for sp in allowed_species]):  # Made non-fatal
                return {}
        else:
            allowed_species = all_species
        seq_ids = {}
        for species in allowed_species:
            if species not in all_germlines["V"][chain_type]:
                continue  # Previously bug.
            for gene, germline_sequence in all_germlines["V"][chain_type][species].items():
                seq_ids[(species, gene)] = get_identity(
                    state_sequence, germline_sequence)
        genes['v_gene'][0] = max(seq_ids, key=lambda x: seq_ids[x])
        genes['v_gene'][1] = seq_ids[genes['v_gene'][0]]

        
        species = genes['v_gene'][0][0]
        if chain_type in all_germlines["J"]:
            if species in all_germlines["J"][chain_type]:
                seq_ids = {}
                for gene, germline_sequence in all_germlines["J"][chain_type][species].items():
                    seq_ids[(species, gene)] = get_identity(
                        state_sequence, germline_sequence)
                genes['j_gene'][0] = max(seq_ids, key=lambda x: seq_ids[x])
                genes['j_gene'][1] = seq_ids[genes['j_gene'][0]]

    return genes


def check_for_j(sequences, alignments, scheme):

    for i in range(len(sequences)):
        # Check the alignment for J region
        if len(alignments[i][1]) == 1:  # Only do for single domain chains.

           
            ali = alignments[i][1][0]

           
            last_state = ali[-1][0][0]
            last_si = ali[-1][1]
            if last_state < 120:  
               
                if last_si + 30 < len(sequences[i][1]):
                    
                    cys_si = dict(ali).get((104, 'm'), None)
                    if cys_si is not None:  # 104 found.

                       
                        cys_ai = ali.index(((104, 'm'), cys_si))

                       
                        _, re_states, re_details = run_hmmer([(sequences[i][0], sequences[i][1][cys_si+1:])],
                                                             bit_score_threshold=10)[0]


                        if re_states and re_states[0][-1][0][0] >= 126 and re_states[0][0][0][0] <= 117:

                           

                            vRegion = ali[:cys_ai+1]
                            jRegion = [(state, index+cys_si+1)
                                       for state, index in re_states[0] if state[0] >= 117]
                            cdrRegion = []
                            next = 105
                            for si in range(cys_si+1, jRegion[0][1]):
                                if next >= 116:
                                    cdrRegion.append(((116, 'i'), si))
                                else:
                                    cdrRegion.append(((next, 'm'), si))
                                    next += 1

                           
                            alignments[i][1][0] = vRegion + cdrRegion + jRegion
                            alignments[i][2][0]['query_end'] = jRegion[-1][1] + 1


def abn(sequences, scheme="imgt", database="ALL", output=False, outfile=None, csv=False, allow=set(["H", "K", "L", "A", "B", "G", "D"]),
           hmmerpath="", ncpu=None, assign_germline=False, allowed_species=['human', 'mouse'], bit_score_threshold=80):
    
    try:
        scheme = scheme_short_to_long[scheme.lower()]
    except KeyError:
        raise AssertionError(
            "Unrecognised or unimplemented scheme: %s" % scheme)

   
    if csv:
        assert outfile, 'If csv output is True then an outfile must be specified'
        _path, _ = os.path.split(outfile)
        assert (not _path) or os.path.exists(
            _path), 'Output directory %s does not exist' % _path

   
    alignments = run_hmmer(sequences, hmm_database=database, hmmerpath=hmmerpath,
                           ncpu=ncpu, bit_score_threshold=bit_score_threshold, hmmer_species=allowed_species)
   
    check_for_j(sequences, alignments, scheme)
    
    numbered, alignment_details, hit_tables = number_sequences_from_alignment(sequences, alignments, scheme=scheme, allow=allow,
                                                                              assign_germline=assign_germline,
                                                                              allowed_species=allowed_species)

   
    if output:
        if csv:
            csv_output(sequences, numbered, details, outfile)
        else:
            outto, close = sys.stdout, False
            if outfile:
                outto, close = open(outfile, 'w'), True
            abn_output(numbered, sequences, alignment_details, outto)
            if close:
                outto.close()

    return numbered, alignment_details, hit_tables



def run_abn(seq, ncpu=1, **kwargs):
    
    if isinstance(seq, list) or isinstance(seq, tuple):
        assert all(len(
            _) == 2 for _ in seq), "If list or tuple supplied as input format must be [ ('ID1','seq1'), ('ID2', 'seq2'), ... ]"
        sequences = seq
    elif os.path.isfile(seq):  # Fasta file.
        # Read the sequences. All are read into memory currently...
        sequences = read_fasta(seq)
        ncpu = int(max(1, ncpu))
    elif isinstance(seq, str):  # Single sequence
        validate_sequence(seq)
        ncpu = 1
        sequences = [["Input sequence", seq]]

    # Handle the arguments to abn.
    output = kwargs.get('output', False)
    outfile = kwargs.get('outfile', False)
    csv = kwargs.get('csv', False)
    if csv:  # Check output arguments before doing work.
        assert outfile, 'If csv output is True then an outfile must be specified'
        _path, _ = os.path.split(outfile)
        assert (not _path) or os.path.exists(
            _path), 'Output directory %s does not exist' % _path

    
    kwargs['ncpu'] = 1
    kwargs['output'] = False  # Overide and write the compiled results here.

    abn_partial = partial(abn, **kwargs)
    chunksize = math.ceil(float(len(sequences))/ncpu)

    
    if ncpu > 1:
        pool = Pool(ncpu)
        results = pool.map_async(
            abn_partial, grouper(chunksize, sequences)).get()
        pool.close()
    else:
        results = list(map(abn_partial, grouper(chunksize, sequences)))

   
    numbered = sum((_[0] for _ in results), [])
    alignment_details = sum((_[1] for _ in results), [])
    hit_tables = sum((_[2] for _ in results), [])


    if output:
        if csv:
            csv_output(sequences, numbered, alignment_details, outfile)
        else:
            outto, close = sys.stdout, False
            if outfile:
                outto, close = open(outfile, 'w'), True
            abn_output(numbered, sequences, alignment_details, outto)
            if close:
                outto.close()

    # Return the results
    return sequences, numbered, alignment_details, hit_tables


# Wrapper function for simple sequence in numbering and chain type out behaviour.
def number(sequence, scheme="imgt", database="ALL", allow=set(["H", "K", "L", "A", "B", "G", "D"]), allowed_species=['human', 'mouse']):
    

    try:
        validate_sequence(sequence)
        scheme = scheme_short_to_long[scheme.lower()]
    except KeyError:
        raise AssertionError(
            "Unrecognised to unimplemented scheme: %s" % scheme)

    # Length check. abn can number fragments of chains well. Encourage full domain numbering.
    if len(sequence) < 70:
        return False, False

    try:
        if not allowed_species:
            numbered, alignment_details, _ = abn(
                [("sequence_0", sequence)], scheme=scheme, database=database, output=False, allow=allow)
        else:
            numbered, alignment_details, _ = abn(
                [("sequence_0", sequence)], scheme=scheme, database=database, output=False, allow=allow, allowed_species=allowed_species)
    except AssertionError:  # Catch where the user has tried to number a TCR with an antibody scheme
        return False, False

    # We return the numbering list and the chain type where kappa and lambda chains are both "L" for light
    if numbered[0]:
        return numbered[0][0][0], chain_type_to_class[alignment_details[0][0]["chain_type"]]
    else:
        return False, False


if __name__ == "__main__":
    # Test and example useage of the abn function.
    sequences = [("12e8:H", "QHLEQSGGGAGGGLVKPGGSLELCCKASGFTFSSYYMCWVRQAPGKGLEWIGCIYAGSSGSRGNTYYGSWVNGRFTLSRDIDQSTGCLQLNSLTVADTAMHYCARLSWXXYWYSGWGTSGAQAP"),
                 ("12e8:L", "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASV"),
                 ("scfv:A", "DIQMTQSPSSLSASVGDRVTITCRTSGNIHNYLTWYQQKPGKAPQLLIYNAKTLADGVPSRFSGSGSGTQFTLTISSLQPEDFANYYCQHFWSLPFTFGQGTKVEIKRTGGGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFDFSRYDMSWVRQAPGKRLEWVAYISSGGGSTYFPDTVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCARQNKKLTWFDYWGQGTLVTVSSHHHHHH"),
                 ("lysozyme:A", "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL")]

    results = abn(sequences, scheme="chothia", output=True)
    numbering, alignment_details, hit_tables = results

    expect_one_VH_domain_numbering, expect_one_VL_domain_numbering, expect_VH_then_VL_numbering, expect_None = numbering
    assert len(expect_one_VH_domain_numbering) == 1
    assert len(expect_one_VL_domain_numbering) == 1
    assert len(expect_VH_then_VL_numbering) == 2
    assert expect_None == None
