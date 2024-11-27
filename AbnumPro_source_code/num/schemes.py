

# Alphabet used for insertion (last (-1th) is a blank space for no insertion)
alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ", "KK", "LL", "MM", "NN", "OO", "PP", "QQ", "RR", "SS", "TT", "UU", "VV", "WW", "XX", "YY", "ZZ", " "]

# Blosum62 matrix. Used in some annotation methods to recognise pre-defined motifs
blosum62 = {('B', 'N'): 3, ('W', 'L'): -2, ('G', 'G'): 6, ('X', 'S'): 0, ('X', 'D'): -1, ('K', 'G'): -2, ('S', 'E'): 0, ('X', 'M'): -1, ('Y', 'E'): -2, ('W', 'R'): -3, ('I', 'R'): -3, ('X', 'Z'): -1, ('H', 'E'): 0, ('V', 'M'): 1, ('N', 'R'): 0, ('I', 'D'): -3, ('F', 'D'): -3, ('W', 'C'): -2, ('N', 'A'): -2, ('W', 'Q'): -2, ('L', 'Q'): -2, ('S', 'N'): 1, ('Z', 'K'): 1, ('V', 'N'): -3, ('Q', 'N'): 0, ('M', 'K'): -1, ('V', 'H'): -3, ('G', 'E'): -2, ('S', 'L'): -2, ('P', 'R'): -2, ('D', 'A'): -2, ('S', 'C'): -1, ('E', 'D'): 2, ('Y', 'G'): -3, ('W', 'P'): -4, ('X', 'X'): -1, ('Z', 'L'): -3, ('Q', 'A'): -1, ('V', 'Y'): -1, ('W', 'A'): -3, ('G', 'D'): -1, ('X', 'P'): -2, ('K', 'D'): -1, ('T', 'N'): 0, ('Y', 'F'): 3, ('W', 'W'): 11, ('Z', 'M'): -1, ('L', 'D'): -4, ('M', 'R'): -1, ('Y', 'K'): -2, ('F', 'E'): -3, ('M', 'E'): -2, ('S', 'S'): 4, ('X', 'C'): -2, ('Y', 'L'): -1, ('H', 'R'): 0, ('P', 'P'): 7, ('K', 'C'): -3, ('S', 'A'): 1, ('P', 'I'): -3, ('Q', 'Q'): 5, ('L', 'I'): 2, ('P', 'F'): -4, ('B', 'A'): -2, ('Z', 'N'): 0, ('M', 'Q'): 0, ('V', 'I'): 3, ('Q', 'C'): -3, ('I', 'H'): -3, ('Z', 'D'): 1, ('Z', 'P'): -1, ('Y', 'W'): 2, ('T', 'G'): -2, ('B', 'P'): -2, ('P', 'A'): -1, ('C', 'D'): -3, ('Y', 'H'): 2, ('X', 'V'): -1, ('B', 'B'): 4, ('Z', 'F'): -3, ('M', 'L'): 2, ('F', 'G'): -3, ('S', 'M'): -1, ('M', 'G'): -3, ('Z', 'Q'): 3, ('S', 'Q'): 0, ('X', 'A'): 0, ('V', 'T'): 0, ('W', 'F'): 1, ('S', 'H'): -1, ('X', 'N'): -1, ('B', 'Q'): 0, ('K', 'A'): -1, ('I', 'Q'): -3, ('X', 'W'): -2, ('N', 'N'): 6, ('W', 'T'): -2, ('P', 'D'): -1, ('B', 'C'): -3, ('I', 'C'): -1, ('V', 'K'): -2, ('X', 'Y'): -1, ('K', 'R'): 2, ('Z', 'R'): 0, ('W', 'E'): -3, ('T', 'E'): -1, ('B', 'R'): -1, ('L', 'R'): -2, ('Q', 'R'): 1, ('X', 'F'): -1, ('T', 'S'): 1, ('B', 'D'): 4, ('Z', 'A'): -1, ('M', 'N'): -2, ('V', 'D'): -3, ('F', 'A'): -2, ('X', 'E'): -1, ('F', 'H'): -1, ('M', 'A'): -1, ('K', 'Q'): 1, ('Z', 'S'): 0, ('X', 'G'): -1, ('V', 'V'): 4, ('W', 'D'): -4, ('X', 'H'): -1, ('S', 'F'): -2, ('X', 'L'): -1, ('B', 'S'): 0, ('S', 'G'): 0, ('P', 'M'): -2, ('Y', 'M'): -1, ('H', 'D'): -1, ('B', 'E'): 1, ('Z', 'B'): 1, ('I', 'E'): -3, ('V', 'E'): -2, ('X', 'T'): 0, ('X', 'R'): -1, ('R', 'R'): 5, ('Z', 'T'): -1, ('Y', 'D'): -3, ('V', 'W'): -3, ('F', 'L'): 0, ('T', 'C'): -1, ('X', 'Q'): -1, ('B', 'T'): -1, ('K', 'N'): 0, ('T', 'H'): -2, ('Y', 'I'): -1, ('F', 'Q'): -3, ('T', 'I'): -1, ('T', 'Q'): -1, ('P', 'L'): -3, ('R', 'A'): -1, ('B', 'F'): -3, ('Z', 'C'): -3, ('M', 'H'): -2, ('V', 'F'): -1, ('F', 'C'): -2, ('L', 'L'): 4, ('M', 'C'): -1, ('C', 'R'): -3, ('D', 'D'): 6, ('E', 'R'): 0, ('V', 'P'): -2, ('S', 'D'): 0, ('E', 'E'): 5, ('W', 'G'): -2, ('P', 'C'): -3, ('F', 'R'): -3, ('B', 'G'): -1, ('C', 'C'): 9, ('I', 'G'): -4, ('V', 'G'): -3, ('W', 'K'): -3, ('G', 'N'): 0, ('I', 'N'): -3, ('Z', 'V'): -2, ('A', 'A'): 4, ('V', 'Q'): -2, ('F', 'K'): -3, ('T', 'A'): 0, ('B', 'V'): -3, ('K', 'L'): -2, ('L', 'N'): -3, ('Y', 'N'): -2, ('F', 'F'): 6, ('L', 'G'): -4, ('B', 'H'): 0, ('Z', 'E'): 4, ('Q', 'D'): 0, ('X', 'B'): -1, ('Z', 'W'): -3, ('S', 'K'): 0, ('X', 'K'): -1, ('V', 'R'): -3, ('K', 'E'): 1, ('I', 'A'): -1, ('P', 'H'): -2, ('B', 'W'): -4, ('K', 'K'): 5, ('H', 'C'): -3, ('E', 'N'): 0, ('Y', 'Q'): -1, ('H', 'H'): 8, ('B', 'I'): -3, ('C', 'A'): 0, ('I', 'I'): 4, ('V', 'A'): 0, ('W', 'I'): -3, ('T', 'F'): -2, ('V', 'S'): -2, ('T', 'T'): 5, ('F', 'M'): 0, ('L', 'E'): -3, ('M', 'M'): 5, ('Z', 'G'): -2, ('D', 'R'): -2, ('M', 'D'): -3, ('W', 'H'): -2, ('G', 'C'): -3, ('S', 'R'): -1, ('S', 'I'): -2, ('P', 'Q'): -1, ('Y', 'A'): -2, ('X', 'I'): -1, ('E', 'A'): -1, ('B', 'Y'): -3, ('K', 'I'): -3, ('H', 'A'): -2, ('P', 'G'): -2, ('F', 'N'): -3, ('H', 'N'): 1, ('B', 'K'): 0, ('V', 'C'): -1, ('T', 'L'): -1, ('P', 'K'): -1, ('W', 'S'): -3, ('T', 'D'): -1, ('T', 'M'): -1, ('P', 'N'): -2, ('K', 'H'): -1, ('T', 'R'): -1, ('Y', 'R'): -2, ('L', 'C'): -1, ('B', 'L'): -4, ('Z', 'Y'): -2, ('W', 'N'): -4, ('G', 'A'): 0, ('S', 'P'): -1, ('E', 'Q'): 2, ('C', 'N'): -3, ('H', 'Q'): 0, ('D', 'N'): 1, ('Y', 'C'): -2, ('L', 'H'): -3, ('E', 'C'): -4, ('Z', 'H'): 0, ('H', 'G'): -2, ('P', 'E'): -1, ('Y', 'S'): -2, ('G', 'R'): -2, ('B', 'M'): -3, ('Z', 'Z'): 4, ('W', 'M'): -1, ('Y', 'T'): -2, ('Y', 'P'): -3, ('Y', 'Y'): 7, ('T', 'K'): -1, ('Z', 'I'): -3, ('T', 'P'): -1, ('V', 'L'): 1, ('F', 'I'): 0, ('G', 'Q'): -2, ('L', 'A'): -1, ('M', 'I'): 1}


def smooth_insertions(state_vector):
    
    enforced_patterns = [ [(25,'m'),(26,'m'),( 27,'m'),( 28,'i')],
                          [(38,'i'),(38,'m'),(39,'m'),(40,'m')],
                          [(54,'m'),(55,'m'),(56,'m'),(57,'i')],
                          [(65,'i'),(65,'m'),(66,'m'),(67,'m')],
                          [(103,'m'),(104,'m'),(105,'m'),(106,'i')],
                          [(117,'i'),(117,'m'),(118,'m'),(119,'m')] ]

    

    state_buffer = [] 
    sv = []
    for (state_id, state_type ), si in state_vector:
        if state_id < 23: # Everything before the cysteine at 23.
            state_buffer.append( ((state_id, state_type ), si) )
            reg = -1   
        elif 25 <= state_id < 28: # Add to the buffer 
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 0
        elif 37 < state_id <= 40: # Add to the buffer 
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 1
        elif 54 <= state_id < 57: # Add to the buffer
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 2
        elif 64 < state_id <= 67: # Add to the buffer
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 3
        elif 103 <= state_id < 106: # Add to the buffer
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 4
        elif 116 < state_id <= 119: # Add to the buffer
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 5
        elif len(state_buffer) != 0: # Add the buffer and reset

           
            nins = sum( 1 for s in state_buffer if s[0][1] == 'i' ) 


            if nins > 0: # We have insertions

                if reg == -1: # FW1, only adjust if there are the same or more N terminal deletions than insertions
                    nt_dels = state_buffer[0][0][0] - 1 # Missing states
                    for (_id, _type ), _si in state_buffer: # Explicit deletion states.
                        if _type == 'd' or _si == None:
                            nt_dels +=1 
                        else: # First residue found
                            break
                    if nt_dels >= nins: # More n terminal deletions than insertions found. Likely misalignment.
                        
                       
                        new_states = [ s for s, _ in state_buffer if s[1] == 'm'] 
                        _first = new_states[0][0]

                        state_buffer = [ s for s in state_buffer if s[0][1] != 'd' ]

                        
                        _add = len( state_buffer ) - len( new_states ) 
                        assert _add >= 0, 'Implementation logic error' # Should be adding a positive number of positions
                        new_states = [ (_,'m') for _ in range( _first - _add, _first ) ] + new_states
                        assert len(new_states)==len(state_buffer), 'Implementation logic error' # Should have the same length

                        
                        for i in range( len(state_buffer ) ):
                            sv.append( ( new_states[i], state_buffer[i][1]) )
                    else:
                        sv += state_buffer 
                else:
                    
                    state_buffer = [ s for s in state_buffer if s[0][1] != 'd' ]
        
                    
                    if reg % 2: 
                        new_states = [enforced_patterns[reg][0]]*max( 0, len(state_buffer)-3) + enforced_patterns[reg][ max( 4-len(state_buffer), 1):]
                    else: 
                        new_states = enforced_patterns[reg][:3] + [enforced_patterns[reg][2]]*max( 0, len(state_buffer)-3)
                    
                    for i in range( len(state_buffer ) ):
                        sv.append( ( new_states[i], state_buffer[i][1]) )
                                
            else: 
                sv += state_buffer


            sv.append( ((state_id, state_type ), si) )


            state_buffer = [] 
            
        else: 
            sv.append( ((state_id, state_type ), si) )


    return sv


# General function to give annotations for regions that have direct mappings onto the hmm alignment (imgt states)
def _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions):
    

    state_vector = smooth_insertions( state_vector )
    
    _regions = [ [] for _ in range(n_regions) ]
    
   
    insertion = -1
    previous_state_id = 1
    previous_state_type = 'd'
    start_index, end_index  = None, None
    
    region = None


    for (state_id, state_type ), si in state_vector:
       
       
        if state_type != "i" or region is None: # BUG_FIX - JD 9/4/15 - do not allow a new region to start as an insertion.
            region = region_index_dict[region_string[state_id-1]] 

       
  
        if state_type == "m": # It is a match
            
            
            if state_string[state_id-1]=="I":
                if previous_state_type != 'd': 
                    insertion +=1 
                rels[region] -= 1 
            else: 
                insertion = -1 
            
                    
            _regions[region].append( ( (state_id + rels[region], alphabet[insertion] ), sequence[si]  ) )
            previous_state_id = state_id 
            if start_index is None:
                start_index = si
            end_index = si

            previous_state_type = state_type
            
        elif state_type == "i": 
            insertion +=1 
            

            _regions[region].append( ( (previous_state_id + rels[region], alphabet[insertion]), sequence[si]  ) )
            if start_index is None:
                start_index = si
            end_index = si

            previous_state_type = state_type

        else: 
            previous_state_type = state_type
            
            
            if state_string[state_id-1]=="I": 
                rels[region] -= 1 
                continue 
            
            insertion = -1 
            previous_state_id = state_id 


        
        if insertion >= 25 and region in exclude_deletions:
            insertion = 0 
        
        assert insertion < 25, "Too many insertions for numbering scheme to handle" # We ran out of letters.
            
    return _regions, start_index, end_index


def number_imgt(state_vector, sequence):
    
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
    
    region_string = '11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777'

    region_index_dict = {
                         "1":0,
                         "2":1,
                         "3":2,
                         "4":3,
                         "5":4,
                         "6":5,
                         "7":6
                         }
    

    rels              =  {0:0, 
                          1:0,
                          2:0,
                          3:0,
                          4:0,
                          5:0,
                          6:0,
                          7:0
                          }
    
    n_regions = 7

    exclude_deletions = [1,3,5]    
    
    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    

    _numbering = [ _regions[0], # Fw1
                   [],          # CDR1
                   _regions[2], # Fw2
                   [],          # CDR2
                   _regions[4], # Fw3
                   [],          # CDR3
                   _regions[6], # Fw4

                 ]

 
    cdr1seq    = "".join([ x[1] for x in _regions[1] if x[1] != "-" ])
    cdr1length = len(cdr1seq) 
    si = 0
    prev_state = 26
    for ann in get_imgt_cdr(cdr1length, 12, 27, 39):
        if not ann:
            _numbering[1].append( ((prev_state+1, ' '), '-') )
            prev_state += 1
        else:
            _numbering[1].append( (ann, cdr1seq[si]) )
            prev_state = ann[0]
            si += 1

 
    cdr2seq    = "".join([ x[1] for x in _regions[3] if x[1] != "-" ])
    cdr2length = len(cdr2seq)
    si = 0
    prev_state = 55
    for ann in get_imgt_cdr(cdr2length, 10, 56, 66):
        if not ann:
            _numbering[3].append( ((prev_state+1, ' '), '-') )
            prev_state += 1
        else:
            _numbering[3].append( (ann, cdr2seq[si]) )
            prev_state = ann[0]
            si += 1

    cdr3seq    = "".join([ x[1] for x in _regions[5] if x[1] != "-" ])
    cdr3length = len(cdr3seq)
    if cdr3length > 117: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    si = 0
    previous_state_id = 104
    for ann in get_imgt_cdr(cdr3length, 13, 105, 118):
        if ann is None:
            _numbering[5].append( ((previous_state_id+1, " "), "-"   ) )
            previous_state_id+=1
        else:
            _numbering[5].append( (ann, cdr3seq[si] ) )
            previous_state_id = ann[0]
            si+=1
  
 
    return gap_missing( _numbering ), startindex, endindex
  
def get_imgt_cdr(length, maxlength, start, end):

    annotations = [ None for _ in range(max(length, maxlength)) ]

    if length == 0:
        return annotations
    elif length == 1:
        annotations[0] = (start, ' ')
        return annotations

    front, back = 0, -1


    az = alphabet[:-1]
    za = az[::-1]

    for i in range(min(length, maxlength)):
        if i % 2:
            annotations[back] = (end + back, " ")
            back -= 1
        else:
            annotations[front] = (start + front, " ")
            front += 1


    centrepoint = [ i for i,v in enumerate(annotations) if v == None ]
    if not centrepoint:
        return annotations

    centre_left  = annotations[min(centrepoint)-1][0] # Get the index right before the first None
    centre_right = annotations[max(centrepoint)+1][0] # Get the index right after  the first None


    if not maxlength % 2:
        frontfactor, backfactor = maxlength//2, maxlength//2

    else:
        frontfactor, backfactor = (maxlength//2)+1, maxlength//2

    for i in range(max(0, length-maxlength)):
        if not i % 2:
            annotations[back] = (centre_right, za[back + backfactor])
            back -= 1
        else:
            annotations[front] = (centre_left, az[front - frontfactor])
            front += 1
    
    return annotations


def number_aho(state_vector, sequence, chain_type):

    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
   
    region_string =  'BBBBBBBBBBCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFHHHHHHHHHHHHHHHHIIIIIIIIIIIIIJJJJJJJJJJJJJKKKKKKKKKKK'
#                     1         2             3               4              5                   7               8            9            10


    region_index_dict = dict( list(zip( "ABCDEFGHIJK", list(range(11)) )) )
    
   
    rels              =  {0:0, 
                         1:0,
                         2:0,
                         3:0,
                         4:2,
                         5:2,
                         6:2,
                         7:2,
                         8:2,
                         9:2,
                         10:21}

    n_regions = 11
    
    exclude_deletions = [1,3,4,5,7,9]    
    
    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    

    _numbering = [ _regions[0], _regions[1], _regions[2],[], _regions[4], [], _regions[6], [], _regions[8],_regions[9],_regions[10] ]

    length = len( _regions[1] )
    if length > 0:
        start = _regions[1][0][0][0] 
        stretch_len = 10 - (start -1)
        if length > stretch_len: # Insertions are present. Place on 8
            annotations = [ (_," ") for _ in range(start,9) ] + [ (8,alphabet[_]) for _ in range( length - stretch_len ) ] + [(9," "),(10," ")]
        else:
            ordered_deletions = [(8," ")] + [(_," ") for _ in range(start, 11) if _ != 8]
            annotations = sorted( ordered_deletions[max(stretch_len-length, 0):] )
        _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in range(length) ]

    _L = 28,36,35,37,34,38,27,29,33,39,32,40,26,30,25,31,41,42
 
    _K = 28,27,36,35,37,34,38,33,39,32,40,29,26,30,25,31,41,42

    _H = 28,36,35,37,34,38,27,33,39,32,40,29,26,30,25,31,41,42 
 
    _A = 28,36,35,37,34,38,33,39,27,32,40,29,26,30,25,31,41,42    

    _B = 28,36,35,37,34,38,33,39,27,32,40,29,26,30,25,31,41,42

    _D = 28,36,35,37,34,38,27,33,39,32,40,29,26,30,25,31,41,42
                              
    _G = 28,36,35,37,34,38,27,33,39,32,40,29,26,30,25,31,41,42

    ordered_deletions = { 'L':_L,'K':_K, 'H':_H, 'A':_A, 'B':_B, 'D':_D, 'G':_G } 

    length = len( _regions[3] )

    annotations = [ (i, ' ') for i in sorted( ordered_deletions[chain_type][ max(18-length, 0): ] ) ]

    insertions = max( length-18 , 0 ) 
    if insertions > 26: 
        return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    elif insertions > 0:
 
        insertat = annotations.index( (36, ' ') )+1 # Always 12 
        assert insertat == 12, 'AHo numbering failed'  
        annotations = annotations[:insertat] + [ (36, alphabet[a]) for a in range( insertions ) ] + annotations[insertat:]

    _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in range(length) ]

    
    if chain_type == 'A':
        ordered_deletions = [74,73,63,62,64,61,65,60,66,59,67,58,68,69,70,71,72,75,76,77]
    else:
        ordered_deletions = [63,62,64,61,65,60,66,59,67,58,68,69,70,71,72,73,74,75,76,77]

    length = len(_regions[5])

    annotations = [ (i, ' ') for i in sorted( ordered_deletions[ max(20-length, 0): ] ) ]

    insertions = max( length-20 , 0 ) 
    if insertions > 26: 
        return [], startindex, endindex 
    elif insertions > 0:

        insertat = annotations.index( (63, ' ') )+1 
        assert insertat == 6, 'AHo numbering failed'  
        annotations = annotations[:insertat] + [ (63, alphabet[a]) for a in range( insertions ) ] + annotations[insertat:]

    _numbering[5] = [ (annotations[i], _regions[5][i][1]) for i in range(length) ]


    ordered_deletions = [86,85,87,84,88,83,89,82,90,81,91,80,92,79,93,78]
    length=len( _regions[7] )

    annotations = [ (i, ' ') for i in sorted( ordered_deletions[ max(16-length, 0): ] ) ]

    insertions = max( length-16 , 0 ) 
    if insertions > 26: 
        return [], startindex, endindex 
    elif insertions > 0:

        insertat = annotations.index( (85, ' ') )+1 
        assert insertat == 8, 'AHo numbering failed'  
        annotations = annotations[:insertat] + [ (85, alphabet[a]) for a in range( insertions ) ] + annotations[insertat:]

    _numbering[7] = [ (annotations[i], _regions[7][i][1]) for i in range(length) ]


    ordered_deletions = [123,124,122,125,121,126,120,127,119,128,118,129,117,130,116,131,115,132,114,133,113,134,112,135,111,
                         136,110,137,109,138,108,107]
    
    length=len( _regions[9] )

    annotations = [ (i, ' ') for i in sorted( ordered_deletions[ max(32-length, 0): ] ) ]


    insertions = max( length-32 , 0 ) 
    if insertions > 26: 
        return [], startindex, endindex 
    elif insertions > 0:

        insertat = annotations.index( (123, ' ') )+1 
        assert insertat == 17, 'AHo numbering failed'  
        annotations = annotations[:insertat] + [ (123, alphabet[a]) for a in range( insertions ) ] + annotations[insertat:]

    _numbering[9] = [ (annotations[i], _regions[9][i][1]) for i in range(length) ]

    numbering = gap_missing( _numbering )
    if len(numbering) > 0:
        if numbering[-1][0] == (148, ' ') and numbering[-1][1] != '-' and endindex+1 < len(sequence):
            numbering.append( ( (149, ' '), sequence[endindex+1]) )
            endindex +=1

    return numbering, startindex, endindex


def number_chothia_heavy(state_vector, sequence):

    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    

    region_string = '11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    

    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [0,2,4,6] # Don't put deletions in these regions

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    


    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]

    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 

    if insertions:
        start = _regions[0][0][0][0] 
       
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in range(start, 7) ] + [ (6, alphabet[_]) for _ in range(insertions) ] + [(7," "),(8," "),(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in range(length) ]
    else:
        _numbering[0] = _regions[0]

    
    length = len( _regions[2] )
    insertions = max(length - 11, 0) 

    if insertions:
        annotations = [(_, " ") for _ in range(23,32)] + [(31, alphabet[i]) for i in range(insertions) ] + [(32," "),(33," ")]
    else:
        annotations = [(_, " ") for _ in range(23,32)][:length-2] + [(32," "),(33," ")][:length]

    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in range(length) ]
 
    length = len( _regions[4] )

    insertions = max(length - 8, 0) 

    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in range(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]


    length = len( _regions[6] )    
    if length > 36: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="heavy")
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in range(length)  ]

    
    return gap_missing( _numbering ), startindex, endindex                                     


def number_chothia_light(state_vector, sequence):
    
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX'
                    
   
    region_string = '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    

    rels              =  {0:0, 
                         1: 0,
                         2:-6,
                         3:-6,
                         4:-13,
                         5:-16,
                         6:-20,
                         }    

    
    n_regions = 7
    
    exclude_deletions = [1,3,4,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    _numbering = [ _regions[0], [], _regions[2], [], _regions[4], [], _regions[6] ]
    
    

    length = len( _regions[1] )
    insertions = max(length - 11, 0) 
    annotations  =  [(24, " "),(25, " "), (26, " "), (27, " "), (28, " "),(29, " "),(30, " ")][:max(0,length)] 
    annotations += [(30, alphabet[i]) for i in range(insertions) ]
    annotations += [(31, " "),(32, " "),(33, " "),(34, " ")][ abs( min(0,length-11) ):] 
    _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in range(length) ]


    length = len( _regions[3] )
    insertions = max( length - 4, 0 )
    if insertions > 0:
        annotations  = [(51, " "),(52, " ")] + [(52, alphabet[i]) for i in range(insertions) ] + [(53, " "),(54, " ")]
        _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in range(length) ]
    else:
        _numbering[3] = _regions[3]
    

    length = len( _regions[4] )
    insertions = max(length - 34, 0)
    if insertions > 0: 
        annotations = [(i," ") for i in range(55,69)]+[(68, alphabet[i]) for i in range(insertions) ]+[(i," ") for i in range(69,89)]
        _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]
    elif length == 33: 
        annotations = [(i," ") for i in range(55,68)]+[(i," ") for i in range(69,89)]            
        _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]
    else: 
        _numbering[4] = _regions[4]


    length = len( _regions[5] )    

    if length > 35: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="light")
    _numbering[5]  = [ (annotations[i], _regions[5][i][1]) for i in range(length)  ]


    return gap_missing( _numbering ), startindex, endindex    


def number_kabat_heavy(state_vector, sequence):

    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    
    
    region_string = '11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    

    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [2,4,6]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    


    
    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]



    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 
    
    if insertions:
        start = _regions[0][0][0][0] 
        
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in range(start, 7) ] + [ (6, alphabet[_]) for _ in range(insertions) ] + [(7," "),(8," "),(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in range(length) ]
    else:
        _numbering[0] = _regions[0]
    
    

    length = len( _regions[2] )
    insertions = max(0,length - 13)
    annotations = [(_,' ') for _ in range(23, 36)][:length] 
    annotations += [(35, alphabet[i]) for i in range(insertions) ]
    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in range(length) ]
 

    length = len( _regions[4] )

    insertions = max(length - 8, 0) 

    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in range(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]

    length = len( _regions[6] )    
    if length > 36: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="kabat", chain_type="heavy") #  Chothia and Kabat the same here
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in range(length)  ]

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return gap_missing( _numbering ), startindex, endindex            
           
# Light chains    
def number_kabat_light(state_vector, sequence):
    
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX'
                    
    
    region_string = '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777'
    
    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    

    rels              =  {0:0, 
                         1: 0,
                         2:-6,
                         3:-6,
                         4:-13,
                         5:-16,
                         6:-20,
                         }    
    
    n_regions = 7
    
    exclude_deletions = [1,3,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    _numbering = [ _regions[0], [], _regions[2], [], _regions[4], [], _regions[6] ]
    
    

    length = len( _regions[1] )
    insertions = max(length - 11, 0) 
    annotations  =  [(24, " "),(25, " "), (26, " "), (27, " ")][:max(0,length)] 
    annotations += [(27, alphabet[i]) for i in range(insertions) ]
    annotations += [(28, " "),(29, " "),(30, " "),(31, " "),(32, " "),(33, " "),(34, " ")][ abs( min(0,length-11) ):] 
    _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in range(length) ]
  

    length = len( _regions[3] )
    insertions = max( length - 4, 0 )
    if insertions > 0:
        annotations  = [(51, " "),(52, " ")] + [(52, alphabet[i]) for i in range(insertions) ] + [(53, " "),(54, " ")]
        _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in range(length) ]
    else: 
        _numbering[3] = _regions[3]


    length = len( _regions[5] )    

    if length > 35: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="kabat", chain_type="light")
    _numbering[5]  = [ (annotations[i], _regions[5][i][1]) for i in range(length)  ]

    return gap_missing( _numbering ), startindex, endindex    




def number_martin_heavy(state_vector, sequence):
    
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    
  
    region_string = '11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'
    
    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    

    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [2,4,5,6]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    


    
    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]
    
    
    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 

    if insertions:
        start = _regions[0][0][0][0] 
       
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in range(start, 9) ] + [ (8, alphabet[_]) for _ in range(insertions) ] + [(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in range(length) ]
    else:
        _numbering[0] = _regions[0]

    

    length = len( _regions[2] )
    insertions = max(length - 11, 0) 
    if insertions:
        annotations = [(_, " ") for _ in range(23,32)] + [(31, alphabet[i]) for i in range(insertions) ] + [(32," "),(33," ")]
    else:
        annotations = [(_, " ") for _ in range(23,32)][:length-2] + [(32," "),(33," ")][:length]
    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in range(length) ]
 
    length = len( _regions[4] )

    insertions = max(length - 8, 0) 
    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in range(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]

    length = len( _regions[5] )
    insertions = max(length - 35, 0)
    if insertions > 0:
        annotations = [(i,' ') for i in range(58,73)]+[(72, alphabet[i]) for i in range(insertions) ]+[(i,' ') for i in range(73,93)]
        _numbering[5] = [ (annotations[i], _regions[5][i][1]) for i in range(length) ]
    else: 
        _numbering[4] = _regions[4]

     

    length = len( _regions[6] )    
    if length > 36: return [], startindex, endindex
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="heavy")
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in range(length)  ]

    return gap_missing( _numbering ), startindex, endindex                                    


def number_martin_light(state_vector, sequence):
    
    return number_chothia_light(state_vector,sequence)



def number_wolfguy_heavy(state_vector, sequence):
    
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

    
    region_string = '11111111111111111111111111222222222222223333333333333344444444444444444444555555555555555555555555555555666666666666677777777777'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    

    rels              =  {0:100, 
                         1:124,
                         2:160,
                         3:196,
                         4:226,
                         5:244,
                         6:283}    
    
    n_regions = 7
    
    exclude_deletions = [1,3,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    

    _numbering = [ _regions[0], [] , _regions[2], [], _regions[4] , [], _regions[6] ]

    ordered_deletions = [151]
    for p1,p2 in zip( list(range(152,176)), list(range(199, 175,-1))): ordered_deletions += [ p1,p2 ]
    length = len( _regions[1] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[1]  = [ ((annotations[i]," "), _regions[1][i][1]) for i in range(length)  ]
    

    ordered_deletions = [251]
    for p1,p2 in zip( list(range(252,271)), list(range(290, 271,-1))): ordered_deletions += [ p1,p2 ]
    ordered_deletions.append( 271 )
    ordered_deletions = list(range( 299, 290, -1)) + ordered_deletions
    length = len( _regions[3] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[3]  = [ ((annotations[i]," "), _regions[3][i][1]) for i in range(length)  ]

    ordered_deletions = []
    for p1,p2 in zip( list(range(356,374)), list(range(391, 373,-1))): ordered_deletions += [ p1,p2 ]
    ordered_deletions = [ 354, 394, 355, 393, 392 ] + ordered_deletions 
    ordered_deletions = [331,332] + [ 399, 398, 351, 352, 397, 353, 396, 395 ] + ordered_deletions
    length = len( _regions[5] )

    if length > len(ordered_deletions): return [], startindex, endindex 
    annotations = sorted(ordered_deletions[:length])
    _numbering[5]  = [ ((annotations[i]," "), _regions[5][i][1]) for i in range(length)  ]    
  
   
    return sum( _numbering, [] ), startindex, endindex   
            

def number_wolfguy_light(state_vector, sequence):
    
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
    
    region_string = '1111111AAABBBBBBBBBBBBB222222222222222223333333333333334444444444444455555555555666677777777777777777777888888888888899999999999'

    region_index_dict = {"1":0,"A":1,"B":2,"2":3,"3":4,"4":5,"5":6,"6":7,"7":8,"8":9,"9":10}
    

    rels              =  {0:500,
                         1:500,
                         2:500,    
                         3:527,
                         4:560,
                         5:595,
                         6:631,
                         7:630,
                         8:630,                                                  
                         9:646,
                         10:683}    
    
    n_regions = 11
    
    exclude_deletions = [1,3,5,7,9]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    

    _numbering = [ _regions[0], [], _regions[2], [] , _regions[4], [], _regions[6], [], _regions[8], [], _regions[10] ]


   
    length = len(_regions[1] )
    annotations = sorted([ (510,' '), (509, ' '), (508, ' ')][ :length ] + [(508,a) for a in alphabet[:max(0, length-3)]])  
    _numbering[1]  = [ (annotations[i], _regions[1][i][1]) for i in range(length)  ]
    

    length = len(_regions[3] )
    annotations = _get_wolfguy_L1( _regions[3], length)
    _numbering[3]  = [ ((annotations[i]," "), _regions[3][i][1]) for i in range(length)  ]
   
    ordered_deletions = []
    for p1,p2 in zip( list(range(652,673)), list(range(694, 672,-1))): ordered_deletions += [ p2,p1 ]
    ordered_deletions = [651] + list(range( 699, 694, -1)) + ordered_deletions + [673]

    length = len( _regions[5] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[5]  = [ ((annotations[i]," "), _regions[5][i][1]) for i in range(length)  ]



    length = len( _regions[7] )
    insertions = max( 0, length - 4 )
    annotations = [(711, ' '), (712, ' '), (713, ' '), (714, ' ')][:length] + [ (714, a) for a in alphabet[:insertions] ]    
    _numbering[7]  = [ (annotations[i], _regions[7][i][1]) for i in range(length)  ]  
    

    ordered_deletions = []
    for p1,p2 in zip( list(range(751,775)), list(range(799, 775,-1))): ordered_deletions += [ p1,p2 ]
    ordered_deletions.append( 775 )
  
    length = len( _regions[9] )
    if length > len(ordered_deletions): return [], startindex, endindex 
    annotations = sorted(ordered_deletions[:length])
    _numbering[9]  = [ ((annotations[i]," "), _regions[9][i][1]) for i in range(length)  ]  
  
    
    return sum( _numbering, [] ), startindex, endindex  


def _get_wolfguy_L1(seq, length):
    
    L1_sequences = {
    9: [['9',     'XXXXXXXXX', [551, 552, 554, 556, 563, 572, 597, 598, 599]]], 
    10: [['10',   'XXXXXXXXXX', [551, 552, 553, 556, 561, 562, 571, 597, 598, 599]]], 
    11: [['11a',  'RASQDISSYLA', [551, 552, 553, 556, 561, 562, 571, 596, 597, 598, 599]], 
         ['11b',  'GGNNIGSKSVH', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]], 
         ['11b.2','SGDQLPKKYAY', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]]], 
    12: [['12a',  'TLSSQHSTYTIE', [551, 552, 553, 554, 555, 556, 561, 563, 572, 597, 598, 599]], 
         ['12b',  'TASSSVSSSYLH', [551, 552, 553, 556, 561, 562, 571, 595, 596, 597, 598, 599]], 
         ['12c',  'RASQSVxNNYLA', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]], 
         ['12d',  'rSShSIrSrrVh', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]]], 
    13: [['13a',  'SGSSSNIGNNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 597, 598, 599]], 
         ['13b',  'TRSSGSLANYYVQ', [551, 552, 553, 554, 556, 561, 562, 563, 571, 572, 597, 598, 599]]], 
    14: [['14a',  'RSSTGAVTTSNYAN', [551, 552, 553, 554, 555, 561, 562, 563, 564, 571, 572, 597, 598, 599]], 
         ['14b',  'TGTSSDVGGYNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 596, 597, 598, 599]]], 
    15: [['15',   'XXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 594, 595, 596, 597, 598, 599]]], 
    16: [['16',   'XXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 594, 595, 596, 597, 598, 599]]], 
    17: [['17',   'XXXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 584, 594, 595, 596, 597, 598, 599]]]
    }    

    if length in L1_sequences: 
        curr_max = None, -10000
        for canonical in L1_sequences[length]:
            sub_score = 0
            for i in range( length ):
                try:
                    sub_score += blosum62[ (seq[i][1].upper(), canonical[1][i].upper() ) ]
                except KeyError:
                    sub_score += blosum62[ (canonical[1][i].upper(), seq[i][1].upper() ) ]
            if sub_score > curr_max[1]:
                curr_max = canonical, sub_score

        return curr_max[0][2]
    else: 
        ordered_deletions = []
        for p1,p2 in zip( list(range(551,575)), list(range(599, 575,-1))): ordered_deletions += [ p2,p1 ]
        ordered_deletions.append(575)
        return sorted( ordered_deletions[:length] )

def gap_missing( numbering ):
    
    num = [ ((0,' '),'-') ]
    for p, a in sum( numbering, [] ):  
        if p[0] > num[-1][0][0]+1:
            for _i in range( num[-1][0][0]+1, p[0] ):
                num.append( ((_i, ' '), '-' ) )
        num.append( (p,a) )
    return num[1:]



    
def get_cdr3_annotations(length, scheme="imgt", chain_type=""):
    
    az = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" 
    za = "ZYXWVUTSRQPONMLKJIHGFEDCBA"
    
    if scheme=="imgt":
        start, end = 105, 118 # start (inclusive) end (exclusive)
        annotations = [None for _ in range(max(length,13))]
        front = 0
        back  = -1
        assert (length-13) < 50, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        for i in range(min(length,13)):
            if i%2:
                annotations[back] = (end+back, " ")
                back -= 1
            else:
                annotations[front] = (start+front, " ")
                front += 1
        for i in range(max(0,length-13)): # add insertions onto 111 and 112 in turn
            if i%2:
                annotations[back] = (112, za[back+6])
                back-=1
            else:
                annotations[front] = (111, az[front-7])
                front +=1        
        return annotations

    elif scheme in [ "chothia", "kabat"] and chain_type=="heavy": # For chothia and kabat

        insertions = max(length - 10, 0)
        assert insertions < 27, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        ordered_deletions = [ (100, ' '), (99,' '), (98,' '), (97,' '), (96,' '), (95,' '), (101,' '),(102,' '),(94,' '), (93,' ') ]
        annotations = sorted( ordered_deletions[ max(0, 10-length): ] + [ (100,a) for a in az[:insertions ] ] )
        return annotations

    elif scheme in [ "chothia", "kabat"] and chain_type=="light":

        insertions = max(length - 9, 0)
        assert insertions < 27, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        ordered_deletions = [ (95,' '),(94,' '),(93,' '),( 92,' '),(91,' '),(96,' '),(97,' '),(90,' '),(89,' ') ]
        annotations = sorted( ordered_deletions[ max(0, 9-length): ] + [ (95,a) for a in az[:insertions ] ] )
        return annotations

    else:
        raise AssertionError("Unimplemented scheme")

