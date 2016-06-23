import RNAstructure
import alignment_parser as a

def energy(seq, db):
    if '&' in seq:
        return "Energy: " + str(duplex_energy(str(seq),str(db)))
    else:
        return "Energy: " + str(single_energy(str(seq),str(db)))

def single_energy(seq, db):
    print "seq",seq
    print "db", db
    r = RNAstructure.RNA.fromString(seq)
    pairs = a.dot2structure(db)
    for i in r.iterIndices():
        if pairs[i-1] is not None:
            ip = pairs[i-1]+1
            r.SpecifyPair(i, pairs[i-1]+1)
    return r.CalculateFreeEnergy(UseSimpleMBLoopRules=True)

def duplex_energy(seq,db):
    seq1, seq2 = seq.split("&")
    db1, db2 = db.split("&")
    r = RNAstructure.HybridRNA.fromString(seq1,seq2)
    r.FoldDuplex()
    #have to call FoldDuplex because there's a bug in HybridRNA constructure
    #underlying structure is not allocated unless you use some method to fold it:
    r.RemovePairs()

    seq1_open = []
    seq1_pairs = []
    seq2_open = []
    seq2_pairs = []
    intermolecular_pairs = []

    for i,c in enumerate(db1):
        if c == '(':
            seq1_open.append(i)
        elif c == '.':
            continue
        elif c==')':
            seq1_pairs.append((seq1_open.pop(),i))
        else:
            raise ValueError, "invalid value %s in dot-bracket structure"%c
    for i,c in enumerate(db2):
        if c == '(':
            seq2_open.append(i)
        elif c == '.':
            continue
        elif c==')':
            if len(seq2_open)>0:
                seq2_pairs.append((seq2_open.pop(),i))
            else:
                intermolecular_pairs.append((seq1_open.pop(),i))
        else:
            raise ValueError, "invalid value %s in dot-bracket structure"%c

    for i,j in seq1_pairs:
        r.SpecifyPair(i+1,j+1)

    for i,j in seq2_pairs:
        r.SpecifyPair(i+len(seq1)+1+3, j+len(seq1)+1+3)

    for i,j in intermolecular_pairs:
        r.SpecifyPair(i+1, j+len(seq1)+1+3)

    return r.CalculateFreeEnergy(UseSimpleMBLoopRules=True)

