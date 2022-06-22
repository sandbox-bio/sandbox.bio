from collections import Counter
def solution_dna(dna):
    counts = Counter(dna)
    return f"{counts['A']} {counts['C']} {counts['G']} {counts['T']}"




def solution_rna(dna):
    return dna.replace('T', 'U')




def solution_revc(dna):
    COMPLEMENT = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    dna_rev = dna[::-1]
    dna_revc = [ COMPLEMENT[base] for base in dna_rev ]

    print(dna_revc)
    return ''.join(dna_revc)




def solution_fib(n, k):
    if n <= 1:
        return n
    else:
        return solution_fib(n-1, k) + k * solution_fib(n-2, k)




def solution_gc(fasta):
    # Calculate GC content
    gc_content_max = None
    gc_content_max_which = None

    lines = fasta.split("\n")
    lines.append('>done')

    sequence = None
    gc = 0
    n = 0
    for line in lines:
        if line.startswith('>'):
            if sequence is not None:
                if gc_content_max is None or gc / n > gc_content_max:
                    gc_content_max = gc / n
                    gc_content_max_which = sequence
            sequence = line.replace('>', '')
            n = 0
            gc = 0
        else:
            for base in line:
                if base == 'G' or base == 'C':
                    gc += 1
                n += 1

    # Output largest
    return f'{gc_content_max_which}\n{gc_content_max*100}'




def solution_hamm(s, t):
    dist = 0
    for i in range(len(t)):
        dist += 0 if s[i] == t[i] else 1

    return dist




def solution_iprb(k, m, n):
    pop = k + m + n
    prob = (4*(k*(k-1)+2*k*m+2*k*n+m*n)+3*m*(m-1))/(4*pop*(pop-1))
    return prob




def solution_prot(rna):
    rnaCodons = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"", "UAG":"",
    "UGU":"C", "UGC":"C", "UGA":"", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
    
    codons = []
    for i in range(0,len(rna),3):
        codons.append(rna[i:i+3])
    protein = ''.join([rnaCodons[codon] for codon in codons])
    return protein





def solution_subs(s, t):
    x = findMotif(s,t)
    answer = ''
    for y in x:
        answer = answer + str(y) + " "
    return answer

def findMotif(s,t):
    start = 0
    while True:
        start = s.find(t,start)
        if start == -1: return
        yield start + 1
        start += 1



