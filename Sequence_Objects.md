```python
#Start of Pt.1
```


```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("CTAGC")

#Createing Sequence with the above string
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
    
#Each letter is given a positional number
```

    0 C
    1 T
    2 A
    3 G
    4 C



```python
print(len(my_seq))
#Length of Sequence 
```

    5



```python
print(my_seq[0])
#Print inital value in sequence
```

    C



```python
print(my_seq[4])
#Print 5th Value
```

    C



```python
Seq("AAAA").count("AA")
#Counts the non-overlap occuraces of AA
```




    2




```python
my_seq2 = Seq("AAAGCTCCGACGTGCTAGTCAGTCGTCAGTCAGCTTCA")
```


```python
len(my_seq2)
# Shortened print format
```




    38




```python
my_seq2.count("A")
#Count occurances of A
```




    9




```python
100 * (my_seq2.count("G") + my_seq2.count("C")) / len(my_seq2)
#Formula for gc fraction of sequence
```




    52.63157894736842




```python
from Bio.SeqUtils import gc_fraction
```


```python
gc_fraction(my_seq2)
#gc fraction calclation simplified
```




    0.5263157894736842




```python
my_seq2[0::3]
# Start at inital and print every 3rd value
```




    Seq('AGCATTTGGACCC')




```python
my_seq2[2::4]
# Start at 3rd value and print every 4th value
```




    Seq('ACCCTTCCT')




```python
my_seq2[2:7]
#Print the 5th through 7th values
```




    Seq('AGCTC')




```python
my_seq2[::-1]
#Print backwards
```




    Seq('ACTTCGACTGACTGCTGACTGATCGTGCAGCCTCGAAA')




```python
fasta_format_string = ">Name\n%s\n" % my_seq2
# Give name tag to sequence in fasta format
```


```python
print(fasta_format_string)
```

    >Name
    AAAGCTCCGACGTGCTAGTCAGTCGTCAGTCAGCTTCA
    



```python
seq3 = Seq("AGTC")
seq4 = Seq("CTAGA")
```


```python
seq3 + seq4
# Added together first infront
```




    Seq('AGTCCTAGA')




```python
#Start of pt.2
```


```python
contigs = [Seq("TAG"), Seq("ATGCTAG"), Seq("TGAGTATGAGTC")]
#Sequence made of chunks
```


```python
spacer = Seq("N" *10)
# N is notation of an unsure Neucleotide 
```


```python
spacer.join(contigs)
#Join the chunks with spacer sequence
```




    Seq('TAGNNNNNNNNNNATGCTAGNNNNNNNNNNTGAGTATGAGTC')




```python
dna_seq = Seq("AtgCCATGaTCcGCAT")
```


```python
dna_seq = dna_seq.lower()
#Change to Upper Case 
```


```python
dna_seq = dna_seq.upper()
#Change to Upper Case 
```


```python
"AtgC" in dna_seq
# Case Sensitive
```




    False




```python
"ATGC" in dna_seq
```




    True




```python
print(my_seq2)
print(my_seq2.complement())
print(my_seq2.reverse_complement())
# Compliment and reverse complement of a given sequence
```

    AAAGCTCCGACGTGCTAGTCAGTCGTCAGTCAGCTTCA
    TTTCGAGGCTGCACGATCAGTCAGCAGTCAGTCGAAGT
    TGAAGCTGACTGACGACTGACTAGCACGTCGGAGCTTT



```python
protien_seq = Seq("EVRNAK")
protien_seq.complement()
# Compliment is nonesense 
```




    Seq('EBYNTM')




```python
coding_dna = Seq("AGCTAGCGTCAGCTATCGCTACTATCGCAGCAGCTGCTACTACTAGTCGCATCTGATCACGTGCTAGCTGTC")
```


```python
template_dna = coding_dna.reverse_complement()
```


```python
messenger_rna = coding_dna.transcribe()
```


```python
#Raw DNA
print(template_dna)
#Coding DNA (Reverse Complement of Template)
print(coding_dna)
#mRNA (Transcribed Coding DNA)
print(messenger_rna)
```

    GACAGCTAGCACGTGATCAGATGCGACTAGTAGTAGCAGCTGCTGCGATAGTAGCGATAGCTGACGCTAGCT
    AGCTAGCGTCAGCTATCGCTACTATCGCAGCAGCTGCTACTACTAGTCGCATCTGATCACGTGCTAGCTGTC
    AGCUAGCGUCAGCUAUCGCUACUAUCGCAGCAGCUGCUACUACUAGUCGCAUCUGAUCACGUGCUAGCUGUC



```python
#mRNA to Coding DNA
print(messenger_rna.back_transcribe())
#Coding DNA to Raw DNA
print(coding_dna.reverse_complement())
```

    AGCTAGCGTCAGCTATCGCTACTATCGCAGCAGCTGCTACTACTAGTCGCATCTGATCACGTGCTAGCTGTC
    GACAGCTAGCACGTGATCAGATGCGACTAGTAGTAGCAGCTGCTGCGATAGTAGCGATAGCTGACGCTAGCT



```python
messenger_rna.translate()
# * is a stop codon
```




    Seq('S*RQLSLLSQQLLLLVASDHVLAV')




```python
#Start of Pt.3
```


```python
coding_dna.translate(table="Vertebrate Mitochondrial")
coding_dna.translate(table= 2)
#Uses table for mitocondria the full word or 2 
# Also "Bacterial" table
```




    Seq('S*RQLSLLSQQLLLLVASDHVLAV')




```python
coding_dna.translate(to_stop = True)
# Codes to stop
coding_dna.translate(to_stop = True, table =2, stop_symbol = "!")
# Multiple clauses can be used
```




    Seq('S')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table= "Bacterial")
#Translating bacterial DNA
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
gene.translate(table= "Bacterial", to_stop = True, stop_symbol = "!")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
gene.translate(table= "Bacterial", to_stop = True, cds = True)
#Indicates the gene is complete so the inital codons are interpreted as M insted of V
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)

```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
#Request for stop codons for mitochondria
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
seq1 = Seq("ACGT")
```


```python
"ACGT" == seq1
#Does seq1 = ACGT?
```




    True




```python
#Start of pt.4
```


```python
unknown_seq = Seq(None, 10)
#returns information, length notated as 10 BP
```


```python
len(unknown_seq)
```




    10




```python
seq2 = Seq({117512683: "TTGAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
 # Data of Location, BP, and total length 
```


```python
seq2[117512687:117512693]
```




    Seq('AACCTG')




```python
seq2[117512687:]
#Examples of partially defined sequences
```




    Seq({0: 'AACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833286)




```python
seq1 + unknown_seq + seq1
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
#most sequances immutable dissallowing editing, but mutable functions can be imported
```


```python
from Bio.Seq import MutableSeq
```


```python
Mutable_seq = MutableSeq("TTGAAACCTGAATGTGAGAGTCAGTCAAGGATAGT")
```


```python
Mutable_seq[5] = "C"
#Replaces 5th BP with C
```


```python
print(Mutable_seq)
```

    TTGAACCCTGAATGTGAGAGTCAGTCAAGGATAGT



```python
Mutable_seq.remove("T")
#removes first T
```


```python
new_seq = Seq(Mutable_seq)
#Reprotects Seq
```


```python

```


```python

```
