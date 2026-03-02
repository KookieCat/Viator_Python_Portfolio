 
# Sequence Objects
## Start of Pt.1


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





## Start of pt.2


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





## Start of Pt.3



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




## Start of pt.4


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

# Sequence Annotation
## Start of Pt.1


```python
from Bio.SeqRecord import SeqRecord
#Important for non-FASFA FILES
```


```python
#help(SeqRecord)
#Help file with code guides
```


```python
from Bio.Seq import Seq

```


```python
simple_seq = Seq("GATCG")
```


```python
simple_seq_r = SeqRecord(simple_seq)
#simple record tag
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATCG'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "ABC314"
```


```python
simple_seq_r.description = "Example Seq for Bio Python for Reccord code"
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATCG'), id='ABC314', name='<unknown name>', description='Example Seq for Bio Python for Reccord code', dbxrefs=[])




```python
simple_seq_r.seq
#Request just the sequence or any particular tag
```




    Seq('GATCG')




```python
simple_seq_r.description
```




    'Example Seq for Bio Python for Reccord code'




```python
simple_seq_r.annotations["Evidence"] = "Examples of spesific request notes"
```


```python
simple_seq_r.annotations
#Request all annotations
```




    {'Evidence': 'Examples of spesific request notes'}




```python
simple_seq_r.annotations["Evidence"]
#Request spesific annotation
```




    'Examples of spesific request notes'




```python
#File Used https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
reccord = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
reccord
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
reccord.annotations
#Returns empty if the requested data is not attached
```




    {}



## Start of pt.2



```python
#File Used https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
reccord = SeqIO.read("NC_005816.gb.txt", "gb")
```


```python
reccord
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
reccord.seq
#Any imformation can also be printed from gb files (GenBank)
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
reccord.annotations
#Gen Bank files have more data generally
```




    {'molecule_type': 'DNA',
     'topology': 'circular',
     'data_file_division': 'BCT',
     'date': '21-JUL-2008',
     'accessions': ['NC_005816'],
     'sequence_version': 1,
     'gi': '45478711',
     'keywords': [''],
     'source': 'Yersinia pestis biovar Microtus str. 91001',
     'organism': 'Yersinia pestis biovar Microtus str. 91001',
     'taxonomy': ['Bacteria',
      'Proteobacteria',
      'Gammaproteobacteria',
      'Enterobacteriales',
      'Enterobacteriaceae',
      'Yersinia'],
     'references': [Reference(title='Genetics of metabolic variations between Yersinia pestis biovars and the proposal of a new biovar, microtus', ...),
      Reference(title='Complete genome sequence of Yersinia pestis strain 91001, an isolate avirulent to humans', ...),
      Reference(title='Direct Submission', ...),
      Reference(title='Direct Submission', ...)],
     'comment': 'PROVISIONAL REFSEQ: This record has not yet been subject to final\nNCBI review. The reference sequence was derived from AE017046.\nCOMPLETENESS: full length.'}




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
#Tag attached to indicate a start after this point
```


```python
end_pos = SeqFeature.BetweenPosition(9, left = 8, right= 9)
#Tag attached to indicate the location is between 9 and 8
```


```python
start_pos2 = SeqFeature.BeforePosition(10)
#Tag attached to indicate the location is before 10
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
int(my_location.start)
#Request location tag included in end seq
```




    5





##Start of Pt.3


```python
from Bio.SeqRecord import SeqRecord
```


```python
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
                       "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK" 
                       "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
                       "SSAC"),
                   id="gi|14150838|gb|AAK54648.1| AF376133_1", 
                   description="chalcone synthase [Cucumis sativus]"
                  )
```


```python
print(record.format("fasta"))
#Print the sequence in the format of a fasta file
```

    >gi|14150838|gb|AAK54648.1| AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC
    



```python
print(record)
```

    ID: gi|14150838|gb|AAK54648.1| AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')



```python
record.Name = "Test"
#Add name to record
```


```python
record = SeqIO.read("NC_005816.gb.txt", "gb")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
#How many deatures in reccord
```




    41




```python
print(record.features[20])
# Print the 20th feature
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
sub_record = record[4300:4800]
#Creates a sub reccord which preserved the feature number, while adding an addional referance of loaction
```


```python
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
print(sub_record.features[1])
#Only location changed since referance is based on given in creation of sub reccord
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record.annotations["Type"] = "I Dont Know What to put"
```


```python
print(sub_record.annotations)
```

    {'molecule_type': 'DNA', 'Type': 'I Dont Know What to put'}



```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
#Edit to clarify this sequence is now partial
```


```python
print(sub_record.format("gb")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA              UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...




##Start of Pt.4



```python
record = SeqIO.read("NC_005816.gb.txt", "gb")
```


```python
record.dbxrefs
#Project #
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
shifted = record[2000:] + record[:2000]
#This DNA is circular bacterial, this notates it as a loop by adding the begining to the end
```


```python
len(record.features)
```




    41




```python
len(shifted.features)
```




    40




```python
shifted.dbxrefs
# The ref code is removed to stop a cascade from edits
```




    []




```python
shifted.dbxrefs = record.dbxrefs
#Copies the referance, also applicable to other data tags
```


```python
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
      #Print first value as string and the next 3 as integers
```

    NC_005816.1 9609 41 1 13



```python
rc = record.reverse_complement(id ="rc_record")
```


```python
print("%s %i %i %i %i" % (rc.id, len(record), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
# Referance removed due to changes, and annotations removed due to shifting (defalut removals)
```

    rc_record 9609 41 0 0



```python

```

314159
