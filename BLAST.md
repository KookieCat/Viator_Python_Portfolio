```python
from Bio.Blast import NCBIWWW
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
NCBIWWW.email = "allyviator@gmail.com"
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
#Number indicates the id 
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("m_cold.fasta", format = "fasta")
#insted with the file
```


```python
print(record)
```

    ID: gi|8332116|gb|BE037100.1|BE037100
    Name: gi|8332116|gb|BE037100.1|BE037100
    Description: gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
    Number of features: 0
    Seq('CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...TTC')



```python
# * means in process due to processing
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
with open("m_cold.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
#Write result handle then close connection
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("***ALIGHTMENT***")
            print("sequence:", alignment.title)
            print("length", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
#printing an alignment for reaablility, sequence values, lengths, e values, the first 75 of query sequence, the matching sequence
#top shows the most likely match
```

    ***ALIGHTMENT***
    sequence: gi|2618480339|ref|XM_048479995.2| PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    length 1028
    e value: 3.60314e-72
    GAA-CAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGA...
    ||| || ||||||||||| |  ||| |||  ||||| |||| |||||||| |   |||  |||| |  ||||  |...
    GAAGCA-AAAATGGGGAG-G--ATGGAGTTTTTGGCTATGAGAACTGATCCA---GCCACGGCTGACTTGATAAA...
    ***ALIGHTMENT***
    sequence: gi|1227938481|ref|XM_022049453.1| PREDICTED: Carica papaya cold-regulated 413 plasma membrane protein 2-like (LOC110820077), mRNA
    length 1009
    e value: 1.30512e-66
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    ||||||||||||| | |  ||| || ||||| ||||| ||||||||   ||||   ||| || | |||  ||| |...
    AGAAAATGGGGAG-G-ATGGAA-TATTTGGCTATGAAGACTGATCA---GGCCACTGCTGATCTCATCACTTCTG...
    ***ALIGHTMENT***
    sequence: gi|2395983798|ref|XM_006466623.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X3, mRNA
    length 1052
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|1350315641|ref|XM_024180293.1| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X4, mRNA
    length 868
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|1350315634|ref|XM_006425717.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X1, mRNA
    length 952
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|1350315636|ref|XM_006425716.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X2, mRNA
    length 881
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|2395983796|ref|XM_025094967.2| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X1, mRNA
    length 980
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|2395983800|ref|XM_006466626.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X5, mRNA
    length 913
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|2395983799|ref|XM_006466625.3| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X4, mRNA
    length 978
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|1350315638|ref|XM_006425719.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X3, mRNA
    length 893
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|2395983797|ref|XM_006466624.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X2, mRNA
    length 968
    e value: 3.65448e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ***ALIGHTMENT***
    sequence: gi|2526866810|ref|XM_057645500.1| PREDICTED: Actinidia eriantha cold-regulated 413 plasma membrane protein 2-like (LOC130785340), mRNA
    length 1152
    e value: 1.31439e-61
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    ||||| ||| ||| |||| | || ||||| ||||| || |||| |  |||  || | ||| ||||| |||| || ...
    AAAAT-GGG-GAG-AATGGATTATTTGGCGATGAAGACCGATCCAGCGGC--TGCCGAAT-TGATCAATTCGGAC...
    ***ALIGHTMENT***
    sequence: gi|1768569081|ref|XM_031406607.1| PREDICTED: Pistacia vera cold-regulated 413 plasma membrane protein 2-like (LOC116120644), mRNA
    length 982
    e value: 3.68043e-57
    TCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATTACGGGTTT...
    || |||||||||||||||||  | ||  ||| || ||| || ||  | || | ||  ||||||  |   || |||...
    TCTGATATCAATGAGCTTAAGCTTGCTGCAAAGAAGCTTATTAACCACGCAACTAAACTCGGTGGTCTTGGCTTT...
    ***ALIGHTMENT***
    sequence: gi|1954740698|ref|XM_038867092.1| PREDICTED: Tripterygium wilfordii cold-regulated 413 plasma membrane protein 2 (LOC120014952), mRNA
    length 999
    e value: 3.68043e-57
    GAACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGAT...
    ||| ||||||| ||| ||| || | | || ||||| ||||| |||||||    ||  ||||   || |||||   ...
    GAAAAGAAAAT-GGG-GAG-AACGGATTATTTGGCGATGAAGACTGATC---CGGTTGTGGACGATTTGATCAGC...
    ***ALIGHTMENT***
    sequence: gi|2526831378|ref|XM_057642435.1| PREDICTED: Actinidia eriantha cold-regulated 413 plasma membrane protein 2-like (LOC130782944), mRNA
    length 1214
    e value: 3.68043e-57
    GAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA...
    |||||| ||| ||| |||| | || ||||| ||||| || |||| |   | | || | ||| ||||| |||||||...
    GAAAAT-GGG-GAG-AATGGATTATTTGGCGATGAAGACCGATCCAGCTG-C-TGCCGAAT-TGATCAATTCCGA...
    ***ALIGHTMENT***
    sequence: gi|1350280614|ref|XM_024170292.1| PREDICTED: Morus notabilis cold-regulated 413 plasma membrane protein 2 (LOC21394987), mRNA
    length 1020
    e value: 2.21505e-54
    TACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCA...
    || |||||||||||||| || | |   |||  |||| || ||||  |||| |||||||| ||||| || || || ...
    TATTTGGCCATGAAAACGGACCCA---GCCACGGCTGATTTGATAAATTCTGATATCAACGAGCTGAAGATCGCT...
    ***ALIGHTMENT***
    sequence: gi|2350562715|ref|XM_052350285.1| PREDICTED: Diospyros lotus cold-regulated 413 plasma membrane protein 2-like (LOC127810688), mRNA
    length 1265
    e value: 1.03056e-52
    GGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAAC...
    || || ||||||   |||||||||||| | |||| | |||| |||||||||||||||||||||||||||||||||...
    GGCTTCGGCACTTCCTTCCTCAAATGGATTGCCTCCTTTGCTGCTATTTACTTGTTGATATTGGATCGAACAAAC...
    ***ALIGHTMENT***
    sequence: gi|1585724761|ref|XM_028202722.1| PREDICTED: Camellia sinensis cold-regulated 413 plasma membrane protein 2-like (LOC114262355), mRNA
    length 910
    e value: 3.70656e-52
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAAT-TGGCCGTGGCTAATATGATCGATTCC...
    |||||||||||||  ||||| |||| ||||| ||||| || |||||    | ||    | |||  |||  |||||...
    AGAAAATGGGGAGGAAAATGGAGTATTTGGCAATGAAGACCGATCATCCAGCCCCAACCCAAT-CGATGAATTCC...
    ***ALIGHTMENT***
    sequence: gi|2583747300|ref|XM_059787294.1| PREDICTED: Cornus florida cold-regulated 413 plasma membrane protein 2-like (LOC132285128), mRNA
    length 1126
    e value: 1.33311e-51
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCC-GTGGCTAATATGATCGATTCC...
    ||||||| ||| |||||| | |||| ||||| |||||||||||||    ||||   | | ||| ||||| |||||...
    AGAAAAT-GGG-GAGAAA-GGAGTATTTGGCTATGAAAACTGATC---CGGCCACAGCCGAAT-TGATCAATTCC...
    ***ALIGHTMENT***
    sequence: gi|2118882425|ref|XM_044613294.1| PREDICTED: Mangifera indica cold-regulated 413 plasma membrane protein 2-like (LOC123198583), mRNA
    length 1083
    e value: 4.79473e-51
    CTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAACTGGAGAACCA...
    ||| ||| ||||||||| | |||| | |||| ||||||||||||||||||||||||||||||||||||||||| |...
    CTC-TTTTCTCAAATGGGTTGCCTCCTTTGCTGCTATTTACTTGTTGATATTGGATCGAACAAACTGGAGAACAA...
    ***ALIGHTMENT***
    sequence: gi|1220094463|ref|XM_021978417.1| PREDICTED: Prunus avium cold-regulated 413 plasma membrane protein 2-like (LOC110773902), mRNA
    length 896
    e value: 2.88569e-48
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    ||||  |||| || |||||||| || || || ||| | | ||  |||||||| ||||| | || ||| |||   |...
    TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAAAGATGCCACTAAGCTTGGTGGGT...
    ***ALIGHTMENT***
    sequence: gi|2537875835|ref|XM_021822463.2| PREDICTED: Hevea brasiliensis cold-regulated 413 plasma membrane protein 2 (LOC110663212), transcript variant X1, mRNA
    length 940
    e value: 2.88569e-48
    CGGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAA...
    ||| ||||| |||  ||| ||||||||| | ||||   |||| || || || ||| |||||||||||||||||||...
    CGGCTTTGGTACTTCTTTTCTCAAATGGGTTGCCTCTTTTGCTGCCATATATTTGCTGATATTGGATCGAACAAA...
    ***ALIGHTMENT***
    sequence: gi|1220047144|ref|XM_021953815.1| PREDICTED: Prunus avium cold-regulated 413 plasma membrane protein 2-like (LOC110753022), mRNA
    length 894
    e value: 2.88569e-48
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    ||||  |||| || |||||||| || || || ||| | | ||  |||||||| ||||| | || ||| |||   |...
    TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAAAGATGCCACTAAGCTTGGTGGGT...
    ***ALIGHTMENT***
    sequence: gi|2118854713|ref|XM_044609109.1| PREDICTED: Mangifera indica cold-regulated 413 plasma membrane protein 2-like (LOC123195424), mRNA
    length 742
    e value: 1.74907e-40
    GGGTTTGGCA-CTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAA...
    || ||||| | ||| ||| ||| ||||||| |||| | |||| |||||||| |||||||||||||||||||||||...
    GGCTTTGGAATCTC-TTTTCTCCAATGGCTTGCCTCCTTTGCTGCTATTTATTTGTTGATATTGGATCGAACAAA...
    ***ALIGHTMENT***
    sequence: gi|2537875838|ref|XM_021822464.2| PREDICTED: Hevea brasiliensis cold-regulated 413 plasma membrane protein 2 (LOC110663212), transcript variant X2, mRNA
    length 924
    e value: 4.93235e-31
    CGGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAA...
    ||| ||||| |||  ||| ||||||||| | ||||   |||| || || || ||| |||||||||||||||||||...
    CGGCTTTGGTACTTCTTTTCTCAAATGGGTTGCCTCTTTTGCTGCCATATATTTGCTGATATTGGATCGAACAAA...
    ***ALIGHTMENT***
    sequence: gi|2739574736|gb|CP157410.1| Ziziphus jujuba isolate TW-2024b chromosome 08
    length 27596186
    e value: 1.08288e-17
    TTACTTGTTGATATTGGATCGAACAAACTGGAGAACCAACATGCTCACGTCACTTTTAGTCCCTTACATATTCCT...
    ||| |||||||||||||||||||||||||||||||| || ||||| || ||||||||||| ||||| || ||| |...
    TTATTTGTTGATATTGGATCGAACAAACTGGAGAACAAATATGCTTACTTCACTTTTAGTTCCTTATATCTTCTT...



```python

```
