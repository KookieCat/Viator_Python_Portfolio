```python
from Bio.Blast import NCBIWWW
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIXML
```


```python
NCBIWWW.email = "allyviator@gmail.com"
```


```python
record = SeqIO.read("DDX41.fasta", format = "fasta")
```


```python
print(record)
```

    ID: 5
    Name: 5
    Description: 5 dna:chromosome chromosome:GRCh38:5:177510969:177517579:-1
    Number of features: 0
    Seq('CTCTGTTAGAGCATTTACAGTAATAACCGTTTACATTTATACAGGGCCTCCCAC...GCA')



```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq, entrez_query = "Pan troglodytes[Organism]")
#For Ease I found the spesification for Chimps
```


```python
with open("DDX41.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
result_handle = open("DDX41.fasta")
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
            print("e value:", hsp.expect)
```

    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 1.14073e-171
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 3.72906e-165
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 7.70656e-117
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 1.39584e-75
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 1.39584e-75
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 2.07161e-73
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 1.07312e-70
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 4.27368e-63
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 1.49166e-62
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 7.727e-60
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 1.59406e-49
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 1.94196e-48
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 2.21578e-41
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 2.36789e-28
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 2.36789e-28
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 8.26474e-28
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 2.88468e-27
    ***ALIGHTMENT***
    sequence: gi|2697806900|ref|XM_063810170.1| PREDICTED: Pan troglodytes DEAD-box helicase 41 (DDX41), transcript variant X1, mRNA
    e value: 9.43006e-21
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 5.90912e-169
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 1.39584e-75
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 1.39584e-75
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 2.07161e-73
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 4.56303e-69
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 4.27368e-63
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 1.49166e-62
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 7.727e-60
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 1.59406e-49
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 1.94196e-48
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 2.21578e-41
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 2.36789e-28
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 2.36789e-28
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 8.26474e-28
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 2.88468e-27
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 9.43006e-21
    ***ALIGHTMENT***
    sequence: gi|525344024|ref|NM_001280484.1| Pan troglodytes DEAD-box helicase 41 (DDX41), mRNA >gi|343961728|dbj|AK305460.1| Pan troglodytes mRNA for probable ATP-dependent RNA helicase DDX41, complete cds, clone: PtsC-32-5_C06
    e value: 1.49562e-05
    ***ALIGHTMENT***
    sequence: gi|2697815253|ref|XM_063811727.1| PREDICTED: Pan troglodytes docking protein 3 (DOK3), transcript variant X3, mRNA
    e value: 5.94051e-55
    ***ALIGHTMENT***
    sequence: gi|2697815255|ref|XM_063811728.1| PREDICTED: Pan troglodytes docking protein 3 (DOK3), transcript variant X4, mRNA
    e value: 4.5711e-31
    ***ALIGHTMENT***
    sequence: gi|116256684|gb|AC145336.19| Pan troglodytes clone rp43-33e5, complete sequence
    e value: 0.000182204
    ***ALIGHTMENT***
    sequence: gi|109689856|gb|AC159041.2| Pan troglodytes BAC clone CH251-35C17 from chromosome unknown, complete sequence
    e value: 0.000182204
    ***ALIGHTMENT***
    sequence: gi|63004019|gb|AC145339.34| Pan troglodytes clone rp43-58e4, complete sequence
    e value: 0.000182204
    ***ALIGHTMENT***
    sequence: gi|126032540|gb|AC195476.2| Pan troglodytes BAC clone CH251-508P9 from chromosome x, complete sequence
    e value: 0.0022197
    ***ALIGHTMENT***
    sequence: gi|145699148|gb|AC202738.1| Pan troglodytes chromosome X clone CH251-574J3 map human ortholog p11.3, complete sequence
    e value: 0.0022197
    ***ALIGHTMENT***
    sequence: gi|116734893|gb|AC190426.2| Pan troglodytes chromosome X clone CH251-156N9 map human ortholog p11.3, complete sequence
    e value: 0.0022197
    ***ALIGHTMENT***
    sequence: gi|1039301578|gb|AC270861.1| Pan troglodytes chromosome 16 clone CH251-963D22, complete sequence
    e value: 0.0077475
    ***ALIGHTMENT***
    sequence: gi|1039301548|gb|AC270831.1| Pan troglodytes chromosome 16 clone CH251-172I21, complete sequence
    e value: 0.0077475
    ***ALIGHTMENT***
    sequence: gi|1603830877|gb|AC278947.1| Pan troglodytes chromosome 16 clone CH251-570O22, complete sequence
    e value: 0.0077475
    ***ALIGHTMENT***
    sequence: gi|1039301532|gb|AC270820.1| Pan troglodytes chromosome 16 clone CH251-497M24, complete sequence
    e value: 0.0077475
    ***ALIGHTMENT***
    sequence: gi|1039679336|gb|AC275224.1| Pan troglodytes chromosome 16 clone CH251-593G6, complete sequence
    e value: 0.0077475
    ***ALIGHTMENT***
    sequence: gi|53828902|gb|AC146099.3| Pan troglodytes BAC clone RP43-5L2 from chromosome 7, complete sequence
    e value: 0.0077475
    ***ALIGHTMENT***
    sequence: gi|2697780376|ref|XR_010154037.1| PREDICTED: Pan troglodytes uncharacterized LOC134809139 (LOC134809139), ncRNA
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|1031987448|gb|AC270728.1| Pan troglodytes chromosome 10 clone CH251-261D7, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|168693799|gb|AC198629.5| Pan troglodytes BAC clone CH251-378F1 from chromosome 7, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|6815214|dbj|AB037522.1| Pan troglodytes gene for natriuretic protein, partial cds
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|1031987305|gb|AC270585.1| Pan troglodytes chromosome 10 clone CH251-370D23, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|1031987264|gb|AC270544.1| Pan troglodytes chromosome 10 clone CH251-74F11, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|146335935|gb|AC189685.3| Pan troglodytes BAC clone CH251-647I11 from chromosome 7, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|148539829|gb|AC195407.2| Pan troglodytes BAC clone CH251-540P22 from chromosome 9, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|2697768642|ref|XR_010153774.1| PREDICTED: Pan troglodytes uncharacterized LOC134809088 (LOC134809088), transcript variant X2, ncRNA
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|118566294|gb|AC160941.2| Pan troglodytes BAC clone CH251-430N24 from chromosome 9, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|117962419|gb|AC187127.3| Pan troglodytes BAC clone CH251-121F11 from chromosome unknown, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|1935598432|gb|AC280428.1| Pan troglodytes BAC CH251-370F10 from chromosome, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|1603830817|gb|AC278887.1| Pan troglodytes chromosome 15 clone CH251-37E05, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|2697718022|ref|XM_063790861.1| PREDICTED: Pan troglodytes eukaryotic translation initiation factor 4E family member 2 (EIF4E2), transcript variant X6, mRNA
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|110796965|gb|AC187851.2| Pan troglodytes chromosome X clone RP43-41J12 map human ortholog p11.23, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|95147499|gb|AC183964.3| Pan troglodytes BAC clone CH251-654L7 from chromosome 12, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|2697718026|ref|XR_010149925.1| PREDICTED: Pan troglodytes eukaryotic translation initiation factor 4E family member 2 (EIF4E2), transcript variant X9, misc_RNA
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|2697768641|ref|XR_010153771.1| PREDICTED: Pan troglodytes uncharacterized LOC134809088 (LOC134809088), transcript variant X1, ncRNA
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|2697718024|ref|XM_054679527.2| PREDICTED: Pan troglodytes eukaryotic translation initiation factor 4E family member 2 (EIF4E2), transcript variant X7, mRNA
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|2697718021|ref|XR_010149923.1| PREDICTED: Pan troglodytes eukaryotic translation initiation factor 4E family member 2 (EIF4E2), transcript variant X3, misc_RNA
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|149893946|gb|AC193446.3| Pan troglodytes BAC clone CH251-699J14 from chromosome 9, complete sequence
    e value: 0.0270414
    ***ALIGHTMENT***
    sequence: gi|102142635|gb|AC183688.2| Pan troglodytes BAC clone CH251-6F11 from chromosome 1, complete sequence
    e value: 0.0270414



```python
#The E-Value when compared to Chimpanzee is 1.14073e-171
```


```python

```
