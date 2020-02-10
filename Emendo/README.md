Attached is a file containing guide RNAs of lengths 20bp-23bp, and their location in the human DNA hg38. <br/>
The location is specified in terms of chromosome (e.g. chr7), and starting position (e.g. 14256990). <br/>
The program (written in python or R) must receive the attached file as input, <br/>
and must generate a single fasta file containing amplicons for each guide RNA in the format: <br/>
 
\>unique_id_1 <br/>
amplicon_sequence_1 <br/>
...
\>unique_id_n <br/>
amplicon_sequence_n

The amplicon for each guide RNA must be extracted from the human genome hg38 automatically,  <br/>
with 200bp upstream and 200bp downstream the guide RNA. <br/>
Therefore, the total length of each amplicon sequence will be between 420bp – 423bp, depending on the guide RNA length. <br/>


<p align="center">
  <img width="460" height="300" src="guide_RNA.png">
</p>
