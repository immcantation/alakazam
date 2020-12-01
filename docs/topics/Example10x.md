**Example10x** - *Small example 10x Genomics Ig V(D)J sequences from CD19+ B cells isolated from PBMCs of a healthy 
human donor. Down-sampled from data provided by 10x Genomics under a Creative Commons Attribute license,
and processed with their Cell Ranger pipeline.*

Description
--------------------

Small example 10x Genomics Ig V(D)J sequences from CD19+ B cells isolated from PBMCs of a healthy 
human donor. Down-sampled from data provided by 10x Genomics under a Creative Commons Attribute license,
and processed with their Cell Ranger pipeline.


Usage
--------------------
```
Example10x
```




Format
-------------------

A data.frame with the following AIRR style columns:

+ `sequence_id`:                Sequence identifier
+ `sequence_alignment`:         IMGT-gapped observed sequence.
+ `germline_alignment`:         IMGT-gapped germline sequence.
+ `v_call`:                     V region allele assignments.
+ `d_call`:                     D region allele assignments.
+ `j_call`:                     J region allele assignments.
+ `c_call`:                     Isotype (C region) assignment.
+ `junction`:                   Junction region sequence.
+ `junction_length`:            Length of the junction region in nucleotides.
+ `np1_length`:                 Combined length of the N and P regions proximal
to the V region.
+ `np2_length`:                 Combined length of the N and P regions proximal
to the J region.
+ `umi_count`:                  Number of unique molecular identifies atttributed to sequence.
+ `cell_id`:                    Cell identifier.
+ `locus`:                      Genomic locus of sequence.



References
-------------------


1. Data source: https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_cd19_b
1. License: https://creativecommons.org/licenses/by/4.0/










