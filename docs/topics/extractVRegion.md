**extractVRegion** - *Extracts FWRs and CDRs from IMGT-gapped sequences*

Description
--------------------

`extractVRegion` extracts the framework and complementarity determining regions of 
the V segment for IMGT-gapped immunoglobulin (Ig) nucleotide sequences according to the 
IMGT numbering scheme.


Usage
--------------------
```
extractVRegion(sequences, region = c("FWR1", "CDR1", "FWR2", "CDR2",
"FWR3"))
```

Arguments
-------------------

sequences
:   character vector of IMGT-gapped nucleotide sequences.

region
:   string defining the region(s) of the V segment to extract. 
May be a single region or multiple regions (as a vector) from
`c("FWR1", "CDR1", "FWR2", "CDR2" ,"FWR3")`.  By default, all
regions will be returned.




Value
-------------------

If only one region is specified in the `region` argument, a character 
vector of the extracted sub-sequences will be returned. If multiple regions 
are specified, then a character matrix will be returned with columns 
corresponding to the specified regions and a row for each entry in 
`sequences`.


References
-------------------


1. Lefranc M-P, et al. IMGT unique numbering for immunoglobulin and T cell 
receptor variable domains and Ig superfamily V-like domains.
Dev Comp Immunol. 2003 27(1):55-77.




Examples
-------------------

```R
# Assign example clone
clone <- subset(ExampleDb, CLONE == 3138)

```

**Error in eval(e, x, parent.frame())**: object 'CLONE' not found
```R

# Get all regions
extractVRegion(clone$SEQUENCE_IMGT)

```

**Error in substr(sequences, IMGT_REGIONS[[x]][1], IMGT_REGIONS[[x]][2])**: object 'clone' not found
```R

# Get single region
extractVRegion(clone$SEQUENCE_IMGT, "FWR1")

```

**Error in substr(sequences, IMGT_REGIONS[[region]][1], IMGT_REGIONS[[region]][2])**: object 'clone' not found
```R

# Get all CDRs
extractVRegion(clone$SEQUENCE_IMGT, c("CDR1", "CDR2"))

```

**Error in substr(sequences, IMGT_REGIONS[[x]][1], IMGT_REGIONS[[x]][2])**: object 'clone' not found
```R

# Get all FWRs
extractVRegion(clone$SEQUENCE_IMGT, c("FWR1", "FWR2", "FWR3"))
```

**Error in substr(sequences, IMGT_REGIONS[[x]][1], IMGT_REGIONS[[x]][2])**: object 'clone' not found

See also
-------------------

IMGT-gapped region boundaries are defined in [IMGT_REGIONS](IMGT_REGIONS.md).






