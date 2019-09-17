**ChangeoClone-class** - *S4 class defining a clone*

Description
--------------------

`ChangeoClone` defines a common data structure for perform lineage recontruction
from Change-O data.






Slots
-------------------



`data`
:   data.frame containing sequences and annotations. Contains the
columns `SEQUENCE_ID` and `SEQUENCE`, as well as any additional 
sequence-specific annotation columns.

`clone`
:   string defining the clone identifier.

`germline`
:   string containing the germline sequence for the clone.

`v_gene`
:   string defining the V segment gene call.

`j_gene`
:   string defining the J segment gene call.

`junc_len`
:   numeric junction length (nucleotide count).




See also
-------------------

See [makeChangeoClone](makeChangeoClone.md) and [buildPhylipLineage](buildPhylipLineage.md) for use.






