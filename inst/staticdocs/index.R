sd_section("Overview",
           "Package overview",
           c("alakazam"))

sd_section("File I/O",
           "Load and write Change-O DB files",
           c("readChangeoDb", "writeChangeoDb"))

sd_section("Sequence cleaning",
           "Sequence cleaning",
           c("maskSeqEnds", "maskSeqGaps", "collapseDuplicates"))

sd_section("Lineage reconstruction",
           "Build Ig lineages",
           c("makeChangeoClone", "buildPhylipLineage", "ChangeoClone-class"))

sd_section("Diversity analysis",
           "Diversity analysis",
           c("countClones", "estimateAbundance", "rarefyDiversity", "testDiversity",
             "calcCoverage", "calcDiversity", "plotAbundance", "plotDiversityCurve", 
             "DiversityCurve-class", "DiversityTest-class"))

sd_section("Ig and TCR sequence annotation",
           "Ig and TCR sequence annotation",
           c("countGenes", "extractVRegion", "getSegment"))

sd_section("Sequence distance calculation",
           "Sequence distance calculation",
           c("getSeqDistance", "getSeqMatrix", "testSeqEqual", 
             "getAAMatrix", "getDNAMatrix"))

sd_section("Amino acid propertes",
           "Amino acid propertes",
           c("translateDNA", "aminoAcidProperties", "countPatterns", "isValidAASeq",
             "aliphatic", "bulk", "charge", "gravy", "polar"))

sd_section("Data and constants",
           "Data and constants",
           c("IUPAC_CODES", "ABBREV_AA", "DEFAULT_COLORS", "IMGT_REGIONS", 
             "ExampleTrees"))