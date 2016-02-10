sd_section("File I/O",
           "Load and write Change-O DB files",
           c("readChangeoDb", "writeChangeoDb"))

sd_section("Sequence cleaning",
           "Sequence cleaning",
           c("maskSeqEnds", "maskSeqGaps", "collapseDuplicates"))

sd_section("Lineage reconstruction",
           "Build Ig lineages",
           c("makeChangeoClone", "buildPhylipLineage"))

sd_section("Diversity analysis",
           "Diversity analysis",
           c("countClones", "estimateAbundance", "rarefyDiversity", "testDiversity",
             "plotAbundance", "plotDiversityCurve"))

sd_section("Ig and TCR sequence annotation",
           "Ig and TCR sequence annotation",
           c("countGenes", "extractVRegion", "getSegment"))

sd_section("Sequence distance calculation",
           "Sequence distance calculation",
           c("getSeqDistance", "getSeqMatrix", "testSeqEqual"))

sd_section("Amino acid propertes",
           "Amino acid propertes",
           c("translateDNA", "aminoAcidProperties", "countPatterns"))
