# load(file.path("tests", "data-tests", "junctionAlignment.rda"))
load(file.path("..", "data-tests", "junctionAlignment.rda"))
load(file.path("..", "data-tests", "junctionAlignment_changeo.rda"))

germline_db <- c(
    "IGHV3-30*04"="CAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTAGCTATGCTATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGAT......GGAAGTAATAAATACTACGCAGACTCCGTGAAG...GGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAGAGA",
    "IGHD6-19*01"="GGGTATAGCAGTGGCTGGTAC",
    "IGHJ6*02"="ATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA",
    "IGHV4-39*01"="CAGCTGCAGCTGCAGGAGTCGGGCCCA...GGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGC......AGTAGTAGTTACTACTGGGGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGAGTATCTATTATAGT.........GGGAGCACCTACTACAACCCGTCCCTCAAG...AGTCGAGTCACCATATCCGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCAGACACGGCTGTGTATTACTGTGCGAGACA",
    "IGHD3-10*01"="GTATTACTATGGTTCGGGGAGTTATTATAAC",
    "IGHJ4*02"="ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
    "IGHV3-48*01"="GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGT......AGTAGTACCATATACTACGCAGACTCTGTGAAG...GGCCGATTCACCATCTCCAGAGACAATGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGA",
    "IGHD3-16*02"="GTATTATGATTACGTTTGGGGGAGTTATCGTTATACC",
    "IGHJ1*01"="GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG",
    "IGLV1-47*01"="CAGTCTGTGCTGACTCAGCCACCCTCA...GCGTCTGGGACCCCCGGGCAGAGGGTCACCATCTCTTGTTCTGGAAGCAGCTCCAACATC............GGAAGTAATTATGTATACTGGTACCAGCAGCTCCCAGGAACGGCCCCCAAACTCCTCATCTATAGGAAT.....................AATCAGCGGCCCTCAGGGGTCCCT...GACCGATTCTCTGGCTCCAAG......TCTGGCACCTCAGCCTCCCTGGCCATCAGTGGGCTCCGGTCCGAGGATGAGGCTGATTATTACTGTGCAGCATGGGATGACAGCCTGAGTGGTCC",
    "IGLJ1*01"="TTATGTCTTCGGAACTGGGACCAAGGTCACCGTCCTAG",
    "IGHV3-11*03"="CAGGTGCAGCTGTTGGAGTCTGGGGGA...GGCTTGGTCAAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTGACTACTACATGAGCTGGATCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGT......AGTAGTTACACAAACTACGCAGACTCTGTGAAG...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGTGCGAGA",
    "IGHV2-10*01"="CAGGTCACCTTGAAGGAGTCTGGTCCT...GCACTGGTGAAACCCACACAGACCCTCATGCTGACCTGCACCTTCTCTGGGTTCTCACTCAGC......ACTTCTGGAATGGGTGTGGGTTAGATCTGTCAGCCCTCAGCAAAGGCCCTGGAGTGGCTTGCACACATTTATTAGAAT.........GATAATAAATACTACAGCCCATCTCTGAAG...AGTAGGCTCATTATCTCCAAGGACACCTCCAAGAATGAAGTGGTTCTAACAGTGATCAACATGGACATTGTGGACACAGCCACACATTACTGTGCAAGGAGAC",
    "IGHD2-21*02"="AGCATATTGTGGTGGTGACTGCTATTCC",
    "IGHV3-23*01"="GAGGTGCAGCTGTTGGAGTCTGGGGGA...GGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTT............AGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAGTGGTAGT......GGTGGTAGCACATACTACGCAGACTCCGTGAAG...GGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAAGA",
    "IGHD1-26*01"="GGTATAGTGGGAGCTACTAC",
    "IGHV3-64*01"="GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTAGCTATGCTATGCACTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAATATGTTTCAGCTATTAGTAGTAAT......GGGGGTAGCACATATTATGCAAACTCTGTGAAG...GGCAGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTTCAAATGGGCAGCCTGAGAGCTGAGGACATGGCTGTGTATTACTGTGCGAGAGA",
    "IGHD5-12*01"="GTGGATATAGTGGCTACGATTAC",
    "IGHJ4*01"="ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG",
    "IGHV3-49*03"="GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTACAGCCAGGGCGGTCCCTGAGACTCTCCTGTACAGCTTCTGGATTCACCTTT............GGTGATTATGCTATGAGCTGGTTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTAGGTTTCATTAGAAGCAAAGCTTATGGTGGGACAACAGAATACGCCGCGTCTGTGAAA...GGCAGATTCACCATCTCAAGAGATGATTCCAAAAGCATCGCCTATCTGCAAATGAACAGCCTGAAAACCGAGGACACAGCCGTGTATTACTGTACTAGAGA",
    "IGLV3-10*01"="TCCTATGAGCTGACACAGCCACCCTCG...GTGTCAGTGTCCCCAGGACAAACGGCCAGGATCACCTGCTCTGGAGATGCATTGCCA..................AAAAAATATGCTTATTGGTACCAGCAGAAGTCAGGCCAGGCCCCTGTGCTGGTCATCTATGAGGAC.....................AGCAAACGACCCTCCGGGATCCCT...GAGAGATTCTCTGGCTCCAGC......TCAGGGACAATGGCCACCTTGACTATCAGTGGGGCCCAGGTGGAGGATGAAGCTGACTACTACTGTTACTCAACAGACAGCAGTGGTAATCATAG",
    "IGLJ2*01"="TGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTAG",
    "IGLV2-23*01"="CAGTCTGCCCTGACTCAGCCTGCCTCC...GTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCACTGGAACCAGCAGTGATGTTGGG.........AGTTATAACCTTGTCTCCTGGTACCAACAGCACCCAGGCAAAGCCCCCAAACTCATGATTTATGAGGGC.....................AGTAAGCGGCCCTCAGGGGTTTCT...AATCGCTTCTCTGGCTCCAAG......TCTGGCAACACGGCCTCCCTGACAATCTCTGGGCTCCAGGCTGAGGACGAGGCTGATTATTACTGCTGCTCATATGCAGGTAGTAGCACTTTAC",
    "IGLJ3*02"="TTGGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTAG",
    "IGKV2D-30*01"="GATGTTGTGATGACTCAGTCTCCACTCTCCCTGCCCGTCACCCTTGGACAGCCGGCCTCCATCTCCTGCAGGTCTAGTCAAAGCCTCGTATACAGT...GATGGAAACACCTACTTGAATTGGTTTCAGCAGAGGCCAGGCCAATCTCCAAGGCGCCTAATTTATAAGGTT.....................TCTAACTGGGACTCTGGGGTCCCA...GACAGATTCAGCGGCAGTGGG......TCAGGCACTGATTTCACACTGAAAATCAGCAGGGTGGAGGCTGAGGATGTTGGGGTTTATTACTGCATGCAAGGTACACACTGGCCTCC",
    "IGKJ1*01"="GTGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAAC",
    "IGKV2-30*01"="GATGTTGTGATGACTCAGTCTCCACTCTCCCTGCCCGTCACCCTTGGACAGCCGGCCTCCATCTCCTGCAGGTCTAGTCAAAGCCTCGTATACAGT...GATGGAAACACCTACTTGAATTGGTTTCAGCAGAGGCCAGGCCAATCTCCAAGGCGCCTAATTTATAAGGTT.....................TCTAACCGGGACTCTGGGGTCCCA...GACAGATTCAGCGGCAGTGGG......TCAGGCACTGATTTCACACTGAAAATCAGCAGGGTGGAGGCTGAGGATGTTGGGGTTTATTACTGCATGCAAGGTACACACTGGCCTCC"
)

test_that("junctionAlignment: counts deleted and in cdr3 germline nucleotides", {
    
    expect_warning(aln_db <- junctionAlignment(db, germline_db ))
    
    expect_equal(aln_db[['e3v_length']][1],0)
    expect_equal(aln_db[['e5d_length']][1],12)
    expect_equal(aln_db[['e3d_length']][1],1)
    expect_equal(aln_db[['e5j_length']][1],6)
    expect_equal(aln_db[['v_cdr3_length']][1],8)
    expect_equal(aln_db[['j_cdr3_length']][1],23)
    
    # plotJunctionAlignment(db[1,],germline_db)$plot
    
    # Has insertion 
    # IMGT:
    # Nucleotide insertions have been detected and automatically removed for this analysis
    # from nt position in user submitted sequence: 253
    # nb of inserted nt: 9
    # causing frameshift: no
    #
    # And deletion
    # from nt position in user submitted sequence: 181
    # nb of deleted nt	: 6
    # causing frameshift: no
    expect_equal(aln_db[['e3v_length']][2],2)
    expect_equal(aln_db[['e5d_length']][2],9)
    expect_equal(aln_db[['e3d_length']][2],5)
    expect_equal(aln_db[['e5j_length']][2],4)
    expect_equal(aln_db[['v_cdr3_length']][2],6)
    expect_equal(aln_db[['j_cdr3_length']][2],10)
    
    # plotJunctionAlignment(db[2,],germline_db)$plot 
    
    
    # Partial: NA junction
    expect_equal(aln_db[['e3v_length']][3],240)
    expect_equal(aln_db[['e5d_length']][3],13)
    expect_equal(aln_db[['e3d_length']][3],9)
    expect_equal(aln_db[['e5j_length']][3],19)
    expect_true(is.na(aln_db[['v_cdr3_length']][3]))
    expect_true(is.na(aln_db[['j_cdr3_length']][3]))
    
    # plotJunctionAlignment(db[3,],germline_db)$plot 
    
    
    # Partial short junction
    # Junction is one nt ("A")
    # IMGT V-quest: no results
    # Should MakeDb.py --partial fail sequences like this?
    # Is v_germline_start correct?
    expect_equal(aln_db[['e3v_length']][4],128)
    expect_equal(aln_db[['e5d_length']][4],26)
    expect_equal(aln_db[['e3d_length']][4],0)
    expect_equal(aln_db[['e5j_length']][4],24)
    expect_true(is.na(aln_db[['v_cdr3_length']][4]))
    expect_true(is.na(aln_db[['j_cdr3_length']][4]))
    
    # plotJunctionAlignment(db[4,],germline_db)$plot 
    
    
    # Partial
    # Nucleotide deletions have been detected in the CDR3
    # Check makeRegionDefinition works with this sequence
    expect_equal(aln_db[['e3v_length']][5],3)
    expect_true(is.na(aln_db[['e5d_length']][5]))
    expect_true(is.na(aln_db[['e3d_length']][5]))
    expect_equal(aln_db[['e5j_length']][5],3)
    expect_equal(aln_db[['v_cdr3_length']][5],24)
    expect_equal(aln_db[['j_cdr3_length']][5],4)
    

    # plotJunctionAlignment(db[5,],germline_db)$plot
    
    # IMGT, for these 2 sequences, IMGT/HighV-QUEST does not find a junction,
    # but IgBLAST does
    #
    # plotJunctionAlignment(db[6,],germline_db)$plot #IgBLAST
    # plotJunctionAlignment(db[9,],germline_db)$plot #IMGT
    # Can't create germline for AGTCCTTAGCCAACAAT_imgt (9), CreateGermlines gives
    # error "Germline sequence differs in length from input sequence by 22 characters."
    expect_equal(aln_db[['e3v_length']][6],2)
    expect_equal(aln_db[['e5d_length']][6],6)
    expect_equal(aln_db[['e3d_length']][6],3)
    expect_equal(aln_db[['e5j_length']][6],12)
    expect_equal(aln_db[['v_cdr3_length']][6],6)
    expect_equal(aln_db[['j_cdr3_length']][6],18)
    
    expect_equal(aln_db[['e3v_length']][9],2)
    expect_true(is.na(aln_db[['e5d_length']][9]))
    expect_true(is.na(aln_db[['e3d_length']][9]))
    expect_equal(aln_db[['e5j_length']][9],0)
    expect_true(is.na(aln_db[['v_cdr3_length']][9]))
    expect_true(is.na(aln_db[['j_cdr3_length']][9]))
    
    # For this nonproductive seq, the junction is all np1 nucleotides (22nt)
    # plotJunctionAlignment(db[7,],germline_db)$plot #IgBLAST
    # plotJunctionAlignment(db[10,],germline_db)$plot #IMGT
    expect_equal(aln_db[['e3v_length']][7],32)
    expect_true(is.na(aln_db[['e5d_length']][7]))
    expect_true(is.na(aln_db[['e3d_length']][7]))
    expect_equal(aln_db[['e5j_length']][7],7)
    expect_equal(aln_db[['v_cdr3_length']][7],0)
    expect_equal(aln_db[['j_cdr3_length']][7],0)
    
    
    # db$sequence_alignment[c(8,11)]
    # > db[c(8,11), "junction"]
    # A tibble: 2 x 1
    # junction                         
    # <chr>                            
    # 1 TGCTCCTCATACAGTAGCACTTATTGGGTGTTC
    # 2 TGCTCCTCATACAGTAGCACTTATTGGGTGTTC
    #  db$sequence[8]
    # [1] "GGGGGTCACAAGAGGCAGCGCTCTCGGGACGTCTCCACCATGGCCTGGGCTCTGCTGCTCCTCACTCTCCTCACTCAGGACACAGGGTCCTGGGCCCAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCACTGGAACCAGCAGTGATTTTGGGAATTATAATTTTGTCTCCTGGTATCAACAACACCCAGACAAAGCCCCCAAACTCATTATTTATGAGGGCAGTAAGCGGCCCTCAGGGGTTTCTAATCACTTCTCTGTCTCCAAGTCTGGCAACACGGCCTCCCTGACAATCTCTGGGCTCCAGGCTGACGACGAGGCTGATTATTACTGCTCCTCATACAGTAGCACTTATTGGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTGGGTCAGCCCAAGGCTGCCCCCTCG"
    #                                                                                                                                                                                                                                                                                                                                                                                             |- position 376
    # This sequence has IMGT gaps in the junction as shown in sequence_alignment, but only when aligned with IMGT.
    # The alignment gaps are not reported in junction
    # No alignment gaps in sequence_alignment when aligned with IgBLAST
    # IMGT/V-QUEST reports a 6nt indel in the input sequence that is located
    # in the junction. This is from position 376 in the input sequence.
    # The 6nt in the reference v germline are gcaggt. Can't create a germline, CreateGermlines.py 
    # gives the error "Germline sequence differs in length from input sequence by 6 characters."
    # "CAGTCTGCCCTGACTCAGCCTGCCTCC...GTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCACTGGAACCAGCAGTGATTTTGGG.........AATTATAATTTTGTCTCCTGGTATCAACAACACCCAGACAAAGCCCCCAAACTCATTATTTATGAGGGC.....................AGTAAGCGGCCCTCAGGGGTTTCT...AATCACTTCTCTGTCTCCAAG......TCTGGCAACACGGCCTCCCTGACAATCTCTGGGCTCCAGGCTGACGACGAGGCTGATTATTACTGCTCCTCATACAGTAGCACTTATTGGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCT"        
    # "CAGTCTGCCCTGACTCAGCCTGCCTCC...GTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCACTGGAACCAGCAGTGATTTTGGG.........AATTATAATTTTGTCTCCTGGTATCAACAACACCCAGACAAAGCCCCCAAACTCATTATTTATGAGGGC.....................AGTAAGCGGCCCTCAGGGGTTTCT...AATCACTTCTCTGTCTCCAAG......TCTGGCAACACGGCCTCCCTGACAATCTCTGGGCTCCAGGCTGACGACGAGGCTGATTATTACTGCTCCTCATAC......AGTAGCACTTATTGGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTGG"
    # plotJunctionAlignment(db[8,],germline_db)$plot
    # plotJunctionAlignment(db[11,],germline_db)$plot
    
    # This sequence has two . ('..') in the "junction" (#13, IMGT alignment)
    # and junction is too short
    # db$junction_length[12:13]
    # [1] 2 0
    # plotJunctionAlignment(db[12,],germline_db)$plot
    # plotJunctionAlignment(db[13,],germline_db)$plot
    
    expect_equal(aln_db[['e3v_length']][12],48)
    expect_true(is.na(aln_db[['e5d_length']][12]))
    expect_true(is.na(aln_db[['e3d_length']][12]))
    expect_equal(aln_db[['e5j_length']][12],3)
    # Junction is too short, no cdr3
    expect_true(is.na(aln_db[['v_cdr3_length']][12]))
    expect_true(is.na(aln_db[['j_cdr3_length']][12]))
    
    # No junction, expect all NA
    expect_equal(aln_db[['e3v_length']][13],0)
    expect_true(is.na(aln_db[['e5d_length']][13]))
    expect_true(is.na(aln_db[['e3d_length']][13]))
    expect_equal(aln_db[['e5j_length']][13],0)
    expect_true(is.na(aln_db[['v_cdr3_length']][13]))
    expect_true(is.na(aln_db[['j_cdr3_length']][13])) ## NOTE: alignment is off, due to changeo issue 178

    db_changeo <- junctionAlignment(db_changeo, 
                          germline_db,sequence_alignment="SEQUENCE_IMGT", 
                          v_call="V_CALL", d_call="D_CALL", j_call="J_CALL",
                          v_germline_start="V_GERM_START_IMGT", v_germline_end="V_GERM_END",
                          d_germline_start="D_GERM_START", d_germline_end="D_GERM_END",
                          j_germline_start="J_GERM_START", j_germline_end="J_GERM_END",
                          np1_length="NP1_LENGTH",np2_length="NP2_LENGTH",
                          junction="JUNCTION", junction_length="JUNCTION_LENGTH")
    
    expect_equal(db_changeo[['e3v_length']][2],3)
    expect_equal(db_changeo[['e5d_length']][2],6)
    expect_equal(db_changeo[['e3d_length']][2],7)
    expect_equal(db_changeo[['e5j_length']][2],9)
    expect_equal(db_changeo[['v_cdr3_length']][2],5)
    expect_equal(db_changeo[['j_cdr3_length']][2],19)
    
})
                                      