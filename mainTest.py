import Bio as b
import primerDesign as p
import sequenceAnalysis as sA



mySeqAnalysis = sA.FastAreader()

header = "THIS IS A TEST FILE"
vector = "ACTATACTATACATGGGTATGGGATGAGAGAGACACATATACACATA"
target = "ATGATCATAGACTGATGATTA"
enzyme1 ="EcoRI"
enzyme2 = "NotI"



myPrimerDesign = p.primerDesign(header, vector, target, enzyme1, enzyme2, 'ATG', "TAA")







