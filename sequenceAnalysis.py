#sequenceAnalysis
#Kyle Johnson (kyajohns)
#Partners : none



class ProteinParam:
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):
        """
        Initalizes the class with parameter of a protein string.  Creates a dictionary of the total amount of amino acids using the provided string.
        Also calculates the total number of amino acids within the dictionary for usage in later methods.
        """
        self.myProtein = protein  # Get the protein
        self.myProtein = self.myProtein.upper()  # Ensure every character from the string, if applicable, is converted to uppercase

        # Establish an amino acid composition dictionary (Numbers are integers)
        self.aaComposit = {
            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0,
        }
        # Putting Values into the dictionary from protein string
        for aa in self.myProtein:
            if aa in self.aaComposit:
                self.aaComposit[aa] += 1

        # Counting every amino acid held within the dicitonary (Would sum() work here?)
        self.aaCounter = 0
        for aa in self.aaComposit:
            self.aaCounter += self.aaComposit.get(aa)

    def aaCount(self):
        """
        Returns the aa count for the program
        """
        return (self.aaCounter)

    def pI(self, precision=2):
        """
        Takes charge and returns the pI of the associated amino acid sequence
        """
        High = 14  # Upper bound of pH
        low = 0  # Lower bound of pH
        midPoint = 0
        while ((High - low) > 10 ** (-precision)):
            midPoint = ((High + low) / 2)
            midCharge = self._charge_(midPoint)
            if midCharge > 0:
                low = midPoint
            else:
                High = midPoint

        return midPoint

    def aaComposition(self):
        """
        Returns the aaComposit dictionary
        """
        return self.aaComposit

    def _charge_(self, pH):
        """
        Calculates the total charge with the usage of an inputted pH
        """

        negCountAA = 0
        posCountAA = 0
        totalCharge = 0
        # Calculating the negative charge
        for aa in ProteinParam.aa2chargeNeg:
            negCountAA += self.aaComposit.get(aa) * ((10 ** pH) / (10 ** ProteinParam.aa2chargeNeg.get(aa) + 10 ** pH))

        # Calculating the positive charge
        for aa in ProteinParam.aa2chargePos:
            posCountAA += self.aaComposit.get(aa) * ((10 ** ProteinParam.aa2chargePos.get(aa)) / (
                        10 ** ProteinParam.aa2chargePos.get(aa) + 10 ** pH))

        # Generating the total charge, accounting for N and C terminus
        negCharge = negCountAA + ((10 ** pH) / (10 ** pH + 10 ** ProteinParam.aaCterm))  # Total of the negative charges
        posCharge = posCountAA + ((10 ** ProteinParam.aaNterm) / (
                    10 ** ProteinParam.aaNterm + 10 ** pH))  # Total of the positive charges
        totalCharge = posCharge - negCharge
        # print(cTerm, nTerm, negCountAA, posCountAA, totalCharge)
        return totalCharge

    def molarExtinction(self):
        """
        Calculates the molar extinction of the amino acid sequence using a pre-built table of requiste values
        """
        extinction = 0
        for aa in self.aaComposit:  # Scan through all of the amino acids in the dictionary
            if aa in ProteinParam.aa2abs280:  # Check to see if the valid amino acids appear
                extinction += self.aaComposit.get(aa) * ProteinParam.aa2abs280.get(
                    aa)  # Generate the exctinction number.
            else:
                pass
        return extinction

    def massExtinction(self):
        """
        Calculates mass Extinction through usage of molecular weight
        """
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):
        """
        Calculates the molecular weight of the amino acid sequence
        """
        totalWeight = 0
        for aa in self.aaComposit:
            totalWeight += self.aaComposit.get(aa) * (ProteinParam.aa2mw.get(aa) - ProteinParam.mwH2O)
        return (totalWeight + ProteinParam.mwH2O)





########################################################################################################################

class NucParams:
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }

    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}  # DNA Codon Table

    # Modifying pre-existing tables in order to generate a codon table for RNA and DNA

    # DNA Codon Composit Table
    # dnaCodonComp = {key.replace('U','T'): value for key, value in rnaCodonComp.items()}

    def __init__(self, inSeq = ' '):  # I don't think I even need a preliminary string sequence to start off, I can probably do without
        self.nucComposit = {n: 0 for n in "ACGTUN"}
        self.aaComposit = {aa: 0 for aa in "AGMSCHNTDIPVEKQWFLRY-"}
        self.codonComp = {
            # RNA Codon Composition Table, may need to rewrite this if it give me an error due to
            # U
            'UUU': 0, 'UCU': 0, 'UAU': 0, 'UGU': 0,  # UxU
            'UUC': 0, 'UCC': 0, 'UAC': 0, 'UGC': 0,  # UxC
            'UUA': 0, 'UCA': 0, 'UAA': 0, 'UGA': 0,  # UxA
            'UUG': 0, 'UCG': 0, 'UAG': 0, 'UGG': 0,  # UxG
            # C
            'CUU': 0, 'CCU': 0, 'CAU': 0, 'CGU': 0,  # CxU
            'CUC': 0, 'CCC': 0, 'CAC': 0, 'CGC': 0,  # CxC
            'CUA': 0, 'CCA': 0, 'CAA': 0, 'CGA': 0,  # CxA
            'CUG': 0, 'CCG': 0, 'CAG': 0, 'CGG': 0,  # CxG
            # A
            'AUU': 0, 'ACU': 0, 'AAU': 0, 'AGU': 0,  # AxU
            'AUC': 0, 'ACC': 0, 'AAC': 0, 'AGC': 0,  # AxC
            'AUA': 0, 'ACA': 0, 'AAA': 0, 'AGA': 0,  # AxA
            'AUG': 0, 'ACG': 0, 'AAG': 0, 'AGG': 0,  # AxG
            # G
            'GUU': 0, 'GCU': 0, 'GAU': 0, 'GGU': 0,  # GxU
            'GUC': 0, 'GCC': 0, 'GAC': 0, 'GGC': 0,  # GxC
            'GUA': 0, 'GCA': 0, 'GAA': 0, 'GGA': 0,  # GxA
            'GUG': 0, 'GCG': 0, 'GAG': 0, 'GGG': 0  # GxG
        }
        self.addSequence(inSeq)

    def addSequence(self, inSeq):
        '''
        Takes a sequence and places relevant data into the proper dictionaries.
        '''
        for nuc in inSeq:  # Scan through the sequence, no changes made to it
            if nuc in self.nucComposit:  # if the character is a part of nucleotide composit, add one of said character to dictionary
                self.nucComposit[nuc] += 1
            else:
                pass

        inSeqU = inSeq.replace("T", "U")

        for pos in range(0, len(inSeqU), 3):  # The primary loop of this method. pos moves in steps of three in order to collect each codon within the sequence
            tempHolder = inSeq[pos:pos + 3]
            if tempHolder in NucParams.rnaCodonTable:  # Check to see if the three characters represent a codon, place into proper dictionaries
                self.codonComp[tempHolder] += 1
                self.aaComposit[NucParams.rnaCodonTable[tempHolder]] += 1
            else:
                pass

    def aaComposition(self):
        '''
        returns the aaComposit dictionary
        '''
        return self.aaComposit

    def nucComposition(self):
        '''
        returns the nucComposit dictionary
        '''
        return self.nucComposit

    def codonComposition(self):
        '''
        returns codonComposit dictionary
        '''
        return self.codonComp

    def nucCount(self):
        '''
        returns a count of the total number of nucleotides
        '''
        return sum(self.nucComposit.values())






#########################################################################################################################

import sys

class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence



class ORFFinder():
    '''      '''



    def __init__(self, seq = " "):
        self.startCodons = ["ATG"]
        self.stopCodons = ["TGA", "TAG", "TAA"]
        self.ORF = []
        self.Orflen = []
        self.startPositions = []
        self.topORF(seq)
        self.bottomORF(seq)
        print(self.reverseComp(seq))


    def reverseComp(self, seq):
        '''reverseComp is used to generate a reverse complement from an inputted sequence.  This is accomplished through the usage of a dicitionary in order to avoid going over any prior replacements.'''

        complement = {"A":"T", "T":"A", "C":"G", "G":"C"}
        return "".join([complement[base] for base in seq[::-1]]) #Counting seq in reverse order

    def topORF (self, seq):
        '''   '''



        for frame in range (0,3):
            self.startPositions = [1]
            for nuc in range (frame, len(seq), 3):

                codon = seq[nuc:nuc+3]

                if codon in self.startCodons:
                    self.startPositions.append(nuc+1)
                    #print(self.startPositions)

                elif codon in self.stopCodons:
                    stopPos = nuc+3
                    for i in self.startPositions:
                        self.Orflen.append(stopPos+1-i)
                    self.ORF.append((frame+1, self.startPositions, stopPos, self.Orflen))
                    self.startPositions = []
                    self.Orflen = []

    def bottomORF (self, seq):
        '''    '''
        reverse = self.reverseComp(seq)
        for frame in range (0,3):
            self.startPositions = [1]
            for nuc in range (frame, len(reverse), 3):

                codon = reverse[nuc:nuc+3]

                if codon in self.startCodons:
                    self.startPositions.append(nuc+1)

                elif codon in self.stopCodons:
                    stopPos = nuc+3
                    for i in self.startPositions:
                        self.Orflen.append(stopPos+1-i)
                    self.ORF.append(((-(frame+1)), self.startPositions, stopPos, self.Orflen))
                    self.startPositions = []
                    self.Orflen = []


    def ORF_List(self):
        return self.ORF
