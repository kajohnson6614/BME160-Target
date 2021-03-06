#! /usr/bin/python3
#Britney Hernandez, Kyle Johnson
#primerDesign


import Bio
import re
import sequenceAnalysis as sa
from Bio import Restriction as r
from Bio.SeqUtils import MeltingTemp as mt
import sys




class primerDesign():
    '''Class takes in a fastA file with header and sequence of target.  Also takes in an enzyme 1 and an enzyme 2.  Has the possibility of taking in an alternative start and stop codon.
        Class then checks to determine if enzymes exist within built in dictionary and if target contains cut sites.  Should those conditions be met the class displays a message to the user and promptly exits.
        Class then determines enzyme order in relation to the built in vector (will be changed in the future to user inputted vector).  Class then determines a spliced version of the vector and finally creates the forward and reverse primers.
        '''


    def __init__(self, head, target='', enzyme1='', enzyme2='', startCodon = '', stopCodon= ''):


        ###########################################################################################
        #Save the vector, target, and enzymes
        self.header = head
        self.vector = 'GAAAACCTGTATTTTCAGGGCGCCATGGATCCGGAATTCAAAGGCCTACGTCGACGAGCTCAACTAGTGCGGCCGCACTCGAGCACCACCACCACCACCACTGAGATC'
        self.target = target
        self.reverseCompTarget = self.buildReverseComp(self.target)
        self.restrictionEnzymeDict = {'AccI': ['CCG', '>90', '20', '2.1', '37'], 'AflIII': ['CC', '>90', '2', '3.1', '37'], 'AscI': ['A', '>90', '2', '1.1, 2.1, 3.1 have equal activity', '37'], 'AvaI': ['CC', '>90', '2', '2.1', '37'],
                                      'BamHI': ['CG', '>90', '2', '3.1', '37'], 'BglII': ['GA', '>90', '20', '3.1', '37'], 'BssHII': ['TTG', '>90', '20', '1.1, 2.1, 3.1 have equal activity', '50'], 'BstXI': ['CTGCAGAA', '>90', '20', '3.1', '37'],
                                      'ClaI': ['CC', '>90', '20', '2.1 and 3.1 have equal activity', '37'], 'EcoRI': ['G', '>90', '2', '3.1', '37'], 'HaeIII': ['GG', '>90', '2', '2.1', '37'], 'HindIII': ['CCC', '75', '20', '2.1', '37'], 
                                      'KpnI': ['GG', '>90', '2', '1.1', '37'], 'MluI': ['CG', '50', '20', '3.1', '37'], 'NcoI': ['CATG', '75', '20', '1.1, 2.1, 3.1 have equal activity', '37'], 'NdeI': ['GGAATTC', '>90', '20', '1.1 and 2.1 have equal acgtivity', '37'],
                                      'NheI': ['CTA', '50', '20', '1.1 and 2.1 have equal activity', '37'], 'NotI': ['AAGGAAAAAA', '>90', '20', '3.1', '37'], 'NsiI': ['CCA', '>90', '2', '3.1', '37'], 'PacI': ['CC', '>90', '20', '1.1', '37'],
                                      'PmeI': ['AGCTTT', '>90', '20', '1.1', '37'], 'PstI': ['AA', '>90', '2', '3.1', '37'], 'PvuI': ['AT', '25', '20', '3.1', '37'], 'SacI': ['C', '10', '2', '1.1', '37'],
                                      'SacII': ['TCC', '>90', '20', '2.1', '37'], 'SalI': ['ACGC', '75', '20', '3.1', '37'], 'ScaI': ['AAA', '75', '2', '1.1 and 2.1 have equal activity', '37'], 'SmaI': ['TCC', '>90', '2', '1.1, 2.1, 3.1 have equal activity', '25'],
                                      'SpeI': ['G', '>90', '20', '2.1', '37'], 'SphI': ['ACAT', '50', '20', '1.1 and 2.1 have equal activity', '37'], 'StuI': ['A', '>90', '2', '2.1', '37'], 'XbaI': ['GC', '>90', '2', '2.1', '37'],
                                      'XhoI': ['CCG', '75', '20', '2.1 and 3.1 have equal activity', '37'], 'XmaI': ['CCC', '>90', '20', '2.1', '37']}
       
        
        if enzyme1 in self.restrictionEnzymeDict: #check the restriction enzymes to the dictionary
            pass
        else:
            print("Error: restriction enzyme not found within program parameters. Please select another restriction enzyme.  A list is available within the readme file of this program.")
            sys.exit()
            
        if enzyme2 in self.restrictionEnzymeDict:
            pass
        else:
            print("Error: restriction enzyme not found within program parameters. Please select another restriction enzyme.  A list is available within the readme file of this program.")
            sys.exit()
            
            
            
        self.enzyme1 = None
        self.enzyme2 = None
        
        self.findEnzymeOrder(enzyme1, enzyme2)
        
        if self.checkTarget() is True: #check to make sure that no cut sites exist within the target
            print("Error: target sequence contains a restriction cut site(s).  Please ensure your target sequence does not contain a cut site within it.")
            sys.exit()
        else:
            pass
        
        
        self.enz1Length = len(self.enzyme1)
        self.enz2Length = len(self.enzyme2)
        
        self.targetStartIndex = None
        self.targetStopIndex = None

        self.revStopCodonList = ['TTA', 'TCA', 'CTA']


        #Other variables
        self.minimumLen = 48 #Required minimum for primer design.
        self.startCodon = self.assignStartCodon(startCodon) #We should make this a user selected input from the main body.
        self.stopCodon = self.assignStopCodons(stopCodon)
        self.fwdFrameCorrection = ''
        self.revFrameCorrection = ''
        self.forwardPrimer = None
        self.reversePrimer = None
        self.finalFwdPrimer = None
        self.finalRevPrimer = None
        self.tempOfFwd = None
        self.tempOfRev = None




        ###########################################################################################
       



        ###########################################################################################

        # Method calls
        self.modVector = self.spliceVector(self.vector, self.enzyme1, self.enzyme2, self.target)
        self.vectorCodons = self.convertAminoAcid()
        self.buildPrimers()
        






    def findEnzymeOrder(self, enzyme1, enzyme2) :
        '''Find order of enzymes within vector sequence.'''
        # Obtain restriction enzyme recognition sequence
        enz1 = r.Restriction_Dictionary.rest_dict[enzyme1]["site"]
        enz2 = r.Restriction_Dictionary.rest_dict[enzyme2]["site"]
        
        # Get index value of sites within vector sequence
        enz1Index = self.vector.find(enz1)
        enz2Index = self.vector.find(enz2)
        
        # Compare index values. The smaller index value is the first enzyme encountered, the biggest index value is the second one encountered.
        if enz1Index > enz2Index :
            self.enzyme1 = enzyme2
            self.enzyme2 = enzyme1
        else :
            self.enzyme1 = enzyme1
            self.enzyme2 = enzyme2
            
        return(self.enzyme1, self.enzyme2)


    #Functions here will be checking functions only
    def checkTarget(self):
        '''Check the target sequence for any restriction enzyme cut sites that may exist within it.  For usage in the main body loop'''

        enz1 = (r.Restriction_Dictionary.rest_dict[self.enzyme1]["site"]) #Get the cut sites using biopython's built in system for restriction enzymes
        enz2 = (r.Restriction_Dictionary.rest_dict[self.enzyme2]["site"]) #Get the cut sites using biopython's built in system for restriction enzymes

        check1 = self.target.find(enz1) #Check for
        check2 = self.target.find(enz2)

        if check1 is -1 or check2 is -1:
            return False
        else:
            return True

        
    def checkTargetLength(self):
        ''' Check the target's length to see if it meets the minimum requirement for a target sequence.  For use in the main body loop'''

        #Return true if the length of the target sequence is less than the required minimum
        if len(self.target) < self.minimumLen:
            return True
        else:
            return False




    #Functions after this point are exclusively action functions

    def assignStartCodon(self,start):#These may be removed with a different version in which the user picks from a drop down menu

        if start is "":
            startCodon = ["ATG"]
        else:
            startCodon = start

        return startCodon


    def assignStopCodons(self,stop): #These may be removed with a different version in which the user picks from a drop down menu

        if stop is "":
            stopCodon = ["TAA"]
        else:
            stopCodon = [stop]
        return stopCodon



    def spliceVector(self, vector, enzyme1, enzyme2, target): #Should we replace the enzyme1 and enzyme2 and just use the self versions?

        #Get the restriction sites for both of the enzymes
        enz1 = (r.Restriction_Dictionary.rest_dict[self.enzyme1]["site"])
        enz2 = (r.Restriction_Dictionary.rest_dict[self.enzyme2]["site"])


        enzymeOneIndex = self.vector.find(enz1) #Find the starting index of the enzyme in relation to the vector.
        enzymeTwoIndex = self.vector.find(enz2) #Find the starting index of the enzyme in relation to the vector

        #Splice the vector into two parts

        #Splicing the first half, from vector index 0 to startIndex+len(enz1)
        modVectorBegin = vector[0:enzymeOneIndex+self.enz1Length]
        #Splicing the second half, from the stop site index to the end of the vector itself.
        modVectorEnd = vector[enzymeTwoIndex:]
        #Save the start and stop indexes for future usage
        self.targetStartIndex = len(modVectorBegin)
        self.targetStopIndex = len(modVectorBegin) + len(self.target)-1

        #Concatenate the vector with the target
        finalModVector = modVectorBegin+target+modVectorEnd
        return finalModVector


    def buildReverseComp(self,seq):

        #Dictionary used for the reverse complement of a given sequence.
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join([complement[base] for base in seq[::-1]]) #Return a new String object that is the reverse Complement
    

        

    def buildPrimers(self):
        '''Build the primers by first setting initial values of the nucleotide number and the optimal temperature.  Go through a while loop
        with the conditional being the optimal temperature.  If the calculated melting point temperature of the primer is too high, decrease the number of nucleotides.
        If it is too low, increase the number of nucleotides.  Do this for both of the primers.  Once the melting temperature is within range of the optimal, exit the loop.
        Then check the frames of the primers in relation to the vector.  Depending on what is found, add nucleotides to the primer to bring it into the correct frame.
        Send primers back to the class.'''
        
        recommendedNuc = 24 #Number of nucleotides used from target sequence. This is expected to change once called inside the two for loops that will be created.
        optimalTempFound = False #boolean control variable for the while loops within the method.
        
        
        #Set up the forward primer
        self.forwardPrimer = self.target[0:recommendedNuc+1]

        
        #Set up the reverse primer
        self.reversePrimer = self.reverseCompTarget[0:recommendedNuc+1]
        if self.reverseCompTarget[0:3] not in self.revStopCodonList:
            self.reversePrimer = 'TTA' + self.reverseCompTarget[0:recommendedNuc+1]
        else :
            pass
        
        #Temperature Variables will be written here for future usage within the loops
        #Will involve the usage of biopython temperature nearest neighbor method class
        self.tempOfFwd = round(mt.Tm_NN(self.forwardPrimer),4)
        self.tempOfRev = round(mt.Tm_NN(self.reversePrimer), 4)
        
        while optimalTempFound is False :
            if self.tempOfFwd >= 60.0 :
                recommendedNuc -= 3
                self.forwardPrimer = self.target[0:recommendedNuc+1]
                self.tempOfFwd = round(mt.Tm_NN(self.forwardPrimer), 4)
            elif self.tempOfFwd <= 54.00 :
                recommendedNuc += 3
                self.forwardPrimer = self.target[0:recommendedNuc+1]
                self.tempOfFwd = round(mt.Tm_NN(self.forwardPrimer), 4)
            else :
                optimalTempFound = True
                
        optimalTempFound = False
        recommendedNuc = 24
        
        while optimalTempFound is False :
            if self.tempOfRev >= 60.00 :
                recommendedNuc -= 3
                self.reversePrimer = self.reverseCompTarget[0:recommendedNuc+1]
                self.tempOfRev = round(mt.Tm_NN(self.reversePrimer), 4)
            elif self.tempOfRev <= 54.00 :
                recommendedNuc += 3
                self.reversePrimer = self.reverseCompTarget[0:recommendedNuc+1]
                self.tempOfRev = round(mt.Tm_NN(self.reversePrimer), 4)
            else :
                optimalTempFound = True
                
        # Concatenate primers from target to restriction enzyme recognition sites, spacer nucleotides, and effeciency nucleotides.
        # Will return as: effeciency_nucleotides enzyme_recognition_site spacer_nucleotides primer_from_target
        self.checkTargetFrame()
        
        self.finalFwdPrimer = self.restrictionEnzymeDict[self.enzyme1][0] + str(r.Restriction_Dictionary.rest_dict[self.enzyme1]["site"]) + self.fwdFrameCorrection + self.forwardPrimer
        self.finalRevPrimer = self.restrictionEnzymeDict[self.enzyme2][0] + str(r.Restriction_Dictionary.rest_dict[self.enzyme2]["site"]) + self.revFrameCorrection + self.reversePrimer
        
        return(self.finalFwdPrimer, self.finalRevPrimer, self.tempOfFwd, self.tempOfRev)
    
        
        
        
        
###### Might delete; might change to convert final target and "vector" sequence into aa sequence    
    def convertAminoAcid(self):
        
        '''Convert the vector sequence into amino acid and find where the tev exists within it'''
        vectorCodonSeq = list()
        
        for nuc in range(0, len(self.vector),3):
            codon = self.vector[nuc:nuc+3]
            if codon in sa.NucParams.dnaCodonTable:
                vectorCodonSeq.append(sa.NucParams.dnaCodonTable[codon])
            else:
                pass
            
        return vectorCodonSeq
    

    def checkTargetFrame(self) :
        '''Find frame of target to make sure it is still in the same frame after cut site.'''
        
        # Check if index of target in vector is in frame (a multiple of 3). If not, add appropriate amount of nucleotides to keep it in frame
        if (self.targetStartIndex-1)%3 == 1:
            # Need to add 2 nuc to end of restriction site in primer
            self.fwdFrameCorrection = 'GA'
            
        if (self.targetStartIndex-1)%3 == 2: 
            # Need to add 1 nuc to end of restriction site in primer
            self.fwdFrameCorrection = 'G'
        else: 
            # Need to add NO nuc to end of restriction site in primer
            self.fwdFrameCorrection = ''
        
        #Check if length of second restriction site is a multiple of 3
        if self.enz2Length%3 == 1:
            self.revFrameCorrection = 'G'
        
        if self.enz2Length%3 == 2:
            self.revFrameCorrection = 'GA'
        
        else:
            self.revFrameCorrection = ''
            
        return (self.fwdFrameCorrection, self.revFrameCorrection)












