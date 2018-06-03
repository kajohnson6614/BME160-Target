
import Bio
from Bio import Restriction as r
from Bio.SeqUtils import MeltingTemp as mt

#Use this version Britney.



class primerDesign():



    def __init__(self, head, vector='', target='', enzyme1='', enzyme2='', startCodon = '', stopCodon= ''):

        #We need a check to make sure everything has been submitted into the class properly.
        #We also need a check to ensure that every target is of sufficient length in order to be an effective target (Completed)
        #We are going to have some try except checks in the future that's for damn sure




        ###########################################################################################
        #Save the vector, target, and enzymes
        self.header = head
        self.vector = vector
        self.target = target
        self.reverseCompTarget = self.buildReverseComp(self.target)
        self.enzyme1 = enzyme1
        self.enzyme2 = enzyme2
        
        self.targetStartIndex = None
        self.targetStopIndex = None

        self.revStopCodonList = ['TTA', 'TCA', 'CTA']


        #Other variables
        self.minimumLen = 48 #Required minimum for primer design.
        self.startCodon = self.assignStartCodon(startCodon) #We should make this a user selected input from the main body.
        self.stopCodon = self.assignStopCodons(stopCodon)




        ###########################################################################################
        #Dictionaries will appear here

        self.restrictionEnzymeDict = {'AccI': ['CCG', '>90', 20], 'AflIII': ['CC', '>90', 2], 'AscI': ['A', '>90', 2], 'AvaI': ['CC', '>90', 2],
                                      'BamHI': ['CG', '>90',2], 'BglII': ['GA', '>90', 20], 'BssHII': ['TTG', '>90', 20], 'BstXI': ['CTGCAGAA', '>90', 20],
                                      'ClaI': ['CC', '>90', 20], 'EcoRI': ['G', '>90', 2], 'HaeIII': ['GG', '>90', 2], 'HindIII': ['CCC', 75, 20], 
                                      'KpnI': ['GG', '>90', 2], 'MluI': ['CG', 50, 20], 'NcoI': ['CATG', 75, 20], 'NdeI': ['GGAATTC', '>90', 20],
                                      'NheI': ['CTA', 50, 20], 'NotI': ['AAGGAAAAAA', '>90', 20], 'NsiI': ['CCA', '>90', 2], 'PacI': ['CC', '>90', 20],
                                      'PmeI': ['AGCTTT', '>90', 20], 'PstI': ['AA', '>90', 2], 'PvuI': ['AT', 25, 20], 'SacI': ['C', 10, 2],
                                      'SacII': ['TCC', '>90', 20], 'SalI': ['ACGC', 75, 20], 'ScaI': ['AAA', 75, 2], 'SmaI': ['TCC', '>90', 2],
                                      'SpeI': ['G', '>90', 20], 'SphI': ['ACAT', 50, 20], 'StuI': ['A', '>90', 2], 'XbaI': ['GC', '>90', 2],
                                      'XhoI': ['CCG', 75, 20], 'XmaI': ['CCC', '>90', 20]}  # Dictionary with nucleotides to add to primer with percent effeciencies. [nucleotides for beginning, nucleotides for end, percent1, hour1, percent2, hour2]




        ###########################################################################################

        #Method calls will be here

        #self.lengthCheck = self.checkTargetLength()
        #self.targetSiteCheck = self.checkTarget()
        self.modVector = self.spliceVector(self.vector, self.enzyme1, self.enzyme2, self.target)









    #Functions here will be checking functions only
    def checkTarget(self):
        '''Check the target sequence for any restriction enzyme cut sites that may exist within it.  For usage in the main body loop'''

        enz1 = (r.Restriction_Dictionary.rest_dict[self.enzyme1]["site"]) #Get the cut sites using biopython's built in system for restriction enzymes
        enz2 = (r.Restriction_Dictionary.rest_dict[self.enzyme2]["site"]) #Get the cut sites using biopython's built in system for restriction enzymes

        check1 = self.target.find(enz1) #Check for
        check2 = self.target.find(enz2)

        if check1 or check2 is -1:
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

    def checkForStart(self): #May get rid of this.
        ''' Check for the first start codon found within the target.  Mark its location'''
        position = self.target[0:len(self.target)+1].find(self.startCodon)
        return position

    def checkForStop(self):
        '''Check for the first stop codon found within the target.  Mark its location'''
        position = self.reverseCompTarget[0:len(self.reverseCompTarget)].find(self.stopCodon)
        return position











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
        print(enzymeOneIndex)
        print(enzymeTwoIndex)


        #Check to determine which of the enzymes is going to be the starting enzyme.
        if enzymeOneIndex > enzymeTwoIndex:
            startIndex = enzymeTwoIndex
            stopIndex = enzymeOneIndex
            beginLength = len(enz2)
            endLength = len(enz1)
        else:
            startIndex = enzymeOneIndex
            stopIndex = enzymeTwoIndex
            beginLength = len(enz1)
            endLength = len(enz2)


        #Splice the vector into two parts

        #Splicing the first half, from vector index 0 to startIndex+len(enz1)
        modVectorBegin = vector[0:startIndex+beginLength]
        #Splicing the second half, from the stop site index to the end of the vector itself.
        modVectorEnd = vector[stopIndex:]
        #Save the start and stop indexes for future usage
        self.targetStartIndex = modVectorBegin 
        self.targetStopIndex = modVectorEnd

        #Concatenate the vector with the target
        finalModVector = modVectorBegin+target+modVectorEnd
        return finalModVector


    def buildReverseComp(self,seq):

        #Dictionary used for the reverse complement of a given sequence.
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join([complement[base] for base in seq[::-1]]) #Return a new String object that is the reverse Complement
    

        

    def buildPrimers(self):
        
        recommendedNuc = 24 #Number of nucleotides used from target sequence. This is expected to change once called inside the two for loops that will be created.
        optimalTempFound = False #boolean control variable for the while loops within the method.
        
        

        #Set up the forward primer
        forwardPrimer = self.target[0:recommendedNuc+1]
        if forwardPrimer.find('ATG') is not 0 :
            forwardPrimer = 'ATG' + self.target[0:recommendedNuc+1] # Change 'ATG' to self.startCodon
            return(forwardPrimer)
        else :
            pass
        
        #Set up the reverse primer
        reversePrimer = self.reverseCompTarget[0:recommendedNuc+1]
        if self.reverseCompTarget[0:3] not in self.revStopCodonList:
            reversePrimer = 'TTA' + self.reverseCompTarget[0:recommendedNuc+1]
            return(reversePrimer)
        else :
            pass
        
        #Temperature Variables will be written here for future usage within the loops
        #Will involve the usage of biopython classes and functions
        tempOfFwd = round(mt.Tm_NN(forwardPrimer),4)
        tempOfRev = round(mt.Tm_NN(reversePrimer), 4)
        
        while optimalTempFound is False :
            if tempOfFwd >= 60.0 :
                recommendedNuc -= 3
                forwardPrimer = self.target[0:recommendedNuc+1]
                tempOfFwd = round(mt.Tm_NN(forwardPrimer), 4)
            elif tempOfFwd <= 54.00 :
                recommendedNuc += 3
                forwardPrimer = self.target[0:recommendedNuc+1]
                tempOfFwd = round(mt.Tm_NN(forwardPrimer), 4)
            else :
                optimalTempFound = True
                
        optimalTempFound = False
        recommendedNuc = 24
        
        while optimalTempFound is False :
            if tempOfRev >= 60.00 :
                recommendedNuc -= 3
                reversePrimer = self.reverseCompTarget[0:recommendedNuc+1]
                tempOfRev = round(mt.Tm_NN(reversePrimer), 4)
            elif tempOfRev <= 54.00 :
                recommendedNuc += 3
                reversePrimer = self.reverseCompTarget[0:recommendedNuc+1]
                tempOfRev = round(mt.Tm_NN(reversePrimer), 4)
            else :
                optimalTempFound = True
        return([forwardPrimer, reversePrimer])

        












