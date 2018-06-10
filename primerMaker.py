
#Kyle Johnson, Britney Hernandez (kyajohns, )
#primerMaker.py



import CommandLine
import sequenceAnalysis
import primerDesign
import Bio
import sys



def main():


    '''Main function for the primer design program.  Imports primerDesign, Bio, CommandLine, and sequenceAnalysis.
    Main function takes in inputs and boolean checks from the command line.  These include the target sequence, restriction
    enzymes one and two, any changes to start and stop codons, and a verbosity check which will enable printing either to
    standard output or a output file.  Output would appear as the following after a successful run of main:


    ############################################################
    
    FastA Header
    
    Forward Primer
    Primer Sequence
    {} nucleotides were added to give {} efficiency after {} hours.
    Buffer {} for digestion at {} degrees.
    
    
    Melting Temperature Forward
    GC Content Percentage Forward
    
    Reverse Primer
    Primer Sequence
    {} nucleotides were added to give {} efficiency after {} hours.
    Buffer {} for digestion at {} degrees.
    
    Melting Temperature Reverse
    GC Content Percentage Reverse
    
    
    ############################################################


    Should errors occur such as improper target sequence or (Fill in the check conditions we can think of here),
    a message will be displayed indicating the potential problem to the user:

    Error:
    This program has detected that (Situation).  Please correct your (Situation) and try again.
    
    After the message is displayed, the program will exit and return back to the terminal line.

    '''
###################################################################################################
    
    #Main method variables 
    cl = CommandLine.Command_Line()
    
    gcForward = None
    gcReverse = None

    #Gather the restriction enzymes
    restrictionEnzyme1 = cl.args.enzymeOne
    restrictionEnzyme2 = cl.args.enzymeTwo

    #Gather the start and stop codon
    startCodon = cl.args.start
    stopCodon = cl.args.stop

    #Verbosity Boolean
    verb = cl.args.verbosity
    
    #File Name
    targetFile = cl.args.target
    
    #Marker Number
    markerNumber = 60
    #Degree Symbol
    degree = "\u00b0"
    
    
    ###################################################################################################
    #Check to see if all of the required elements are in place.  If any value is at none, terminate the program with a message
    if targetFile is None:
        print("No Target Sequence inputted.  Please retry with a proper input file in the appropriate location.  See -h for command line help.")
        sys.exit()
    else:
        pass
    
    if restrictionEnzyme1 is None:
        print("No Enzyme One inputted.  Please retry with an Enzyme One in the appropriate location.  See -h for command line help.")
        sys.exit()
    else:
        pass
    
    if restrictionEnzyme2 is None:
        print("No Enzyme Two inputted.  Please retry with an Enzyme Two in the appropriate location.  See -h for command line help.")
        sys.exit()
    else:
        pass
    
    
        



    
    
    
    
        


###################################################################################################
        
    fastA = sequenceAnalysis.FastAreader(targetFile)#Read from the file collected by CommandLine class
    for head, seq in fastA.readFasta(): #Ideally there should be only one fastA to read given a run.

        #If we decide otherwise, place all of the class method calls within the loop
        createdPrimer = primerDesign.primerDesign(head,seq, restrictionEnzyme1, restrictionEnzyme2, startCodon, stopCodon)
        
        
    #createdPrimer.buildPrimers() #Build the primers using the built in buildPrimers method.  Results are stored in the class object
    
    nucForward = sequenceAnalysis.NucParams(str(createdPrimer.forwardPrimer))
    nucReverse = sequenceAnalysis.NucParams(str(createdPrimer.reversePrimer))
    
    gcForward = (nucForward.nucComposit["G"]+nucForward.nucComposit["C"])/nucForward.nucCount()
    gcReverse = (nucReverse.nucComposit["G"]+nucReverse.nucComposit["C"])/nucReverse.nucCount()
    
###################################################################################################    
    #Printing Section
    #Print either to an output file or std out depending upon verbosity condition
    
    if verb is True: #Verbosity mode output.  If enabled, writes to a file instead of std out.
        with open("PrimerOut.txt", "w") as p:
            p.write("#"*markerNumber+"\n\n")
            p.write(createdPrimer.header+"\n\n")
            p.write("Forward Primer\n")
            p.write(createdPrimer.finalFwdPrimer+"\n")
            p.write("\n")
            p.write("'{0}' nucleotides were added to give {1} efficiency after {2} hours.\n".format(createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][0],createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][1], createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][2]))
            p.write("Buffer {} for digestion at {} Degrees.\n\n".format(createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][3], createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][4]))
            p.write("Melting Temperature = "+str(createdPrimer.tempOfFwd)+degree+"C"+"\n")
            p.write("\n")
            p.write(str.format("{0:.4f}", gcForward)+" % GC Content\n")
            p.write("\n")
            p.write("Reverse Primer\n")
            p.write(createdPrimer.finalRevPrimer+"\n")
            p.write("\n")
            p.write("'{0}' nucleotides were added to give {1} efficiency after {2} hours.\n".format(createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][0],createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][1], createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][2]))
            p.write("Buffer {} for digestion at {} Degrees.\n\n".format(createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][3], createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][4]))
            p.write("Melting Temperature = "+str(createdPrimer.tempOfRev)+degree+"C"+"\n")
            p.write("\n")
            p.write(str.format("{0:.4f}",gcReverse)+" % GC Content\n")
            p.write("#"*markerNumber)
            
    else: #Print to std out instead of to a file
        print("#"*markerNumber)
        print(createdPrimer.header)
        print()
        print("Forward Primer")
        print(createdPrimer.finalFwdPrimer)
        print()
        print("'{0}' nucleotides were added to give {1} efficiency after {2} hours.\n".format(createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][0],createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][1], createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][2]))
        print("Buffer {} for digestion at {} Degrees.\n".format(createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][3], createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme1][4]))
        print("Melting Temperature = "+str(createdPrimer.tempOfFwd)+"C")
        print()
        print(str.format("{0:.4f}", gcForward)+" % GC Content\n")
        print()
        print("Reverse Primer")
        print(createdPrimer.finalRevPrimer)
        print()
        print("'{0}' nucleotides were added to give {1} efficiency after {2} hours.\n".format(createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][0],createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][1], createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][2]))
        print("Buffer {} for digestion at {} Degrees.\n".format(createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][3], createdPrimer.restrictionEnzymeDict[createdPrimer.enzyme2][4]))
        print()
        print("Melting Temperature = "+ str(createdPrimer.tempOfRev)+"C")
        print()
        print(str.format("{0:.4f}",gcReverse)+" % GC Content\n")
        print()
        print("#"*markerNumber)



###################################################################################################

main()
