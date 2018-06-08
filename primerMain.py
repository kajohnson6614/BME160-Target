#Env
#Kyle Johnson, Britney Hernandez (kyajohns, )
#primerMain.py



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


    Put output format in here



    Should errors occur such as improper target sequence or (Fill in the check conditions we can think of here),
    a message will be displayed indicating the potential problem to the user:

    Error:
    This program has detected that (Place example situtation here).  Other warning messages can be placed around here detailing
    what the user needs to do in order to avoid further issues.

    '''


    #We are going to need a large series of checks in order to ensure no errors pop up.  Keep tabs about this point.


    cl = CommandLine.Command_Line()

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


    #Check to see if all of the required elements are in place.  If any value is at none, terminate the program with a message
    #if cl.args.target or restrictionEnzyme1 or restrictionEnzyme2 is None:
     #   print("Error:Please ensure that all required variables have been inputted into the command line: \n","These required inputs are: Target Sequence, Restriction Enzyme One, Restriction Enzyme Two\n")
      #  sys.exit
    #else:
    #    pass
    
    
    
    
        



    fastA = sequenceAnalysis.FastAreader(targetFile)#Read from the file collected by CommandLine class
    for head, seq in fastA.readFasta(): #Ideally there should be only one fastA to read given a run.

        #If we decide otherwise, place all of the class method calls within the loop
        createdPrimer = primerDesign.primerDesign(head,seq, restrictionEnzyme1, restrictionEnzyme2, startCodon, stopCodon)
        
        
    createdPrimer.buildPrimers() #Build the primers using the built in buildPrimers method.  Results are stored in the class object
    
    #Checks on the primers will be conducted here.  All output should be to standard output instead of to a file
    print(verb)
    
    if verb is True: #Verbosity mode output.  If enabled, writes to a file instead of std out.
        with open("PrimerOut.txt", "w") as p:
            p.write(createdPrimer.target)
            p.write("\n")
            p.write(createdPrimer.finalFwdPrimer+"\n")
            p.write(str(createdPrimer.tempOfFwd))
            p.write("\n")
            p.write(createdPrimer.finalRevPrimer+"\n")
            p.write(str(createdPrimer.tempOfRev))
            
    else: #Print to std out instead of to a file
        print(createdPrimer.target)
        print()
        print(createdPrimer.finalFwdPrimer)
        print(createdPrimer.finalRevPrimer)
        print(str(createdPrimer.tempOfFwd))
        print(str(createdPrimer.tempOfRev))


    #Here is where we shall put all of the checks and balances for the main function.

    #Check if the restriction enzymes are within the dictionary:  Fail if they are not

    #Check if the target sequence is too short to be of effective usage:  Fail if the minimum is not met

    #Check if the



    #Final printing stage
    #Check the value of the verbosity input
    #If the value is True, Send the output to a specific txt file that can be accessible to the user
    # with open("primerOutput.txt", "w") as p:
    #   Do the work required to give out a decent amount of information.  Especially with the primers


    #Else if the value of verbosity is false (by default), print the information out to standard output.
    #This will involve print statements in order to accomplish this goal.



main()