# CommandLine.py
# Kyle Johnson, Britney Hernandez
#

###################################################################################################


class Command_Line():



    def __init__(self):
        '''Takes in inputs from the command line and allows usage in other objects once called upon.
            Class takes in the inputs: target, enzymeOne, enzymeTwo, start, stop, verbosity
            start, stop, and verbosity  have default values associated and are therefore not necessarily required for the class
            All other inputs are required for program runtime.'''
        import argparse

        self.parser = argparse.ArgumentParser(description= "Primer Design Argument Parser.  This program takes in a target fastA file, two restriction enzymes, and optional alternate start and stop codons",
                                              epilog="What we would want this thing to say at the end of our help section.  For now bananas.", add_help=True, prefix_chars= "-", usage='%(prog)s [options]')


        self.parser.add_argument("--target", "-t", help="Takes in a target fastA file as input for file",action= 'store')

        self.parser.add_argument("--enzymeOne", "-e1", help="Takes in the first desired restriction enzyme", type=str, action='store')
        self.parser.add_argument("--enzymeTwo", "-e2", help="Takes in the last desired restriction enzyme", type=str, action='store')

        self.parser.add_argument("--start","-st", help="Takes in a desired start codon.  Default is set to ATG", type=str, const="ATG", nargs='?',action='store') #These would be incredibly easy to convert to a list should the desire arise.
        self.parser.add_argument("--stop","-sp", help="Takes in a desired stop codon.  Default is set to TAA", type=str, const="TAA", nargs='?', action='store')
        self.parser.add_argument("--verbosity", "-v", help="Turns on verbosity mode, determining if output is to standard output or to a txt file.", action='store_true', default=False)


        self.args = self.parser.parse_args() #Insantiates a list of objects in order to call upon





###################################################################################################