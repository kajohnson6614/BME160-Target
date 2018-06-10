# BME160-Target
primerMaker Version 1.0
Kyle Johnson and Britney Hernandez


######################################################################################################################################

This program is designed to take in a user target sequence fastA file along with two restriction enzymes and generate the appropriate forward and reverse primers eligble for PCR based plasmid cloning.  The program requires the latest version of python which as of this writing is python 3.6.  Python is available either directly or through Anaconda.  Links are provided below for both versions.

Python 3.6
https://www.python.org/downloads/ (select the version that is compatible for your operating system)

Anaconda
https://www.anaconda.com/download/ (Latest version as 6/10/18 is 5.2)
(To find where python is installed with Anaconda, go to the anaconda console once installed and type "where python".  This will give you the directory location of python)


This program will also require the usage of a command line interface.  On windows the prefered shell is git bash as this provides
a unix styled command line for windows machines.  In order to run the program, the shell must have a path set up to the latest version of python.  Instructions on how to set a path up for python 3.6 on windows machines is provided below.

Path Installation Tutorial
https://docs.python.org/3/using/windows.html#excursus-setting-environment-variables

Furthermore, this program will required that biopython is installed to your version of python.  An installation guide for both versions of python (direct or through Anaconda) is provided.

https://biopython.org/wiki/Download (Standard python)
https://biopython.org/wiki/Packages (Anaconda Package)

To run the program, three inputs are required for successful activation of the program.  An example command line is shown below:

python primerMaker.py --target "File Name" --enzymeOne "Enzyme One" --enzymeTwo "Enzyme Two"

or

python primerMaker.py -t "File Name" -e1 "Enzyme One" -e2 "Enzyme Two"

Please note that the restriction enzymes must be spelled with correct spelling.  For instance, EcoRI is an acceptable input but ecori or ecoR1 is not.  This is expected to change in later versions.

The command line also has optional commands.  These commands include -h to display help options, --start and --stop to input custom start and stop codons (this feature has not been implemented as of version 1.0), and -v for a verbosity mode that changes output from either standard output (default) to a written text file.





######################################################################################################################################






