#!/usr/bin/python

# for commandline options:
from optparse import OptionParser, OptionGroup

from matrices import *

class Sequence:
    """Stores a sequence object"""
    
    def __init__(self, Label="", Sequence="" ):
        """Initialize a new Sequence object

        Label -- identifier of sequence (text)
        Sequence -- sequence string in single-letter alphabet
        """
        self.Label       = Label
        self.Sequence    = Sequence

    # this makes that you can do 'print sequence' and get nice output:
    def __str__(self):
        """Return string representation of a Sequence object"""
        # newline-delimited values of all the attributes
        return ">%s\n%s" % (self.Label, self.Sequence)


def readSequences(lines):
    """Return Sequences object
    
    lines -- list of lines or any object that behaves like it
    
    This routine parses a fasta file and returns a list of Sequence objects
    containing the sequences with label and sequence data set
    """
    seqs = []
    label = None
    seq_lines = []
    for line in lines:
        line = line.strip()      # strip off white space
        if not line:             # skip empty lines
            continue
        if line.startswith(';'): # ignore comment lines
            continue
        # check for start of next sequence:
        if line.startswith('>'): # label line
            # first, store the previous sequence if we had one:
            if seq_lines:
                seqs.append(Sequence(label, ''.join(seq_lines)))
                seq_lines = []
            # get the label (name) for the next sequence
            label = line[1:].strip()
        else:
            # collect all lines with sequence information for this sequence:
            seq_lines.append(line)
    # take care of the last sequence in the file
    seqs.append(Sequence(label, ''.join(seq_lines)))
    return seqs


def check_args(parser, options, args):
    
    if not options.fasta:
        # check if we have an option left (to be used as input filename):
        if args:
            options.fasta = args.pop()
        else:
            print "Need at least an input file (fasta)\n"
            parser.print_help()
            print "\nERROR: no input file given"
            exit(-1)

    # check alignment type:
    align_options = [options.align_local, options.align_global, options.align_semiglobal]
    # check if at least one alignment option was true, else choose global
    if align_options.count(True)==0:
        print "No alignment type given, using Global"
        options.align_global=True
    # check if not more than one alignment option was true, else error and exit 
    if align_options.count(True)>1:
        print "ERROR: multiple alignment types chosen"
        exit(-1)

    # check for any leftover command line arguments:
    if len(args):
        warning("ignoring additional arguments "+str(args))
    
    return options
    
    
def parse_commandline():
    usage = "%prog <fasta> [options]"
    version = "0.0a"
    description = \
        "%prog aligns two sequences."
    epilog = \
        "Copyright (c) 2011 K. Anton Feenstra -- "\
        "feenstra@few.vu.nl -- www.few.vu.nl/~feenstra"
    parser = OptionParser(usage=usage, description=description,
                          version="%prog "+version, epilog=epilog)

    # sequence/alignment options:
    parser.add_option("-f", "--fasta",  dest="fasta", metavar="<file>",
                     help="input alignment file (fasta)")
    parser.set_defaults(fasta=None)
    
    parser.add_option("-e", "",  dest="exchange_matrix", 
                     help="Exchange matrix: pam250, blosum62 or identity (%default%)")
    parser.set_defaults(exchange_matrix="pam250")
    
    parser.add_option("-l", "",  dest="align_local",  action="store_true",
                     help="align local")
    parser.set_defaults(align_local=False)
    
    parser.add_option("-g", "",  dest="align_global", action="store_true",
                     help="align global")
    parser.set_defaults(align_global=False)
   
    parser.add_option("-v", "", dest="print_Matrix",action="store_true",
                     help="print Matrix")
    parser.set_defaults(print_Matrix=False)
    
    parser.add_option("-s", "",  dest="align_semiglobal", action="store_true",
                     help="align semi global")
    parser.set_defaults(align_semiglobal=False)
    
    parser.add_option("-p", "",  dest="gap_penalty", type="int",
                     help="Gap penalty (%default%)")
    parser.set_defaults(gap_penalty=2)
    
    # get the options:
    (options, args) = parser.parse_args()

    #check if everything is as expected
    options = check_args(parser, options, args)

    # clean up (recommended):
    del(parser)
    
    return options, args

# main function:
def main():
    # get command line options
    options, args = parse_commandline()
    
    #print the value of print_Matrix to check wether the value is true or false
    print "Options.print_Matrix ", options.print_Matrix
    # set substitution matrix:
    if options.exchange_matrix == "pam250":
        exchangeMatrix = pam250
    elif options.exchange_matrix == "blosum62":
        exchangeMatrix = blosum62
    elif options.exchange_matrix == "identity":
        exchangeMatrix = identity
    else:
        print "unknown exchange matrix", options.exchange_matrix
        exit(-1)
    #if the print_Matrix is True then keep the value true
    if options.print_Matrix == True:
        print_Matrix = True
    #if the value of the print_Matrix is not true than put the value False to it.
    else:
        print_Matrix = False
    # read sequences from fasta file, and catch error reading file
    try:
        sequences = readSequences(open(options.fasta))
    except IOError:
        print "ERROR: cannot open or read fasta input file:", fastafile
        exit(-1)

    for seq in sequences:
        print seq
    
    # call alignment routine(s):
    print "print_Matrix=", print_Matrix

    if options.align_global:
        # put the result of the function do_global_alginment to the variable score_matrix
        score_matrix = do_global_alignment(sequences, exchangeMatrix, options.gap_penalty,options.print_Matrix)
        # call the function give_gobal_traceback to show the traceback of the local alignment 
        give_global_traceback(sequences,exchangeMatrix,score_matrix,options.gap_penalty)      
    elif options.align_semiglobal:
        score_matrix = do_semiglobal_alignment(sequences,exchangeMatrix, options.gap_penalty,options.align_semiglobal)
        give_semiglobal_traceback(sequences,exchangeMatrix,score_matrix,options.gap_penalty)
    elif options.align_local:
        score_matrix = do_local_alignment(sequences,exchangeMatrix, options.gap_penalty,options.align_local)
        give_local_traceback(sequences,exchangeMatrix,score_matrix,options.gap_penalty)

    else:
        print "BUG! this should not happen."
        exit(-1)
if __name__ == "__main__":
    main()

# last line
