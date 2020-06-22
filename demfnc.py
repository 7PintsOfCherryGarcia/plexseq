###############################################################################
# Julan Regalado - julian.regalado@tuebingen.mpg.de
#
#
# GNU GENERAL PUBLIC LICENSE
# Version 2, June 1991
#
#
# Demultiplex fastq files acording to 1) barcode sequences present on reads.
#
###############################################################################

import sys
import time
import gzip
from optparse import OptionParser
import os

"""
The following funtions provide the main functionality to plexseq.py
Errors are logged to stderr
"""
def error(error_string):
    """
    Log error message to stderr and exit with error code
    ARGS:
        INPUT - error_string(str), string to output
                err_code(int), exit status
        RETURN  - None
    """
    sys.stderr.write(error_string)
    usage()
    sys.exit()


def usage():
    """
    Prints usage to stderr
    ARGS:
        INPUT - None
        RETURN - None
    """
    sys.stderr.write("\nUsage:\npython plexseq.py -r <read_file> -i <index"\
                     "file> -o outdir [--r2 <read2_file> --ri <interleaved_read_file> --nostrict(defoult off)]"\
                     "\nOr:\ncat seqfile | python demultiplex.py -i <indexfil"\
                     "e> -o outdir --s [--nostrict(defoult off)]\n")


def getoptions():
    """
    Parse sys.argv[1:] to get script options
    ARGS:
        INPUT - argline, sys.argv[1:](default)
        RETURN - checked_params, dictionary of correctly used parameters
    """
    # set of options to parse, single dashed parameters are required, doubled
    # dashed parameters are optional

    # set parser
    parser = OptionParser()

    # set options, description in help string
    parser.add_option("-i","--index_file",action="store",type="string",\
                      dest="index_filename",
                      help="barcode file with barcode names and sequences")

    parser.add_option("-r","--read_file",action="store",type="string",\
                      dest="read_filename",
                      help="read file to demultiplex, specify --r2 if paired")

    parser.add_option("-o","--outdir",action="store",type="string",\
                      dest="outdir",help="output directory")

    parser.add_option("--r2",action="store",type="string",\
                      dest="read2_filename",
                      help="read2 filename for paired reads")

    parser.add_option("--ri",action="store",type="string",\
                      dest="int_filename",
                      help="interleaved read filename if r1 and r2 in same fi"\
                           "le, incompatible with -r --r2, and --s")

    parser.add_option("--s",action="store_true",default=False,\
                      dest="stdin",
                      help="demultiplex reads from STDIN, incopatible with -r"\
                           " and --r2")
    parser.add_option("--nostrict",action = "store_true",dest="no_strict",default=False,\
                      help="allows barcode sequences to mismatch" \
                      " read sequence (False)")

    # parse options in sys.argv[1:]
    (options, args) = parser.parse_args()

    # parameters should be set properly and not all parameter combinations
    # are valid.

    # index file requiered
    if not options.index_filename:
        error("ERROR: incorrect command line input\nINDEX_FILENAME missing"\
              "\n")

    # read file or data from STDIN or interleaved file required
    if not options.read_filename:
        if not options.stdin:
            if not options.int_filename:
                error("ERROR: incorrect command line input\nREAD_FILENAME "\
                      " missing, and neither --s nor --ri options set\n")

    # output directory required (for now only in command line, output is still
    # written to current directory)
    if not options.outdir:
        error("ERROR: incorrect command line input\nOUTDIR missing"\
              "\n")


    # if STDIN data stream is set, no -r or --r2 otions can be specified
    if (options.read_filename or options.read2_filename) and options.stdin:
        error("ERROR: incorrect command line input\n-r(--r2) and --s optio"\
              "ns are mutially exclusive\n")

    # same with interleaved fasqfiledemfnc.strictd format
    checked_params = checkparams(options)

    return checked_params


def checkparams(options):
    """
    Checks for valid parameters, exits otherwise
    ARGS:
        INPUT - paramlist, list of parameters
        OUTPUT - returnlist, list of checked parameters
    """
    returndict = {}

    # set strict flag
    if options.no_strict:
         sys.stdout.write("\nWARNING: --nostrict flag set, barcodes will be "\
                         "matched to read sequences with mismatches. Although"\
                         "matches have to be unamiguous, this may misgroup a "\
                         "read depending on the nature of the barcodes provid"\
                         "ed\n")
         returndict['no_strict'] = True
    else:
        returndict['no_strict'] = False

    # open index file and add file handle
    try:
        indexfile = open(options.index_filename,'r')
        returndict['index_file'] = indexfile
    # error if file missing or not readable
    except IOError:
        error("ERROR: filename error\nUnable to open \""+\
              str(options.index_filename)+"\" file not found or no reading pe"\
              "rmissions\n")

    # if --s flag, check for data in STDIN
    if options.stdin:
        if not sys.stdin.isatty():
            fastqfile = sys.stdin
            returndict['stdin'] = sys.stdin
        else:
            error("ERROR: STDIN error\n--s option specified but no data co"\
                  "mming from STDIN\n")

    # if STDIN is false, check for files
    elif options.int_filename:
        if not sys.stdin.isatty():
            error("ERROR: STDIN error\n--s option not set but data is comm"\
                  "ing from STDIN\n")

        # check that file exists or exit
        if os.path.isfile(options.int_filename):
            returndict['int_fqfile'] = options.int_filename
        else:
             error("ERROR[2]: filename error\n\""+\
                  str(options.int_filename)+"\" does not exist\n")
    else:
        # check if file exists or exit
        if os.path.isfile(options.read_filename):
            # --s and -r options incompatible
            if not sys.stdin.isatty():
                error("ERROR: STDIN error\n--s option not set but data is "\
                      "comming from STDIN\n")
            else:
                returndict['fastqfile'] = options.read_filename
            # same
            if options.read2_filename:
                if os.path.isfile(options.read2_filename):
                    returndict['fastq2filename'] = options.read2_filename
                else:
                   error("ERROR: filename error\n\""+\
                  str(options.read2_filename)+"\" does not exist\n")
        else:
             error("ERROR: filename error\n\""+\
                  str(options.read_filename)+"\" does not exist\n")

    # Create output directory if it does not exist
    if not os.path.exists(options.outdir):
        try:
            os.makedirs(options.outdir)
            sys.stdout.write("Writting to: "+options.outdir+"\n")
            returndict['outdir'] = options.outdir
        except OSError:
            error("ERROR: filename error\nCould not create"+options.outdir+\
                  "directory, check for correct writing permissions.")
    else:
       returndict['outdir'] = options.outdir

    return returndict


def file_type(filename):
    """
    Determines if "filename" is compressed, and returns compression format.
    Checks for the first bytes in the file to see if they match the standard
    first bytes of typical compression formats. This function can determine if
    "filename" is gzip, bgzip or zip compressed.

    NOTE*** A file may have the corresponding starting bytes without actually
    being compressed in the determined format.
    ARGS:
        INPUT - filename, file name to determine format
        RETURN - None
    """
    # Dictionary of "madic numbers" == starting bytes : file type
    magic_dict = { "\x1f\x8b\x08": "gz",
                   "\x42\x5a\x68": "bz2",
                   "\x50\x4b\x03\x04": "zip"
                 }
    # Number of bytes to read at beggining of file
    max_len = max(len(x) for x in magic_dict)

    with open(filename) as f:
        # Read from start of file max_len bytes
        file_start = f.read(max_len)
    # for every pair of starting bytes and file type check if there is a match
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            f.close()
            # Return match
            return filetype
    f.close()
    return "None"


def readindex(filehandle,outdir,gzipflag):
    """
    Reads a tab separated barcode file with the following format:
    BARCODE_NAME(id)	SEQUENCE	<SEQUENCE2>
    First and second columns are mandatory. If a second sequence is provided,
    the script will also require the second sequence to match the read.
    ARGS:
        INPUT - filehandle: open file handle pointing to the index file
                gzipflag: flag indicating if output files should write in
                          gzip format.
    """
    checkfirst = False
    maxindexlen = 0
    # barcode ids are stored in a dictionary where barcode sequences are keys
    indexdict = dict()
    # loop over index file lines
    for line in filehandle:
        # ignore lines beginning with '#' (comments)
        if line[0] == '#' or line == '\n':
            pass
        else:
            indexlist = line.rstrip().split('\t')
            if len(indexlist) < 2 or len(indexlist) > 3:
                error("ERROR: format error:\nIndex file seems to be of inc"\
                      "orrect format.\nIndex file must be a tab separated fil"\
                      "e:\nBARCODE_NAME(id)i\tSEQUENCE\t<SEQUENCE2>")

            if not checkfirst:
                # if a second barcode sequence is specified, warn user and
                # continue
                if len(indexlist) == 3:
                    sys.stderr.write("\nWARNING: More than two barcodes have bee"\
                                     "n specified\ndemultiplex.py will look f"\
                                     "or index1 in read1 5' end and index2 in"\
                                     "read2 5' end\n")
                    checkfirst = 3
            # store first and second entry,barcode name and sequence
            name = indexlist[0]
            index = indexlist[1]
            if len(index) > maxindexlen:
                maxindexlen = len(index)
            if checkfirst == 3:
                try:
                    index2 = indexlist[2]
                except IndexError:
                    error("ERROR: format error\nAll barcodes must have 2 bar"\
                          "code sequences. Barcode "+str(name)+" only has 1\n")
            else:
                index2 = 'None'
            # store in dictinary with adapter sequence as key.
            if gzipflag: # not implemented right now
                pass
                #indexdict[index] =[name,gzip.open("BCDE"+name+'.fastq.gz','wb',6)]
            else:
                if index not in indexdict.keys():
                    filename = outdir + '/' + name
                    indexdict[index] =[ [ name,
                                          open(filename+'.fq','w'),
                                          index2,
                                          ("",0)
                                      ] ]
                # in case barcode sequences are repeated and no second sequence
                # is present to distinguish them, warn user and store both
                # barcodes to same file
                elif index2 == 'None':
                    sys.stderr.write("\nWARNING: Barcode "+index+" is repeated "\
                                     "and no second barcode was provided. All"\
                                     " sequences with "+index+" barcode will "\
                                     "be saved to same file.\n")
                elif index2 in indexdict[index][0]:
                    sys.stderr.write("\nWARNING: Barcode pair "+index+'\t'+\
                                     index2+" is repeated. All sequences with"\
                                     " these barcodes will be saved to same"
                                     " file "+indexdict[index][0][0]+"\n")
                else:
                    filename = outdir + '/' + name
                    indexdict[index].append([name,open(filename+'.fq','w'),index2,("",0)])
            """
            NOTE* indexdict stores a list of lists as values to each key. Each
            individual list correspond to a diferent barcode. Two barcodes
            share the same key (len(indexdict[key])>1) when their sequences are
            identical and thus, the second sequence is used to distiguish
            between them.
            Each barcode list has 4 elements in the following order:
                indexdict[key][list_number][0] - barcodeID
                indexdict[key][list_number][0] - filehandle
                indexdict[key][list_number][0] - second barcode sequence
                indexdict[key][list_number][0] - tuple buffer (not implemented)
            """

    if gzipflag:
        pass
        #noexactmatch = gzip.open("NOEXACTMATCH.fq.gz",'wb',6)
    else:
        noexactmatch = open(outdir+'/'+"NOEXACTMATCH.fq",'a')

    return [indexdict,noexactmatch]


def parsefq(filehandle,file2handle=False,pairedflag=False,stacklim=500):
    """
    read fastq file and parse data into a list of lists where each list
    contains the data of one read (header, sequence quality string)
    this function parses "stacklim" reads at a time. This is done in order for
    small memory machines to be able to demultiplex large datasets.
    Paired reads are atomatically detected as well as format errors
    ARGS:
        INUPUT - filehandle, file containing reads to parse
                 file2handle, file with paired reads to parse
                 pairedflag, a priory information of paired reads
                 stacklim, number of reads to parse a a time
        RETURNS - (returnstack,bool), list of parsed reads and flag indicating
                                      if there are still more reads to parse
    """
    # list of read lists to return
    returnstack = []
    # if reads should be paired check thet they are correctly paired
    if pairedflag:
        # get read
        read1 = getread(filehandle)
        # if pared read file exists, get read from there
        if file2handle:
            read2 = getread(file2handle)
        # else get next read in file
        else:
            read2 = getread(filehandle)
        # check that reads are correctly paired
        if ispaired(read1,read2):
            returnstack.append([read1,read2])
        # log errors
        else:
            if file2handle:
                error("ERROR: file error\nread files provided for pared re"\
                      "ad data (-r and --r2), but reads in these files are no"\
                      "t properly paired\n")
            # auto detect paired information if data comes from STDIN
            elif not sys.stdin.isatty():
                #error("\nERROR: reads provided from STDIN are not correctly"\
                #      " paired, error in:\n"+str(read1[0])+str(read2[0])+"\n")
                sys.stdout.write("reads unpaired\n")
                pairedflag = False
                returnstack.append(read1)
                if read2 != 'EOF':
                    returnstack.append(read2)
            else:
                error("ERROR: file error\ninterleaved read file provided ("\
                      "--ri), but reads are not properly paired\n")
    # Once paired info has been detected, parse rest of the reads accordingly
    i = 1
    while i != stacklim:
        # get paires until there are no more reads to parse
        if pairedflag:
            read1 = getread(filehandle)
            if read1 == 'EOF':
                return (returnstack,False)
                #break
            if file2handle:
                read2 = getread(file2handle)
            else:
                read2 = getread(filehandle)
            if read2 == 'EOF':
                return (returnstack,False)
            # check for correctly paired reads through whole data set
            if ispaired(read1,read2):
                returnstack.append([read1,read2])
                i += 1

            # treat reads if unpaired if not properly paired
            else:
                error("ERROR: Paired reads were set but reads:\n"\
                                 ""+read1[0]+"\nand\n"+read2[0]+"\nare not paire"\
                                 "d\n")
        # unpaired reads parsing
        else:
            read1 = getread(filehandle)
            if read1 == 'EOF':
                break
            else:
                returnstack.append(read1)
                i += 1

    return (returnstack,True)


def getread(filehandle):
    """
    parse a single read from file with each part of the read as elements of a
    list
    ARGS:
        INPUT - filehandle, file with read data
        RETURN - readlist, read information
    """
    # list to return
    readlist = []
    # read 4 non-empty lines corresponding to a single read
    for i in xrange(0,4):
        # ignore empty lines
        readstr = readnoempty(filehandle)
        # return End Of File when no more data to parse
        if readstr == 'EOF':
            return readstr
        # remove end of line new line (chomp)
        readstr = readstr.rstrip()
        # Add read element to readlist
        readlist.append(readstr)
    # loosly check for correct formating
    if (readlist[0][0] != '@') or (readlist[2][0] != '+'):
        error("ERROR: format error\ndata provided not in fasta format")
    return readlist


def readnoempty(filehandle):
    """
    read line from file ignoring empty lines (only new line) and exit if end
    of file has been reached
    ARGS:
        INPUT - filehandle, file to read data from
        RETURN - line form file
    """
    # begin assuming file is not empty
    noempty = True
    # as long as there is data in file
    while noempty:
        # get line from file
        string = filehandle.readline()
        # if string is empty, end of file has been reached
        if string == '':
            return 'EOF'
        # ignore empty lines (only newline)
        elif string == '\n':
             pass
        # break loop and return line
        else:
             noempty = False
    return string


def ispaired(read1,read2):
    """
    check that two reads are properly paired (shared header identifier except
    for mate flag)
    ARGS:
        INPUT - read1|read2, reads to compare
        RETURN - bool, True for paired, False for not paired
    """
    # get header ID
    pairid = (read1[0].split(' ')[0], read2[0].split(' ')[0])
    # If not equal, not paired
    if pairid[0] != pairid[1]:
        return False
    # paired
    else:
        return True

def writefq(filehandle, header, seq, qual,\
            headerp = False, seqp = False, qualp = False):
    """
    write to properly formated fastq file
    ARGS:
        INPUT - filehandle, file to write
                header, self explanatory
                seq, read sequence
                qual, quality string for sequence
        RETURN - None
    """
    # writting for paired reads
    if headerp:
         line = header+'\n'+seq+'\n'+'+'+'\n'+qual+'\n'+\
                headerp+'\n'+seqp+'\n'+'+'+'\n'+qualp+'\n'
         filehandle.write(line)
    # writting for unpaired reads
    else:
        line = header+'\n'+seq+'\n'+'+'+'\n'+qual+'\n'
        filehandle.write(line)
    return None



def demultiplex(readstack,indexdict,readnum,pairflag = False):
    """
    gets read lists from readstack and writes them to a fastq file specific for
    a barcode. If no barcode matches, reads are wriyten to NOEXACTMATCH.fq
    ARGS:
        INPUT - readstack, list of read lists to demultiplex
                indexdict, dictionary of barcode information
                readnum, number of reads demultiplexed so far
        RETURN - (nomatch,readnum), reads with no barcode match, number of
                                    reads demultiplexed
    """
    # string of non demultiplexed reads in fastq format
    global nummismatched
    nomatch = ''
    while True:
        try:
            # get read list
            readlist = readstack.pop()
            # determine if it's paired or unpaired
            # unpaired
            if len(readlist) == 4:
                r1_head = readlist[0]
                r1_seq = readlist[1]
                r1_qual = readlist[3]
                read = [r1_head, r1_seq, r1_qual]
            # paired
            elif len(readlist) == 2:
                pairflag = True
                r1_head = readlist[0][0]
                r1_seq = readlist[0][1]
                r1_qual = readlist[0][3]
                r2_head = readlist[1][0]
                r2_seq = readlist[1][1]
                r2_qual = readlist[1][3]
                read = [r1_head, r1_seq, r1_qual, r2_head, r2_seq, r2_qual]
            # something very bad
            else:
                print "Something happended bery bad"
                print readlist
                #sys.exit(1)

            # look for exact match between read and barcode sequence
            try:
                barcode = indexdict[r1_seq[0:9]]
                nomatch += arrange(barcode, read, pairflag)
            except KeyError:
                #try to find barcode with mismatches
                #print "mismatched barcode"
                if  nostrict:
                    barcodes = hammingS(r1_seq[0:9],indexdict)
                    if len(barcodes) == 1:
                        nomatch += arrange(indexdict[barcodes[0]],read,pairflag)
                        nummismatched += 1
                    else:
                        if pairflag:
                            line = r1_head+'\n'+r1_seq+'\n'+'+'+'\n'+r1_qual+'\n'+\
                            r2_head+'\n'+r2_seq+'\n'+'+'+'\n'+r2_qual+'\n'
                        # unpaired
                        else:
                            line = r1_head+'\n'+r1_seq+'\n'+'+'+'\n'+r1_qual+'\n'
                        nomatch = nomatch + line
                else:
                # for paired
                    if pairflag:
                        line = r1_head+'\n'+r1_seq+'\n'+'+'+'\n'+r1_qual+'\n'+\
                               r2_head+'\n'+r2_seq+'\n'+'+'+'\n'+r2_qual+'\n'
                    # unpaired
                    else:
                        line = r1_head+'\n'+r1_seq+'\n'+'+'+'\n'+r1_qual+'\n'
                    nomatch = nomatch + line
            # increase read number
            readnum += 1
        # end loop if no more sequences are available for demultiplexing
        except IndexError:
            break
    return (nomatch,readnum)

def arrange(barcode, read, pairflag):
    """
    Set pairing information and write to file
    In light of rewritten code, this function has been created
    !! Further analysis is needed!!!
    """
    global numdem
    nomatch = ''
    r1_head = read[0]
    r1_seq = read[1]
    r1_qual = read[2]
    if pairflag:
        r2_head = read[3]
        r2_seq = read[4]
        r2_qual = read[5]

    # double barcode indexfile
    if len(barcode) > 1:
        breakflag = 0
        for code in barcode:
            if r2_seq[0:9] == code[2]:
                writefq(code[1], r1_head, r1_seq, r1_qual,\
                    r2_head, r2_seq, r2_qual)
                numdem += 1
                breakflag = 1
                break
        if breakflag == 1:
            pass
        else:
            # If second barcode does not match any input, write to NOMATCH file
            line = r1_head+'\n'+r1_seq+'\n'+'+'+'\n'+r1_qual+'\n'+\
                   r2_head+'\n'+r2_seq+'\n'+'+'+'\n'+r2_qual+'\n'
            nomatch += line
            return nomatch

    # only one barcode sequence provided
    # if no second sequence provided, write to file
    elif barcode[0][2] == 'None':
        # either paired
        if pairflag:
            writefq(barcode[0][1], r1_head, r1_seq, r1_qual,\
                r2_head, r2_seq, r2_qual)
            numdem += 1
        # unpaired
        else:
            writefq(barcode[0][1], r1_head, r1_seq, r1_qual)
            numdem += 1

    # when two barcode sequences are provided, second sequence
    # must also match
    else:
        # for paired
        if pairflag:
            if r2_seq[0:9] == barcode[0][2]:
                writefq(barcode[0][1], r1_head, r1_seq, r1_qual,\
                        r2_head, r2_seq, r2_qual)
                numdem += 1
            else:
                line = r1_head+'\n'+r1_seq+'\n'+'+'+'\n'+r1_qual+'\n'+\
                       r2_head+'\n'+r2_seq+'\n'+'+'+'\n'+r2_qual+'\n'
                nomatch += line
                return nomatch
        # if unpaired and two barcodes are provided, second barcode must be
        # present in read1 3' end
        elif r1_seq[-9:] == barcode[0][2]:
            writefq(barcode[0][1], r1_head, r1_seq, r1_qual,\
                    r2_head, r2_seq, r2_qua)
            numdem += 1

        # if second barcode sequence does not match, discard
        else:
            line = r1_head+'\n'+r1_seq+'\n'+'+'+'\n'+r1_qual+'\n'
            nomatch = nomatch + line
            return nomatch
    return nomatch


# Compute hamming dustance between sequence index and all other indexes
# NAIVE
def hammingS(index,indexdict):
    # list of matching distances
    matchlist = []
    # Loop over indexes
    for keyindex in indexdict.keys():
        includeflag = True
        # Check for same size for computation of hamming distance
        if len(index) == len(keyindex):
            hammings = 0
            for i in xrange(0,len(index)):
                if index[i] != keyindex[i]:
                    hammings += 1
                if hammings > 2: # hamming distance tolerance
                    includeflag = False
                    break
        else:
            print "Strings\n",index,keyindex,"of different lenght"
            sys.exit()

        if includeflag:
            matchlist.append(keyindex)
        else:
            pass
    # Compare base by base until finished of hamming distance exceeded
    return matchlist


def analizeinput(datainput,indexdict,noexactmatch,datainput2=False,pairedflag=True):
    """
    from input filehandle, parse fq sequences and demultiplex them
    """
    # Total read counter
    readnum = 0
    # Flag to keep track of further to demultiplex sequences
    stackflag = True
    while stackflag:
        # parse fastq file according to paired information
        if datainput2:
            readstack = parsefq( datainput,datainput2,pairedflag=True)
        else:
            readstack = parsefq( datainput, pairedflag=pairedflag )
        # if all reads have been parsed, turn off stackflag
        if not readstack[1]:
            stackflag = False
        #print len(readstack[0])
        # demultiplex reads and return non matched sequences
        nomatch, readnum = demultiplex( readstack[0], indexdict, readnum )
        noexactmatch.write(nomatch)
        # log reads parsed
        if not ( readnum % 1000 ):
            sys.stdout.write( str( readnum ) + " reads parsed" )
            sys.stdout.flush()
            sys.stdout.write( '\r' )
    return readnum
