#!/usr/bin/python
###############################################################################
# Julan Regalado - julian.regalado@tuebingen.mpg.de
#
#
# GNU GENERAL PUBLIC LICENSE
# Version 2, June 1991
#
#
# Demultiplex fastq files
#
###############################################################################

# modules
import sys
import time
import gzip
import demfnc # functions for this script present in demfnc.py

# main program
def main():
    logfile = open("plexSeq.log",'w')
    start = time.time() # log time
    gzipflag = False    # output in gzip compressed (not implemented so far)
    stackflag = True    # data for demultiplexing exists
    demfnc.numdem = 0   # demultiplexed read counters
    demfnc.nummismatched = 0 # undemultiplexed read counter
    demfnc.usage()


    # get options
    optdict = demfnc.getoptions()
    sys.stdout.write("\nParameters:\n")
    logfile.write("\nParameters:\n")
    for opt in optdict.keys():
        sys.stdout.write( opt + ":\t" + str(optdict[opt]) + "\n" )
        logfile.write(opt + ":\t" + str(optdict[opt]) + "\n")
    sys.stdout.write('\n')
    logfile.write('\n')
    time.sleep(1)


    # change strict flag according to parameters
    if optdict['no_strict']:
        demfnc.nostrict = True
    else:
        demfnc.nostrict = False


    # read index file and open noexactmatch file handler
    indexdict, noexactmatch = demfnc.readindex( optdict['index_file'], \
                                               optdict['outdir'], gzipflag )


    # get data input, parse and demultiplex TODO: detect if reads from stdin are
    # paired in sequence
    if 'stdin' in optdict.keys(): # check if data comes from STDIN
        # stackflag remains true as long as there is data to demultiplex
        nomatch,readnum = demfnc.analizeinput(optdict['stdin'],indexdict,noexactmatch,pairedflag=True)

    # check if two read files are provided (paired)
    elif 'fastq2filename' in optdict.keys():
        # check if data is gzip compressed
        # for file 1
        ftype = demfnc.file_type(optdict['fastqfile'])
        if ftype == 'gz':
            fastqfile = gzip.open(optdict['fastqfile'],'rb')
        else:
            fastqfile = open(optdict['fastqfile'],'r')
        # for file 2
        ftype = demfnc.file_type(optdict['fastq2filename'])
        if ftype == 'gz':
            fastq2file = gzip.open(optdict['fastq2filename'],'rb')
        else:
            fastq2file = open(optdict['fastq2filename'],'r')
        #demultiplex
        readnum = demfnc.analizeinput(fastqfile,indexdict,noexactmatch,fastq2file)
        #noexactmatch.write( nomatch )

    # check if data comes from an interleaved fastq file (paired)
    elif 'int_fqfile' in optdict.keys():
        # check for compression
        ftype = demfnc.file_type(optdict['int_fqfile'])
        if ftype == 'gz':
            intfqfile = gzip.open(optdict['int_fqfile'],'rb')
        else:
            intfqfile = open(optdict['int_fqfile'],'r')

        #demultiplex
        nomatch,readnum = demfnc.analizeinput(intfqfile,indexdict,noexactmatch,stackflag=True)


    # check for single file (unpaired)
    else:
        # check for compression
        ftype = demfnc.file_type(optdict['fastqfile'])
        if ftype == 'gz':
            fastqfile = gzip.open(optdict['fastqfile'],'rb')
        else:
            fastqfile = open(optdict['fastqfile'],'r')

        # demultiplex
        nomatch,readnum = demfnc.analizeinput(fastqfile,indexdict,noexactmatch,pairedflag=False)


    # log time
    end = time.time()
    total = (end - start) / float(60)
    pctdemplexed = demfnc.numdem / float(readnum)
    pctmismatched = demfnc.nummismatched / float(demfnc.numdem)

    sys.stdout.write(str(readnum)+" total reads.")
    sys.stdout.write("\nFinished in: "+str(total)+" minutes\n")
    sys.stdout.write(str(demfnc.numdem)+" reads demultiplexed, "+str(pctdemplexed*100)+"%\n")
    sys.stdout.write("Of these, "+str(demfnc.nummismatched)+" reads with mismatches, "+str(pctmismatched*100)+"%\n")

    logfile.write(str(readnum)+" total reads.\n"+"\nFinished in: "+str(total)+\
                  " minutes\n"+str(demfnc.numdem)+" reads demultiplexed, " +\
                  "Of these, " + str(pctdemplexed*100)+"%\n" +\
                  str(demfnc.nummismatched)+" reads with mismatches, " +\
                  str(pctmismatched*100)+"%\n")
    sys.exit(0)



if __name__ == "__main__":
    sys.stdout.write("\nseqplex\nSequencing read demultiplexer\nv-1.0\n\n")
    main()
