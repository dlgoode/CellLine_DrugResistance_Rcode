### Useful python code for running through directories and performing operations:

import os
import sys
import glob

pwd = os.path.dirname(os.path.abspath(__file__))


print "Current dir is "+pwd

for arg in sys.argv:
        print "Arg is "+arg

 F_out = os.path.join(FREECoutdir,dirname)
        #       print os.path.exists( F_out )

        if os.path.exists( F_out ) is False:
                print "Making an output dir for "+dirname
                os.system( "mkdir %s" %( F_out ) )


indir = "".join([basedir,dirname,"/Bam"])

        print indir;

        ### Run QDNAseq
        if os.path.exists( indir ) is True:
 ### Get the bamfile name
                for bamfile in glob.glob(os.path.join(indir, '*.bam')):
                        print bamfile
                        infile = bamfile
                #"FREEC_",dirname,"=`",
                FREEC_CMD = "".join([qsub_base," -N FREEC_",dirname,FREEC_cmd,infile,FREEC_window,F_\
out,"'"])
                print FREEC_CMD

#               os.system( FREEC_CMD )

                Freec_outfile = infile.split("/")
                print len(Freec_outfile)
                outdir = "".join([F_out,"/",Freec_outfile[len(Freec_outfile)-1],"/"])

                                ### " -W depend=afterok:$FREEC_",dirname,
                ADD_PVALS_CMD = "".join([qsub_base,assess_sig_cmd,outdir,Freec_outfile[len(Freec_outfile)-1],"_CNVs ",outdir,Freec_outfile[len(Freec_outfile)-1],"_ratio.txt'"])
                print ADD_PVALS_CMD

                os.system( ADD_PVALS_CMD )


 files=glob.glob(os.path.join('*/Bam/*.bam'))

>>> files[0].split("/")
['FD02699715', 'Bam', 'CBCRRANXX.sorted.bam']
>>> files[0].split("/")[0]
'FD02699715'

>>> "".join([infile_name[0],"_sorted.bam"])
'FD02699715_sorted.bam'

outfile="".join([infile_name[0],"_sorted.bam"])


mydir="".join([infile_name[0],"/",infile_name[1]])

mvcmd="".join(["mv ",files[0]," ",mydir,"/",outfile])