#!/usr/bin/env python3

## The Following code has been adapted from the RSeQC infer_experiment function.
## RSeQC is distributed under GNU General Public License (GPLv3)
## 
## This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
##
## Contact
## Liguo Wang: wangliguo78@gmail.com
## Shengqin Wang: wzsqwang@gmail.com
## Wei Li: superliwei@gmail.com
## 
## Reference
## Wang, L., Wang, S., & Li, W. (2012). RSeQC: quality control of RNA-seq experiments. Bioinformatics (Oxford, England), 28(16), 2184–2185. http://doi.org/10.1093/bioinformatics/bts356
## Wang, L., Nie, J., Sicotte, H., Li, Y., Eckel-Passow, J. E., Dasari, S., et al. (2016). Measure transcript integrity using RNA-seq data. BMC Bioinformatics, 17(1), 1–16. http://doi.org/10.1186/s12859-016-0922-z
## 

import os,sys
if sys.version_info[0] != 3:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
	sys.exit()	

import subprocess
import pkg_resources

required = {'pysam','RSeQC'} 
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    # implement pip as a subprocess:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install',*missing])

import pysam
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
from time import strftime

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="4.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('class.log','a')
	print(mesg, file=sys.stderr)
	print(mesg, file=LOG)


def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Input alignment file in SAM or BAM format")
	parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed fomat.")
	parser.add_option("-s","--sample-size",action="store",type="int",dest="sample_size",default=200000, help="Number of reads sampled from SAM/BAM file. default=%default")	
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be considered as \"uniquely mapped\". default=%default")
	parser.add_option("-o","--out",action="store",type="string",dest="output_file",default="infer_result",help=" default=%infer_result")

	(options,args)=parser.parse_args()

	if not (options.input_file and options.refgene_bed):
		parser.print_help()
		print('\n\n' + __doc__, file=sys.stderr)
		sys.exit(0)
	for f in (options.input_file,options.refgene_bed):
		if not os.path.exists(f):
			print('\n\n' + f + " does NOT exists." + '\n', file=sys.stderr)
			sys.exit(0)
	if options.sample_size <1000:
		print("Warn: Sample Size too small to give a accurate estimation", file=sys.stderr)
	pysam.index(options.input_file)
	obj = SAM.ParseBAM(options.input_file)
	(protocol,sp1,sp2,other)=obj.configure_experiment(refbed=options.refgene_bed, sample_size = options.sample_size, q_cut = options.map_qual)
	if other <0: other=0.0
	file_object  = open(options.output_file + ".txt", "w")
	if protocol == "PairEnd":
		file_object.write("This is Paired End Data\n")
		file_object.write("Fraction of reads failed to determine: %.4f" % other + "\n")
		file_object.write("Fraction of reads explained by \"1++,1--,2+-,2-+\": %.4f" % sp1 + "\n")
		file_object.write("Fraction of reads explained by \"1+-,1-+,2++,2--\": %.4f" % sp2 + "\n")
		if sp1 > 2 * sp2:
			file_object.write("\nExperiment is likely \"1++,1--,2+-,2-+\" (HTSeq.count --forward)\n")
		if sp2 > 2 * sp1:
			file_object.write("\nExperiment is likely \"1+-,1-+,2++,2--\" (HTSeq.count --reverse)\n")

	elif protocol == "SingleEnd":
		file_object.write("This is Single End Data\n")
		file_object.write("Fraction of reads failed to determine: %.4f" % other + "\n")
		file_object.write("Fraction of reads explained by \"++,--\": %.4f" % sp1 + "\n")
		file_object.write("Fraction of reads explained by \"+-,-+\": %.4f" % sp2 + "\n")
		if sp1 > 2 * sp2:
			file_object.write("\nExperiment is likely \"++,--\" (HTSeq.count --forward)\n")
		if sp2 > 2 * sp1:
			file_object.write("\nExperiment is likely \"+-,-+\" (HTSeq.count --reverse)\n")
		
	else:
		file_object.write("Unknown Data type\n")
	#print mesg
	file_object.close()

if __name__ == '__main__':
	main()
