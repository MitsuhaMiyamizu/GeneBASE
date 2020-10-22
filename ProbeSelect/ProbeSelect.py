#Exon Array Probe Selection
#Xing, Y., Kapur, K., and Wong, W.H. 2006. Probe selection and expression index computation of Affymetrix Exon Arrays. PLoS ONE in press.
#For questions and comments, please contact yi-xing@uiowa.edu or whwong@stanford.edu

import re
import sys
import os
import getopt

data_dir=sys.argv[1]
os.cur_dir=data_dir
array_list=os.listdir(os.cur_dir)
n_array=len(array_list)
FileHandlers=[]
array_name=[]

write_to_dir='out'

#read command line options
opts=sys.argv[2:]
try:
    optlist,args=getopt.getopt(opts,'t:r:')
except getopt.GetoptError:
    print "Usage: ProbeSelect.py DATA_DIR -t [T/F] -r [T/F]"
    sys.exit(2)
for o,a in optlist:
    if o=="-t":
        trim_zero=a
    if o=="-r":
        run_r=a

for Array in array_list:
    array_name.append(Array[:-4])
    data_file=data_dir+'/'+Array
    data_ifile=open(data_file,"r")
    FileHandlers.append(data_ifile)

#read one line from all files
def readline_all(FileHandlers):
    n_file=len(FileHandlers)
    lines=[]
    for i in range(n_file):
        line=FileHandlers[i].readline()
        line=line[:-1]
        lines.append(line)
    return lines

current_transcript_cluster_id='-1'

#skip the first two lines
readline_all(FileHandlers)
readline_all(FileHandlers)

while(1):
    lines=readline_all(FileHandlers)

    line=lines[0]
    line_data=re.split('\t',line)

    #Stop at transcript clusters with only extended and full probes
    if line_data[0]=='#ExtendedTranscriptClusters':
        ofile.close()
        break
    
    transcript_cluster_id=line_data[0]
    probeset_id=line_data[1]
    probe_level=line_data[2]
    probe_id=line_data[3]
    probe_intensity=line_data[5]

    if (transcript_cluster_id!=current_transcript_cluster_id):
        if (current_transcript_cluster_id!='-1'):
            ofile.close()
        
        current_transcript_cluster_id=transcript_cluster_id
    
        ofile=open(write_to_dir+'/'+str(transcript_cluster_id)+'.txt','w')
        print "Processing "+str(transcript_cluster_id)+'...'
        ofile.write(transcript_cluster_id)
        for Array in array_list:
            ofile.write("\t"+Array)
        ofile.write("\n")
    ofile.write(str(probe_id)+"|"+str(probeset_id)+"|"+str(probe_level))


    for line in lines:
        line_data=re.split("\t",line)
        probe_intensity=line_data[5]

        if (trim_zero!='F'):
            if float(probe_intensity)<=0:
                probe_intensity=0
        ofile.write("\t"+str(probe_intensity))
    ofile.write("\n")

if (run_r!='F'):
    import rpy
    rpy.r.source("ProbeSelect.R")




        



