#!/usr/bin/python

import os
import sys
import getopt

genome_size = 'hg19.chrom.size'
script = 'script'

def opt(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:",["ifile=","ofile=","name="])
    except getopt.GetoptError:
        print 'call_domains.py -i <input directory> -o <output directory> -n <name>'
        sys.exit(2)
    global cell,hic_dir,output_DIR
    if len(opts) != 3:
        print "Missing arguments!"
        print 'call_domains.py -i <input directory> -o <output directory> -n <name>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'call_domains.py -i <input directory> -o <output directory> -n <name>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            hic_dir = arg
        elif opt in ("-o", "--ofile"):
            output_DIR = arg
        elif opt in ("-n","--name"):
            cell = arg

def make_DI():
    j=1
    while j<=23:
        print(j)
        if j==23:
            chr='chrX'
        else:
            chr='chr'+str(j)
        os.system('python '+ script + '/juicer2DI.py -i ' + hic_dir + '/' + cell + '.observed.KR.' + chr + '.BP.10000.txt.gz -g ' + genome_size + ' -b 10Kb -w 2Mb -c ' + chr + ' -o '+ output_DIR + '/DI.' + chr + '.10000.2000000')
        j+=1

def make_hmm_result():
	i=1
	while i<=23:
		print i
		if i==23:
			chr='chrX'
		else:
			chr='chr'+str(i)
		input1=open(output_DIR+'/DI.'+chr+'.10000.2000000.output','r')
		input2=open(output_DIR+'/DI.'+chr+'.10000.2000000','r')
		all_input1=input1.readlines()
		all_input2=input2.readlines()
		j=0
		state_3_list=[] # set positive
		state_2_list=[] # set negative
		state_1_list=[] # set neutral
		while j<len(all_input2):
			if float(all_input1[j+2].split()[1])==3:
				state_3_list.append(float(all_input2[j].split()[3]))
			if float(all_input1[j+2].split()[1])==2:
				state_2_list.append(float(all_input2[j].split()[3]))
			if float(all_input1[j+2].split()[1])==1:
				state_1_list.append(float(all_input2[j].split()[3]))
			j+=1

		myList=[sum(state_3_list),sum(state_2_list),sum(state_1_list)]
		idx=[j[0] for j in sorted(enumerate(myList), key=lambda x:x[1])]
		print myList, idx
		map_dic={}
		map_dic[[3,2,1][idx[0]]]=1
		map_dic[[3,2,1][idx[1]]]=2
		map_dic[[3,2,1][idx[2]]]=3
		print map_dic
		order=[[3,2,1][idx[0]],[3,2,1][idx[1]],[3,2,1][idx[2]]]
		print order
		output1=open(output_DIR+'/hmm_7colfile.chr'+str(i),'w')
		j=0
		while j<len(all_input2):
			each1=all_input2[j].split()
			output1.write(each1[0]+'\t'+each1[1]+'\t'+each1[2]+'\t')
			each2=all_input1[j+2].split()
			output1.write(str(map_dic[int(float(each2[1]))])+'\t'+each2[order[0]+1]+'\t'+each2[order[1]+1]+'\t'+each2[order[2]+1]+'\n')
			j+=1
		output1.close()
		input1.close()
		input2.close()
		i+=1

if __name__ == "__main__":
    opt(sys.argv[1:])
    make_DI()
    os.system('matlab -nodisplay -nojvm < HMM_calls.m > '+output_DIR+'/dumpfile')
    make_hmm_result()

    chr_list=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr23']

    for i in chr_list:
        print i
        os.system('perl '+script+'/hmm_probablity_correcter.pl '+output_DIR+'/hmm_7colfile.'+i+' 2 0.99 10000 | perl '+script+'/hmm-state_caller.pl '+genome_size+' '+i+' | perl  '+script+'/hmm-state_domains.pl > '+output_DIR+'/finaldomainclls.'+i)

    input1=open(output_DIR+'/finaldomainclls.'+i,'r')
    all_input1=input1.readlines()
    input1.close()

    output1=open(output_DIR+'/finaldomainclls.'+i,'w')
    for j in all_input1:
        if len(j.split())==3:
            output1.write(j)
    output1.close()
