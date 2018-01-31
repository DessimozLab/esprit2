#!/bin/bash

usage(){
cat << EOF
$0 - find split genes using Esprit 2.0

$0 [options] gene_prefix gene_family_folder orthoxml

This program computes split genes in a given genome using the methodoly 
described in http://doi.org/...

The program currently requires a SGE (Sun Grid Engine) environment 
to run properly. This might still change in the future.

The following options and parameters can be set.

Parameters:
-----------

gene_prefix     a prefix of the gene identifiers that indicates possible
                split genes, i.e. a prefix of the genes in the genome 
                you want to search for split genes. This argument is 
                required.
                Example: "Traes_1DS" will search for split genes in the 
                     wheat 1DS chomosome (all ids start with this prefix)

gene_family_folder  
                A folder that contains one fasta file per gene family. The
                files in there should be sequencially nummered and having 
                a HOG prefix, e.g. HOG00001.fa
                If it is not provided, it defaults to ./HOGFasta

orthoxml        the file that contains the HOGs for all genomes in the
                analysis. This file should match the data in the gene 
                family folder.
                If not set, it defaults to ./HierarchicalGroups.orthoxml


Options:
--------

-l              minimum length of candidate genes, defaults to 50 if not set.
                Shorter genes will be excluded.

-o              maximum overlap in the alignment, defaults to 0.1 if not set.
                This refers to the maximum possible overlap that two fragments
                might have in the alignment to still be considered a possible 
                split gene.

-s              min size of the gene family including the fragements, defaults
                to 5 if not set. This means, a family must have at least 3 
                reference sequences and 2 fragments.
                
-b              number of bootstrap replicates, defaults to 100 if not set.

-c              threshold for collapsing, defaults to 0.95 if not set. This 
                refers to the number of bootstrap replicates that need to
                contain a certain split in order to not be collapsed.

-a              significance of the likelihood ratio test, defaults to 0.01.

-g              filename of gff input file. IDs of gene features need to match
                those of the input sequences. By default the merged output
                gff file will be named same as input + '.merged'. Use the -w 
                to overwrite.

-w              filename of gff output file. Only affects if -g is specified.

EOF
}

min_len="50"
ovlp_len="0.1"
family_size="5"
n_samples="100"
col_thr="95"
lrt_sign="0.01"
gff=""
gff_out=""
while getopts "hvl:o:s:b:c:a:g:w:" opt; do
    case $opt in
    h) usage
       exit 0
       ;;
    v) echo "Esprit 2.0"
       exit 0
       ;;
    l) min_len="$OPTARG"
       if ! [[ $min_len =~ ^[1-9][0-9]*$ ]] ; then
          echo "-$opt requires positive integer argument"
          usage
          exit 1
       fi
       ;;
    o) ovlp_len="$OPTARG"
       ;;
    s) family_size="$OPTARG"
       ;;
    b) n_samples="$OPTARG"
       ;;
    c) col_thr="$OPTARG"
       ;;
    a) lrt_sign="$OPTARG"
       ;;
    g) gff="$OPTARG"
       ;;
    w) gff_out="$OPTARG"
       ;;
    esac
done
shift $((OPTIND-1))

unique_str="$1"
gene_fam_dir="${2:-./HOGFasta}"
orthoxml="${3:-./HierarchicalGroups.orthoxml}"

if [ -z "$unique_str" ] || [ ! -f "$orthoxml" ] || [ ! -d "$gene_fam_dir" ] ; then
    >&2 usage 
    >&2 echo "invalid parameters: \"$unique_str\" \"$gene_fam_dir\" \"$orthoxml\""
    exit 1
fi


source ./load_env
pip install -r ./requirements.txt
if [ "$?" != "0" ] ; then
    >&2 echo "failed to install python dependencies."
    exit 1
fi

# check that mafft, phyml and fasttree are installed and in PATH
allok="true"
for prog in mafft FastTree ; do
    if [ ! $(which ${prog} 2>/dev/null) ] ; then
        >&2 echo "Could not detect \"$prog\" in PATH. Please make sure \"$prog\" is "
        >&2 echo "installed and accessible from the PATH."
        allok="false"
    fi
done
if [ "$allok" != "true" ] ; then
    exit 1
fi


#go to the working directory, put HierarchicalGroups.orthoxml there
#get OMA <-> Ensembl mapping
cat organise_files.txt > mapping.sh
cat >> mapping.sh << EOF
grep "protId=" $orthoxml | grep $unique_str | awk -F\" '{ for(i=2; i<=NF; i=i+2){ a = a"\""\$i"\""",\t";} {print a; a="";}}' | awk '{\$2="";print}'|tr -d '"'| tr -d ',' > mapping.txt
EOF

qsub -N mapping mapping.sh 

#get HOGs of interest

cat > orthoxmlquery.py << EOF
#Author: Adrian Altenhoff

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future.builtins import str
from future import standard_library
standard_library.install_hooks()


class ElementError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return str(self.msg)


class OrthoXMLQuery(object):
    """Helper class with predefined queries on an orthoxml tree."""

    ns = {"ns0": "http://orthoXML.org/2011/"}   # xml namespace

    @classmethod
    def getToplevelOrthologGroups(cls, root):
        """returns a list with the toplevel orthologGroup elements
        of the given root element."""
        xquery = ".//{{{ns0}}}groups/{{{ns0}}}orthologGroup".format(**cls.ns)
        return root.findall(xquery)

    @classmethod
    def getTaxRangeNodes(cls, root, recursively=True):
        xPrefix = ".//" if recursively else "./"
        xquery = '{}{{{}}}property[@name="TaxRange"]'.format(xPrefix,
                                                             cls.ns['ns0'])
        return root.findall(xquery)

    @classmethod
    def getGeneRefNodes(cls, root, recursively=True):
        iterfn = root.iter if recursively else root.iterchildren
        iterator = iterfn('{{{}}}geneRef'.format(cls.ns['ns0']))
        return list(iterator)

    @classmethod
    def getGeneFromId(cls, id_, root):
        xquery = ".*//{{{}}}gene[@id='{}']".format(cls.ns['ns0'], id_)
        genes = root.findall(xquery)
        if len(genes) > 1:
            raise ElementError('several gene nodes with id {} '
                               'exist'.format(id_))
        gene = genes[0] if len(genes)>0 else None
        return gene

    @classmethod
    def getGroupsAtLevel(cls, level, root):
        """returns a list with the orthologGroup elements which have a
        TaxRange property equals to the requested level."""
        xquery = (".//{{{0}}}property[@name='TaxRange'][@value='{1}']/..".
                  format(cls.ns['ns0'], level))
        return root.findall(xquery)

    @classmethod
    def getSubNodes(cls, targetNode, root, recursively=True):
        """method which returns a list of all (if recursively
        is set to true) or only the direct children nodes
        having 'targetNode' as their tagname.
        The namespace is automatically added to the tagname."""
        xPrefix = ".//" if recursively else "./"
        xquery = "{}{{{}}}{}".format(xPrefix, cls.ns['ns0'], targetNode)
        return root.findall(xquery)

    @classmethod
    def is_geneRef_node(cls, element):
        """check whether a given element is an instance of a geneRef
        element."""
        return element.tag == '{{{ns0}}}geneRef'.format(**cls.ns)

    @classmethod
    def getLevels(cls, element):
        """returns a list of the TaxRange levels associated to the
        passed orthologGroup element. If the element does not have
        any TaxRange property tags associated, an empty list is
        returned."""
        propTags = cls.getSubNodes("property", element, recursively=False)
        res = [t.get('value') for t in propTags if t.get('name') == 'TaxRange']
        return res

    @classmethod
    def getInputGenes(cls, root, species=None):
        """returns a list of all gene elements in the orthoxml inside
        <species><database> tags, i.e. the list of genes prior to running
        OMA-HOGS. Optionally filtered by species."""
        filter_ = ('[@name="{}"]'.format(species)
                   if species is not None else '')
        if filter_ > '':
            xquery = ('/ns:orthoXML/ns:species{}/ns:database/'
                      'ns:genes//ns:gene'.format(filter_))
        else:
            xquery = '//ns:gene'
        return root.xpath(xquery, namespaces={'ns': cls.ns['ns0']})

    @classmethod
    def getGroupedGenes(cls, root, species=None):
        """ returns a list of all geneRef elements inside <group> tags, i.e.
        the list of genes clustered into families after running OMA-HOGS.
        Optionally filtered by species."""
        filter_ = ('[@name="TaxRange"and@value="{}"]'.format(species)
                   if species is not None else '')
        if filter_ > '':
            xquery = ('/ns:orthoXML/ns:groups/ns:orthologGroup//ns:property{}/'
                      'following-sibling::ns:geneRef'.format(filter_))
        else:
            xquery = '//ns:geneRef'
        return root.xpath(xquery, namespaces={'ns': cls.ns['ns0']})
EOF



cat > get_candidates.py << EOF
# -*- coding: utf-8 -*- 

import lxml.etree as etree
import sys, os
import orthoxmlquery as oq
from Bio import SeqIO
import requests, random, subprocess
from xml.dom import minidom
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Phylo.Applications import PhymlCommandline 
import glob
import os

def read_fasta(loc): #adapted for empty lines in HOGs.fasta
    seqs = dict()
    tmp = open(loc, 'r')
    lines = tmp.readlines()
    key = lines[0].split('>')[-1].split()[-2]
    value = ''
    for i in range(1, len(lines)):
        if len(lines[i]) == 1:
            seqs[key] = value
        elif lines[i].split()[0][0] == '>':
            key = lines[i].split('>')[-1].split()[-2]
            value = ''
        else:
            value += lines[i].split()[0]
    return seqs

def write_dictionary(d, out):
    '''Writes a dictionary (d) to an output file (out) in the following format:
    > key
    d[key]
    >key
    d[key]
    ...'''
    f = open(out, 'w')
    for key in d.keys():
        f.write('>' + key + '\n')
        f.write(d[key] + '\n')
    f.close()
    return


def process_stdout(stdout):
    '''Parse MSA from Mafft stdout (stdout),
    save all genes to a dictionary (msa),
    returns the dictionary (msa)'''
    stdout_process = stdout.split(">")
    msa = dict()
    for i in range(len(stdout_process)):
        if len(stdout_process[i]) > 0:
            key = stdout_process[i].split('\n')[0].split()[0]
            value = ''.join(stdout_process[i].split('\n')[1:])
            msa[key] = value
    return msa


def check_overlap_length(msa, id1, id2):
    ''' '''
    seq1 = msa[id1].split("-")
    lengths1 = [len(i) for i in seq1]
    non01 = [i for i,v in enumerate(lengths1) if v > 0]
    k11 = non01[0] #index of the 1st residual in seq1
    k12 = len(msa[id1])- (len(lengths1)-non01[-1]) #index of the last residual in seq1
    seq2 = msa[id2].split("-")
    lengths2 = [len(i) for i in seq2]
    non02 = [i for i,v in enumerate(lengths2) if v > 0]
    k21 = non02[0] #index of the 1st residual in seq2
    k22 = len(msa[id2])- (len(lengths2)-non02[-1]) #index of the last residual in seq2
    s1 = list(range(k11, k12 + 1))
    s2 = list(range(k21, k22 + 1))
    intrsc = list(set(s1) & set(s2))
    return k11, k12, k21, k22, sorted(intrsc), len(intrsc)/float(len(s1)), len(intrsc)/float(len(s2))


def edit_overlapping_seqs(msa, id1, id2, k11, k21, o_start, o_end):
    new_msa = dict()
    for key in msa:
        new_msa[key] = msa[key]
    h = (o_start - 1) + (o_end - o_start + 1)//2 + ((o_end - o_start + 1)%2) * random.randint(0,1)
    if k11 < k21:
        left = id1
        right = id2
    else:
        left = id2
        right = id1
    new_msa[left] = new_msa[left][0 : h + 1] + 'X' * (o_end - h) + new_msa[left][o_end + 1 : ]    
    new_msa[right] = new_msa[right][0 : o_start] + 'X' * (h - o_start + 1) + new_msa[right][h + 1 :  ]
    return new_msa, h + 1



def concatenate_overlapping_seqs(msa, id1, id2, o_start, o_end):
    new_msa = dict()
    for key in msa:
        new_msa[key] = msa[key]
    seq = ''
    for i in range(0, o_start):
        if new_msa[id1][i] != '-':
            seq += new_msa[id1][i]
        elif new_msa[id2] != '-':
            seq += new_msa[id2][i]
        else:
            seq += '-'
    for i in range(o_start, o_end + 1):
        if new_msa[id1][i] != 'X':
            seq += new_msa[id1][i]
        else:
            seq += new_msa[id2][i]
    for i in range(o_end + 1, len(new_msa[id1])):
        if new_msa[id1][i] != '-':
            seq += new_msa[id1][i]
        elif new_msa[id2] != '-':
            seq += new_msa[id2][i]
        else:
            seq += '-'        
    new_msa[id1 + '_' + id2] = seq
    new_msa.pop(id1)
    new_msa.pop(id2)
    return new_msa


def concatenate_nonoverlapping_seqs(msa, id1, id2):
    new_msa = dict()
    for key in msa:
        new_msa[key] = msa[key]
    seq = ''
    for i in range(len(new_msa[id1])):
        if new_msa[id1][i] != '-':
            seq += new_msa[id1][i]
        elif new_msa[id2] != '-':
            seq += new_msa[id2][i]
        else:
            seq += '-'
    new_msa[id1 + '_' + id2] = seq
    new_msa.pop(id1)
    new_msa.pop(id2)
    return new_msa


loc_mapping = os.getcwd() + '/mapping.txt'
mapping = dict()
target_genes = []
tmp = open(loc_mapping, 'r')
lines = tmp.readlines()
for i in range(len(lines)):
    temp = lines[i].split()
    mapping[temp[0]] = " ".join(temp[1:])
    target_genes.append(temp[0])
tmp.close()

XML = "$orthoxml"
doc = etree.parse(XML)
root_XML = doc.getroot()

topLevel_OG = oq.OrthoXMLQuery.getToplevelOrthologGroups(root_XML)

HOGs = dict()
candidates_number = dict()
for OG in topLevel_OG:
    OG_genes = [ e.get("id") for e in OG.getiterator() if e.tag == "{http://orthoXML.org/2011/}geneRef" ]
    if len(set(OG_genes).intersection(target_genes)) >= 2:
        HOGs[OG.get("id")] = OG_genes
        candidates_number[OG.get("id")] = list(set(OG_genes).intersection(target_genes))


candidates_id = dict()
for key in candidates_number:
    cand = []
    for i in range(len(candidates_number[key])):
        if ';' in mapping[candidates_number[key][i]]:
            cand.append(mapping[candidates_number[key][i]].split(';')[-1][1:])
        else:
            cand.append(mapping[candidates_number[key][i]])
    candidates_id[key] = cand


positions = list()
cuts = list()
lengths = list()


for key in candidates_id:
    hog_loc = "$gene_fam_dir" +'/HOG' + key + '.fa'
    #hog gene sequences
    hog_seqs = read_fasta(hog_loc)
    #write fasta 
    w_temp = os.getcwd() + '/hog_temp/HOG' + key + '.fa'
    write_dictionary(hog_seqs, w_temp)
    #hog alignment
    mafft_cline = MafftCommandline(input=w_temp)
    msa_fasta_hog = os.getcwd() + '/hog_aln/HOG' + key + '.aln' #MSA in FASTA format
    stdout, stderr = mafft_cline()
    with open(msa_fasta_hog, 'w') as handle:
        handle.write(stdout)
    #
    msa = process_stdout(stdout)
    if len(msa) >= `echo $family_size`:
        for i in range(len(candidates_id[key]) - 1):
            for j in range(i + 1, len(candidates_id[key])):
                lengths.append((key, candidates_id[key][i], len(hog_seqs[candidates_id[key][i]]), candidates_id[key][j], len(hog_seqs[candidates_id[key][j]])))
                if len(hog_seqs[candidates_id[key][i]]) >= `echo $min_len` and len(hog_seqs[candidates_id[key][j]]) >= `echo $min_len`:
                    k11, k12, k21, k22, intrsc, o1, o2 = check_overlap_length(msa, candidates_id[key][i], candidates_id[key][j])
                    if o1 < `echo $ovlp_len` and o2 < `echo $ovlp_len`:
                        name = candidates_id[key][i] + '_' + candidates_id[key][j]
                        msa_fasta = os.getcwd() + '/aln/' + name + '.aln'
                        msa_phy = os.getcwd() + '/phy/' + name + '.phy' #MSA in phylip format
                        out_c = os.getcwd() + '/aln_c/' + name + '_c.aln'
                        out_c_phy = os.getcwd() + '/phy_c/' + name + '_c.phy'
                    #convert msa to phylip
                    #concatenate
                        if len(intrsc) == 0:                        
                            write_dictionary(msa, msa_fasta)                        
                            AlignIO.convert(msa_fasta, "fasta", msa_phy, "phylip-relaxed")
                            msa_fasta_con = concatenate_nonoverlapping_seqs(msa, candidates_id[key][i], candidates_id[key][j])
                            positions.append((key, candidates_id[key][i], candidates_id[key][j], k11, k12, k21, k22, -1, -1, o1, o2))
                            n = random.randint(min(k12, k22) + 1, max(k11, k21))
                            if k21 < k11:
                                cuts.append((key, candidates_id[key][j], candidates_id[key][i], n))
                            else:
                                cuts.append((key, candidates_id[key][i], candidates_id[key][j], n))                        
                            write_dictionary(msa_fasta_con, out_c)                       
                            AlignIO.convert(out_c, "fasta", out_c_phy, "phylip-relaxed")
                        else:
                            msa_edited, n = edit_overlapping_seqs(msa, candidates_id[key][i], candidates_id[key][j], k11, k21, intrsc[0], intrsc[-1])
                            if k21 < k11:
                                cuts.append((key, candidates_id[key][j], candidates_id[key][i], n))
                            else:
                                cuts.append((key, candidates_id[key][i], candidates_id[key][j], n))
                            write_dictionary(msa_edited, msa_fasta)      
                            AlignIO.convert(msa_fasta, "fasta", msa_phy, "phylip-relaxed")                    
                            msa_fasta_con = concatenate_overlapping_seqs(msa_edited, candidates_id[key][i], candidates_id[key][j], intrsc[0], intrsc[-1])
                            positions.append((key, candidates_id[key][i], candidates_id[key][j], k11, k12, k21, k22, intrsc[0], intrsc[-1], o1, o2))
                            write_dictionary(msa_fasta_con, out_c)
                            AlignIO.convert(out_c, "fasta", out_c_phy, "phylip-relaxed")


hogs_overlap_positions = os.getcwd() + '/alignment_positions.txt' 
f = open(hogs_overlap_positions, 'w')
for i in range(len(positions)):
    f.write(str(positions[i][0]) + '\t' + str(positions[i][1]) + '\t' + str(positions[i][2]) + '\t' + str(positions[i][3]) + '\t' + str(positions[i][4]) + '\t' + str(positions[i][5]) + '\t' + 
        str(positions[i][6]) + '\t' + str(positions[i][7]) + '\t' + str(positions[i][8]) + '\t' + str(positions[i][9])  + '\t' + str(positions[i][10]) + '\n')
f.close()

hogs_length_sequences = os.getcwd() + '/sequence_lengths.txt'
f = open(hogs_length_sequences, 'w')
for i in range(len(lengths)):
    f.write(lengths[i][0] + '\t' + lengths[i][1] + '\t' + str(lengths[i][2]) + '\t' + lengths[i][3] + '\t' + str(lengths[i][4]) + '\n')
f.close()

hogs_cuts = os.getcwd() + '/cuts.txt'
f = open(hogs_cuts, 'w')
for i in range(len(cuts)):
    f.write(cuts[i][0] + '\t' + cuts[i][1] + '\t' + cuts[i][2] + '\t' + str(cuts[i][3]) + '\n')
f.close()
EOF


cat organise_files.txt > get_candidates.sh
cat >> get_candidates.sh << EOF

mkdir -p hog_temp hog_aln aln \
         phy aln_c phy_c

source ./load_env
python get_candidates.py

rm -rf hog_temp
EOF

qsub -sync y -N get_cand -hold_jid mapping get_candidates.sh

#computing n-1 trees

cd phy_c

for f in *.phy
do
    cat ../one_tree.txt > $f.sh
    cat >> $f.sh << EOF
    source ../load_env
    FastTree $f > ${f}_tree.txt
EOF
done


i=1
for f in *.phy.sh
do
    qsub -N job_n_1_$i -cwd $f
    ((i++))
done

cd ..

cat organise_files.txt > splitting_trees.sh
cat >> splitting_trees.sh << EOF
mkdir -p n_1_trees
find phy_c/ -name "*_tree.txt" -exec mv {} n_1_trees \; &
wait

source ./load_env
python splitting_trees.py n_1_trees/
EOF


cat > splitting_trees.py << EOF
import random
import sys, glob, os

def findName(t, id_f):
    id_s = t.index(id_f)
    for i in range(id_s, len(t)):
        if t[i]==':':
            return t[id_s:i]

def splitBr_difid(t, id_old, id1, id2):
    to_repl = '(' + id1 + ':0,' + id2 + ':0' + ')' + '0.5'
    return t.replace(id_old, to_repl)

folder = os.getcwd() + '/' + sys.argv[1]
files = glob.glob(folder + '*tree.txt')
for file in files:
    temp = open(file, 'r')
    t = temp.readlines()[0].split()[0]
    id_old = file.split('/')[-1].split('_c.')[0]
    id1 = id_old.split('_' + '`echo $unique_str`')[0] 
    id2 = '`echo $unique_str`' + id_old.split('_' + '`echo $unique_str`')[1]
    t_s = splitBr_difid(t, id_old, id1, id2)
    output_name = file.split('.txt')[0] + '_s.txt'
    f = open(output_name, 'w')
    f.write(t_s)
    f.close()
EOF


qsub -N s_trees -hold_jid "job_n_1_*" splitting_trees.sh

cat organise_files.txt > organise_n_1.sh

cat >> organise_n_1.sh << EOF
mkdir n_1_trees_s
find n_1_trees/ -name "*_s.txt" -exec mv {} n_1_trees_s \; &
wait
tar -zcvf n_1_trees.tar.gz n_1_trees --remove-files
tar -zcvf n_1_trees_s.tar.gz n_1_trees_s --remove-files

mkdir n_1_res 
find phy_c/ -name "job*" -exec mv {} n_1_res \; &
wait
tar -zcvf n_1_res.tar.gz n_1_res --remove-files

tar -zcvf phy_c.tar.gz phy_c --remove-files
tar -zcvf hog_aln.tar.gz hog_aln --remove-files
EOF


qsub -N o_n_1 -hold_jid s_trees organise_n_1.sh


#computing n trees - without input topology

cd phy

for f in *.phy
do
    cat ../one_tree.txt > $f.notop.sh
    cat >> $f.notop.sh << EOF
    source ../load_env
    FastTree $f > ${f}_tree_notop.txt
EOF
done

i=1
for f in *.phy.notop.sh
do
    qsub -N job_n_notop_$i -cwd $f
    ((i++))
done

cd ..

cat organise_files.txt > organise_n_notop.sh

cat >> organise_n_notop.sh << EOF
mkdir n_trees_notop
find phy/ -name "*_tree_notop.txt" -exec mv {} n_trees_notop \; &
wait
tar -zcvf n_trees_notop.tar.gz n_trees_notop --remove-files

mkdir n_notop_res 
find phy/ -name "job_n_notop*" -exec mv {} n_notop_res \; &
wait
tar -zcvf n_notop_res.tar.gz n_notop_res --remove-files
EOF

qsub -N o_n_not -hold_jid "job_n_notop_*" organise_n_notop.sh

#computing n trees - with input topology

#put input trees to the folder
cat organise_files.txt > moving_trees.sh
cat >> moving_trees.sh <<EOF

tar -zxvf n_1_trees_s.tar.gz
find n_1_trees_s/ -name "*.txt" -exec mv {} phy \; &
wait
rmdir n_1_trees_s
EOF

qsub -N move_tr -hold_jid o_n_1 moving_trees.sh

cd phy

for f in *.phy
do
    cat ../one_tree.txt > $f.top.sh
    cat >> $f.top.sh << EOF
    source ../load_env
    FastTree -intree  ${f%.phy}_c.phy_tree_s.txt $f > ${f}_tree_top.txt
EOF
done

i=1
for f in *.top.sh
do
    qsub -N job_n_top_$i -hold_jid move_tr -cwd $f
    ((i++))
done

cd ..

cat organise_files.txt > organise_n_top.sh
cat >> organise_n_top.sh << EOF
mkdir n_top_res 
find phy/ -name "job_n_top*" -exec mv {} n_top_res \; &
wait
tar -zcvf n_top_res.tar.gz n_top_res --remove-files
EOF

qsub -N o_n_top -hold_jid "job_n_top_*" organise_n_top.sh

cat organise_files.txt > organise_n.sh
cat >> organise_n.sh << EOF
tar -zcvf phy.tar.gz phy --remove-files
tar -zcvf aln.tar.gz aln --remove-files
EOF

qsub -N o_n -hold_jid "o_n_*" organise_n.sh

#get bootstrap replicates

cat bootstrap.txt > bootstrap.sh
cat >> bootstrap.sh << EOF

source ./load_env
python bootstrap.py --loc_folder='aln_c/' --num_bootstraps=$n_samples -j 1 --verbose=2
EOF


cat > bootstrap.py << EOF
#Author: Kevin Gori

import sys, os
from ruffus import transform, subdivide, cmdline, suffix, formatter, mkdir
from Bio import AlignIO
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from abc import ABCMeta, abstractmethod, abstractproperty
from locale import getpreferredencoding
import shlex
from subprocess import PIPE, Popen
import threading
from collections import defaultdict
import random, glob

IS_PY3 = sys.version_info[0] == 3
POSIX = 'posix' in sys.builtin_module_names
DEFAULT_ENCODING = getpreferredencoding() or "UTF-8"

if IS_PY3:
    from queue import Queue
else:
    from Queue import Queue

class AbstractWrapper(object):
    """
    Abstract Base Class:
    Run an external program as a subprocess with non-blocking collection of stdout.
    This class uses subprocess.Popen to handle the system call, and threading.Thread
    to allow non-blocking reads of stderr and stdout.

    Specific program wrappers inherit from this. They should implement two methods:
    1: _default_exe and 2: _set_help.

    Example:
    class SomeProgram(AbstractInternal):
        @property
        def _default_exe(self):
            return 'some_program'           # The usual name of the program binary

        def _set_help(self):
            self('--help', wait=True)       # calls 'some_program --help'
            self._help = self.get_stdout()  # puts the help output into self._help

        #### Possible Extra Requirement ####
        @property
        def _hyphen_policy(self):
            return 1       # If some_program arguments take a single leading hyphen
                           # then set this to 1 (i.e. to the number of leading hyphens)
                           # This only affects calling some_program using named arguments

    That's it!

    See blog post: http://www.zultron.com/2012/06/python-subprocess-example-running-a-background-subprocess-with-non-blocking-output-processing/
    Also, StackOverflow: https://stackoverflow.com/questions/375427/non-blocking-read-on-a-subprocess-pipe-in-python
    Keyword handling borrowed from sh: https://amoffat.github.io/sh

    """
    __metaclass__ = ABCMeta

    def __init__(self, executable=None, verbose=True):
        """
        Sets up the wrapper. A custom executable can be passed in, otherwise
        it will search the PATH.
        :param executable: Path to a custom executable to use
        :return:
        """
        exe = None
        if executable:
            exe = self._search_for_executable(executable)

        if exe is None:
            exe = self._search_for_executable(self._default_exe)
            if not exe:
                raise IOError(executable if executable else self._default_exe)

        self.exe = exe           # The wrapped executable
        self.verbose = verbose   # Controls printing of output
        self.stdout_q = Queue()  # The output-reading threads pipe output into
        self.stderr_q = Queue()  # these queues.
        self.stdout_l = list()   # Queued output gets placed in these lists,
        self.stderr_l = list()   # making it available to the caller
        self.process = None      # This holds the running process once the wrapped program is called
        self._help = None        # A place to hold a help string

    def __repr__(self):
        return '{}(executable=\'{}\')'.format(self.__class__.__name__, self.exe)

    # Private
    @abstractproperty
    def _default_exe(self):
        pass

    @property
    def _hyphen_policy(self):
        """
        Returns 'n', where 'n' is the number of hyphens prepended to commandline flags.
        Used internally when constructing command line arguments from parameters
        passed to __call__ as keywords.
        :return: 2 (as in 2 hyphens, '--')
        """
        return 2

    @abstractmethod
    def _set_help(self):
        pass

    def _log_thread(self, pipe, queue):
        """
        Start a thread logging output from pipe
        """

        # thread function to log subprocess output (LOG is a queue)
        def enqueue_output(out, q):
            for line in iter(out.readline, b''):
                q.put(line.rstrip())
            out.close()

        # start thread
        self.t = threading.Thread(target=enqueue_output,
                                  args=(pipe, queue))
        self.t.daemon = True  # thread dies with the program
        self.t.start()

    def _search_for_executable(self, executable):
        """
        Search for file give in "executable". If it is not found, we try the environment PATH.
        Returns either the absolute path to the found executable, or None if the executable
        couldn't be found.
        """
        if os.path.isfile(executable):
            return os.path.abspath(executable)
        else:
            envpath = os.getenv('PATH')
            if envpath is None:
                return
            for path in envpath.split(os.pathsep):
                exe = os.path.join(path, executable)
                if os.path.isfile(exe):
                    return os.path.abspath(exe)

    # Public
    def __call__(self, cmd=None, wait=False, **flags):
        """
        Spawns the subprocess and the threads used to monitor stdout and stderr without blocking.
        :param cmd: Pass the command line arguments as a string
        :param wait: Block until the process returns
        :param flags: Pass the commandline arguments as a dictionary. Will be appended to any content in cmd.

        :return:
        """
        # Check there is not already a process running
        if self.running():
            self.kill()

        # Wipe any stdout/stderr from previous processes
        self.stderr_l = []
        self.stdout_l = []

        # Assemble command line
        if cmd is None:
            cmd = ''
        if flags:
            cmd = ' '.join([cmd.strip(), _kwargs_to_args(flags, self._hyphen_policy)])
        self.cmd = '{} {}'.format(self.exe, cmd)

        # spawn
        self.process = Popen(shlex.split(self.cmd),
                             shell=False, stdout=PIPE, stderr=PIPE, bufsize=1, close_fds=POSIX)
        if self.verbose:
            print('Launched {} with PID {}'.format(self.exe, self.process.pid))

        # start stdout and stderr logging threads
        self._log_thread(self.process.stdout, self.stdout_q)
        self._log_thread(self.process.stderr, self.stderr_q)

        if wait:
            self.process.wait()

    @property
    def help(self):
        """
        Returns a helpful string, preferably derived from the wrapped program
        :return:
        """
        if self._help is None:
            self._set_help()
        return self._help

    def get_stderr(self, tail=None):
        """
        Returns current total output written to standard error.
        :param tail: Return this number of most-recent lines.
        :return: copy of stderr stream
        """
        while not self.stderr_q.empty():
            self.stderr_l.append(self.stderr_q.get_nowait())
        if tail is None:
            tail = len(self.stderr_l)
        return _py2_and_3_joiner('\n', self.stderr_l[:tail])

    def get_stdout(self, tail=None):
        """
        Returns current total output written to standard output.
        :param tail: Return this number of most-recent lines.
        :return: copy of stdout stream
        """
        while not self.stdout_q.empty():
            self.stdout_l.append(self.stdout_q.get_nowait())
        if tail is None:
            tail = len(self.stdout_l)
        return _py2_and_3_joiner('\n', self.stdout_l[:tail])

    def finished(self):
        """
        Check if the running process is finished. Raises an exception if no process has ever been launched.
        :return: bool
        """
        if self.process is None:
            raise ExternalProcessError('No process has been launched from this instance')
        return self.process.poll() is not None

    def kill(self):
        """
        Kill the running process (if there is one)
        :return: void
        """
        if self.running():
            if self.verbose:
                print('Killing {} with PID {}'.format(self.exe, self.process.pid))
            self.process.kill()

            # Thread *should* tidy up itself, but we do it explicitly
            if self.t.is_alive():
                self.t.join(1)

    def running(self):
        """
        True if there is a running process. False if either no process is associated with this instance,
        or if the associated process has finished.
        :return: bool
        """
        if self.process is None:
            return False
        return self.process.poll() is None



class Mafft(AbstractWrapper):
    @property
    def _default_exe(self):
        return 'mafft'

    def _set_help(self):
        self(help=True, wait=True)
        self._help = self.get_stdout()

def bootstrap(alignment):
    length = alignment.get_alignment_length()
    columns = [alignment[:,n] for n in range(length)]
    sample = _sample_wr(columns, length)
    sequences = [''.join(x) for x in zip(*sample)]
    ids = [seqrec.id for seqrec in alignment]
    alphabets = [seqrec.seq.alphabet for seqrec in alignment]
    return _assemble_msa(sequences, ids, alphabets)

def _sample_wr(population, k):
    _int = int
    _random = random.random
    n = len(population)
    return [population[_int(_random() * n)] for _ in range(k)]


def _assemble_msa(strings, ids, alphabets):
    seqs = []
    for (s, i, a) in zip(strings, ids, alphabets):
        seq = Seq(s, alphabet=a)
        seqs.append(SeqRecord(seq, id=i))
    return MultipleSeqAlignment(seqs)

def concatenate(alignments):
    """
    Concatenates a list of Bio.Align.MultipleSeqAlignment objects.
    If any sequences are missing the are padded with unknown data
    (Bio.Seq.UnknownSeq).
    Returns a single Bio.Align.MultipleSeqAlignment.
    Limitations: any annotations in the sub-alignments are lost in
    the concatenated alignment.
    """

    # Get the full set of labels (i.e. sequence ids) for all the alignments
    all_labels = set(seq.id for aln in alignments for seq in aln)

    # Make a dictionary to store info as we go along
    # (defaultdict is convenient -- asking for a missing key gives back an empty list)
    tmp = defaultdict(list)

    # Assume all alignments have same alphabet
    alphabet = alignments[0]._alphabet

    for aln in alignments:
        length = aln.get_alignment_length()

        # check if any labels are missing in the current alignment
        these_labels = set(rec.id for rec in aln)
        missing = all_labels - these_labels

        # if any are missing, create unknown data of the right length,
        # stuff the string representation into the tmp dict
        for label in missing:
            new_seq = UnknownSeq(length, alphabet=alphabet)
            tmp[label].append(str(new_seq))

        # else stuff the string representation into the tmp dict
        for rec in aln:
            tmp[rec.id].append(str(rec.seq))

    # Stitch all the substrings together using join (most efficient way),
    # and build the Biopython data structures Seq, SeqRecord and MultipleSeqAlignment
    return MultipleSeqAlignment(SeqRecord(Seq(''.join(v), alphabet=alphabet), id=k)
               for (k,v) in tmp.items())


parser = cmdline.get_argparse(description='Sequences->Alignments->Bootstraps,Trees->Bootstrapped Tree')
parser.add_argument('-lf', '--loc_folder', type=str, help='Folder with input files')
parser.add_argument('-b', '--num_bootstraps', type=int, default=100, help='Number of replicates to bootstrap')
options = parser.parse_args()

NUM_BTSTPS = options.num_bootstraps
loc = options.loc_folder
#  standard python logger which can be synchronised across concurrent Ruffus tasks
logger, logger_mutex = cmdline.setup_logging (__name__, options.log_file, options.verbose)

sequence_files = glob.glob(os.getcwd() + '/' + loc + '*.aln')

@transform(sequence_files, suffix(".aln"), ".phy")

def convert_to_phy(input_file, output_file):
    """
    PIPELINE STAGE 'convert_to_phy'
    ===============================

    Does:   Converts fasta formatted alignment to phylip-relaxed
    Input:  <filename>.aln
    Output: <filename>.phy

    """
    AlignIO.convert(input_file, "fasta", output_file, "phylip-relaxed")

@mkdir(convert_to_phy,
       formatter(),
       "{path[0]}/bootstrap")

@subdivide(convert_to_phy,
           formatter(),
           "{path[0]}/bootstrap/{basename[0]}.*.boot",
           "{path[0]}/bootstrap/{basename[0]}")

def generate_bootstraps(input_file, output_files, output_stub):
    """
    PIPELINE STAGE 'generate_bootstraps'
    ====================================

    Does:   Generates bootstrap replicates
    Input:  <filename>.phy = Phylip format alignments from "convert_to_phy"
    Output: <filename>.<number>.boot

    Description:
    Take the output of stage "convert_to_phy" and generate 10 bootstrap replicates.
    Write these to <filename>.<number>.boot, in phylip-relaxed format.
    """
    aln = AlignIO.read(input_file, 'phylip-relaxed')
    for i in range(NUM_BTSTPS):
        boot = bootstrap(aln)

        # {:0>3} is shorthand for "make the string 3 chars long by padding
        # with 0 to the left -- so 1 -> 001, 23 -> 023, etc.
        output_filename = '{0}.{1:0>{2}}.boot'.format(output_stub, i+1, len(str(NUM_BTSTPS)))

        with open(output_filename, 'w') as handle:
            AlignIO.write(boot, handle, 'fasta')

def summ(infiles, summfile):
    with open(summfile, 'w') as outp:
        for f in infiles:
            with open(f) as inp:
                out.write(f.read())

cmdline.run(options)
EOF

qsub -N bs -cwd bootstrap.sh

#convert to phy
cat organise_files.txt > boot_phy.sh
cat >> boot_phy.sh << EOF

source ./load_env
python boot_phy.py aln_c/bootstrap/
EOF

cat > boot_phy.py << EOF
import sys
import requests, random, subprocess
from xml.dom import minidom
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Phylo.Applications import PhymlCommandline 
import glob
import os

loc = sys.argv[1]
boot_files = glob.glob(os.getcwd() + '/' + loc + '*.boot')

for file in boot_files:
    msa_phy = loc + file.split('/')[-1] + '.phy' #MSA in phylip format
    AlignIO.convert(file, "fasta", msa_phy, "phylip-relaxed")
EOF

qsub -N bs_phy -hold_jid bs boot_phy.sh

#organise n_1_b_phy

cat organise_files.txt > organise_br_n_1_phy.sh
cat >> organise_br_n_1_phy.sh << EOF

mkdir bootstrap_phy
find aln_c/bootstrap/ -name "*.phy" -exec mv {} bootstrap_phy \; &
wait

tar -zcvf bootstrap_phy.tar.gz bootstrap_phy --remove-files
EOF

qsub -N ob -hold_jid bs_phy organise_br_n_1_phy.sh 

#split replicates

cat organise_files.txt > split_aln.sh
cat >> split_aln.sh << EOF

source ./load_env
python split_aln.py cuts.txt aln_c/bootstrap/

cat > split_aln.py << EOF
import sys
from Bio import AlignIO
import glob, os
import logging
logger = logging.getLogger("split_aln")
logging.basicConfig(level=logging.DEBUG, format="%(asctime)-15s - %(levelname)s: %(message)s")

def find_key(d, s1, s2):
    for key in d.keys():
        if s1 in key[0] and s2 in key[1]:
            return key

def find_key_full(d, s):
        for key in d.keys():
                if s in key:
                        return key

def read_msa(loc):
    msa = dict()
    tmp = open(loc, 'r')
    lines = tmp.readlines()
    key = lines[0].split()[0][1:]
    value = ''
    for i in range(1, len(lines)):
        if lines[i].split()[0][0] == '>':
            msa[key] = value
            key = lines[i].split()[0][1:]
            value = ''
        else:
            value += lines[i].split()[0]
    msa[key] = value
    return msa

def split_sequence(msa, id_full, id1, id2, n):
    logger.info('splitting sequence {} at {} into {} and {}'.format(id_full, n, id1, id2))
    msa[id1] = msa[id_full][0:n] + '-' * (len(msa[id_full]) - n)
    msa[id2] = '-' * n + msa[id_full][n:]
    msa.pop(id_full)
    return msa

def write_dictionary(d, out):
    f = open(out, 'w')
    for key in d.keys():
        f.write('>' + key + '\n')
        f.write(d[key] + '\n')
    f.close()
    return

cuts = dict()
fn = os.path.join(os.getcwd(), sys.argv[1])
logger.info('loading cuts from '+fn)
with open(fn, 'r') as cf:
    for line in cf:
        temp  = line.split()
        cuts[temp[1], temp[2]] = int(temp[3])
        logger.debug('cuts: {} / {}: {}'.format(temp[1], temp[2], temp[3])) 


boot_dir = os.getcwd() + '/' + sys.argv[2]
files = glob.glob(boot_dir + "*.boot")

for i in range(len(files)):
    i1 = files[i].split('/')[-1].split('_c')[0].split('_' + '`echo $unique_str`')[0]
    i2 = '`echo $unique_str`' + files[i].split('/')[-1].split('_c')[0].split('_'  + '`echo $unique_str`')[1]
    c = True
    if find_key(cuts, i1, i2) != None:
        id1, id2 = find_key(cuts, i1, i2)       
    elif find_key(cuts, i2, i1) != None:
        id1, id2 = find_key(cuts, i2, i1)
    else:
        c = False
    logger.debug('checking {}, {} for cuts: {}'.format(i1, i2, c)) 
    if c:   
        n = cuts[id1, id2]
        msa = read_msa(files[i])
        id_full = find_key_full(msa, id1)
        if len(msa[id_full]) > n:
            msa_s = split_sequence(msa, id_full, id1, id2, n)
            out = files[i] + '_s.fa'
            write_dictionary(msa_s, out)
            out_phy = files[i] + '_s.phy'
            AlignIO.convert(out, "fasta", out_phy, "phylip-relaxed")

EOF


qsub -N bs_s -hold_jid ob split_aln.sh

#organise the rest
cat organise_files.txt > organise_b_rest.sh
cat >> organise_b_rest.sh << EOF

cd aln_c
mkdir bootstrap_s_aln
mkdir bootstrap_s_phy
find bootstrap/ -name "*_s.fa" -exec mv {} bootstrap_s_aln \; &
wait
find bootstrap/ -name "*_s.phy" -exec mv {} bootstrap_s_phy \; &
wait
mv bootstrap_s_aln ../
mv bootstrap_s_phy ../

mv bootstrap bootstrap_aln
mv bootstrap_aln ../

cd ../

tar -zcvf bootstrap_aln.tar.gz bootstrap_aln --remove-files
tar -zcvf bootstrap_s_aln.tar.gz bootstrap_s_aln --remove-files
tar -zcvf bootstrap_s_phy.tar.gz bootstrap_s_phy --remove-files
tar -zcvf aln_c.tar.gz aln_c --remove-files
EOF

qsub -N obr -hold_jid bs_s organise_b_rest.sh

#computing bootstrap n-1 trees

cat organise_files.txt > b_n_1_prep.sh
cat >> b_n_1_prep.sh << EOF

tar -zxvf bootstrap_phy.tar.gz 
rm bootstrap_phy.tar.gz 

cd bootstrap_phy

for f in *.phy; do echo \${f/_c.*.phy/};done | sort| uniq > foldernames.txt

for line in \`cat foldernames.txt\`
do  
    mkdir \$line
    find . -maxdepth 1 -name "\$line*.phy" -exec mv {} \$line/ \; &
    wait
done

for line in \`cat foldernames.txt\`
do
    cat ../bootstrap_trees.txt > \$line.sh
    cat >> \$line.sh << EOA
    source load_env
    cd bootstrap_phy/\$line
    for f in *.phy; do 
        FastTree \\\$f > \\\$f.tree.txt
    done
EOA
done
EOF

qsub -sync y -N b_n_1_prep -hold_jid "ob","obr" b_n_1_prep.sh

i=1
for f in bootstrap_phy/*.sh
do
    qsub -N job_b_n_1_$i -cwd $f
    ((i++))
done


cat organise_files.txt > splitting_trees_b.sh
cat >> splitting_trees_b.sh << EOF

mkdir -p n_1_b_trees
find bootstrap_phy/ -type f -name "*tree.txt" -exec mv {} n_1_b_trees \; &
wait

source ./load_env
python splitting_trees.py n_1_b_trees/
EOF

qsub -N s_trees_b -hold_jid "job_b_n_1_*" splitting_trees_b.sh

cat organise_files.txt > organise_b_n_1.sh
cat >> organise_b_n_1.sh << EOF

mkdir n_1_b_trees_s
find n_1_b_trees/ -type f -name "*_s.txt" -exec mv {} n_1_b_trees_s \; &
wait
tar -zcvf n_1_b_trees.tar.gz n_1_b_trees --remove-files
tar -zcvf n_1_b_trees_s.tar.gz n_1_b_trees_s --remove-files

mkdir n_1_b_res
find . -name "job_b_n_1*" -exec mv {} n_1_b_res \; 2>/dev/null &
wait
tar -zcvf n_1_b_res.tar.gz n_1_b_res --remove-files

tar -zcvf bootstrap_phy.tar.gz bootstrap_phy --remove-files
EOF

qsub -N o_b_n_1 -hold_jid s_trees_b organise_b_n_1.sh

#computing bootstrap n trees

cat organise_files.txt > b_n_prep.sh 
cat >> b_n_prep.sh << EOF

tar -zxvf bootstrap_s_phy.tar.gz 
rm bootstrap_s_phy.tar.gz 
tar -zxvf n_1_b_trees_s.tar.gz
find n_1_b_trees_s/ -name "*.txt" -exec mv {} bootstrap_s_phy \; &
wait
rm -rf n_1_b_trees_s
cd bootstrap_s_phy

for f in *.phy; do echo \${f/_c.*.phy/};done | sort| uniq > foldernames.txt

for line in \`cat foldernames.txt\`
do  
    mkdir \$line
    find . -maxdepth 1 -name "\$line*.phy" -exec mv {} \$line/ \; &
    wait
    find . -maxdepth 1 -name "\$line*.txt" -exec mv {} \$line/ \; &
    wait
done

for line in \`cat foldernames.txt\`
do
    cat ../bootstrap_trees.txt > \$line.notop.sh
    cat >> \$line.notop.sh << EOA
    source load_env
    cd  bootstrap_s_phy/\$line
    for f in *.phy; do 
        FastTree \\\$f > \\\$f.tree_notop.txt
    done
EOA
done

for line in \`cat foldernames.txt\`
do
    cat ../bootstrap_trees.txt > \$line.top.sh
    cat >> \$line.top.sh << EOB
    source load_env
    cd bootstrap_s_phy/\$line
    for f in *.phy; do 
        FastTree -intree \\\${f/_s.phy/}.phy.tree_s.txt \\\$f > \\\$f.tree_top.txt
    done
EOB
done
EOF

qsub -sync y -N b_n_prep -hold_jid obr,o_b_n_1, b_n_prep.sh


i=1
for f in bootstrap_s_phy/*.notop.sh
do
    qsub -N job_b_n_notop_$i -cwd $f
    ((i++))
done

i=1
for f in bootstrap_s_phy/*.top.sh
do
    qsub -N job_b_n_top_$i -cwd $f
    ((i++))
done

cat organise_files.txt > organise_b_n.sh
cat >> organise_b_n.sh << EOF

mkdir n_b_notop_res
find . -name "job_b_n_notop*" -exec mv {} n_b_notop_res \; & 
wait
tar -zcvf n_b_notop_res.tar.gz n_b_notop_res --remove-files

mkdir n_b_top_res
find . -name "job_b_n_top*" -exec mv {} n_b_top_res \; & 
wait
tar -zcvf n_b_top_res.tar.gz n_b_top_res --remove-files

tar -zcvf bootstrap_s_phy.tar.gz bootstrap_s_phy --remove-files
EOF

qsub -N obn -hold_jid "job_b_n_notop_*","job_b_n_top_*" organise_b_n.sh



#LRT
cat organise_files.txt > prepare_lrt.sh
cat >> prepare_lrt.sh << EOF

tar -zxvf n_1_res.tar.gz
tar -zxvf n_notop_res.tar.gz
tar -zxvf n_top_res.tar.gz
tar -zxvf n_1_b_res.tar.gz
tar -zxvf n_b_notop_res.tar.gz
tar -zxvf n_b_top_res.tar.gz

mkdir lrt_summaries
EOF

qsub -N prep_lrt -hold_jid o_n_1,o_n_not,o_n_top,o_b_n_1,obn prepare_lrt.sh

cat organise_files.txt > do_lrt.sh
cat >> do_lrt.sh << EOF

source ./load_env
python lrt.py n_1_res/ n_notop_res/ n_top_res/ n_1_b_res/ n_b_notop_res/ n_b_top_res/
EOF

cat > lrt.py << EOF

import glob, os, sys

def grepLoglk_fasttree_1(file):
    f = open(file, 'r')
    lines = f.readlines()
    for line in lines:
        if 'Optimize all lengths:' in line:
            return float(line.split()[-3])
    return False

def grepLoglk_fasttree_m(file):
    tmp = list()
    f = open(file, 'r')
    lines = f.readlines()
    for line in lines:
        if 'Optimize all lengths:' in line:
            tmp.append(float(line.split()[-3]))
    return tmp  

def findElement(files, s1, s2):
    for i in range(len(files)):
        if s1 in files[i] and s2 in files[i]:
            return files[i]

def findElements(files, s1, s2):
    temp = list()
    for i in range(len(files)):
        if s1 in files[i] and s2 in files[i]:
            temp.append(files[i])
    return temp

def readlkboot_ft(files, s1, s2):
    ls = list()
    fs = findElements(files, s1, s2)
    fs.sort() 
    for f in fs:
        ls = ls + grepLoglk_fasttree_m(f)
    return ls

def rename_files(path):
    for filename in glob.glob(path + '*.o*'):
        if filename.split()[-1].split('/')[-1].startswith('job'):
            f = open(filename, 'r')
            lines = f.readlines()
            new_name = filename.split()[-1].split('job')[0] + lines[1].split()[1].split('_c')[0] + '.res.txt'
            os.rename(filename, new_name)


loc_n_1 = os.getcwd() + '/' + sys.argv[1]
loc_n_notop = os.getcwd() + '/' + sys.argv[2]
loc_n_top = os.getcwd() + '/' + sys.argv[3]
loc_n_1_b = os.getcwd() + '/' + sys.argv[4]
loc_n_b_notop = os.getcwd() + '/' + sys.argv[5] 
loc_n_b_top = os.getcwd() + '/' + sys.argv[6]

rename_files(loc_n_1)
rename_files(loc_n_notop)
rename_files(loc_n_top)
rename_files(loc_n_1_b)
rename_files(loc_n_b_notop)
rename_files(loc_n_b_top)

files_n_1 = glob.glob(loc_n_1 + '*.res.txt')
files_n_notop = glob.glob(loc_n_notop + '*.res.txt')
files_n_top = glob.glob(loc_n_top + '*.res.txt')
files_n_1_b = glob.glob(loc_n_1_b + '*.res.txt')
files_n_b_notop = glob.glob(loc_n_b_notop + '*.res.txt')
files_n_b_top = glob.glob(loc_n_b_top + '*.res.txt')

test_st = dict()
distr = dict()
lrt = dict()

miss_lk = list()
for file in files_n_1:
    key = [file.split('/')[-1].split('_' + '`echo $unique_str`')[0], '`echo $unique_str`' + file.split('/')[-1].split('_' + '`echo $unique_str`')[1].split('.res')[0]]
    v1 = grepLoglk_fasttree_1(file)
    file21 = findElement(files_n_top, key[0], key[1])
    v21 = grepLoglk_fasttree_1(file21)
    file22 = findElement(files_n_notop, key[0], key[1])
    v22 = grepLoglk_fasttree_1(file22)
    b_temp = True
    if v1 == False:
        miss_lk.append((key[0], key[1], 'n-1'))
        b_temp = False
    if v21 == False:
        miss_lk.append((key[0], key[1], 'n top'))
        b_temp = False
    if v22 == False:
        miss_lk.append((key[0], key[1], 'n no top'))
        b_temp = False
    if b_temp == True:
        v2 = max(v21, v22)
        #T= -2ln(H0) + 2ln(H1) = -2v1 + 2v2
        test_st[key[0], key[1]] = [v1, v21, v22, -2*v1 + 2*v2]
    boot1 = readlkboot_ft(files_n_1_b, key[0], key[1])
    boot21 = readlkboot_ft(files_n_b_top, key[0], key[1])
    boot22 = readlkboot_ft(files_n_b_notop, key[0], key[1])
    if len(boot1) < `echo $n_samples`:
        miss_lk.append((key[0], key[1], 'boot n-1'))
    if len(boot21) < `echo $n_samples`:
        miss_lk.append((key[0], key[1], 'boot n top'))
    if len(boot22) < `echo $n_samples`:
        miss_lk.append((key[0], key[1], 'boot n no top'))
    if len(boot1) == `echo $n_samples` and len(boot21) == `echo $n_samples` and len(boot22) == `echo $n_samples` and b_temp == True:
        test_st_b = list()
        for i in range(len(boot1)):
            test_st_b.append(-2*boot1[i] + 2*max(boot21[i], boot22[i]))
        distr[key[0], key[1]] = [boot1, boot21, boot22, test_st_b]
        p = (sum(i >= test_st[key[0], key[1]][3] for i in test_st_b)+1)/(float(len(test_st_b))+1)
        lrt[key[0], key[1]] = [test_st[key[0], key[1]][3], p]
        file_name = key[0] + '_' + key[1] + '_lrt.txt'
        f = open(os.getcwd() + '/lrt_summaries/'+ file_name, 'w')
        f.write('Log-likelihood n-1; Log-likelihood n without input topology; Log-likelihood n with input topology; Test statistic \n')
        f.write(str(v1) + '\t' + str(v22) + '\t' + str(v21) + '\t' + str(test_st[key[0], key[1]][3]) + '\n')
        f.write('Bootstrap replicates \n')
        for i in range(len(boot1)):
            f.write(str(boot1[i]) + '\t' + str(boot22[i]) + '\t' + str(boot21[i]) + '\t' + str(test_st_b[i]) + '\n')
        f.close()

f = open(os.getcwd() + '/' + 'lrt_summary.txt', 'w')
f.write('Gene1; Gene2; Test statistic; p-value \n')
for key in lrt.keys():
    f.write(key[0] + '\t' + key[1] + '\t' + str(lrt[key[0], key[1]][0])  + '\t' + str(lrt[key[0], key[1]][1]) + '\n')
f.close()

f = open(os.getcwd() + '/' + 'missing_lk.txt', 'w')
f.write('Gene1; Gene2; Where \n')
for i in miss_lk:
    f.write(i[0] + '\t' + i[1] + '\t' + i[2] + '\n')
f.close()

EOF

qsub -N lrt -hold_jid prep_lrt do_lrt.sh

cat organise_files.txt > organise_lrt.sh
cat >> organise_lrt.sh << EOF

tar -zcvf lrt_summaries.tar.gz lrt_summaries --remove-files

tar -zcvf n_1_res.tar.gz n_1_res --remove-files
tar -zcvf n_notop_res.tar.gz n_notop_res --remove-files
tar -zcvf n_top_res.tar.gz n_top_res --remove-files
tar -zcvf n_1_b_res.tar.gz n_1_b_res --remove-files
tar -zcvf n_b_notop_res.tar.gz n_b_notop_res --remove-files
tar -zcvf n_b_top_res.tar.gz n_b_top_res --remove-files
EOF

qsub -N o_lrt -hold_jid lrt organise_lrt.sh

cat > collapse.py << EOF
#!/usr/bin/env python
# Author: Kevin Gori
# Date: 18 Aug 2014
from __future__ import print_function
import dendropy as dpy
import os
import sys

class SupportValueError(Exception):
    pass

def collapse(tree, threshold, keep_lengths=True, support_key=None, length_threshold=0):
    t = dpy.Tree(tree)
    to_collapse = []
    for node in t.postorder_node_iter():
        if node.is_leaf():
            if node.edge_length < length_threshold:
                node.edge_length = 0
            continue
        if node is t.seed_node:
            continue
        try:
            if support_key:
                support = float(node.annotations.get_value(support_key))
                node.label = support
            else:
                support = float(node.label)
        except TypeError as e:
            support=2.0
            #raise SupportValueError('Inner node with length {} has no support value'.format(node.edge_length), e)
        #except ValueError, e:
            #raise SupportValueError('Inner node with length {} has a non-numeric support value {}'.format(node.edge_length), e)
        if support < threshold:
            to_collapse.append(node.edge)
        elif node.edge_length < length_threshold:
            to_collapse.append(node.edge)

    for edge in to_collapse:
        if keep_lengths:
            for child in edge.head_node.child_nodes():
                child.edge.length += edge.length
        edge.collapse()
    return t

def parse_args():
    """ Parses command line arguments """
    import argparse
    help_msgs = dict(tree='REQUIRED: path to the tree file',
                     threshold='OPTIONAL: support value threshold - node with support below this value will be collapsed. DEFAULT=%(default)s',
                     length_threshold='OPTIONAL: edge length threshold - edges shorter than this value will be collapsed. DEFAULT=%(default)f',
                     format='OPTIONAL: file format of the input tree - can be newick or nexus. DEFAULT=%(default)s',
                     keep_lengths='OPTIONAL: if set, add the length of the branch being removed to the child branches, so that root-to-tip distances are preserved. DEFAULT=%(default)s',
                     outfile='OPTIONAL: write the result to disk. DEFAULT=print to terminal',
                     support_key=('OPTIONAL: the key that retrieves the support value, if the support value'
                                  ' is given in some extended newick format (i.e. as key-value pair inside a comment)'),
                     output_format=('OPTIONAL: file format of output tree - newick or nexus. DEFAULT: same as input.'
                                    ' CAVEAT: nexus format may be incompatible with FigTree'))
    parser = argparse.ArgumentParser()
    parser.add_argument('tree', type=str, help=help_msgs['tree'])
    parser.add_argument('--threshold', type=float, default=0.5, help=help_msgs['threshold'])
    parser.add_argument('--length_threshold', type=float, default=0.0001, help=help_msgs['length_threshold'])
    parser.add_argument('--format', type=str, default='newick',
                        choices=['nexus', 'newick'], help=help_msgs['format'])
    parser.add_argument('--output_format', type=str, help=help_msgs['output_format'])
    parser.add_argument('--keep_lengths', action='store_true', help=help_msgs['keep_lengths'])
    parser.add_argument('--outfile', type=str, help=help_msgs['outfile'])
    parser.add_argument('--test', action='store_true', help=argparse.SUPPRESS) # just used as a quick debug test - not relevant for user
    parser.add_argument('--support_key', type=str, help=help_msgs['support_key'])
    return parser.parse_args()

def read_tree_file(tree_file, tree_format):
    filename = tree_file
    if not os.path.exists(filename):
        raise IOError('File not found: {}'.format(filename))
    try:
        t = dpy.Tree.get_from_path(filename, tree_format, extract_comment_metadata=True)
    except dpy.utility.error.DataParseError:
        file_format, = ({'nexus', 'newick'} - {tree_format})
        t = dpy.Tree.get_from_path(filename, file_format, extract_comment_metadata=True)
    return t

def test():
    n = ('((A:1,(B:0.5,C:0.5)80:0.1)30:1,'
         '(D:1.2,(E:0.7,F:0.7)20:0.5)60:0.8)100:0;')
    t = dpy.Tree.get_from_string(n, 'newick')
    t.print_plot(plot_metric='length')
    fh = sys.stdout
    zero = collapse(t, 25)
    zero.print_plot(plot_metric='length')
    zero.write(fh, 'newick', suppress_rooting=True)
    one = collapse(t, 50)
    one.print_plot(plot_metric='length')
    one.write(fh, 'newick', suppress_rooting=True)
    two = collapse(t, 65)
    two.print_plot(plot_metric='length')
    two.write(fh, 'newick', suppress_rooting=True)
    three = collapse(t, 65, False)
    three.print_plot(plot_metric='length')
    three.write(fh, 'newick', suppress_rooting=True)
    four = collapse(t, 25, True, length_threshold=0.11)
    four.print_plot(plot_metric='length')
    four.write(fh, 'newick', suppress_rooting=True)

def get_file_handle(outfile):
    if outfile == 'stdout' or outfile is None:
        return sys.stdout
    else:
        return open(outfile, 'w')

def main():
    args = parse_args()
    if args.output_format is None:
        args.output_format = args.format
    if args.test:
        return test()
    tree = read_tree_file(args.tree, args.format)
    collapsed_tree = collapse(tree, args.threshold, args.keep_lengths, args.support_key, args.length_threshold)
    f = get_file_handle(args.outfile)
    collapsed_tree.write(file = f, schema="newick")
    

if __name__ == '__main__':
    sys.exit(main())
EOF

cat organise_files.txt > collapse.sh
cat >> collapse.sh << EOF

mkdir collapsed_$col_thr
tar -zxvf n_trees_notop.tar.gz
source ./load_env

cd n_trees_notop/
for f in *phy_tree_notop.txt
do
    python ../collapse.py \$f --threshold 0.$col_thr > ../collapsed_$col_thr/\$f
done
EOF

qsub -N coll -hold_jid o_n_not collapse.sh

cat > check_sisters.py << EOF
import os
from ete3 import Tree
import sys

def isSister(t, id1, id2):
    sisters = t.search_nodes(name=id2)[0].get_sisters()
    for i in range(len(sisters)):
        if id1 in sisters[i]:
            return True
    return False

#get a list of files and genes
files = list()
genes = list()
for file in os.listdir(os.getcwd() + '/n_trees_notop/'):
    if file.endswith('phy_tree_notop.txt'):
        files.append(file)
        genes.append((file.split('.')[0].split('_' + '`echo $unique_str`')[0], '`echo $unique_str`' + file.split('.')[0].split('_' + '`echo $unique_str`')[1]))

out = os.getcwd() + '/collapsing_results.txt'

res = dict()

for i in range(len(genes)):
    t_before = Tree(os.getcwd() + '/n_trees_notop/' + files[i])
    fragment1 = [f for f in t_before.get_leaf_names() if genes[i][0] in f]
    fragment2 = [f for f in t_before.get_leaf_names() if genes[i][1] in f]
    if len(fragment1)!=1 or len(fragment2)!=1:
        print('Problem: ' + files[i])

    sister_before = isSister(t_before, fragment1[0], fragment2[0]) and isSister(t_before, fragment2[0], fragment1[0])
    t_after = Tree(os.getcwd() + '/collapsed_' + '`echo $col_thr`' + '/' + files[i])
    sister_after = isSister(t_after, fragment1[0], fragment2[0]) and isSister(t_after, fragment2[0], fragment1[0])

    res[fragment1[0], fragment2[0]] = [sister_before, sister_after]

out_file = open(out, 'w')
out_file.write('gene1; gene2; before collapsing; after collapsing with threshold 0.' + '`echo $col_thr`' + '\n')
for key in res:
    out_file.write(key[0] + '\t' + key[1] + '\t' + str(res[key][0]) + '\t' + str(res[key][1]) + '\n')
out_file.close()
EOF

cat organise_files.txt > check_sisters.sh
cat >> check_sisters.sh << EOF

source ./load_env
python check_sisters.py

rm -rf n_trees_notop/
tar -zcvf collapsed_$col_thr.tar.gz collapsed_$col_thr
rm -rf collapsed_$col_thr
EOF

qsub -N check_sis -hold_jid coll check_sisters.sh
 
cat > predict.py << EOF
#http://stackoverflow.com/questions/21646703/grouping-elements-from-2-tuples-recursively

import sys, os
import logging

logger = logging.getLogger(__name__)

def find_clusters( tuples ):
    # clusterlist contains at each position either a set
    # representing an actual cluster, or an int referring
    # to another cluster that has eaten this one here.
    # the cluster id is its position within this list
    clusterlist=[]
    # clustermap maps an element to the id of the containing
    # cluster within clusterlist
    clustermap = {}

    # here we find the current cluster id for elem, by following the
    # chain within clusterlist, and replace that entire chain
    # with the new cluster id n.   We return the old cluster id.
    def set_cluster_id( elem, n ):
        if elem not in clustermap:
            return None
        k = clustermap[elem]
        # clusters may be replaced by references to other clusters,
        # we follow that chain
        while k < n and isinstance( clusterlist[k], int ):
            k1 = clusterlist[k]
            # this is optional, we make the chain shorter
            # by making this entry point directly to the current cluster
            clusterlist[k] = n
            k = k1
        return k

    for t in tuples:
        # for each tuple we create a new cluster
        thiscluster = set(t)
        n = len( clusterlist ) # the id of thiscluster
        for x in t:
            # we absorb existing clusters into the new one
            # if there is overlap
            k = set_cluster_id(x, n)
            if k is not None and k != n:
                thiscluster.update( clusterlist[k] )
                # we replace the existing cluster
                # with a reference to the new one
                clusterlist[k] = n 
            clustermap[x] = n
        clusterlist.append(thiscluster)

    return [ tuple(x) for x in clusterlist if isinstance( x, set ) ]

def check_cl_ovlp(cluster, positions):
    for i in range(len(cluster)-1):
        r1 = list(range(positions[cluster[i]][0], positions[cluster[i]][1] + 1))
        for j in range(i+1, len(cluster)):
            r2 = list(range(positions[cluster[j]][0], positions[cluster[j]][1] + 1))
            intrsc = list(set(r1) & set(r2))
            if len(intrsc)/float(len(r1)) >= `echo $ovlp_len` or len(intrsc)/float(len(r2)) >= `echo $ovlp_len`:
                return True 
    return False


collapsing_res = dict()
tmp = open(os.getcwd() + '/collapsing_results.txt', 'r')
lines = tmp.readlines()
for i in range(1,len(lines)):
    collapsing_res[lines[i].split()[0], lines[i].split()[1]] = lines[i].split()[3]
tmp.close()

lrt_res = dict()
tmp = open(os.getcwd() + '/lrt_summary.txt', 'r')
lines = tmp.readlines()
for i in range(1,len(lines)):
    lrt_res[lines[i].split()[0], lines[i].split()[1]] = float(lines[i].split()[3])
tmp.close()


res_temp = list()

for key in collapsing_res:
    if collapsing_res[key] == 'True' and lrt_res[key] > float('`echo $lrt_sign`'):
        res_temp.append(key)


ambiguities = [[0]*len(res_temp) for x in range(len(res_temp))]

for i in range(len(res_temp)-1):
    for j in range(i+1, len(res_temp)):
        if res_temp[i][0] in res_temp[j] or res_temp[i][1] in res_temp[j]:
            ambiguities[i][j] = 1
            ambiguities[j][i] = 1

unamb = list()
amb_tmp = list()
for i in range(len(res_temp)):
    if sum(ambiguities[i]) == 0:
        unamb.append([res_temp[i][0], res_temp[i][1]])
    else:
        amb_tmp.append([res_temp[i][0], res_temp[i][1]])


positions = dict()
tmp = open(os.getcwd() + '/alignment_positions.txt', 'r')
lines = tmp.readlines()
for i in range(len(lines)):
    positions[lines[i].split()[1]] = [int(lines[i].split()[3]), int(lines[i].split()[4])]
    positions[lines[i].split()[2]] = [int(lines[i].split()[5]), int(lines[i].split()[6])]    
tmp.close()

clusters = find_clusters(amb_tmp)

nonovlp_cl = list()
for cluster in clusters:
    if not check_cl_ovlp(cluster, positions):
        nonovlp_cl.append(cluster)

nonovlp_cl = [item for sublist in nonovlp_cl for item in sublist]

amb = list()
for pair in amb_tmp:
    if pair[0] in nonovlp_cl or pair[1] in nonovlp_cl:
        unamb.append(pair)
    else:
        amb.append(pair)

unamb_out = open(os.getcwd() + '/predictions_unambiguous.txt', 'w')
unamb_out.write('Collapsing with threshold 0.' + '`echo $col_thr`' + ', LRT significance level ' + '`echo $lrt_sign`' + '\n')
for i in range(len(unamb)):
    unamb_out.write(unamb[i][0] + '\t' + unamb[i][1] + '\n')
unamb_out.close()

amb_out = open(os.getcwd() + '/predictions_ambiguous.txt', 'w')
amb_out.write('Collapsing with threshold 0.' + '`echo $col_thr`' + ', LRT significance level ' + '`echo $lrt_sign`' + '\n')
for i in range(len(amb)):
    amb_out.write(amb[i][0] + '\t' + amb[i][1] + '\n')
amb_out.close()

if len('$gff') > 0:
    import gff
    try:
        merger = gff.GFFOperations('$gff')
    except Exception:
        logger.exception('input gff file is invalid')
        return

    for pair in unamb:
        try:
            merger.merge_genes(pair)
        except MergeError as e:
            logger.warning('cannot merge {} gene pair in gff file: {}', pair, str(e))
        except KeyError as e:
            logger.error('invalid gene for gff file: {}', str(e))
    outfile = '$gff_out' if len('$gff_out') > 0 else '$gff'+'.merged'
    merger.write(outfile)

print('all finished. Bye!')
EOF

cat organise_files.txt > predict.sh
cat >> predict.sh << EOF
source ./load_env
python predict.py 
EOF

qsub -N predict -hold_jid lrt,check_sis predict.sh

