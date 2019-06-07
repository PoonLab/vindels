from xml.etree.ElementTree import ElementTree as Tree
from xml.etree.ElementTree import Element as Node
from xml.etree.ElementTree import tostring
import sys
import os

import argparse

def convert_fasta (handle):
    result = []
    h = None
    sequence = ''
    for line in handle:
        if line.startswith('$'): # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h,sequence])
                sequence = ''   # reset
            h = line.strip('>#\n')
        else:
            sequence += line.strip('\n')
    if h is not None:
        result.append([h,sequence]) # handle last entry
    return result

def parse(template_file, fasta_file, stem, outfile, time_unit='days', nreps=1):
    # import sequences from NEXUS file
    handle = open(fasta_file, 'rU')
    alignment = convert_fasta(handle.readlines())
    handle.close()
    
    template = Tree()

    t_root = template.parse(template_file)

    # reset TAXA and ALIGNMENT blocks
    t_taxa = template.findall('taxa')[0]
    id_tx = t_taxa.get("id")
    # t_taxa._children = []  # deprecated since Python 2.7
    t_taxa.clear()

    t_taxa.set('id',id_tx)

    t_aln = template.find('alignment')
    alntype = t_aln.get("dataType")
    id_aln = t_aln.get("id")
    # t_aln._children = []
    t_aln.clear()

    t_aln.set('dataType', alntype)
    t_aln.set('id',id_aln)




    #t_aln = alntype
    all_dates = []

    # update blocks
    for id, seq in alignment:
        # TAXON
        date_val = float(id.split('_')[-1])
        if date_val not in all_dates:
            all_dates.append(date_val)
            
        try:
            date = Node('date', {'units':time_unit, 'direction':'forwards', 'value':str(date_val)})
        except:
            raise
        
        taxon = Node('taxon', {'id':id})
        taxon.append(date)
        t_taxa.append(taxon)
    
        # SEQUENCE
        seqtag = Node('sequence', {})
        staxon = Node('taxon', {'idref':id})
        staxon.tail = '\n\t\t\t'+seq.upper()  # mimic formatting in BEAST XML
        seqtag.append(staxon)
        t_aln.append(seqtag)

    #print 'all_dates = ' + str(all_dates)
    
    if len(all_dates) < 2:
        print ('ERROR: only one sample date in file')

        sys.exit()
    
    # find log elements
    t_mcmc = template.find('mcmc')
    log_elements = t_mcmc.findall('log')
    log_tree_elements = t_mcmc.findall('logTree')
    if log_tree_elements is None or not len(log_tree_elements):
        print ("ERROR: Failed to find tree log, exiting without write!")
        sys.exit()    



    # revise log file paths and prior settings
    if nreps > 1:
        for rep in range(nreps):
            t_mcmc.set("operatorAnalysis", stem + ".%d.ops" % rep)
            of = outfile.replace('.xml', '.%d.xml' % rep)
            for log in log_elements:
                if log.get('id') == 'fileLog':
                    log.set('fileName', stem+'.%d.log'%rep)
                elif log.get("id") == "stateLog":
                    log.set('fileName', stem+'.%d.states.log'%rep)
                    ancestral_trait = log.find("ancestralTrait")
                    ancestral_trait.set("name", os.path.basename(stem))
                    break

            for log_tree_element in log_tree_elements:
                if log_tree_element.get('id') == "treeFileLog":
                    log_tree_element.set('fileName', stem+'.%d.time.trees' % rep)
                elif log_tree_element.get('id') == "substTreeFileLog":
                    log_tree_element.set('fileName', stem+'.%d.subst.trees' % rep)
            template.write(of)
    else:
        t_mcmc.set("operatorAnalysis", stem + ".ops")
        for log in log_elements:
            if log.get('id') == 'fileLog':
                log.set('fileName', stem+'.log')
            elif log.get("id") == "stateLog":
                log.set('fileName', stem+'.states.log')
                #ancestralTrait name="patient3794.realn.841_1140" traitName="states">
                ancestral_trait = log.find("ancestralTrait")
                ancestral_trait.set("name", os.path.basename(stem))
                break
        # revise tree file path
        for log_tree_element in log_tree_elements:
            if log_tree_element.get('id') == "treeFileLog":
                log_tree_element.set('fileName', stem+'.time.trees')
            elif log_tree_element.get('id') == "substTreeFileLog":
                log_tree_element.set('fileName', stem+'.subst.trees')

        template.write(outfile)


def main():
    parser = argparse.ArgumentParser(description="Insert contents of FASTA file into BEAST XML template")
    parser.add_argument('fasta', help='<input> FASTA containing sequence data.  Sequence headers must be'
                        'underscore delimited with tip date in last position.')
    parser.add_argument('template', help='<input> BEAST XML to use as template')
    parser.add_argument('-stem', default='', help='<optional> file path to write log and tree files to')
    parser.add_argument('out', help='<output> file to write new XML')
    parser.add_argument('-time', choices=['days', 'years'], help='<option> time units')
    parser.add_argument('-nreps', type=int, default=1, help='Number of replicate XMLs to generate.')
    args = parser.parse_args()
    
    if not args.out.endswith('.xml'):
        print ('<out> must end with .xml')
        sys.exit()
    
    parse(template_file=args.template, fasta_file=args.fasta, stem=args.stem, outfile=args.out, nreps=args.nreps)


if __name__ == '__main__':
    main()
