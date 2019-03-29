#!/bin/bash

#####################################################################################################################
### Command line argument should be: `fast[a|q] filename(s)`, N.B. filesnames must not contain spaces!            ###
###     Example: adapter_detection_script.sh *.fastq > ouput.txt                                                  ###
#####################################################################################################################

input_file="$@"
threads=13

# NexteraPE-PE.fa	(Source: trimmomatics_0.36_adapters_list)
NexteraPE_PE='>PrefixNX/1:AGATGTGTATAAGAGACAG
>PrefixNX/2:AGATGTGTATAAGAGACAG
>Trans1:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc:CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc:CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'

# TruSeq2-PE.fa		(Source: trimmomatic_0.36_adapters_list)
TruSeq2_PE='>PrefixPE/1:AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2:CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer1:AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR_Primer1_rc:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_Primer2:CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer2_rc:AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>FlowCell1:TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
>FlowCell2:TTTTTTTTTTCAAGCAGAAGACGGCATACGA'

# TruSeq2-SE.fa		(Source: trimmomatic_0.36_adapters_list)
TruSeq2_SE='>TruSeq2_SE:AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_f:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>TruSeq2_PE_r:AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG'

# TruSeq3-PE.fa		(Source: trimmomatic_0.36_adapters_list)
TruSeq3_PE='>PrefixPE/1:TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'

# TruSeq3-PE-2.fa	(Source: trimmomatic_0.36_adapters_list)
TruSeq3_PE_2='>PrefixPE/1:TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1:TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

# TruSeq3-SE.fa		(Source: trimmomatic_0.36_adapters_list)
TruSeq3_SE='>TruSeq3_IndexedAdapter:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3_UniversalAdapter:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'

#####

function count_adapter()
{
	library_prep_name=${1}
	adapter=${2}
	adapter_name=${2/:**/}
	adapter_seq=${2/**:/}
	input_file=${3}
	echo -e "${input_file}\t${library_prep_name}\t${adapter_name}\t${adapter_seq}\t$(grep -c "${adapter_seq}" ${input_file})"
}

export -f count_adapter

echo -e "Input_filename\tLibrary_prep_name\tAdapter_name\tAdapter_seq\tOccurence_count"
parallel -kj ${threads} count_adapter ::: "NexteraPE_PE" ::: ${NexteraPE_PE} ::: ${input_file}
parallel -kj ${threads} count_adapter ::: "TruSeq2_SE" ::: ${TruSeq2_SE} ::: ${input_file}
parallel -kj ${threads} count_adapter ::: "TruSeq2_PE" ::: ${TruSeq2_PE} ::: ${input_file}
parallel -kj ${threads} count_adapter ::: "TruSeq3_SE" ::: ${TruSeq3_SE} ::: ${input_file}
parallel -kj ${threads} count_adapter ::: "TruSeq3_PE" ::: ${TruSeq3_PE} ::: ${input_file}
parallel -kj ${threads} count_adapter ::: "TruSeq3_PE_2" ::: ${TruSeq3_PE_2} ::: ${input_file}
