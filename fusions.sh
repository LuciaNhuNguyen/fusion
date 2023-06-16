#!/bin/bash

# This script converts Arriba fusion calls into a format that is compatible with cBioPortal.

set -e -o pipefail

function print_usage {
	echo Usage: $(basename "$0") "--fusion-files FILE_CONTAINING_FUSION_FILES --output OUTPUT_FILE --annotation GTF_FILE --valid-gene-symbols FILE [--ncbi-build GRCh37|GRCh38]"
	echo "The fusion file must have two tab-separated columns (no header line): SAMPLE_ID	FUSION_FILE"
	exit 1
}

# define default parameters
NCBI_BUILD="GRCh37"

# check parameters
while [ $# -gt 0 ]; do
	PARAMETER="$1"
	shift
	case "$PARAMETER" in
		--fusion-files) FUSION_FILES="$1";;
		--output) OUTPUT_FILE="$1";;
		--annotation) ANNOTATION="$1";;
		--help|-h) print_usage;;
		--ncbi-build)
			NCBI_BUILD="$1"
			if ! [[ "$NCBI_BUILD" =~ ^GRCh3[78]$ ]]; then
				echo "Invalid NCBI build: $NCBI_BUILD" 1>&2
				exit 1
			fi
		;;
		--valid-gene-symbols) VALID_GENE_SYMBOLS="$1";;
		*) echo "Invalid parameter: $PARAMETER" 1>&2 && print_usage;;
	esac
	shift
done
for PARAMETER in "$FUSION_FILES" "$OUTPUT_FILE" "$ANNOTATION" "$VALID_GENE_SYMBOLS"; do
	if [ -z "$PARAMETER" ]; then
		echo "Missing mandatory parameter" 1>&2
		print_usage
	fi
done

# make sure input files exist
(
	echo "$FUSION_FILES"
	echo "$ANNOTATION"
	echo "$VALID_GENE_SYMBOLS"
	cut -f2 "$FUSION_FILES" 2> /dev/null || true
	
) | while read FILE; do
	if [ ! -e "$FILE" ]; then
		echo "File not found: $FILE" 1>&2
		exit 1
	fi
done

# convert GTF to BED
ANNOTATION_BED=$(
	tr -d '"' < "$ANNOTATION" |
	sort -k3,3r | # sort to make sure transcripts come before exons
	awk -F '\t' -v OFS='\t' -v VALID_GENE_SYMBOLS="$VALID_GENE_SYMBOLS" '

		# load gene symbols recognized by cBioPortal
		FILENAME==VALID_GENE_SYMBOLS { valid_gene_symbols[$0]=1 }

		# parse GTF attributes
		function get_gtf_attribute(name) {
			if (index($9, name)) {
				result=$9
				sub(".*" name " ", "", result)
				sub(";.*", "", result)
				return(result)
			} else
				return(".")
		}

		# load transcript sizes/status and compute priority
		FILENAME=="/dev/stdin" && $3=="transcript" {
			transcript_id=get_gtf_attribute("transcript_id")
			transcript_priority[transcript_id]=$5-$4+1 # transcript length
			if (/tag appris_principal_1;/) { transcript_priority[transcript_id]+=120000000 }
			else if (/tag appris_principal_2;/) { transcript_priority[transcript_id]+=110000000 }
			else if (/tag appris_principal_3;/) { transcript_priority[transcript_id]+=100000000 }
			else if (/tag appris_principal_4;/) { transcript_priority[transcript_id]+=90000000 }
			else if (/tag appris_principal_5;/) { transcript_priority[transcript_id]+=80000000 }
                        else if (/tag appris_principal;/) { transcript_priority[transcript_id]+=70000000 }
                        else if (/tag appris_candidate_longest;/) { transcript_priority[transcript_id]+=60000000 }
                        else if (/tag appris_candidate;/) { transcript_priority[transcript_id]+=50000000 }
                        else if (/tag appris_alternative_1;/) { transcript_priority[transcript_id]+=40000000 }
                        else if (/tag appris_alternative_2;/) { transcript_priority[transcript_id]+=30000000 }
                        else if (/tag appris_alternative;/) { transcript_priority[transcript_id]+=20000000 }
                        else if (/tag CCDS;/) { transcript_priority[transcript_id]+=10000000 }
		}

		# convert GTF to BED
		FILENAME=="/dev/stdin" && $3=="exon" {
			gene_name=get_gtf_attribute("gene_name")
			gene_id=get_gtf_attribute("gene_id")
			transcript_id=get_gtf_attribute("transcript_id")
			exon_number=get_gtf_attribute("exon_number")
			if (valid_gene_symbols[gene_name])
				print $1, ($4-1), $5, gene_name, gene_id, transcript_id, transcript_priority[transcript_id], exon_number
		}

	' "$VALID_GENE_SYMBOLS" /dev/stdin |
	sort -k1,1V -k2,2n -k3,3n
)

# reannotate breakpoints in case the annotation by Arriba is incomplete or not compatible with cBioPortal in terms of gene symbols
REANNOTATED_BREAKPOINTS=$(
	for ORIENTATION in upstream downstream; do
		cut -f2 "$FUSION_FILES" |
		xargs -d $'\n' cat |
		awk -F '\t' -v ORIENTATION=$ORIENTATION '
			# get column headers
			/^#/ { for (i=1; i<=NF; i++) col[$i]=i }
			# find breakpoints with given orientation (upstream/downstream)
			!(/^#/) && $col["direction1"]!=ORIENTATION { print $col["breakpoint1"], $col["direction1"] }
			!(/^#/) && $col["direction2"]!=ORIENTATION { print $col["breakpoint2"], $col["direction2"] }
		' |
		awk '!duplicate[$0]++ {
			# convert breakpoint to BED format
			chrom=$1; sub(/:[0-9]*$/, "", chrom)
			pos=$1; sub(/.*:/, "", pos)
			print chrom "\t" (pos-1) "\t" pos "\t" $2
		}' |
		bedtools closest -D ref -i${ORIENTATION:0:1} -a /dev/stdin -b <(cat <<<"$ANNOTATION_BED")
	done
)

# concatenate fusion files and prefix each line with sample ID
FUSIONS=$(
	cat "$FUSION_FILES" | while read LINE; do
		SAMPLE_ID=$(cut -f1 <<<"$LINE")
		FUSION_FILE=$(cut -f2 <<<"$LINE")
		awk -v SAMPLE_ID="$SAMPLE_ID" '
			FNR==1 { sub(/^#/, "sample_id\t", $0); print }
			FNR>1 { print SAMPLE_ID "\t" $0 }
		' "$FUSION_FILE"
	done
)

# print column headers of output file
echo "Sample_ID	Site1_Hugo_Symbol	Site1_Entrez_Gene_Id	Site1_Ensembl_Transcript_Id	Site1_Exon	Site1_Chromosome	Site1_Position	Site1_Description	Site2_Hugo_Symbol	Site2_Entrez_Gene_Id	Site2_Ensembl_Transcript_Id	Site2_Exon	Site2_Chromosome	Site2_Position	Site2_Description	Site2_Effect_On_Frame	NCBI_Build	DNA_Support	RNA_Support	Normal_Read_Count	Tumor_Read_Count	Normal_Variant_Count	Tumor_Variant_Count	Normal_Paired_End_Read_Count	Tumor_Paired_End_Read_Count	Normal_Split_Read_Count	Tumor_Split_Read_Count	Annotation	Breakpoint_Type	Connection_Type	Event_Info	Class	Length	Comments	External_Annotation" > "$OUTPUT_FILE"

# convert Arriba format to cBioPortal format
cat <<<"$REANNOTATED_BREAKPOINTS" |
awk -F '\t' -v OFS='\t' -v NCBI_BUILD="$NCBI_BUILD" -v REANNOTATED_BREAKPOINTS="/dev/stdin" -v VALID_GENE_SYMBOLS="$VALID_GENE_SYMBOLS" '

	# load reannotation of breakpoints by bedtools
	FILENAME==REANNOTATED_BREAKPOINTS {
		annotation_breakpoint[FNR]=$1 ":" $3
		annotation_direction[FNR]=$4
		annotation_chrom[FNR]=$5
		annotation_start[FNR]=$6
		annotation_end[FNR]=$7
		annotation_gene_name[FNR]=$8
		annotation_gene_id[FNR]=$9
		annotation_transcript_id[FNR]=$10
		annotation_transcript_prio[FNR]=$11
		annotation_exon_number[FNR]=$12
		annotation_index[annotation_breakpoint[FNR] " " annotation_direction[FNR]]=annotation_index[annotation_breakpoint[FNR] " " annotation_direction[FNR]] " " FNR
	}

	# load gene symbols recognized by cBioPortal
	FILENAME==VALID_GENE_SYMBOLS { valid_gene_symbols[$0]=1 }

	# get column headers of fusion file
	FILENAME!=REANNOTATED_BREAKPOINTS && FILENAME!=VALID_GENE_SYMBOLS && /^sample_id\t/ {
		for (i=1; i<=NF; i++) col[$i]=i
		fusion_number=0
	}

	# iterate through fusion file
	FILENAME!=REANNOTATED_BREAKPOINTS && FILENAME!=VALID_GENE_SYMBOLS && !(/^sample_id\t/) {

		fusion_number++

		# find reannotation of breakpoint
		split("", reannotated_breakpoint)
		for (brk=1; brk<=2; brk++) {
			split(annotation_index[$col["breakpoint" brk] " " $col["direction" brk]], annotation_hits, " ")
			for (a in annotation_hits) {
				i=annotation_hits[a]
				if (i!="") {
					# check if annotation record matches given breakpoint
					if (!valid_gene_symbols[$col["gene" brk]] ||
					     match($col["gene" brk], "(^|,)" annotation_gene_name[i] "($|[(])") &&
					     (!("gene_id" brk in col) || $col["gene_id" brk]==annotation_gene_id[i] || $col["gene_id" brk]==".") &&
					     (!("transcript_id" brk in col) || $col["transcript_id" brk]==annotation_transcript_id[i] || $col["transcript_id" brk]==".")) {
						# if a breakpoint is annotated with multiple records, pick the one with the longest transcript
						if (!(brk in reannotated_breakpoint) ||
						    annotation_exon_number[reannotated_breakpoint[brk]]~/^(|[.])$/ ||
						    annotation_exon_number[i]!~/^(|[.])$/ && annotation_transcript_prio[i] > annotation_transcript_prio[reannotated_breakpoint[brk]]) {
							reannotated_breakpoint[brk]=i
						}
					}
				}
			}
		}

		Sample_ID=$col["sample_id"]

		Site1_Hugo_Symbol=($col["site1"]!="intergenic" && valid_gene_symbols[$col["gene1"]]) ? $col["gene1"] : annotation_gene_name[reannotated_breakpoint[1]]

		Site1_Entrez_Gene_Id="" # is annotated by cBioPortal

		Site1_Ensembl_Transcript_Id=("transcript_id1" in col && $col["transcript_id1"]!=".") ? $col["transcript_id1"] : annotation_transcript_id[reannotated_breakpoint[1]]

		Site1_Exon=annotation_exon_number[reannotated_breakpoint[1]]

		Site1_Chromosome=$col["breakpoint1"]
		sub(/:[0-9]*$/, "", Site1_Chromosome)

		Site1_Position=$col["breakpoint1"]
		sub(/.*:/, "", Site1_Position)

		Site1_Description=$col["site1"]

		Site2_Hugo_Symbol=($col["site2"]!="intergenic" && valid_gene_symbols[$col["gene2"]]) ? $col["gene2"] : annotation_gene_name[reannotated_breakpoint[2]]

		Site2_Entrez_Gene_Id="" # is annotated by cBioPortal

		Site2_Ensembl_Transcript_Id=("transcript_id2" in col && $col["transcript_id2"]!=".") ? $col["transcript_id2"] : annotation_transcript_id[reannotated_breakpoint[2]]

		Site2_Exon=annotation_exon_number[reannotated_breakpoint[2]]

		Site2_Chromosome=$col["breakpoint2"]
		sub(/:[0-9]*$/, "", Site2_Chromosome)

		Site2_Position=$col["breakpoint2"]
		sub(/.*:/, "", Site2_Position)

		Site2_Description=$col["site2"]

		Site2_Effect_On_Frame=("reading_frame" in col && $col["reading_frame"]!=".") ? $col["reading_frame"] : ""
		sub("in-frame", "IN_FRAME", Site2_Effect_On_Frame)
		sub("out-of-frame|stop-codon", "FRAMESHIFT", Site2_Effect_On_Frame)

		NCBI_Build=NCBI_BUILD

		DNA_Support="no"

		RNA_Support="yes"

		Normal_Read_Count=""

		Tumor_Read_Count=""

		Normal_Variant_count=""

		Tumor_Variant_Count=$col["split_reads1"]+$col["split_reads2"]+$col["discordant_mates"]

		Normal_Paired_End_Read_Count=""

		Tumor_Paired_End_Read_Count=$col["discordant_mates"]

		Normal_Split_Read_Count=""

		Tumor_Split_Read_Count=$col["split_reads1"]+$col["split_reads2"]

		Annotation=Site1_Hugo_Symbol "--" Site2_Hugo_Symbol " fusion"

		Breakpoint_Type=($col["split_reads1"]+$col["split_reads2"] > 0) ? "PRECISE" : "IMPRECISE"

		Connection_Type=($col["type"]~/5.-5/) ? "5to5" : (($col["type"]~/3.-3/) ? "3to3" : "5to3")

		Event_Info="confidence=" $col["confidence"]

		Class=$col["type"]
		sub(/.*translocation.*/, "TRANSLOCATION", Class)
		sub(/.*duplication.*/, "DUPLICATION", Class)
		sub(/.*inversion.*/, "INVERSION", Class)
		sub(/.*deletion.*/, "DELETION", Class)

		Length=(Class=="TRANSLOCATION") ? "" : Site1_Position - Site2_Position
		if (Class!="TRANSLOCATION" && Length<0) Length *= -1

		Comments=""

		External_Annotation=("tags" in col && $col["tags"]!=".") ? $col["tags"] : ""

		# cBioPortal requires that either all transcripts & exons be set or none
		if (Site1_Ensembl_Transcript_Id=="" || Site2_Ensembl_Transcript_Id=="" || Site1_Exon=="" || Site2_Exon=="") {
			Site1_Ensembl_Transcript_Id=Site2_Ensembl_Transcript_Id="NA"
			Site1_Exon=Site2_Exon=-1
		}

		if (Site1_Hugo_Symbol!="" && Site2_Hugo_Symbol!="") {
			print Sample_ID,Site1_Hugo_Symbol,Site1_Entrez_Gene_Id,Site1_Ensembl_Transcript_Id,Site1_Exon,Site1_Chromosome,Site1_Position,Site1_Description,Site2_Hugo_Symbol,Site2_Entrez_Gene_Id,Site2_Ensembl_Transcript_Id,Site2_Exon,Site2_Chromosome,Site2_Position,Site2_Description,Site2_Effect_On_Frame,NCBI_Build,DNA_Support,RNA_Support,Normal_Read_Count,Tumor_Read_Count,Normal_Variant_count,Tumor_Variant_Count,Normal_Paired_End_Read_Count,Tumor_Paired_End_Read_Count,Normal_Split_Read_Count,Tumor_Split_Read_Count,Annotation,Breakpoint_Type,Connection_Type,Event_Info,Class,Length,Comments,External_Annotation
		} else {
			print "Skipping fusion #" fusion_number " of sample " Sample_ID ", because one of the breakpoints cannot be mapped to a gene symbol recognized by cBioPortal" > "/dev/stderr"
		}
	}
' /dev/stdin "$VALID_GENE_SYMBOLS" <(cat <<<"$FUSIONS") |

# empty fields containing only a dot
awk -F '\t' -v OFS='\t' '{
	for (col=1; col<=NF; col++)
		if ($col == ".")
			$col = ""
	print
}' >> "$OUTPUT_FILE"

