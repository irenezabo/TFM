
#!/bin/bash

# Irene modification Nov 2024. Original script: github Fabian Jetzinger https://github.com/FabianJetzinger/sort_gtf
# The script sorts the lines of a GTF based on chromosome, position, feature...
# Modified to ensure feature order (Bmo GTF was problematic for DISDOM analysis)


# Initialize variables for input and output files and patterns
input=""
output=""
gene_id_pattern='gene_id "([^"]+)"'  # Default gene_id pattern
transcript_id_pattern='transcript_id "([^"]+)"'  # Default transcript_id pattern
exon_num_pattern='exon_number "([^"]+)"'

# Parse command-line options
while getopts ":i:o:g:t:e:-:" opt; do
    case "${opt}" in
        i)
            input="${OPTARG}"
            ;;
        o)
            output="${OPTARG}"
            ;;
        g)
            gene_id_pattern="${OPTARG}"
            ;;
        t)
            transcript_id_pattern="${OPTARG}"
            ;;
        e)
            exon_num_pattern="${OPTARG}"
            ;;
        *)
            echo "Usage: $0 [-i|< input_file] [-o|> output_file] [-g gene_id_pattern] [-t transcript_id_pattern] [-e exon_num_pattern]"
            echo "  -i   Specify the input GTF file."
            echo "  -o   Specify the output GTF file."
            echo "  -g   Specify the gene ID matching pattern."
            echo "  -t   Specify the transcript ID matching pattern."
	    echo "  -e   Specify the exon number matching pattern."
            exit 1
            ;;
    esac
done

# Check if no input file provided, then read from stdin
if [[ -z "$input" ]]; then
    input="/dev/stdin"
fi

# Check if no output file provided, then write to stdout
if [[ -z "$output" ]]; then
    output="/dev/stdout"
fi

echo $gene_id_pattern
echo $transcript_id_pattern
echo $exon_num_pattern

awk -F'\t' '{
    # Initial sorting of input by start position of entry, keeping chromosome order intact
    # should ensure that we find the correct gene/transcript start position in the first entry for every gene_id/transcript_id

    # remove commented lines
    if($1 ~ /^#/) next;

    if(!seen[$1]) {
        seen[$1] = ++chrom_count
    }
    chrom_id=seen[$1]
    print $0 "\t" chrom_id
}' "$input" | sort -t$'\t' -k10,10n -k4,4n | \
awk -F'\t' -v gene_pattern="$gene_id_pattern" -v transcript_pattern="$transcript_id_pattern" -v exon_pattern="$exon_num_pattern" '{
    # Main processing

    # Clear the match arrays at the start of each record processing
    split("", gene_match)
    split("", transcript_match)
    split("", exon_match)
    
    # match gene_id
    match($9, gene_pattern, gene_match)
    gene_id=gene_match[1]

    # match transcript_id
    match($9, transcript_pattern, transcript_match)
    transcript_id=transcript_match[1]
    
    # match exon_num
    match($9, exon_pattern, exon_match)
    exon_num=exon_match[1]
    
    # Keep track of gene and transcript start position
    # This is done to handle overlapping genes/transcripts; so e.g. transcript 1 and all its exons will appear before transcript 2 and all its exons, 
    # rather than having e.g. transcript 1, exon 1.1, transcript 2, exon 2.1, exon 1.2, exon 2.2
    gene_start_pos = $4
    transcript_start_pos = $4
    chrom_id=$10
    if(gene_id != "") {
        chrom_gene = chrom_id "_" gene_id
        if(!gene_start_pos_map[chrom_gene]) {
            gene_start_pos_map[chrom_gene] = $4
        }
        gene_start_pos = gene_start_pos_map[chrom_gene]

        if(transcript_id != "") {
            chrom_gene_transcript = chrom_id "_" gene_id "_" transcript_id
            if(!transcript_start_pos_map[chrom_gene_transcript]) {
                transcript_start_pos_map[chrom_gene_transcript] = $4
            }
            transcript_start_pos = transcript_start_pos_map[chrom_gene_transcript]
        }
    }


    # determine entry type so gene > transcript > exon > CDS > start/stop > others
  
    if($3=="gene") {
        transcript_id = "0"
        entry_type = "0"
	exon_num = "0"
    } else if($3=="transcript") {
        entry_type = "1"
	exon_num = "0"
    } else if($3=="exon") {
        entry_type = "2"
    } else if($3=="CDS") {
        entry_type = "3"
    } else if($3=="start_codon") {
        entry_type = "4"
    } else if($3=="stop_codon") {
        entry_type = "4"
    } else {
        # for other entries: if they have a gene id, put them after their gene; if they also have a transcript id, put them after their transcript and after other defined features.
        if(gene_id == "") gene_id = "0"
        if(transcript_id == "") transcript_id = "00"
	if(exon_num == "") exon_num = "000"
        entry_type = "5"
    }
    print $0 "\t" gene_start_pos "\t" transcript_start_pos "\t" gene_id "\t" transcript_id "\t" exon_num "\t" entry_type "\t" $3
}' | sort -t$'\t' -k10,10n -k11,11n -k12,12n -k13,13 -k14,14 -k15,15n -k16,16n -k17,17 |  cut -f1-9  > "$output"

# explanation:  1. -k10,10n: sort chromosomes
#               2. -k11,11n: sort all entries within a chromosome by gene start position (if no gene_id, use the entry's own start position)
#                   --> this means, if e.g. gene 1 and 2 overlap, all entries associated with gene 1 should come before gene 2, even if e.g. the start position of the last exon of gene 1 is later than the start of gene 2
#               3. -k12,12n: sort all entries within the same gene start position (either same gene or two genes with same start) by transcript start position (or its own start position if no transcript_id)
#                   --> same logic but on transcript level; tie between gene entry and its transcripts with the same start position
#               4. -k13,13: sort by gene_id
#                   --> for genes with the same start position, sort them (and all entries belonging to them) by their ID alphanumerically
#               5. -k14,14: sort by transcript_id
#                   --> same logic but on transcript level; tie between: gene entry + first transcript + all its exons
#               ###6. -k4,4n: sort by start position of the entry itself (IRENE MODIFIED)
#               6. -k15,15n: sort by exon number
#		7. -k16,16n: sort by kind of entry, where gene > transcript > others > exon > CDS > start/stop > other
#                   --> puts entries with remaining ties in order so gene is on top, then its first transcript, then any other entries associated with the transcript, then its first exon
#               8. -k17,17: Resolve any remaining ties between other entry types by alphanumerically sorting the name of the entry type
