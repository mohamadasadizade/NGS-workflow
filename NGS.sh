#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <user_id> <job_id> [-ignrqc]"
    exit 1
fi

user_id="$1"
job_id="$2"
ignore_qc=false 

if [ "$#" -eq 3 ] && [ "$3" == "-ignrqc" ]; then
    ignore_qc=true
fi

input_fastq="/root/NextFlow_project/test/${user_id}/${job_id}/inputfiles/TEST_R1.fastq.gz"
input_vcf="/root/NextFlow_project/test/${user_id}/${job_id}/inputfiles/input.vcf"
annotation_db="/root/NextFlow_project/test/db/clinvar.vcf"
ref_genome="/root/NextFlow_project/test/hg38.fa"
snpEff="/root/snpEff/SnpEff.jar"
snpSift="/root/snpEff/SnpSift.jar"
cores=8

qc_dir="${user_id}/${job_id}/qc"
alignment_dir="${user_id}/${job_id}/alignment"
vcf_dir="${user_id}/${job_id}/vcf"
annotation_dir="${user_id}/${job_id}/annotation"
filtered_dir="${user_id}/${job_id}/filtered"
parsed_dir="${user_id}/${job_id}/parsed"

if [ ! -f "$input_fastq" ] && [ ! -f "$input_vcf" ]; then
    echo "Error: Neither VCF nor FASTQ input files found. Please provide at least one of them."
    exit 1
elif [ -f "$input_fastq" ] && [ -f "$input_vcf" ]; then
    echo "Both FASTQ and VCF files found. Starting pipeline from Step 1 (FastQC)..."
    start_step=1
elif [ -f "$input_fastq" ]; then
    echo "Only FASTQ file found. Starting pipeline from Step 1 (FastQC)..."
    start_step=1
elif [ -f "$input_vcf" ]; then
    echo "Only VCF file found. Starting pipeline from Step 4 (SnpEff Annotation)..."
    start_step=4
fi

# 1. FastQC
if [ "$start_step" -le 1 ] && [ "$ignore_qc" = false ] && [ ! -f "${qc_dir}/run.log" ]; then
    echo "Running FastQC..."
    mkdir -p $qc_dir
    fastqc -o $qc_dir $input_fastq
    echo 'FastQC completed successfully.' >> ${qc_dir}/run.log
else
    if [ "$ignore_qc" = true ]; then
        echo "Skipping FastQC due to -ignrqc flag."
    else
        echo "FastQC already completed or not required. Skipping..."
    fi
fi

# 2. BWA MEM and Samtools sort
if [ "$start_step" -le 1 ] && [ ! -f "${alignment_dir}/output.bam" ]; then
    echo "Running BWA MEM and Samtools sort..."
    mkdir -p $alignment_dir
    bwa mem -t $cores $ref_genome $input_fastq | samtools sort -@$cores -o ${alignment_dir}/output.bam
    echo 'BWA MEM and Samtools sort completed successfully.' >> ${alignment_dir}/run.log
else
    echo "BWA MEM and Samtools sort already completed or not required. Skipping..."
fi

# 3. FreeBayes
if [ "$start_step" -le 1 ] && [ ! -f "${vcf_dir}/output.vcf" ]; then
    echo "Running FreeBayes..."
    mkdir -p $vcf_dir
    freebayes -f $ref_genome ${alignment_dir}/output.bam > ${vcf_dir}/output.vcf
    echo 'FreeBayes completed successfully.' >> ${vcf_dir}/run.log
else
    echo "FreeBayes already completed or not required. Skipping..."
fi

# 4. SnpEff Annotation
if [ "$start_step" -le 4 ] && [ ! -f "${annotation_dir}/output.annotated.vcf" ]; then
    echo "Running SnpEff annotation..."
    mkdir -p $annotation_dir

    # Determine which VCF to use: the one from Step 3 or the input VCF
    if [ -f "${vcf_dir}/output.vcf" ]; then
        input_vcf_file="${vcf_dir}/output.vcf"
    else
        input_vcf_file="$input_vcf"  # Use provided input VCF
    fi

    java -Xmx8g -jar $snpEff -v -stats ${annotation_dir}/VarStats.html GRCh37.75 $input_vcf_file > ${annotation_dir}/output.annotated.vcf
    echo 'SnpEff annotation completed successfully.' >> ${annotation_dir}/run.log
else
    echo "SnpEff annotation already completed or not required. Skipping..."
fi


# 5. Variant Filtering with SnpSift
if [ "$start_step" -le 4 ] && [ ! -f "${filtered_dir}/Final.vcf" ]; then
    echo "Running Variant Filtering..."
    mkdir -p $filtered_dir
    java -jar $snpSift filter "((ANN[*] =~ '.*\\|.*\\|HIGH\\|.*') | (ANN[*] =~ '.*\\|.*\\|MODERATE\\|.*'))" ${annotation_dir}/output.annotated.vcf > ${filtered_dir}/output.filtered.vcf

    java -jar $snpSift Annotate $annotation_db ${filtered_dir}/output.filtered.vcf > ${filtered_dir}/output.filtered.clinvar.vcf
    java -jar $snpSift Filter "(exists CLNDBN) & ((ANN[*].EFFECT has 'stop_gained')|(ANN[*].EFFECT has 'frameshift_variant')|(ANN[*].EFFECT has 'stop_lost'))" ${filtered_dir}/output.filtered.clinvar.vcf > ${filtered_dir}/Final.vcf
    echo 'Variant filtering completed successfully.' >> ${filtered_dir}/run.log
else
    echo "Variant filtering already completed or not required. Skipping..."
fi

# 6. Variant Parsing and Checking
if [ "$start_step" -le 4 ] && [ ! -f "${parsed_dir}/parsed.txt" ]; then
    echo "Running Variant Parsing and Checking..."
    mkdir -p $parsed_dir
    python3 parse_and_check_variants.py ${filtered_dir}/Final.vcf ${parsed_dir}/parsed.txt
    parse_exit_code=$?

    if [ $parse_exit_code -eq 1 ]; then
        echo "No variants found, running alternative filtering..."
        
        java -jar $snpSift Filter "(exists CLNDBN) & ((ANN[*].EFFECT has 'missense_variant')|(ANN[*].EFFECT has 'inframe_deletion'))" ${filtered_dir}/output.filtered.clinvar.vcf > ${filtered_dir}/Final_alt.vcf
        echo 'Alternative filtering completed successfully.' >> ${filtered_dir}/run_alt.log

        python3 parse_and_check_variants.py ${filtered_dir}/Final_alt.vcf ${parsed_dir}/parsed_alt.txt
        echo 'Final variant parsing with alternative filter completed successfully!' >> ${parsed_dir}/run_alt.log

    elif [ $parse_exit_code -eq 0 ]; then
        echo "Variants found, no need for alternative filtering."
    fi

else
    echo "Variant parsing and checking already completed or not required. Skipping..."
fi
