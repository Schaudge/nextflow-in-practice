VERSION="1.0"

log.info "===================================================="
log.info "YunYing NGS Analysis Nextflow Pipeline (v${VERSION})"
log.info "Process begins at ${workflow.start}!"
log.info "===================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run panel_somatic_gatk4.nf -c nextflow_basic.conf -profile docker --SampleId *** --PanelName *** --fq1 R1.fq.gz --fq2 R2.fq.gz --ctl_fq1 X_R1.fq.gz --ctl_fq2 X_R2.fq.gz --bed target_bed_file"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --SampleId         STRING             Sample name (dafault: fq1 basename)"
  log.info "    --PanelName        STRING             Panel name, used for exported gene filter (required)"
  log.info "    --fq1              FILE               tumor fastq(.gz) file path for read R1 (required)"
  log.info "    --fq2              FILE               tumor fastq(.gz) file path for read R2 (required)"
  log.info "    --ctl_fq1          FILE               control fastq(.gz) file path for read R1 (required for paired sample)"
  log.info "    --ctl_fq2          FILE               control fastq(.gz) file path for read R2 (required for paired sample)"
  log.info "    --bed              FILE               the (target) capture bed file (absolute) path (required)"
  log.info "    --st_bed           FILE               the bed file (absolute) path used for region statistics (optional)"
  log.info "    --sv_bed           FILE               the bed file (absolute) path used for sv analysis (optional)"
  log.info "    --baseline         STRING             ONCOCNV control clade flag (optional, for cnv is required)"
  log.info "    --msi_db           STRING             msi database subdir in related to directory /yunying/ref/... (optional, for msi is required)"
  log.info "    --msi_depth        STRING             msi locus depth threshold (optional, for msi is required)"
  log.info "    --output           DIR                Output directory (default: /yunying/data/product/result/dev_test/)"
  log.info "    --RG               STRING             Read group tag (dafault: fq1 basename)"
  log.info " "
  log.info "===================================================================="
  exit 1
}


tumor_fq1 = file("$params.fq1")
tumor_fq2 = file("$params.fq2")
params.ctl_fq1 = ""
params.ctl_fq2 = ""
ctl_fq1 = params.ctl_fq1 ? file("$params.ctl_fq1") : ""
ctl_fq2 = params.ctl_fq2 ? file("$params.ctl_fq2") : ""
bed_target = Channel.value("$params.bed")
params.output = "/yunying/data/product/result/dev_test"
params.SampleId = tumor_fq1.name.split('\\_')[0]
params.RG = tumor_fq1.name.split('\\_')[0]
params.high_recall = false
params.threads = 16
params.remove_dup = true
params.recal = 'elprep'
params.st_bed = ''
params.sv_bed = ''
params.baseline = ''
params.msi_db = ''
params.codinglength = ''
params.queue_name = "dev"
params.update_db = "no"


def runID = "uuid"
def raw_unfilter_dir = "${params.output}/unfiltered_raw_vars"
def filtered_review_dir = "${params.output}/"
parent_top_dir = file("${params.output}").parent
def chip_batch_id = ("${parent_top_dir}" - ~/.*\\//)
def unified_bam_dir = "${parent_top_dir}/bam"
st_bed_ch = params.st_bed ? Channel.value("$params.st_bed") : Channel.value("$params.bed")
sv_bed_ch = params.sv_bed ? Channel.value("$params.sv_bed") : Channel.value("$params.bed")
sv_process_mark = params.sv_bed =~ /null.bed$/ ? false : true


memory_scale1 = Math.ceil(tumor_fq1.size()/1024/1024/1024 + tumor_fq2.size()/1024/1024/1024)
memory_scale2 = Math.ceil(ctl_fq1.size()/1024/1024/1024 + ctl_fq2.size()/1024/1024/1024)
memory_scale1 = params.high_recall ? memory_scale1 * 2 : memory_scale1
memory_scale2 = params.high_recall ? memory_scale2 * 2 : memory_scale2
println "The memory scale factor is ${memory_scale1}, ---, ${memory_scale2} * 15G for memory configuration!"
split_region_count = params.threads/2

log.info "Sample ${params.SampleId} was running ..."
log.info "And the result files will be put into directory ${filtered_review_dir}."

// alignment database and gatk bundle database
def genome_dir = "/yunying/ref/human/b37/"
def ref_fa = "${genome_dir}b37_Homo_sapiens_assembly19.fasta"
def ref_elfa = "${genome_dir}b37_Homo_sapiens_assembly19.elfasta"
def ref_fa_dict = "${genome_dir}b37_Homo_sapiens_assembly19.dict"
def gatk_db_dir = "${genome_dir}gatk"
def oncocnv_db_dir = "${genome_dir}oncocnv"
def snpeff_db_dir = "${genome_dir}snpeff"
def hotspots_mutscan_db = file("${genome_dir}mutscan/hotspots.csv")
def core_hotspots_conf = file("${genome_dir}variation_hotspots/core_hotspots_v1.0.vcf")
def core_hotspots_bed = file("${genome_dir}variation_hotspots/core_hotspots_v1.0.bed")
def complex_hotspot_bed = file("${genome_dir}variation_hotspots/complex_hotspot_mutation.bed")
def chemotherapy_rs_vcf = file("${genome_dir}variation_hotspots/chemotherapy_rs_v1.0.vcf")
def chemotherapy_rs_bed = file("${genome_dir}variation_hotspots/chemotherapy_rs_v1.0.bed")
def radiotherapy_rs_vcf = file("${genome_dir}variation_hotspots/radiotherapy_rs_v1.0.vcf")
def chemo_radio_rs_vcf = file("${genome_dir}variation_hotspots/chemotherapy_radiotherapy_rs_v1.0.vcf")
def chemo_radio_rs_bed = file("${genome_dir}variation_hotspots/chemotherapy_radiotherapy_rs_v1.0.bed")
def clinvar_vcf_conf = file("${genome_dir}snpsift/clinvar_20210828.vcf.gz")
def variation_hotspots_conf = file("${genome_dir}variation_hotspots/yunying_hotspots_GRCh37_v1.0.txt")
def msi_database_dir = "${genome_dir}/msi/${params.msi_db}"

// all configuration including qc, alignment, and gatk process configuration file (module ---> software/toolkit ---> conf)
def prefix_conf_dir = "/yunying/codes/product/module"
def preprocess_conf = file("${prefix_conf_dir}/software/preprocess/fastp/conf/fastp_standard.conf")
def mutscan_conf = file("${prefix_conf_dir}/software/somatic/Mutscan/conf/mutscan_filter.conf")
def align_conf = file("${prefix_conf_dir}/toolkit/align_sambam_vcf/conf/align_standard.conf")
def gatk_toolkit_dir = "${prefix_conf_dir}/toolkit/gatk/"
def recal_conf = file("${gatk_toolkit_dir}refinement/conf/refinement_standard.conf")
def sv_svaba_conf = file("${prefix_conf_dir}/software/sv/svaba/conf/somatic_sv_paired.conf")
def sv_hotspots_conf = file("${genome_dir}anno/sv_gene_extend_region.txt")
def sv_partner_conf = file("${genome_dir}anno/sv_hotspots_partner.txt")
def cnv_oncocnv_conf = file("${prefix_conf_dir}/software/cnv/oncocnv/conf/somatic_cnv_tumor.conf")
somatic_calling_conf = params.high_recall ? "somatic_high_depth_cfdna.conf" : "somatic_tumor_standard.conf"
trim_in_conf = params.high_recall ? "trim_include_conservative.conf" : "trim_include_standard.conf"
trim_ex_conf = params.high_recall ? "trim_exclude_conservative2.conf" : "trim_exclude_standard.conf"
def mutect2_conf = file("${gatk_toolkit_dir}somatic/conf/${somatic_calling_conf}")
def secondary_annotate_conf = file("${gatk_toolkit_dir}annotation/conf/secondary_bam_annotate.conf")
def germline_calling_conf = file("${gatk_toolkit_dir}germline/conf/haplotype_caller.conf")
def anno_snpeff_conf = file("${prefix_conf_dir}/software/annotation/SnpEff/conf/snpeff_canon_only.conf")
def somatic_filter_conf = file("${gatk_toolkit_dir}filter/conf/somatic_high_recall.conf")
def vcf_trim_in_conf = file("${prefix_conf_dir}/toolkit/align_sambam_vcf/conf/${trim_in_conf}")
def vcf_trim_ex_conf = file("${prefix_conf_dir}/toolkit/align_sambam_vcf/conf/${trim_ex_conf}")
def perbase_conf = file("${prefix_conf_dir}/software/calling/perbase/conf/perbase_standard.conf")


process preprocess_fastp_tumor_standard {
     // use fastp to trim low quality bases from reads and do adapters trimming
     container '10.168.2.67:5000/bromberglab/fastp:latest'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*qc.json'

     cpus params.threads
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     input:
         file tumor_fq1
         file tumor_fq2

     output:
         path "${tumor_fastq_r1}" into trimmed_tumor_fq_r1_ch
         path "${tumor_fastq_r2}" into trimmed_tumor_fq_r2_ch
         path "${fastp_qc_json}" into fastq_qc_json

     script:
         tumor_fastq_r1 = "${params.SampleId}_tumor_r1_trimmed.fastq"
         tumor_fastq_r2 = "${params.SampleId}_tumor_r2_trimmed.fastq"
         fastp_qc_json = "${params.SampleId}_tumor_fastp_qc.json"
     """
     cat ${preprocess_conf} | xargs \
     fastp -i ${tumor_fq1} -I ${tumor_fq2} -w ${task.cpus} -j ${fastp_qc_json} -o ${tumor_fastq_r1} -O ${tumor_fastq_r2}
     """
}


process preprocess_fastp_control_standard {
     // use fastp to trim low quality bases from reads and do adapters trimming
     container '10.168.2.67:5000/bromberglab/fastp:latest'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*qc.json'

     cpus params.threads
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3
     
     when:
         params.ctl_fq1 != ""

     input:
         file ctl_fq1
         file ctl_fq2

     output:
         path "${ctl_fastq_r1}" into trimmed_ctl_fq_r1_ch
         path "${ctl_fastq_r2}" into trimmed_ctl_fq_r2_ch
         path "${fastp_qc_json}" into ctl_fastq_qc_json

     script:
         ctl_fastq_r1 = "${params.SampleId}_ctl_r1_trimmed.fastq"
         ctl_fastq_r2 = "${params.SampleId}_ctl_r2_trimmed.fastq"
         fastp_qc_json = "${params.SampleId}_ctl_fastp_qc.json"
     """
     cat ${preprocess_conf} | xargs \
     fastp -i ${ctl_fq1} -I ${ctl_fq2} -w ${task.cpus} -j ${fastp_qc_json} -o ${ctl_fastq_r1} -O ${ctl_fastq_r2}
     """
}

process tumor_align_bwa_sort {
     //container 'cmopipeline/bwa-samtools-gatk4-sambamba-samblaster:latest'
     container '10.168.2.67:5000/align_sambam_vcf:1.0'

     cpus params.threads
     memory { 20.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     input:
         path('trimmed_fastq_r1') from trimmed_tumor_fq_r1_ch
         path('trimmed_fastq_r2') from trimmed_tumor_fq_r2_ch

     output:
         path("${sorted_bam}") into tumor_sorted_bam_ch

     script:
         sorted_bam = "${params.SampleId}_tumor_sorted.bam"
     """
     cat ${align_conf} | xargs \
     bwa mem -R '@RG\\tID:${params.SampleId}\\tSM:${params.SampleId}\\tLB:YUNYING\\tPL:Illumina' -t ${task.cpus} ${ref_fa} \
     trimmed_fastq_r1 trimmed_fastq_r2 | sambamba view \
     -l 0 -f bam -S -h /dev/stdin |sambamba sort -m 5G -t 8 -o ${sorted_bam} /dev/stdin

     ### alternatively, we would use samtools for sort bam file, but some restrictions occurs!
     ### argument "-Y" (or "-C") used by bwa does not support samtools sort directly, so we omit it (!!!)
     ### and some read with problemed coordinate sorted!
     # cat ${align_conf} | xargs \
     # bwa mem -R '@RG\\tID:${params.SampleId}\\tSM:${params.SampleId}\\tLB:YUNYING\\tPL:Illumina\\tPU:chip' \
     # -M -t ${task.cpus} ${ref_fa} trimmed_fastq_r1 \
     # trimmed_fastq_r2 | samtools sort -m 2G -t ${task.cpus} -O bam -o ${sorted_bam} /dev/stdin
     # samtools index ${sorted_bam}

     """
}


process control_align_bwa_sort {
     //container 'cmopipeline/bwa-samtools-gatk4-sambamba-samblaster:latest'
     container '10.168.2.67:5000/align_sambam_vcf:1.0'

     cpus params.threads
     memory { 20.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3
     
     when:
          params.ctl_fq1 != ""

     input:
         path('trimmed_fastq_r1') from trimmed_ctl_fq_r1_ch
         path('trimmed_fastq_r2') from trimmed_ctl_fq_r2_ch

     output:
         path("${sorted_bam}") into ctl_sorted_bam_ch

     script:
         sorted_bam = "${params.SampleId}_ctl_sorted.bam"
     """
     cat ${align_conf} | xargs \
     bwa mem -R '@RG\\tID:${params.SampleId}_X\\tSM:${params.SampleId}_X\\tLB:YUNYING\\tPL:Illumina' -t ${task.cpus} ${ref_fa} \
     trimmed_fastq_r1 trimmed_fastq_r2 | sambamba view \
     -l 0 -f bam -S -h /dev/stdin |sambamba sort -m 5G -t 8 -o ${sorted_bam} /dev/stdin

     ### alternatively, we would use samtools for sort bam file, but some restrictions occurs!
     ### argument "-Y" (or "-C") used by bwa does not support samtools sort directly, so we omit it (!!!)
     ### and some read with problemed coordinate sorted!
     # cat ${align_conf} | xargs \
     # bwa mem -R '@RG\\tID:${params.SampleId}\\tSM:${params.SampleId}\\tLB:YUNYING\\tPL:Illumina\\tPU:chip' \
     # -M -t ${task.cpus} ${ref_fa} trimmed_fastq_r1 \
     # trimmed_fastq_r2 | samtools sort -m 2G -t ${task.cpus} -O bam -o ${sorted_bam} /dev/stdin
     # samtools index ${sorted_bam}

     """
}

process refinement_gatk_elprep_tumor_standard {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'

     cpus params.threads/2
     memory { 20.GB * memory_scale1 * task.attempt }
     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5

     input:
          path "sorted.bam" from tumor_sorted_bam_ch
          path "target.bed" from bed_target

     output:
          path "${recal_bam}" into tumor_recal_bam_ch1, tumor_recal_bam_ch2, tumor_recal_bam_ch3, tumor_recal_bam_ch4, tumor_recal_bam_ch5, tumor_recal_bam_ch6, tumor_recal_bam_ch7, tumor_recal_bam_ch8, tumor_recal_bam_ch9, tumor_recal_bam_ch10
          path "${recal_bai}" into tumor_recal_bai_ch1, tumor_recal_bai_ch2, tumor_recal_bai_ch3, tumor_recal_bai_ch4, tumor_recal_bai_ch5, tumor_recal_bai_ch6, tumor_recal_bai_ch7, tumor_recal_bai_ch8, tumor_recal_bai_ch9, tumor_recal_bai_ch10

     script:
          recal_bam = "${params.SampleId}_tumor_recal.bam"
          recal_bai = "${params.SampleId}_tumor_recal.bai"
     if (params.remove_dup == true && params.recal == 'elprep')
        """
        elprep filter sorted.bam ${recal_bam} --mark-duplicates --mark-optical-duplicates output.metrics \
               --sorting-order keep --bqsr output.recal --bqsr-reference ${ref_elfa} --known-sites ${gatk_db_dir}/dbsnp_138.b37.elsites,${gatk_db_dir}/1000G_phase1.indels.b37.elsites,${gatk_db_dir}/Mills_and_1000G_gold_standard.indels.b37.elsites --nr-of-threads ${task.cpus}
        gatk BuildBamIndex -I ${recal_bam}
        """
     else if (params.remove_dup == false && params.recal == 'elprep')
        """
        elprep filter sorted.bam ${recal_bam} --sorting-order keep --bqsr output.recal \
               --bqsr-reference ${ref_elfa} --known-sites ${gatk_db_dir}/dbsnp_138.b37.elsites,${gatk_db_dir}/1000G_phase1.indels.b37.elsites,${gatk_db_dir}/Mills_and_1000G_gold_standard.indels.b37.elsites --nr-of-threads ${task.cpus}
        gatk BuildBamIndex -I ${recal_bam}
        """
     else if (params.remove_dup == false && params.recal == 'gatk')
        """
        cat ${recal_conf} | xargs \
        gatk BaseRecalibrator -I sorted.bam -R ${ref_fa} --known-sites ${gatk_db_dir}/dbsnp_138.b37.vcf \
             --known-sites ${gatk_db_dir}/1000G_phase1.indels.b37.vcf --known-sites ${gatk_db_dir}/Mills_and_1000G_gold_standard.indels.b37.vcf \
             -L target.bed -O recal_table

        gatk ApplyBQSR -R ${ref_fa} -I sorted.bam --bqsr-recal-file recal_table --create-output-bam-index true -O ${recal_bam}

        """
     else
        error "We does support other options for refinement process !!! >-<"
}


process refinement_gatk_elprep_control_standard {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'

     cpus params.threads/2
     memory { 10.GB * memory_scale2 * task.attempt }
     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5
     
     when:
          params.ctl_fq1 != ""

     input:
          path 'sorted.bam' from ctl_sorted_bam_ch
          path 'target.bed' from bed_target

     output:
          path "${recal_bam}" into ctl_recal_bam_ch1, ctl_recal_bam_ch2, ctl_recal_bam_ch4
          path "${recal_bai}" into ctl_recal_bai_ch1, ctl_recal_bai_ch2, ctl_recal_bai_ch4

     script:
          recal_bam = "${params.SampleId}_ctl_recal.bam"
          recal_bai = "${params.SampleId}_ctl_recal.bai"
     if (params.remove_dup == true && params.recal == 'elprep')
        """
        elprep filter sorted.bam ${recal_bam} --mark-duplicates --mark-optical-duplicates output.metrics \
               --sorting-order keep --bqsr output.recal --bqsr-reference ${ref_elfa} --known-sites ${gatk_db_dir}/dbsnp_138.b37.elsites,${gatk_db_dir}/1000G_phase1.indels.b37.elsites,${gatk_db_dir}/Mills_and_1000G_gold_standard.indels.b37.elsites --nr-of-threads ${task.cpus}
        gatk BuildBamIndex -I ${recal_bam}
        """
     else if (params.remove_dup == false && params.recal == 'elprep')
        """
        elprep filter sorted.bam ${recal_bam} --sorting-order keep --bqsr output.recal \
               --bqsr-reference ${ref_elfa} --known-sites ${gatk_db_dir}/dbsnp_138.b37.elsites,${gatk_db_dir}/1000G_phase1.indels.b37.elsites,${gatk_db_dir}/Mills_and_1000G_gold_standard.indels.b37.elsites --nr-of-threads ${task.cpus}
        gatk BuildBamIndex -I ${recal_bam}
        """
     else if (params.remove_dup == false && params.recal == 'gatk')
        """
        cat ${recal_conf} | xargs \
        gatk BaseRecalibrator -I sorted.bam -R ${ref_fa} --known-sites ${gatk_db_dir}/dbsnp_138.b37.vcf \
             --known-sites ${gatk_db_dir}/1000G_phase1.indels.b37.vcf --known-sites ${gatk_db_dir}/Mills_and_1000G_gold_standard.indels.b37.vcf \
             -L target.bed -O recal_table

        gatk ApplyBQSR -R ${ref_fa} -I sorted.bam --bqsr-recal-file recal_table --create-output-bam-index true -O ${recal_bam}

        """
     else
        error "We does support other options for refinement process !!! >-<"
}


process sv_svaba_paired {
     container '10.168.2.67:5000/svaba-yy:1.1.2'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus { params.threads/2 }
     memory { 16.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5
     
     when:
          params.ctl_fq1 != "" && sv_process_mark == true

     input:
          path 'tumor_recal.bam'  from tumor_recal_bam_ch2
          path 'tumor_recal.bai' from tumor_recal_bai_ch2
          path 'ctl_recal.bam'  from ctl_recal_bam_ch2
          path 'ctl_recal.bai' from ctl_recal_bai_ch2
          path 'target.bed' from sv_bed_ch

     output:
          path '*.unfiltered.somatic.sv.vcf' into sv_svaba_somatic_ch
          path '*.unfiltered.germline.sv.vcf' into sv_svaba_germline_ch
          path '*.svaba.somatic.indel.vcf' into indel_svaba_somatic_ch
          path '*.svaba.germline.indel.vcf' into indel_svaba_germline_ch

     script:
          prefix_svaba_vcf = "${params.SampleId}"
     """
     cat ${sv_svaba_conf} | xargs \
     svaba run -t tumor_recal.bam -n ctl_recal.bam -p ${task.cpus} -G ${ref_fa} -a `pwd`/${prefix_svaba_vcf} -k target.bed
     """
}

process sv_svaba_tumor_only {
     container '10.168.2.67:5000/svaba-yy:1.1.2'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus { params.threads/2 }
     memory { 16.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5
     
     when:
         params.ctl_fq1 == "" && sv_process_mark == true

     input:
          path 'tumor_recal.bam'  from tumor_recal_bam_ch2
          path 'tumor_recal.bai' from tumor_recal_bai_ch2
          path 'target.bed' from sv_bed_ch

     output:
          path '*.svaba.unfiltered.sv.vcf' into sv_svaba_tumor_ch
          path '*.svaba.unfiltered.indel.vcf' into indel_svaba_tumor_ch
          path '*.svaba.indel.vcf' into indel_svaba_ch

     script:
          prefix_svaba_vcf = "${params.SampleId}"
     """
     cat ${sv_svaba_conf} | xargs \
     svaba run -t tumor_recal.bam -p ${task.cpus} -G ${ref_fa} -a `pwd`/${prefix_svaba_vcf} -k target.bed
     """
}

process cnv_oncocnv_tumor_controls {
     container '10.168.2.67:5000/oncocnv-yy:6.9m'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*.{summary.txt,profile.txt,profile.png}'
     publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.stats.txt'

     cpus 2
     memory { 8.GB }

     errorStrategy { task.exitStatus in 2..140 ? 'retry' : 'ignore' }
     maxRetries 5

     input:
          path "${params.SampleId}_tumor_recal.bam"  from tumor_recal_bam_ch3
          path "${params.SampleId}_tumor_recal.bai" from tumor_recal_bai_ch3
          path 'target.bed' from bed_target

     when:
          params.baseline =~ /^ZB.*/

     output:
          path '*' into cnv_oncocnv_ch

     script:
          prefix_cnv = "${params.SampleId}"
     """
     IFS=','; read -ra clades <<< "${params.baseline}"
     for clade in "\${clades[@]}"
     do
         #get normalized read counts for the tumor samples
         perl /opt/ONCOCNV_getCounts.pl getSampleStats -m Ampli -c ${oncocnv_db_dir}/\${clade}.Control.stats.txt -s ${params.SampleId}_tumor_recal.bam -o `pwd`/${prefix_cnv}.\${clade}.Test.stats.txt 
   
         #process test samples and predict CNA and CNVs:
         cat /opt/processSamples.v6.4.R | R --slave --args `pwd`/${prefix_cnv}.\${clade}.Test.stats.txt ${oncocnv_db_dir}/\${clade}.Control.stats.Processed.txt `pwd`/${prefix_cnv}.\${clade}.Test.output.txt cghseg
 
     done

     """
}

process cnv_oncocnv_tumor_covers {
     container '10.168.2.67:5000/oncocnv-yy:6.9m'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*.{summary.txt,profile.txt,profile.png}'

     cpus 2
     memory { 8.GB }

     errorStrategy { task.exitStatus in 2..140 ? 'retry' : 'ignore' }
     maxRetries 2

     when:
          params.SampleId =~ /Q_g_639/

     input:
          path "${params.SampleId}_tumor_recal.bam"  from tumor_recal_bam_ch9
          path "${params.SampleId}_tumor_recal.bai" from tumor_recal_bai_ch9
          val cover from Channel.from("cnv_cover_5", "cnv_cover_10", "cnv_cover_20")

     output:
          path '*summary.txt' into oncocnv_cover_summary
          path '*profile.*' into oncocnv_cover_profile

     script:
          prefix_cnv = "${params.SampleId}"
     """
     perl /opt/ONCOCNV_getCounts.pl getSampleStats -m Ampli -c ${oncocnv_db_dir}/panel_639_cover/${cover}.Control.stats.txt -s ${params.SampleId}_tumor_recal.bam -i ${cover}. -o `pwd`/${prefix_cnv}.${cover}.Test.stats.txt
     cat /opt/processSamples.v6.4.R | R --slave --args `pwd`/${prefix_cnv}.${cover}.Test.stats.txt ${oncocnv_db_dir}/panel_639_cover/${cover}.Control.stats.Processed.txt `pwd`/${prefix_cnv}.${cover}.Test.output.txt cghseg
     """
}

process refinement_enzyme_correct_control {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'
     publishDir "${unified_bam_dir}", mode: 'copy', pattern: '*.{bam,bai}'
 
     cpus 2
     memory { 8.GB * task.attempt }
 
     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 5
     
     when:
          params.ctl_fq1 != ""
 
     input:
          path 'ctl_recal.bam' from ctl_recal_bam_ch1
          path 'ctl_recal.bai' from ctl_recal_bai_ch1
          path 'target.bed' from bed_target
          path 'chemo_radio_rs.bed' from chemo_radio_rs_bed
 
     output:
          path "${ctl_correct_bam}" into control_refined_bam_ch1, control_refined_bam_ch2
          path "${ctl_correct_bai}" into control_refined_bai_ch1, control_refined_bai_ch2
          path "merged.interval_list" into merged_region_ch
 
     script:
          ctl_correct_bam = "${params.SampleId}_ctl_recal.bam"
          ctl_correct_bai = "${params.SampleId}_ctl_recal.bai"

     if (params.SampleId =~ "PCR_")
        """
        cp -d ctl_recal.bam ${ctl_correct_bam}
        cp -d ctl_recal.bai ${ctl_correct_bai}
        """
     else
        """
        java -jar /gatk/gatk-package-4.1.7.0-local.jar BedToIntervalList -I target.bed -O germline.interval_list -SD ${ref_fa_dict}
        java -jar /gatk/gatk-package-4.1.7.0-local.jar BedToIntervalList -I chemo_radio_rs.bed -O chemo_radio_therapy.interval_list -SD ${ref_fa_dict}
        java -jar /gatk/gatk-package-4.1.7.0-local.jar IntervalListTools --ACTION UNION -I germline.interval_list -I chemo_radio_therapy.interval_list -O merged.interval_list

        recsc -r ${ref_fa} -g merged.interval_list -i ctl_recal.bam -o ${ctl_correct_bam}
        """
}

process germline_gatk4_for_control {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'

     cpus 2
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5
     
     when:
          params.ctl_fq1 != ""

     input:
          path 'ctl_recal.bam' from control_refined_bam_ch2
          path 'ctl_recal.bai' from control_refined_bai_ch2
          path 'merged.interval_list' from merged_region_ch

     output:
          path "${germline_vcf}" into control_germline_ch

     script:
          germline_vcf = "${params.SampleId}_germline_for_control.vcf"
     """
     gatk HaplotypeCaller -I ctl_recal.bam -R ${ref_fa} \
             -D ${gatk_db_dir}/dbsnp_138.b37.vcf --alleles ${chemo_radio_rs_vcf} \
             -L merged.interval_list -stand-call-conf 10 \
             --interval-padding 20 -A AlleleFraction -O ${germline_vcf}
     """
}

process statistic_sinotools_paired {
     container '10.168.2.67:5000/sinotools-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*_qc.txt'
     publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.cov'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3
     
     when:
          params.ctl_fq1 != ""

     input:
          path "${params.SampleId}_tumor_recal.bam"  from tumor_recal_bam_ch4
          path "${params.SampleId}_tumor_recal.bai" from tumor_recal_bai_ch4
          path "${params.SampleId}_ctl_recal.bam"  from ctl_recal_bam_ch4
          path "${params.SampleId}_ctl_recal.bai" from ctl_recal_bai_ch4
          path 'target.bed' from st_bed_ch

     output:
          path "*_qc.txt" into statistic_paired_ch
          path "*.cov" into coverage_paired_ch

     script:
         tumor_statistics = "${params.SampleId}_tumor_qc.txt"
         ctl_statistics = "${params.SampleId}_ctl_qc.txt"
     """
     sinotools bam_qc -i ${params.SampleId}_tumor_recal.bam -bed target.bed -hotspot ${core_hotspots_bed} -o ${tumor_statistics}
     sinotools bam_qc -i ${params.SampleId}_ctl_recal.bam -bed target.bed -hotspot ${core_hotspots_bed} -o ${ctl_statistics}
     """
}

process statistic_sinotools_tumor_only {
     container '10.168.2.67:5000/sinotools-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*_qc.txt'
     publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.cov'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3
     
     when:
          params.ctl_fq1 == ""

     input:
          path "${params.SampleId}_tumor_recal.bam"  from tumor_recal_bam_ch4
          path "${params.SampleId}_tumor_recal.bai" from tumor_recal_bai_ch4
          path "target.bed" from st_bed_ch

     output:
          path "*_qc.txt" into statistic_tumor_ch
          path "*.cov" into coverage_tumor_ch

     script:
         tumor_statistics = "${params.SampleId}_tumor_qc.txt"
     """
     sinotools bam_qc -i ${params.SampleId}_tumor_recal.bam -bed target.bed -hotspot ${core_hotspots_bed} -o ${tumor_statistics}
     """
}


process somatic_gatk_interval_split {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'
     
     cpus 1
     memory { 2.GB * task.attempt }
 
     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3
 
     input:
          path 'target.bed' from bed_target
 
     output:
          path 'split_interval/*' into split_bed_ch
     
     script:
     """
     gatk SplitIntervals -R $ref_fa -L target.bed -scatter ${split_region_count} -O split_interval
     """
}

process refinement_enzyme_correct_tumor {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'
     publishDir "${unified_bam_dir}", mode: 'copy'
 
     cpus 2
     memory { 8.GB * task.attempt }
 
     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 5
 
     input:
           path "tumor_recal.bam" from tumor_recal_bam_ch1
           path "tumor_recal.bai" from tumor_recal_bai_ch1
           path "target.bed" from bed_target
 
     output:
           path "${tumor_correct_bam}" into tumor_correct_bam_ch1, tumor_correct_bam_ch2
           path "${tumor_correct_bai}" into tumor_correct_bai_ch1, tumor_correct_bai_ch2
 
     script:
           tumor_correct_bam = "${params.SampleId}_tumor_recal.bam"
           tumor_correct_bai = "${params.SampleId}_tumor_recal.bai"

     if (params.SampleId =~ "PCR_")
        """
        cp -d tumor_recal.bam ${tumor_correct_bam}
        cp -d tumor_recal.bai ${tumor_correct_bai}
        """
     else
        """
        recsc -r ${ref_fa} -g target.bed -i tumor_recal.bam -o ${tumor_correct_bam}
        """
}

process germline_gatk4_for_tumor {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'

     cpus 2
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5
     
     when:
          params.ctl_fq1 == ""

     input:
          path "tumor_refined.bam" from tumor_correct_bam_ch2
          path "tumor_refined.bai" from tumor_correct_bai_ch2
          path "target.bed" from bed_target
          path "chemo_radio_rs.bed" from chemo_radio_rs_bed

     output:
          path "${tumor_germline_vcf}" into tumor_germline_ch

     script:
          tumor_germline_vcf = "${params.SampleId}_germline_for_tumor.vcf"
     """
     java -jar /gatk/gatk-package-4.2.2.0-local.jar BedToIntervalList -I target.bed -O germline.interval_list -SD ${ref_fa_dict}
     java -jar /gatk/gatk-package-4.2.2.0-local.jar BedToIntervalList -I chemo_radio_rs.bed -O chemo_radio_therapy.interval_list -SD ${ref_fa_dict}
     java -jar /gatk/gatk-package-4.2.2.0-local.jar IntervalListTools --ACTION UNION -I germline.interval_list -I chemo_radio_therapy.interval_list -O merged.interval_list

     gatk HaplotypeCaller -I tumor_refined.bam -R ${ref_fa} \
             -D ${gatk_db_dir}/dbsnp_138.b37.vcf --alleles ${chemo_radio_rs_vcf} \
             -L merged.interval_list -stand-call-conf 10 \
             --interval-padding 20 -A AlleleFraction -O ${tumor_germline_vcf}
     """
}

process somatic_gatk_mutect2_parallel_paired {
    container '10.168.2.67:5000/gatk-yy:4.2.2.1'

    cpus 2
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 123..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    when:
          params.ctl_fq1 != ""

    input:
          path 'tumor_refined.bam' from tumor_correct_bam_ch1
          path 'tumor_refined.bai' from tumor_correct_bai_ch1
          path 'ctl_recal.bam' from control_refined_bam_ch1 
          path 'ctl_recal.bai' from control_refined_bai_ch1
          path 'splited.interval_list' from split_bed_ch.flatten()

    output:
          path "${raw_vcf}" into raw_vcf_paired_ch
          path "${vcf_stats}" into vcf_stat_paired_ch
          path "${f1r2_stat}" into f1r2_stat_paired_ch
          path "${pileup_stat}" into pileup_stat_paired_ch
    
    script:
         raw_vcf = "${params.SampleId}_mutect2_raw.vcf"
         vcf_stats = "${params.SampleId}_mutect2.vcf.stats"
         f1r2_stat = "${params.SampleId}_f1r2.tar.gz"
         pileup_stat = "${params.SampleId}_getpileupsummaries.table"
    """
    cat ${mutect2_conf} | xargs \
    gatk Mutect2 -I tumor_refined.bam -I ctl_recal.bam -normal "${params.SampleId}_X" -R ${ref_fa} -alleles ${core_hotspots_conf} \
                 -germline-resource ${gatk_db_dir}/b37_af-only-gnomad.raw.sites.vcf -pon ${gatk_db_dir}/b37_yunying_pon.vcf \
                 --f1r2-tar-gz ${f1r2_stat} -L splited.interval_list --native-pair-hmm-threads ${task.cpus} -O "${params.SampleId}_mutect2.vcf"

    cat ${secondary_annotate_conf} | xargs \
    gatk VariantAnnotator -I tumor_refined.bam -R ${ref_fa} -V "${params.SampleId}_mutect2.vcf" -O ${raw_vcf}

    gatk GetPileupSummaries -I tumor_refined.bam -V ${gatk_db_dir}/b37_small_exac_common_3.vcf \
                 --interval-set-rule UNION -L splited.interval_list -O ${pileup_stat}

    """
}


process somatic_gatk_mutect2_tumor_only {
    container '10.168.2.67:5000/gatk-yy:4.2.2.1'

    cpus 2
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 123..140 ? 'retry' : 'terminate' }
    maxRetries 5

    when:
          params.ctl_fq1 == ""

    input:
          path 'tumor_refined.bam' from tumor_correct_bam_ch1
          path 'tumor_refined.bai' from tumor_correct_bai_ch1
          path 'splited.interval_list' from split_bed_ch.flatten()

    output:
          path "${raw_vcf}" into raw_vcf_tumor_ch
          path "${vcf_stats}" into vcf_stat_tumor_ch
          path "${f1r2_stat}" into f1r2_stat_tumor_ch
          path "${pileup_stat}" into pileup_stat_tumor_ch
    
    script:
         raw_vcf = "${params.SampleId}_mutect2_raw.vcf"
         vcf_stats = "${params.SampleId}_mutect2.vcf.stats"
         f1r2_stat = "${params.SampleId}_f1r2.tar.gz"
         pileup_stat = "${params.SampleId}_getpileupsummaries.table"
    """
    cat ${mutect2_conf} | xargs \
    gatk Mutect2 -I tumor_refined.bam -R ${ref_fa} -alleles ${core_hotspots_conf} \
                 -germline-resource ${gatk_db_dir}/b37_af-only-gnomad.raw.sites.vcf -pon ${gatk_db_dir}/b37_yunying_pon.vcf \
                 --f1r2-tar-gz ${f1r2_stat} -L splited.interval_list --native-pair-hmm-threads ${task.cpus} -O "${params.SampleId}_mutect2.vcf"

    cat ${secondary_annotate_conf} | xargs \
    gatk VariantAnnotator -I tumor_refined.bam -R ${ref_fa} -V "${params.SampleId}_mutect2.vcf" -O ${raw_vcf}

    gatk GetPileupSummaries -I tumor_refined.bam -V ${gatk_db_dir}/b37_small_exac_common_3.vcf \
                 --interval-set-rule UNION -L splited.interval_list -O ${pileup_stat}

    """
}

process somatic_gatk_filter_mark_paired {
    container '10.168.2.67:5000/gatk-yy:4.2.0.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy'
    // *** this mount volume is used for local executation for docker, not necessary for kubernetes!!!
    containerOptions '-v /home/yunying/:/home/yunying/'

    cpus 2
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    when:
          params.ctl_fq1 != ""

    input:
         val vcf_items from raw_vcf_paired_ch.reduce{ previous, item -> previous + " -I " + item }
         val stats_items from vcf_stat_paired_ch.reduce{ previous, item -> previous + " -stats " + item }
         val f1r2_items from f1r2_stat_paired_ch.reduce{ previous, item -> previous + " -I " + item }
         val pileup_items from pileup_stat_paired_ch.reduce{ previous, item -> previous + " -I " + item }

    output:
         path "${filter_marked_vcf}" into filter_marked_paired_ch
    
    script:
         filter_marked_vcf = "${params.SampleId}_mutect2_filter_marked.vcf"
    """
    gatk MergeVcfs -I ${vcf_items} -O mutect2_merged.vcf
    gatk MergeMutectStats -stats ${stats_items} -O mutect2_merged.vcf.stats
    gatk LearnReadOrientationModel -I ${f1r2_items} -O read-orientation-model.tar.gz
    gatk GatherPileupSummaries -I ${pileup_items} --sequence-dictionary ${ref_fa_dict} -O getpileupsummaries.table
    
    gatk CalculateContamination -I getpileupsummaries.table \
           -tumor-segmentation segments.table \
           -O contamination.table

    cat ${somatic_filter_conf} | xargs \
    gatk FilterMutectCalls -R ${ref_fa} -V mutect2_merged.vcf --tumor-segmentation segments.table \
           --contamination-table contamination.table --ob-priors read-orientation-model.tar.gz \
           -O ${filter_marked_vcf}
    """
}

process somatic_gatk_filter_mark_tumor_only {
    container '10.168.2.67:5000/gatk-yy:4.2.0.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy'
    // *** this mount volume is used for local executation for docker, not necessary for kubernetes!!!
    containerOptions '-v /home/yunying/:/home/yunying/'

    cpus 2
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    when:
          params.ctl_fq1 == ""

    input:
         val vcf_items from raw_vcf_tumor_ch.reduce{ previous, item -> previous + " -I " + item }
         val stats_items from vcf_stat_tumor_ch.reduce{ previous, item -> previous + " -stats " + item }
         val f1r2_items from f1r2_stat_tumor_ch.reduce{ previous, item -> previous + " -I " + item }
         val pileup_items from pileup_stat_tumor_ch.reduce{ previous, item -> previous + " -I " + item }

    output:
         path "${filter_marked_vcf}" into filter_marked_tumor_ch
    
    script:
         filter_marked_vcf = "${params.SampleId}_mutect2_filter_marked.vcf"
    """
    gatk MergeVcfs -I ${vcf_items} -O mutect2_merged.vcf
    gatk MergeMutectStats -stats ${stats_items} -O mutect2_merged.vcf.stats
    gatk LearnReadOrientationModel -I ${f1r2_items} -O read-orientation-model.tar.gz
    gatk GatherPileupSummaries -I ${pileup_items} --sequence-dictionary ${ref_fa_dict} -O getpileupsummaries.table
    
    gatk CalculateContamination -I getpileupsummaries.table \
           -tumor-segmentation segments.table \
           -O contamination.table

    cat ${somatic_filter_conf} | xargs \
    gatk FilterMutectCalls -R ${ref_fa} -V mutect2_merged.vcf --tumor-segmentation segments.table \
           --contamination-table contamination.table --ob-priors read-orientation-model.tar.gz \
           -O ${filter_marked_vcf}
    """
}

process anno_snpeff_hgvs_paired {
    container '10.168.2.67:5000/snpeff-yy:4.4m'

    cpus 4
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 3

    when:
          params.ctl_fq1 != ""
    
    input:
         path 'all_marked.vcf' from filter_marked_paired_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

    output:
         path '*.anno' into haplotype_anno_paired_ch
         path "${prefix_anno_vcf}_hgvs.vcf" into hgvs_anno_paired_ch

    script:
        prefix_anno_vcf = "${params.SampleId}_filter_marked"
    """
    python3 /opt/vcf_multiple_split.py all_marked.vcf all_marked_splitted.vcf
    cat anno_snpeff.conf | xargs \
    java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}.haplotype.anno GRCh37.p13.RefSeq all_marked_splitted.vcf > ${prefix_anno_vcf}_hgvs.vcf
    """
}

process anno_snpsift_database_paired {
     container '10.168.2.67:5000/snpsift-yy:1.0'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
         params.ctl_fq1 != ""
       
     input:
         path 'hgvs.vcf' from hgvs_anno_paired_ch

     output:
         path "${paired_anno_vcf}" into paired_anno_hgvs_db_ch

     script:
         paired_anno_vcf = "${params.SampleId}_filter_marked_hgvs_db.vcf"
     """
     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_vcf_conf} hgvs.vcf > ${paired_anno_vcf}
     """
}

process anno_snpeff_hgvs_tumor_only {
    container '10.168.2.67:5000/snpeff-yy:4.4m'

    cpus 4
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    when:
         params.ctl_fq1 == ""

    input:
         path 'all_marked.vcf' from filter_marked_tumor_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

    output:
         path '*.anno' into haplotype_anno_tumor_ch
         path "${prefix_anno_vcf}_hgvs.vcf" into hgvs_anno_tumor_ch

    script:
         prefix_anno_vcf = "${params.SampleId}_filter_marked"
    """
    python3 /opt/vcf_multiple_split.py all_marked.vcf all_marked_splitted.vcf
    cat anno_snpeff.conf | xargs \
    java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}.haplotype.anno GRCh37.p13.RefSeq all_marked_splitted.vcf > ${prefix_anno_vcf}_hgvs.vcf
    """
}

process anno_snpsift_database_tumor_only {
     container '10.168.2.67:5000/snpsift-yy:1.0'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
         params.ctl_fq1 == ""
       
     input:
         path 'hgvs.vcf' from hgvs_anno_tumor_ch

     output:
         path "${tumor_anno_vcf}" into tumor_anno_hgvs_db_ch

     script:
         tumor_anno_vcf = "${params.SampleId}_filter_marked_hgvs_db.vcf"
     """
     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_vcf_conf} hgvs.vcf > ${tumor_anno_vcf}
     """
}

process mutect2_report_paired {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 != ""

     input:
         path 'tumor_recal.bam'  from tumor_recal_bam_ch8
         path 'tumor_recal.bai' from tumor_recal_bai_ch8
         path 'marked_all_anno.vcf' from paired_anno_hgvs_db_ch
         path 'haplotype.anno' from haplotype_anno_paired_ch

     output:
         path '*.xls' into small_var_report_paired

     script:
         report_output = "${params.SampleId}_review_for_report.xls"
     """
     python3 /opt/vcf2mut_report_review.py marked_all_anno.vcf haplotype.anno ${variation_hotspots_conf} all_sites.mpileup ${report_output}
     """
}

process mutect2_report_tumor_only {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 2
     
     when:
         params.ctl_fq1 == ""

     input:
         path 'tumor_recal.bam'  from tumor_recal_bam_ch8
         path 'tumor_recal.bai' from tumor_recal_bai_ch8
         path 'marked_all_anno.vcf' from tumor_anno_hgvs_db_ch
         path 'haplotype.anno' from haplotype_anno_tumor_ch

     output:
         path '*.xls' into small_var_report

     script:
         report_output = "${params.SampleId}_review_for_report.xls"
     """
     python3 /opt/vcf2mut_report_review.py marked_all_anno.vcf haplotype.anno ${variation_hotspots_conf} all_sites.mpileup ${report_output}
     """
}

process filter_svaba_review_paired {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
     maxRetries 3
     
     when:
         params.ctl_fq1 != ""

     input:
         path 'somatic_sv.vcf' from sv_svaba_somatic_ch
         path 'germline_sv.vcf' from sv_svaba_germline_ch

     output:
         path '*' into sv_paired_report_ch

     script:
         sv_report_output = "${params.SampleId}_sv_report.xls"
     """
     python3 /opt/sv2report_review.py somatic_sv.vcf germline_sv.vcf ${sv_hotspots_conf} ${sv_partner_conf} ${sv_report_output}
     """
}

process filter_svaba_review_tumor_only {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
     maxRetries 3

     when:
         params.ctl_fq1 == ""

     input:
         path 'unfiltered_sv.vcf' from sv_svaba_tumor_ch

     output:
         path '*' into sv_tumor_report_ch

     script:
         sv_report_output = "${params.SampleId}_sv_report.xls"
     """
     python3 /opt/sv2report_review.py unfiltered_sv.vcf ${sv_hotspots_conf} ${sv_partner_conf} ${sv_report_output}
     """
}

process msi_pipline_tumor_controls {
     container '10.168.2.67:5000/msi_pipeline-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*MSI_result*'
     publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*model_*'
	
     cpus 2
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     input:
          path 'recal.bam'  from tumor_recal_bam_ch5
          path 'recal.bai' from tumor_recal_bai_ch5

     when:
          params.msi_db != ""
     
     output:
          path '*xls' into msi_result_ch

     script:
          prefix_msi = "${params.SampleId}"
          set_depth_num = "${params.msi_depth}"
     """
     python3 /opt/integration_model_123_version3.py recal.bam ${prefix_msi} ${prefix_msi} ${msi_database_dir} ${set_depth_num}
     """
}

process perbase_complex_standard {
     container '10.168.2.67:5000/perbase-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 2
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
     maxRetries 2

     input:
         path 'realigned.bam' from tumor_recal_bam_ch10
         path 'realigned.bai' from tumor_recal_bai_ch10

     output:
         path "${perbase_out_txt}" into perbase_out_ch

     script:
         perbase_out_txt = "${params.SampleId}_perbase_complex_tips.txt"
     """
     cat ${perbase_conf} | xargs \
     perbase base-depth realigned.bam -b ${complex_hotspot_bed} -r ${ref_fa} -t ${task.cpus} -o ${perbase_out_txt}
     """
}


process Calc_TMB_tumor_controls {
     container '10.168.2.67:5000/calc2tmb-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
          params.codinglength != "" && params.ctl_fq1 != "" 
       
     input:
         path 'report.xls' from small_var_report_paired

     output:
         path '*' into tmb_report_paired

     script:
         report_muts_output = "${params.SampleId}_review_muts_TMB_for_report.txt"
         report_tmb_output = "${params.SampleId}_review_TMB_for_report.txt"
     """
     python3 /opt/calc_TMB_match_normal.py report.xls ${params.codinglength} ${report_muts_output} ${report_tmb_output}
     """
}
 
process Calc_TMB_tumor_only {
     container '10.168.2.67:5000/calc2tmb-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
         params.codinglength != "" && params.ctl_fq1 == ""
       
     input:
         path 'report.xls' from small_var_report

     output:
         path '*' into tmb_report_ch

     script:
         report_muts_output = "${params.SampleId}_review_muts_TMB_for_report.txt"
         report_tmb_output = "${params.SampleId}_review_TMB_for_report.txt"
     """
     python3 /opt/calc_TMB_tumor_only.py report.xls ${params.codinglength} ${report_muts_output} ${report_tmb_output}
     """
}

process anno_germline_snpeff_control {
    container '10.168.2.67:5000/snpeff-yy:4.4m'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'
    publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*.xls'

    cpus 4
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
          params.ctl_fq1 != ""
    
    input:
         path 'germline_control.vcf' from control_germline_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

    output:
         path "${germline_anno_vcf}" into germline_anno_ch
         path '*.xls' into germline_chemo_report_ch

    script:
        germline_anno_vcf = "${params.SampleId}_germline.vcf"
        germline_chemo_xls = "${params.SampleId}_chemotherapy_rs_germline.xls"
        germline_radio_xls = "${params.SampleId}_radiotherapy_rs.xls"
    """
    python3 /opt/vcf_multiple_split.py germline_control.vcf germline_control_splitted.vcf
    cat anno_snpeff.conf | xargs \
    java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
                     -haplotype-ouput haplotype.anno GRCh37.p13.RefSeq germline_control_splitted.vcf > ${germline_anno_vcf}
    python3 /opt/chemotherapy_to_report.py germline_control.vcf ${chemotherapy_rs_vcf} ${germline_chemo_xls} 
    python3 /opt/chemotherapy_to_report.py germline_control.vcf ${radiotherapy_rs_vcf} ${germline_radio_xls} 
    """
}


process anno_germline_snpeff_tumor {
    container '10.168.2.67:5000/snpeff-yy:4.4m'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'
    publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*.xls'

    cpus 4
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
          params.ctl_fq1 == ""
    
    input:
         path 'germline_tumor.vcf' from tumor_germline_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

    output:
         path "${germline_tumor_anno_vcf}" into tumor_germline_anno_ch, tumor_germline_anno_ch2
         path '*.xls' into tumor_chemo_report_ch

    script:
        germline_tumor_anno_vcf = "${params.SampleId}_germline_tumor.vcf"
        tumor_chemo_xls = "${params.SampleId}_chemotherapy_rs_tumor.xls"
        tumor_radio_xls = "${params.SampleId}_radiotherapy_rs.xls"
    """
    python3 /opt/vcf_multiple_split.py germline_tumor.vcf germline_tumor_splitted.vcf
    cat anno_snpeff.conf | xargs \
    java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
                     -haplotype-ouput haplotype.anno GRCh37.p13.RefSeq germline_tumor_splitted.vcf > ${germline_tumor_anno_vcf}
    python3 /opt/chemotherapy_to_report.py germline_tumor.vcf ${chemotherapy_rs_vcf} ${tumor_chemo_xls}
    python3 /opt/chemotherapy_to_report.py germline_tumor.vcf ${radiotherapy_rs_vcf} ${tumor_radio_xls}
    """
}


process germline_clinvar_filter_paired {
     container '10.168.2.67:5000/snpsift-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
         params.ctl_fq1 != "" 
       
     input:
         path 'germline.vcf' from germline_anno_ch

     output:
         path '*.xls' into germline_clinvar_ch

     script:
         clinvar_anno_vcf = "${params.SampleId}_germline_clinvar.vcf"
         clinvar_anno_xls = "${params.SampleId}_germline_clinvar.xls"
     """
     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_vcf_conf} germline.vcf > ${clinvar_anno_vcf}
     python3 /opt/germline_clinvar_filter.py ${clinvar_anno_vcf} ${clinvar_anno_xls}
     """
}

process germline_clinvar_tumor_only {
     container '10.168.2.67:5000/snpsift-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
         params.ctl_fq1 == ""
       
     input:
         path 'pesudo_germline.vcf' from tumor_germline_anno_ch

     output:
         path '*.xls' into tumor_clinvar_ch

     script:
         clinvar_anno_vcf = "${params.SampleId}_tumor_clinvar.vcf"
         clinvar_anno_xls = "${params.SampleId}_tumor_clinvar.xls"
     """
     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_vcf_conf} pesudo_germline.vcf > ${clinvar_anno_vcf}
     python3 /opt/germline_clinvar_filter.py ${clinvar_anno_vcf} ${clinvar_anno_xls}
     """
}
//
//process tumor_clinvar_pesudo_germline {
//     container '10.168.2.67:5000/snpsift-yy:1.0'
//     publishDir "${filtered_review_dir}", mode: 'copy'
//
//     cpus 2
//     memory { 4.GB * task.attempt }
//
//     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
//     maxRetries 3
//
//     when:
//         params.ctl_fq1 == "" && params.SampleId =~ /_g_38/
//       
//     input:
//         path 'pesudo_germline2.vcf' from tumor_germline_anno_ch2
//
//     output:
//         path '*.xls' into tumor_clinvar_ch2
//
//     script:
//         clinvar_anno_vcf = "${params.SampleId}_tumor_clinvar2.vcf"
//         clinvar_anno_xls = "${params.SampleId}_tumor_clinvar2.xls"
//     """
//     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_vcf_conf} pesudo_germline2.vcf > ${clinvar_anno_vcf}
//     python3 /opt/germline_clinvar_filter.py ${clinvar_anno_vcf} ${clinvar_anno_xls}
//     """
//}

workflow.onComplete {
    def send_msg_success = false
    def send_retry_counts = 0
    while (!send_msg_success && send_retry_counts < 10) {
        def proc = "python3 ${workflow.projectDir}/bin/workflow_rabbitmq_pika.py ${chip_batch_id} ${params.SampleId} ${workflow.success ? 1 : 0} ${params.queue_name}".execute()
        proc.waitFor()
        send_retry_counts += 1
        if (proc.exitValue() == 0)
            send_msg_success = true
        else
            sleep(3000)
        log.info("Send rabbitmq queue with status: ${proc.exitValue()}!")
        def update_process = "python3 ${workflow.projectDir}/bin/update_sample_record.py ${chip_batch_id} ${params.SampleId} ${params.update_db}".execute()
        update_process.waitFor()
        log.info("Update history record with status: ${update_process.exitValue()}!")
    }
    state_flag_file = file("${params.output}/success")
    state_flag_file.text = (workflow.success ? "1\n" : "0\n")
    log.info("Sample (${params.SampleId}) execution status: ${workflow.success ? 'OK' : 'failed'}")
    log.info("Workflow completed at: ${workflow.complete}")
}

