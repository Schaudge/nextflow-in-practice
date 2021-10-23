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
  log.info "nextflow run somatic_vardict_pindel.nf -c nextflow_basic.conf -profile docker --fq1 R1.fq.gz --fq2 R2.fq.gz --ctl_fq1 X_R1.fq.gz --ctl_fq2 X_R2.fq.gz --bed target_bed_file"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --SampleId         STRING             Sample name (dafault: fastq1 basename)"
  log.info "    --PanelName        STRING             Panel name, used for exported gene filter (required)"
  log.info "    --fq1              FILE               tumor fastq(.gz) file for read1 (required)"
  log.info "    --fq2              FILE               tumor fastq(.gz) file for read2 (required)"
  log.info "    --ctl_fq1          FILE               control fastq(.gz) file for read1 (required for paired sample)"
  log.info "    --ctl_fq2          FILE               control fastq(.gz) file for read2 (required for paired sample)"
  log.info "    --bed              FILE               target_bed_file (required)"
  log.info "    --baseline         STRING             ONCOCNV control clade flag (optional, for cnv is required)"
  log.info "    --msi_db           STRING             msi database subdir in related to directory /yunying/ref/... (optional, for msi is required)"
  log.info "    --msi_depth        STRING             msi locus depth threshold (optional, for msi is required)"
  log.info "    --output           DIR                Output directory (default: /yunying/data/product/)"
  log.info "    --RG               STRING             Read group tag (dafault: fastq1 basename)"
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
params.threads = 16
params.remove_dup = true
params.pindel_mark = true
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
pcr_panel_mark = params.SampleId =~ /PCR_/  ? true : false
somatic_process_mark = (params.bed =~ /H.bed$/ || params.bed =~ /MSI.bed$/ || params.bed =~ /M.bed$/) ? false : true
st_bed_ch = params.st_bed ? Channel.value("$params.st_bed") : Channel.value("$params.bed")
sv_bed_ch = params.sv_bed ? Channel.value("$params.sv_bed") : Channel.value("$params.bed")
sv_process_mark = params.sv_bed =~ /null.bed$/ ? false : true
pindel_process_mark = somatic_process_mark ? true : false
germline_process_mark = (params.bed =~ /M.bed$/ || params.bed =~ /MSI.bed$/) ? false : true
dedup_process_mark = (params.SampleId =~ /ct_447$/ || params.SampleId =~ /ct_580$/ || params.SampleId =~ /ct_113$/) ? true : false

memory_scale1 = Math.ceil(tumor_fq1.size()/1024/1024/1024 + tumor_fq2.size()/1024/1024/1024)
memory_scale2 = Math.ceil(ctl_fq1.size()/1024/1024/1024 + ctl_fq2.size()/1024/1024/1024)
println "The memory scale factor is ${memory_scale1}, ---, ${memory_scale2} * 15G for memory configuration!"
split_region_count = params.threads/2

log.info "Sample ${params.SampleId} was running ..."
log.info "And the result files will be put into directory ${filtered_review_dir}."

// alignment database and gatk bundle database
def genome_dir = "/yunying/ref/human/b37/"
def ref_fa = "${genome_dir}b37_Homo_sapiens_assembly19.fasta"
ref_fa = params.SampleId =~ /_M$/ ? "${genome_dir}mgmt/mgmt.fa" : ref_fa
def ref_elfa = "${genome_dir}b37_Homo_sapiens_assembly19.elfasta"
ref_elfa = params.SampleId =~ /_M$/ ? "${genome_dir}mgmt/mgmt.elfa" : ref_elfa
def gatk_db_dir = "${genome_dir}gatk"
def oncocnv_db_dir = "${genome_dir}oncocnv"
def snpeff_db_dir = "${genome_dir}snpeff"
def realignment_bed = "${gatk_db_dir}/core_realignment.bed"
def chemotherapy_rs_vcf = file("${genome_dir}variation_hotspots/chemotherapy_rs_v1.0.vcf")
def chemotherapy_rs_bed = file("${genome_dir}variation_hotspots/chemotherapy_rs_v1.0.bed")
def radiotherapy_rs_vcf = file("${genome_dir}variation_hotspots/radiotherapy_rs_v1.0.vcf")
def chemo_radio_rs_vcf = file("${genome_dir}variation_hotspots/chemotherapy_radiotherapy_rs_v1.0.vcf")
def chemo_radio_rs_bed = file("${genome_dir}variation_hotspots/chemotherapy_radiotherapy_rs_v1.0.bed")
def rs_1p19q_bed = file("${genome_dir}mgmt/1p19q.rs.bed")
def cosmic_database = "${genome_dir}snpsift/CosmicCodingMuts.vcf.gz"
def clinvar_database = "${genome_dir}snpsift/clinvar_20210828.vcf.gz"
def gnomad_database = "${genome_dir}snpsift/b37_af-only-gnomad.raw.sites.vcf"
def dbsnp_database = "${genome_dir}snpsift/dbsnp_138.b37.vcf"
def core_hotspots_bed = "${genome_dir}variation_hotspots/core_hotspots_v1.0.bed"
def core_hotspots_vcf = "${genome_dir}variation_hotspots/core_hotspots_v1.0.vcf"
def snp_hotspots_vcf = "${genome_dir}variation_hotspots/snp_hotspots_v1.0.vcf"
def vardict_hotspots_list = file("${genome_dir}variation_hotspots/core_hotspots_v1.0.list")
def msi_database_dir = "${genome_dir}/msi/${params.msi_db}"
def hg19_ref_fa = "/yunying/ref/human/hg19/hg19.fa"
def rna_fusion_conf = "${genome_dir}/genefuse/cancer.hg19.csv"
def rna_fusion_fa = "${genome_dir}rna_fusion/6gene_RNA_fusion.fa"
def rna_fusion_bed = "/yunying/data/product/bed/yunying_bed/6gene_bed9.bed"


// all configuration including qc, alignment, and gatk process configuration file (module ---> software/toolkit ---> conf)
def prefix_conf_dir = "/yunying/codes/product/module"
def preprocess_conf = file("${prefix_conf_dir}/software/preprocess/fastp/conf/fastp_standard.conf")
def genefuse_standard_conf = "${prefix_conf_dir}/software/sv/genefuse/conf/genefuse_standard.conf"
def pcr_preprocess_conf = file("${prefix_conf_dir}/software/preprocess/fastp/conf/fastp_pcr_amplicon.conf")
def mutscan_conf = file("${prefix_conf_dir}/software/somatic/Mutscan/conf/mutscan_filter.conf")
def align_conf = file("${prefix_conf_dir}/toolkit/align_sambam_vcf/conf/align_standard.conf")
def pcr_align_conf = file("${prefix_conf_dir}/toolkit/align_sambam_vcf/conf/align_pcr_amplicon.conf")
def gatk_toolkit_dir = "${prefix_conf_dir}/toolkit/gatk/"
def recal_conf = file("${gatk_toolkit_dir}refinement/conf/refinement_standard.conf")
def samtools_mpileup_conf = file("${prefix_conf_dir}/toolkit/varscan/conf/samtools_mpileup.conf")
def pcr_mpileup_conf = file("${prefix_conf_dir}/toolkit/varscan/conf/samtools_pcr_mpileup.conf")
def germline_calling_conf = file("${gatk_toolkit_dir}germline/conf/haplotype_caller.conf")
def somatic_varscan_conf = file("${prefix_conf_dir}/toolkit/varscan/conf/somatic_high_depth_ctdna.conf")
def somatic_varscan_paired_conf = file("${prefix_conf_dir}/toolkit/varscan/conf/somatic_paired_ctdna.conf")
def pcr_varscan_conf = file("${prefix_conf_dir}/toolkit/varscan/conf/somatic_pcr_ctdna.conf")
def indels_pindel_conf = file("${prefix_conf_dir}/software/indels/pindel/conf/pindel_tumor_standard.conf")
def indel_hotspots_bed = file("${prefix_conf_dir}/software/indels/pindel/conf/indel_hotspots.bed")
def sv_svaba_conf = file("${prefix_conf_dir}/software/sv/svaba/conf/somatic_sv_high_recall.conf")
def sv_hotspots_conf = file("${genome_dir}anno/sv_gene_extend_region.txt")
def sv_partner_conf = file("${genome_dir}anno/sv_hotspots_partner.txt")
def variation_hotspots_conf = file("${genome_dir}variation_hotspots/yunying_hotspots_GRCh37_v1.0.txt")
def cnv_oncocnv_conf = file("${prefix_conf_dir}/software/cnv/oncocnv/conf/somatic_cnv_tumor.conf")
def anno_snpeff_conf = file("${prefix_conf_dir}/software/annotation/SnpEff/conf/snpeff_canon_only.conf")
def perbase_conf = file("${prefix_conf_dir}/software/calling/perbase/conf/perbase_standard.conf")
def somatic_vardict_conf = "${prefix_conf_dir}/software/calling/VarDict/conf/ctdna_standard.conf"
def secondary_annotate_conf = file("${gatk_toolkit_dir}annotation/conf/secondary_bam_annotate.conf")


process preprocess_fastp_tumor_standard {
     // use fastp to trim low quality bases from reads and do adapters trimming
     container '10.168.2.67:5000/bromberglab/fastp:latest'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*qc.json'

     cpus params.threads
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 1

     input:
         file tumor_fq1
         file tumor_fq2

     output:
         path "${tumor_fastq_r1}" into trimmed_tumor_fq_r1_ch, trimmed_tumor_fq_r1_ch2, trimmed_tumor_fq_r1_ch3
         path "${tumor_fastq_r2}" into trimmed_tumor_fq_r2_ch, trimmed_tumor_fq_r2_ch2, trimmed_tumor_fq_r2_ch3
         path "${fastp_qc_json}" into fastq_qc_json

     script:
         tumor_fastq_r1 = "${params.SampleId}_tumor_r1_trimmed.fastq"
         tumor_fastq_r2 = "${params.SampleId}_tumor_r2_trimmed.fastq"
         fastp_qc_json = "${params.SampleId}_tumor_fastp_qc.json"

     if (params.SampleId =~ "PCR_g_MSI" || params.SampleId =~ "PCR_ct_MSI")
         """
         cp ${tumor_fq1} ${tumor_fastq_r1}
         cp ${tumor_fq2} ${tumor_fastq_r2}
         touch ${fastp_qc_json}
         """
     else if (params.SampleId =~ "PCR_")
         """
         cat ${pcr_preprocess_conf} | xargs \
         fastp -i ${tumor_fq1} -I ${tumor_fq2} -w ${task.cpus} -j ${fastp_qc_json} -o ${tumor_fastq_r1} -O ${tumor_fastq_r2}
         """
     else
         """
         cat ${preprocess_conf} | xargs \
         fastp -i ${tumor_fq1} -I ${tumor_fq2} -w ${task.cpus} -j ${fastp_qc_json} -o ${tumor_fastq_r1} -O ${tumor_fastq_r2}
         """
}


process fusion_rna_call {
     container '10.168.2.67:5000/fusion_rna_py:v1.0'
     publishDir "${unified_bam_dir}", mode: 'copy', pattern: '*sorted.bam'
     publishDir "${unified_bam_dir}", mode: 'copy', pattern: '*sorted.bam.bai'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*_fusion.txt'

     cpus params.threads
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'ignore' }
     maxRetries 1

     when:
         params.SampleId =~ /PCR_g_12$/ || params.SampleId =~ /PCR_ct_12$/

     input:
         file "fastq_r1.fq" from trimmed_tumor_fq_r1_ch3
         file "fastq_r2.fq" from trimmed_tumor_fq_r2_ch3

     output:
         path "${fusion_result_file}" into fusion_result_file_ch

     script:
         fusion_bam = "${params.SampleId}_fusion_sorted.bam"
         fusion_result_file = "${params.SampleId}_rna_fusion.txt"

     """
     set -o pipefail
     cat ${pcr_align_conf} | xargs \
     bwa mem -R '@RG\\tID:${params.SampleId}\\tSM:${params.SampleId}\\tLB:YUNYING\\tPL:Illumina' -t 4 ${rna_fusion_fa} \
     fastq_r1.fq fastq_r2.fq | samtools fixmate -m - - | samtools sort -m 2G -@4 -T /tmp/${params.SampleId} -o ${fusion_bam} -
     samtools index ${fusion_bam}

     python3 /script/fusion_check.py ${rna_fusion_fa} ${fusion_bam} > ${fusion_result_file}
     """
}

process genefuse_rna_tumor_only {
     container '10.168.2.67:5000/genefuse-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*_report.*'

     cpus 2
     memory { 16.GB * task.attempt }

     errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
     maxRetries 2
     
     when:
         params.SampleId =~ /PCR_g_12$/ || params.SampleId =~ /PCR_ct_12$/

     input:
         path 'tumor_r1.fq.gz' from trimmed_tumor_fq_r1_ch2
         path 'tumor_r2.fq.gz' from trimmed_tumor_fq_r2_ch2
         path 'fusion.csv' from file("${rna_fusion_conf}")

     output:
         path '*' into genefuse_rna_ch

     script:
         report_html = "${params.SampleId}_rna_fusion_report.html"
         report_plain_txt = "${params.SampleId}_rna_fusion_report.txt"
     """
     cat ${genefuse_standard_conf} | xargs \
     genefuse -r ${hg19_ref_fa} -t 1 -f fusion.csv -1 tumor_r1.fq.gz -2 tumor_r2.fq.gz -h ${report_html} > ${report_plain_txt}
     """
}

process preprocess_fastp_control_standard {
     // use fastp to trim low quality bases from reads and do adapters trimming
     container '10.168.2.67:5000/bromberglab/fastp:latest'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*qc.json'

     cpus params.threads
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 1
     
     when:
          params.ctl_fq1 != ""

     input:
         file ctl_fq1
         file ctl_fq2

     output:
         path "${ctl_fastq_r1}" into trimmed_ctl_fq_r1_ch, trimmed_ctl_fq_r1_ch2
         path "${ctl_fastq_r2}" into trimmed_ctl_fq_r2_ch, trimmed_ctl_fq_r2_ch2
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
     maxRetries 2

     input:
         path 'trimmed_fastq_r1' from trimmed_tumor_fq_r1_ch
         path 'trimmed_fastq_r2' from trimmed_tumor_fq_r2_ch

     output:
         path "${tumor_sorted_bam}" into tumor_sorted_bam_ch
         path "${tumor_sorted_bai}" into tumor_sorted_bai_ch

     script:
         tumor_sorted_bam = "${params.SampleId}_tumor_sorted.bam"
         tumor_sorted_bai = "${params.SampleId}_tumor_sorted.bam.bai"
     if (params.SampleId =~ "PCR_")
         """
         cat ${pcr_align_conf} | xargs \
         bwa mem -R '@RG\\tID:${params.SampleId}\\tSM:${params.SampleId}\\tLB:YUNYING\\tPL:Illumina' -t ${task.cpus} ${ref_fa} trimmed_fastq_r1 trimmed_fastq_r2 | sambamba view \
         -l 0 -f bam -S -h /dev/stdin |sambamba sort -m 5G -t 10 -o ${tumor_sorted_bam} /dev/stdin

         """
     else
         """
         cat ${align_conf} | xargs \
         bwa mem -R '@RG\\tID:${params.SampleId}\\tSM:${params.SampleId}\\tLB:YUNYING\\tPL:Illumina' -t ${task.cpus} ${ref_fa} trimmed_fastq_r1 trimmed_fastq_r2 | sambamba view \
         -l 0 -f bam -S -h /dev/stdin |sambamba sort -m 5G -t 10 -o ${tumor_sorted_bam} /dev/stdin

         """
}


process control_align_bwa_sort {
     //container 'cmopipeline/bwa-samtools-gatk4-sambamba-samblaster:latest'
     container '10.168.2.67:5000/align_sambam_vcf:1.0'

     cpus params.threads
     memory { 20.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
          params.ctl_fq1 != ""
     
     input:
         path 'trimmed_fastq_r1' from trimmed_ctl_fq_r1_ch
         path 'trimmed_fastq_r2' from trimmed_ctl_fq_r2_ch

     output:
         path "${ctl_sorted_bam}" into ctl_sorted_bam_ch
         path "${ctl_sorted_bai}" into ctl_sorted_bai_ch

     script:
         ctl_sorted_bam = "${params.SampleId}_ctl_sorted.bam"
         ctl_sorted_bai = "${params.SampleId}_ctl_sorted.bam.bai"
     """
     cat ${align_conf} | xargs \
     bwa mem -R '@RG\\tID:${params.SampleId}_X\\tSM:${params.SampleId}_X\\tLB:YUNYING\\tPL:Illumina' -t ${task.cpus} ${ref_fa} \
     trimmed_fastq_r1 trimmed_fastq_r2 | sambamba view \
     -l 0 -f bam -S -h /dev/stdin |sambamba sort -m 5G -t 10 -o ${ctl_sorted_bam} /dev/stdin

     ### alternatively, we would use samtools for sort bam file, but some restrictions occurs!
     ### argument "-Y" (or "-C") used by bwa does not support samtools sort directly, so we omit it (!!!)
     ### and some read with problemed coordinate sorted!
     # cat ${align_conf} | xargs \
     # bwa mem -R '@RG\\tID:${params.SampleId}\\tSM:${params.SampleId}\\tLB:YUNYING\\tPL:Illumina\\tPU:chip' \
     # -M -t ${task.cpus} ${ref_fa} trimmed_fastq_r1 \
     # trimmed_fastq_r2 | samtools sort -m 2G -t ${task.cpus} -O bam -o ${ctl_sorted_bam} /dev/stdin
     # samtools index ${ctl_sorted_bam}

     """
}

process refinement_process_tumor_standard {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'

     cpus params.threads/2
     memory { 20.GB * memory_scale1 * task.attempt }
     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5

     input:
          path 'tumor_sorted.bam' from tumor_sorted_bam_ch
          path 'tumor_sorted.bai' from tumor_sorted_bai_ch
          path 'target.bed' from bed_target

     output:
          path "${tumor_recal_bam}" into tumor_recal_bam_ch1, tumor_recal_bam_ch2, tumor_recal_bam_ch3, tumor_recal_bam_ch4, tumor_recal_bam_ch5, tumor_recal_bam_ch6, tumor_recal_bam_ch7, tumor_recal_bam_ch8, tumor_recal_bam_ch9, tumor_recal_bam_ch10
          path "${tumor_recal_bai}" into tumor_recal_bai_ch1, tumor_recal_bai_ch2, tumor_recal_bai_ch3, tumor_recal_bai_ch4, tumor_recal_bai_ch5, tumor_recal_bai_ch6, tumor_recal_bai_ch7, tumor_recal_bai_ch8, tumor_recal_bai_ch9, tumor_recal_bai_ch10

     script:
          tumor_recal_bam = "${params.SampleId}_tumor_recal.bam"
          tumor_recal_bai = "${params.SampleId}_tumor_recal.bai"

     """
     elprep filter tumor_sorted.bam ${tumor_recal_bam} --sorting-order keep --bqsr output.recal \
            --bqsr-reference ${ref_elfa} --known-sites ${gatk_db_dir}/dbsnp_138.b37.elsites,${gatk_db_dir}/1000G_phase1.indels.b37.elsites,${gatk_db_dir}/Mills_and_1000G_gold_standard.indels.b37.elsites --nr-of-threads ${task.cpus}
     gatk BuildBamIndex -I ${tumor_recal_bam}
     """
}

process dedup_refinement_process_tumor {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'
     publishDir "${unified_bam_dir}", mode: 'copy'

     cpus params.threads
     memory { 20.GB * memory_scale1 * task.attempt }
     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5

     when:
          pcr_panel_mark == false

     input:
          path 'tumor_recal.bam' from tumor_recal_bam_ch8
          path 'tumor_recal.bai' from tumor_recal_bai_ch8
          path 'target.bed' from bed_target

     output:
          path "${tumor_dedup_bam}" into tumor_dedup_bam_ch1, tumor_dedup_bam_ch2
          path "${tumor_dedup_bai}" into tumor_dedup_bai_ch1, tumor_dedup_bai_ch2

     script:
          tumor_dedup_bam = "${params.SampleId}_tumor_dedup.bam"
          tumor_dedup_bai = "${params.SampleId}_tumor_dedup.bai"

     """
     gencore -i tumor_recal.bam -o tumor_deduped.bam -r ${ref_fa} -b target.bed
     gatk SortSam -SO coordinate -I tumor_deduped.bam -O ${tumor_dedup_bam}
     gatk BuildBamIndex -I ${tumor_dedup_bam}
     """
}


process refinement_process_control_standard {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'
     publishDir "${unified_bam_dir}", mode: 'copy'

     cpus params.threads/2
     memory { 15.GB * memory_scale2 * task.attempt }
     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5
     
     when:
          params.ctl_fq1 != ""

     input:
          path 'ctl_sorted.bam' from ctl_sorted_bam_ch
          path 'ctl_sorted.bai' from ctl_sorted_bai_ch
          path 'ctl_target.bed' from bed_target

     output:
          path "${ctl_recal_bam}" into ctl_recal_bam_ch1, ctl_recal_bam_ch2, ctl_recal_bam_ch3, ctl_recal_bam_ch4, ctl_recal_bam_ch5, ctl_recal_bam_ch6, ctl_recal_bam_ch7, ctl_recal_bam_ch10
          path "${ctl_recal_bai}" into ctl_recal_bai_ch1, ctl_recal_bai_ch2, ctl_recal_bai_ch3, ctl_recal_bai_ch4, ctl_recal_bai_ch5, ctl_recal_bai_ch6, ctl_recal_bai_ch7, ctl_recal_bai_ch10

     script:
          ctl_recal_bam = "${params.SampleId}_ctl_recal.bam"
          ctl_recal_bai = "${params.SampleId}_ctl_recal.bai"

     """
     java -Xmx16G -jar /opt/abra2.jar --in ctl_sorted.bam --out ctl_realigned.bam --ref ${ref_fa} --threads 6 --targets ${realignment_bed}
     elprep filter ctl_realigned.bam ${ctl_recal_bam} --mark-duplicates --mark-optical-duplicates output.metrics --sorting-order keep --bqsr output.recal \
            --bqsr-reference ${ref_elfa} --known-sites ${gatk_db_dir}/dbsnp_138.b37.elsites,${gatk_db_dir}/1000G_phase1.indels.b37.elsites,${gatk_db_dir}/Mills_and_1000G_gold_standard.indels.b37.elsites --nr-of-threads ${task.cpus}
     gatk BuildBamIndex -I ${ctl_recal_bam}
     """
}

process dedup_refinement_process_control {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'
     publishDir "${unified_bam_dir}", mode: 'copy'

     cpus params.threads
     memory { 15.GB * memory_scale2 * task.attempt }
     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5
     
     when:
          params.ctl_fq1 != "" && pcr_panel_mark == false

     input:
          path 'ctl_recal.bam' from ctl_recal_bam_ch6
          path 'ctl_recal.bai' from ctl_recal_bai_ch6
          path 'target.bed' from bed_target

     output:
          path "${ctl_dedup_bam}" into ctl_dedup_bam_ch1
          path "${ctl_dedup_bai}" into ctl_dedup_bai_ch1

     script:
          ctl_dedup_bam = "${params.SampleId}_ctl_dedup.bam"
          ctl_dedup_bai = "${params.SampleId}_ctl_dedup.bai"

     """
     gencore -i ctl_recal.bam -o ctl_deduped.bam -r ${ref_fa} -b target.bed
     gatk SortSam -SO coordinate -I ctl_deduped.bam -O ${ctl_dedup_bam}
     gatk BuildBamIndex -I ${ctl_dedup_bam}                 
     """
}

process refinement_enzyme_correct_tumor {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'
     publishDir "${unified_bam_dir}", mode: 'copy'
 
     cpus 2
     memory { 8.GB * task.attempt }
 
     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 2
     
 	 when:
         somatic_process_mark == true
 
     input:
           path 'tumor_recal.bam'  from tumor_recal_bam_ch1
           path 'tumor_recal.bai'  from tumor_recal_bai_ch1
           path "target.bed" from bed_target
 
     output:
           path "${tumor_correct_bam}" into tumor_refined_bam_ch1, tumor_refined_bam_ch2, tumor_refined_bam_ch3
           path "${tumor_correct_bai}" into tumor_refined_bai_ch1, tumor_refined_bai_ch2, tumor_refined_bai_ch3
           // path "${tumor_correct_bam}" into tumor_correct_bam_ch
           // path "${tumor_correct_bai}" into tumor_correct_bai_ch
 
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

// TODO: long time consuming realignment process
// process realignment_abra2_tumor_standard {
//      container '10.168.2.67:5000/gatk-yy:4.2.0.0'
// 
//      cpus params.threads/2
//      memory { 20.GB * memory_scale1 * task.attempt }
//      errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
//      maxRetries 2
// 
//      input:
//           path 'tumor_correct.bam' from tumor_correct_bam_ch
//           path 'tumor_correct.bai' from tumor_correct_bai_ch
// 
//      output:
//           path "${tumor_realign_bam}" into tumor_refined_bam_ch1, tumor_refined_bam_ch2, tumor_refined_bam_ch3
//           path "${tumor_realign_bai}" into tumor_refined_bai_ch1, tumor_refined_bai_ch2, tumor_refined_bai_ch3
// 
//      script:
//           tumor_realign_bam = "${params.SampleId}_tumor_recal.bam"
//           tumor_realign_bai = "${params.SampleId}_tumor_recal.bai"
//      """
//      java -Xmx16G -jar /opt/abra2.jar --in tumor_correct.bam --out ${tumor_realign_bam} --ref ${ref_fa} --threads 6 --targets ${realignment_bed}
//      gatk BuildBamIndex -I ${tumor_realign_bam}
//      """
// }

process vardict_variation_tumor_only {
     container '10.168.2.67:5000/vardict-yy:1.8.2m'

     cpus 6
     memory { 16.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'ignore' }
     maxRetries 2

     when:
         params.ctl_fq1 == ""

     input:
         path "refined.bam" from tumor_refined_bam_ch2
         path "target.bed" from bed_target

     output:
         path "${vardict_out_vcf}" into vardict_out_ch

     script:
         vardict_out_vcf = "${params.SampleId}_raw_vardict.vcf"
     """
     cat ${somatic_vardict_conf} | xargs java -Xmx12g -XX:-UseParallelGC -jar /opt/VarDict-1.8.2.jar -b refined.bam -G ${ref_fa} -th ${task.cpus} -N ${params.SampleId} -hotspot ${vardict_hotspots_list} -c 1 -S 2 -E 3 target.bed > "${params.SampleId}_vardict.var"
     cat "${params.SampleId}_vardict.var" | /opt/teststrandbias.R | /opt/var2vcf_valid.pl -N ${params.SampleId} -A -E -f 0.00085 > raw_vardict_tumor.vcf
     bcftools sort raw_vardict_tumor.vcf | awk -F "\t" '{if (\$4 != "." && \$5 != "." && \$5 !~ /^</) print \$0}' > ${vardict_out_vcf}
     """
}

process secondary_filter_annotation_tumor_only {
    container '10.168.2.67:5000/gatk-yy:4.2.2.1'

    cpus 2
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
    maxRetries 3
    
    when:
        params.ctl_fq1 == ""

    input:
        path 'tumor_deduped.bam' from tumor_dedup_bam_ch2
        path 'tumor_deduped.bai' from tumor_dedup_bai_ch2
        path 'raw_sorted.vcf' from vardict_out_ch

    output:
        path "${secondary_stats_vcf}" into vardict_add_stats_ch
    
    script:
        secondary_stats_vcf = "${params.SampleId}_vardict_stats.vcf"
    """
    cat ${secondary_annotate_conf} | xargs \
    gatk VariantAnnotator -I tumor_deduped.bam -R ${ref_fa} -V raw_sorted.vcf -O ${secondary_stats_vcf}
    """
}

process anno_snpeff_tumor_only2 {
     container '10.168.2.67:5000/snpeff-yy:4.4m'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 120..140 ? 'retry' : 'ignore' }
     maxRetries 2

     when:
         params.ctl_fq1 == ""

     input:
         path "raw_vardict_var.vcf" from vardict_add_stats_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 
       
     output:
         path "${prefix_anno_vcf}_vardict_hgvs.vcf" into tumor_anno_vcf_ch2
         path "${prefix_anno_vcf}.haplotype.anno" into haplotype_anno_tumor_ch

     script:
         prefix_anno_vcf = "${params.SampleId}_raw"

     """
     cat anno_snpeff.conf | xargs \
     java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}.haplotype.anno GRCh37.p13.RefSeq raw_vardict_var.vcf > ${prefix_anno_vcf}_vardict_hgvs.vcf

     """
}

process anno_snpsift_database_tumor_only {
     container '10.168.2.67:5000/snpsift-yy:1.0'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 120..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
         params.ctl_fq1 == ""
       
     input:
         path 'hgvs.vcf' from tumor_anno_vcf_ch2

     output:
         path "${tumor_anno_vcf}" into tumor_hgvs_clinvar_ch

     script:
         tumor_anno_vcf = "${params.SampleId}_filter_marked_hgvs_db.vcf"
     """
     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_database} hgvs.vcf > ${tumor_anno_vcf}
     """
}

process vardict_report_tumor_only {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'ignore' }
     maxRetries 2
     
     when:
         params.ctl_fq1 == ""

     input:
         path 'marked_all_anno.vcf' from tumor_hgvs_clinvar_ch
         path 'haplotype.anno' from haplotype_anno_tumor_ch

     output:
         path '*.xls' into vardict_report_ch

     script:
         report_output = "${params.SampleId}_vardict_for_report.xls"
     """
     python3 /opt/vcf2mut_report_review.py marked_all_anno.vcf haplotype.anno ${variation_hotspots_conf} all_sites.mpileup ${report_output} cfdna
     """
}

process samtools_tumor_mpileup {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 6.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        somatic_process_mark == true

    input:
        path "recal.bam" from tumor_refined_bam_ch1
        path "recal.bai" from tumor_refined_bai_ch1
        path "samtools_mpileup_conf" from samtools_mpileup_conf
        path "pcr_mpileup_conf" from pcr_mpileup_conf
        path "target.bed" from bed_target

    output:
        path "${tumor_mpileup_file}" into tumor_mpileup_ch1, tumor_mpileup_ch2
    
    script:
        tumor_mpileup_file = "${params.SampleId}_T_mpileup"
    if (params.SampleId =~ "PCR_")
        """
        cat pcr_mpileup_conf | xargs \
        samtools mpileup -f ${ref_fa} -l target.bed recal.bam -o ${tumor_mpileup_file}
        """
    else
        """
        awk '{print \$1"\t"\$2-20"\t"\$3+20}' target.bed > target.ext20.bed
        cat samtools_mpileup_conf | xargs \
        samtools mpileup -f ${ref_fa} -l target.ext20.bed recal.bam -o ${tumor_mpileup_file}
        """
}

process samtools_tumor_dedup_mpileup {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 6.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        dedup_process_mark == true && somatic_process_mark == true

    input:
        path "dedup.bam" from tumor_dedup_bam_ch1
        path "dedup.bai" from tumor_dedup_bai_ch1
        path "samtools_mpileup_conf" from samtools_mpileup_conf
        path "target.bed" from bed_target

    output:
        path "${tumor_dedup_mpileup_file}" into tumor_dedup_mpileup_ch1, tumor_dedup_mpileup_ch2        
    
    script:
        tumor_dedup_mpileup_file = "${params.SampleId}_dedup_T_mpileup"

    """
    cat samtools_mpileup_conf | xargs \
    samtools mpileup -f ${ref_fa} -l target.bed dedup.bam -o ${tumor_dedup_mpileup_file}
    """
}

process samtools_control_mpileup {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 6.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        params.ctl_fq1 != ""

    input:
        path "ctl.bam" from ctl_recal_bam_ch1
        path "ctl.bai" from ctl_recal_bai_ch1
        path "samtools_mpileup_conf" from samtools_mpileup_conf
        path "target.bed" from bed_target

    output:
        path "${control_mpileup_file}" into control_mpileup_ch
    
    script:
        control_mpileup_file = "${params.SampleId}_N_mpileup"

    """
    awk '{print \$1"\t"\$2-20"\t"\$3+20}' target.bed > target.ext20.bed
    cat samtools_mpileup_conf | xargs \
    samtools mpileup -f ${ref_fa} -l target.ext20.bed ctl.bam -o ${control_mpileup_file}

    """
}

process samtools_control_dedup_mpileup {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 6.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        params.ctl_fq1 != "" && dedup_process_mark == true

    input:
        path "ctl_dedup.bam" from ctl_dedup_bam_ch1
        path "ctl_dedup.bai" from ctl_dedup_bai_ch1
        path "samtools_mpileup_conf" from samtools_mpileup_conf
        path "target.bed" from bed_target

    output:
        path "${ctl_dedup_mpileup_file}" into ctl_dedup_mpileup_ch
    
    script:
        ctl_dedup_mpileup_file = "${params.SampleId}_dedup_N_mpileup"

    """
    cat samtools_mpileup_conf | xargs \
    samtools mpileup -f ${ref_fa} -l target.bed ctl_dedup.bam -o ${ctl_dedup_mpileup_file}
    """
}


process somatic_varscan_tumor_paired {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 4.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        params.ctl_fq1 != ""

    input:
        path "${params.SampleId}_T_mpileup" from tumor_mpileup_ch1
        path "${params.SampleId}_N_mpileup" from control_mpileup_ch
        path "somatic_varscan_conf" from somatic_varscan_paired_conf

    output:
        path "${raw_snp_vcf}" into raw_snp_vcf_ch1
        path "${raw_indel_vcf}" into raw_indel_vcf_ch1
        path "${germline_snp_vcf}" into germline_snp_vcf_ch1
        path "${germline_indel_vcf}" into germline_indel_vcf_ch1
    
    script:
        raw_snp_vcf = "${params.SampleId}_raw.snp.Somatic.hc.vcf"
        raw_indel_vcf = "${params.SampleId}_raw.indel.Somatic.hc.vcf"
        germline_snp_vcf = "${params.SampleId}_raw.snp.Germline.hc.vcf"
        germline_indel_vcf = "${params.SampleId}_raw.indel.Germline.hc.vcf"

    """
    cat somatic_varscan_conf | xargs \
    java -jar /opt/VarScan.v2.3.11.jar somatic ${params.SampleId}_N_mpileup ${params.SampleId}_T_mpileup ${params.SampleId}_raw --hotspot-vcf ${core_hotspots_vcf}
    java -jar /opt/VarScan.v2.3.11.jar processSomatic ${params.SampleId}_raw.snp.vcf --min-tumor-freq 0.001 --max-normal-freq 0.02 --p-value 0.05 --hotspot-vcf ${snp_hotspots_vcf}
    java -jar /opt/VarScan.v2.3.11.jar processSomatic ${params.SampleId}_raw.indel.vcf --min-tumor-freq 0.001 --max-normal-freq 0.02 --p-value 0.05 --hotspot-vcf ${core_hotspots_vcf}

    """
}

process somatic_vardict_paired_parallel {
     container '10.168.2.67:5000/vardict-yy:1.8.2m'
     
     cpus 6
     memory { 16.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'ignore' }
     maxRetries 2

     when:
         params.ctl_fq1 != "" && somatic_process_mark == true

     input:
         path "refined.bam" from tumor_refined_bam_ch2
         path "control.bam" from ctl_recal_bam_ch7
         path "control.bai" from ctl_recal_bai_ch7
         path "target.bed" from bed_target

     output:
         path "${vardict_out_vcf}" into paired_vardict_ch

     script:
         vardict_out_vcf = "${params.SampleId}_raw_vardict.vcf"
     """
     cat ${somatic_vardict_conf} | xargs java -Xmx12g -XX:-UseParallelGC -jar /opt/VarDict-1.8.2.jar -b "refined.bam|control.bam" -G ${ref_fa} -th ${task.cpus} -N ${params.SampleId} -hotspot ${vardict_hotspots_list} -c 1 -S 2 -E 3 target.bed > "${params.SampleId}_vardict.var"
     cat "${params.SampleId}_vardict.var" | /opt/testsomatic.R | /opt/var2vcf_paired.pl -N "${params.SampleId}|${params.SampleId}_X" -A -f 0.00085 > raw_vardict_tumor.vcf
     bcftools sort raw_vardict_tumor.vcf | awk -F "\t" '{if (\$4 != "." && \$5 != "." && \$5 !~ /^</) print \$0}' > ${vardict_out_vcf}
     """
}

process secondary_filter_annotation_paired {
    container '10.168.2.67:5000/gatk-yy:4.2.2.1'

    cpus 2
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
    maxRetries 3
    
    when:
        params.ctl_fq1 != ""

    input:
        path 'tumor_deduped.bam' from tumor_dedup_bam_ch2
        path 'tumor_deduped.bai' from tumor_dedup_bai_ch2
        path 'paired_raw_sorted.vcf' from paired_vardict_ch

    output:
        path "${secondary_stats_vcf}" into vardict_stats_paired_ch

    script:
        secondary_stats_vcf = "${params.SampleId}_vardict_stats.vcf"
    """
    cat ${secondary_annotate_conf} | xargs \
    gatk VariantAnnotator -I tumor_deduped.bam -R ${ref_fa} -V paired_raw_sorted.vcf -O ${secondary_stats_vcf}
    """
}

process anno_snpeff_somatic_paired {
     container '10.168.2.67:5000/snpeff-yy:4.4m'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 120..140 ? 'retry' : 'ignore' }
     maxRetries 2

     when:
         params.ctl_fq1 != ""

     input:
         path "paired_vardict_var.vcf" from vardict_stats_paired_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

     output:
         path "${prefix_anno_vcf}_vardict_hgvs.vcf" into paired_anno_vcf_ch
         path "${prefix_anno_vcf}.haplotype.anno" into haplotype_pseudo_ch

     script:
         prefix_anno_vcf = "${params.SampleId}_raw"

     """
     cat anno_snpeff.conf | xargs \
     java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}.haplotype.anno GRCh37.p13.RefSeq paired_vardict_var.vcf > ${prefix_anno_vcf}_vardict_hgvs.vcf

     """
}

process anno_snpsift_database_paired {
     container '10.168.2.67:5000/snpsift-yy:1.0'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 120..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
         params.ctl_fq1 != ""
       
     input:
         path "hgvs.vcf" from paired_anno_vcf_ch

     output:
         path "${tumor_anno_vcf}" into paired_hgvs_clinvar_ch

     script:
         tumor_anno_vcf = "${params.SampleId}_filter_marked_hgvs_db.vcf"
     """
     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_database} hgvs.vcf > ${tumor_anno_vcf}
     """
}

process vardict_report_paired {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'ignore' }
     maxRetries 2
     
     when:
         params.ctl_fq1 != ""

     input:
         path "marked_all_anno.vcf" from paired_hgvs_clinvar_ch
         path "haplotype.anno" from haplotype_pseudo_ch

     output:
         path "*.xls" into vardict_paired_report

     script:
         report_output = "${params.SampleId}_vardict_for_report.xls"
     """
     python3 /opt/vcf2mut_report_review.py marked_all_anno.vcf haplotype.anno ${variation_hotspots_conf} all_sites.mpileup ${report_output} cfdna
     """
}

process somatic_varscan_tumor_dedup_paired {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 4.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        params.ctl_fq1 != "" && dedup_process_mark == true

    input:
        path "${params.SampleId}_dedup_T_mpileup" from tumor_dedup_mpileup_ch1
        path "${params.SampleId}_dedup_N_mpileup" from ctl_dedup_mpileup_ch
        path "somatic_varscan_conf" from somatic_varscan_paired_conf


    output:
        path "${dedup_raw_snp_vcf}" into dedup_raw_snp_vcf_ch1
        path "${dedup_raw_indel_vcf}" into dedup_raw_indel_vcf_ch1
        path "${dedup_germline_snp_vcf}" into dedup_germline_snp_vcf_ch1
        path "${dedup_germline_indel_vcf}" into dedup_germline_indel_vcf_ch1
    
    script:
        dedup_raw_snp_vcf = "${params.SampleId}_dedup_raw.snp.Somatic.hc.vcf"
        dedup_raw_indel_vcf = "${params.SampleId}_dedup_raw.indel.Somatic.hc.vcf"
        dedup_germline_snp_vcf = "${params.SampleId}_dedup_raw.snp.Germline.hc.vcf"
        dedup_germline_indel_vcf = "${params.SampleId}_dedup_raw.indel.Germline.hc.vcf"

    """
    cat somatic_varscan_conf | xargs \
    java -jar /opt/VarScan.v2.3.10.jar somatic ${params.SampleId}_dedup_N_mpileup ${params.SampleId}_dedup_T_mpileup ${params.SampleId}_dedup_raw  --hotspot-vcf ${core_hotspots_vcf}
    java -jar /opt/VarScan.v2.3.10.jar processSomatic ${params.SampleId}_dedup_raw.snp.vcf --min-tumor-freq 0.001 --max-normal-freq 0.02 --p-value 0.07  --hotspot-vcf ${core_hotspots_vcf}
    java -jar /opt/VarScan.v2.3.10.jar processSomatic ${params.SampleId}_dedup_raw.indel.vcf --min-tumor-freq 0.001 --max-normal-freq 0.02 --p-value 0.07  --hotspot-vcf ${core_hotspots_vcf} 
    """
}


process concat_snv_indel_paired {
    container '10.168.2.67:5000/bcftools-yy:1.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'

    cpus 4
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
        params.ctl_fq1 != ""
    
    input:
        path 'raw_snp_vcf' from raw_snp_vcf_ch1
        path 'raw_indel_vcf' from raw_indel_vcf_ch1
        path 'germline_hc_snp_vcf' from germline_snp_vcf_ch1
        path 'germline_hc_indel_vcf' from germline_indel_vcf_ch1
    
    output:
        path "${raw_snp_indel_vcf}" into raw_paired_vcf_ch        
        path "${germline_snp_indel_vcf}" into germline_paired_vcf_ch

    script:
        raw_snp_indel_vcf = "${params.SampleId}_raw.snp.indel.Somatic.hc.vcf"
        germline_snp_indel_vcf = "${params.SampleId}_raw.snp.indel.Germline.hc.vcf"
          
    """      
    bgzip -c raw_snp_vcf > raw_snp_bgzip.gz
    bcftools index raw_snp_bgzip.gz
    bgzip -c raw_indel_vcf > raw_indel_bgzip.gz
    bcftools index raw_indel_bgzip.gz
    bcftools concat -a -D raw_snp_bgzip.gz raw_indel_bgzip.gz -o ${raw_snp_indel_vcf}

    bgzip -c germline_hc_snp_vcf > germline_snp_bgzip.gz
    bcftools index germline_snp_bgzip.gz
    bgzip -c germline_hc_indel_vcf > germline_indel_bgzip.gz
    bcftools index germline_indel_bgzip.gz
    bcftools concat -a -D germline_snp_bgzip.gz germline_indel_bgzip.gz -o ${germline_snp_indel_vcf}  
    """     
}

process concat_snv_indel_dedup_paired {
    container '10.168.2.67:5000/bcftools-yy:1.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'

    cpus 4
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
        params.ctl_fq1 != "" && dedup_process_mark == true
    
    input:
        path 'dedup_raw_snp_vcf' from dedup_raw_snp_vcf_ch1
        path 'dedup_raw_indel_vcf' from dedup_raw_indel_vcf_ch1
        path 'dedup_germline_hc_snp_vcf' from dedup_germline_snp_vcf_ch1
        path 'dedup_germline_hc_indel_vcf' from dedup_germline_indel_vcf_ch1
    
    output:
        path "${dedup_raw_snp_indel_vcf}" into dedup_raw_paired_vcf_ch
        path "${dedup_germline_snp_indel_vcf}" into dedup_germline_paired_vcf_ch

    script:
        dedup_raw_snp_indel_vcf = "${params.SampleId}_dedup_raw.snp.indel.Somatic.hc.vcf"
        dedup_germline_snp_indel_vcf = "${params.SampleId}_dedup_raw.snp.indel.Germline.hc.vcf"
          
    """      
    bgzip -c dedup_raw_snp_vcf > dedup_raw_snp_bgzip.gz
    bcftools index dedup_raw_snp_bgzip.gz
    bgzip -c dedup_raw_indel_vcf > dedup_raw_indel_bgzip.gz
    bcftools index dedup_raw_indel_bgzip.gz
    bcftools concat -a -D dedup_raw_snp_bgzip.gz dedup_raw_indel_bgzip.gz -o ${dedup_raw_snp_indel_vcf}

    bgzip -c dedup_germline_hc_snp_vcf > dedup_germline_snp_bgzip.gz
    bcftools index dedup_germline_snp_bgzip.gz
    bgzip -c dedup_germline_hc_indel_vcf > dedup_germline_indel_bgzip.gz
    bcftools index dedup_germline_indel_bgzip.gz
    bcftools concat -a -D dedup_germline_snp_bgzip.gz dedup_germline_indel_bgzip.gz -o ${dedup_germline_snp_indel_vcf}           
    """ 
}          

process somatic_varscan_tumor_snp {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 4.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        params.ctl_fq1 == ""

    input:
        path "${params.SampleId}_mpileup" from tumor_mpileup_ch1
        path "samtools_mpileup_conf" from samtools_mpileup_conf
        path "somatic_varscan_conf" from somatic_varscan_conf
        path "somatic_pcr_varscan_conf" from pcr_varscan_conf

    output:
        path "${raw_snp_vcf}" into raw_snp_vcf_ch2

    script:
        raw_snp_vcf = "${params.SampleId}_raw.snp.vcf"
    if (params.SampleId =~ "PCR_")
        """
        cat somatic_pcr_varscan_conf | xargs \
        java -jar /opt/VarScan.v2.3.10.jar mpileup2snp ${params.SampleId}_mpileup --hotspot-vcf ${core_hotspots_vcf} > ${raw_snp_vcf}
        """
    else
        """
        cat somatic_varscan_conf | xargs \
        java -jar /opt/VarScan.v2.3.10.jar mpileup2snp ${params.SampleId}_mpileup --hotspot-vcf ${core_hotspots_vcf} > ${raw_snp_vcf}   
        """
}


process somatic_varscan_tumor_dedup_snp {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 4.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        params.ctl_fq1 == "" && dedup_process_mark == true

    input:
        path "${params.SampleId}_dedup_mpileup" from tumor_dedup_mpileup_ch1
        path "samtools_mpileup_conf" from samtools_mpileup_conf
        path "somatic_varscan_conf" from somatic_varscan_conf

    output:
       path "${dedup_raw_snp_vcf}" into dedup_raw_snp_vcf_ch2

    script:
       dedup_raw_snp_vcf = "${params.SampleId}_dedup_raw.snp.vcf"
         """
         cat somatic_varscan_conf | xargs \
         java -jar /opt/VarScan.v2.3.10.jar mpileup2snp ${params.SampleId}_dedup_mpileup --hotspot-vcf ${core_hotspots_vcf} > ${dedup_raw_snp_vcf}       
         """
}

process somatic_varscan_tumor_indel {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 4.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        params.ctl_fq1 == ""

    input:
        path "${params.SampleId}_mpileup" from tumor_mpileup_ch2
        path "samtools_mpileup_conf" from samtools_mpileup_conf
        path "somatic_varscan_conf" from somatic_varscan_conf
        path "somatic_pcr_varscan_conf" from pcr_varscan_conf

    output:
        path "${raw_indel_vcf}" into raw_indel_vcf_ch2

    script:
        raw_indel_vcf = "${params.SampleId}_raw.indel.vcf"
    if (params.SampleId =~ "PCR_")
        """
        cat somatic_pcr_varscan_conf | xargs \
        java -jar /opt/VarScan.v2.3.10.jar mpileup2indel ${params.SampleId}_mpileup --hotspot-vcf ${core_hotspots_vcf} > ${raw_indel_vcf}
        """
    else
        """
        cat somatic_varscan_conf | xargs \
        java -jar /opt/VarScan.v2.3.10.jar mpileup2indel ${params.SampleId}_mpileup > ${raw_indel_vcf}
        """
}

process somatic_varscan_tumor_dedup_indel {
    container '10.168.2.67:5000/varscan-yy:1.0'

    cpus 2
    memory { 4.GB * task.attempt }
    time { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
        params.ctl_fq1 == "" && dedup_process_mark == true

    input:
        path "${params.SampleId}_dedup_mpileup" from tumor_dedup_mpileup_ch2
        path "samtools_mpileup_conf" from samtools_mpileup_conf
        path "somatic_varscan_conf" from somatic_varscan_conf

    output:
       path "${dedup_raw_indel_vcf}" into dedup_raw_indel_vcf_ch2

    script:
       dedup_raw_indel_vcf = "${params.SampleId}_dedup_raw.indel.vcf"

    """
    cat somatic_varscan_conf | xargs \
    java -jar /opt/VarScan.v2.3.10.jar mpileup2indel ${params.SampleId}_dedup_mpileup > ${dedup_raw_indel_vcf}      
    """
}


process concat_snv_indel_tumor {
    container '10.168.2.67:5000/bcftools-yy:1.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'

    cpus 4
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
        params.ctl_fq1 == ""
    
    input:
        path 'raw_snp_vcf' from raw_snp_vcf_ch2
        path 'raw_indel_vcf' from raw_indel_vcf_ch2
    
    output:
        path "${raw_snp_indel_vcf}" into raw_tumor_vcf_ch     
        path "${germline_snp_indel_vcf}" into germline_tumor_vcf_ch

    script:
        raw_snp_indel_vcf = "${params.SampleId}_raw.snp.indel.Somatic.hc.vcf"
        germline_snp_indel_vcf = "${params.SampleId}_raw.snp.indel.Germline.hc.vcf"
          
    """
    bgzip -c raw_snp_vcf > raw_snp_bgzip.gz
    bcftools index raw_snp_bgzip.gz
    bgzip -c raw_indel_vcf > raw_indel_bgzip.gz
    bcftools index raw_indel_bgzip.gz
    bcftools concat -a -D raw_snp_bgzip.gz raw_indel_bgzip.gz -o ${raw_snp_indel_vcf}    
    
    # Filter >= 20.0% 
    java -jar /opt/VarScan.v2.3.10.jar filter ${raw_snp_indel_vcf} --min-var-freq 0.20 --output-file ${germline_snp_indel_vcf}
    """
}

process concat_snv_indel_dedup_tumor {
    container '10.168.2.67:5000/bcftools-yy:1.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'

    cpus 4
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
        params.ctl_fq1 == "" && dedup_process_mark == true
    
    input:
        path 'dedup_raw_snp_vcf' from dedup_raw_snp_vcf_ch2
        path 'dedup_raw_indel_vcf' from dedup_raw_indel_vcf_ch2
    
    output:
        path "${dedup_raw_snp_indel_vcf}" into dedup_raw_tumor_vcf_ch
        path "${dedup_germline_snp_indel_vcf}" into dedup_germline_tumor_vcf_ch


    script:
        dedup_raw_snp_indel_vcf = "${params.SampleId}_dedup_raw.snp.indel.Somatic.hc.vcf"
        dedup_germline_snp_indel_vcf = "${params.SampleId}_dedup_raw.snp.indel.Germline.hc.vcf"
          
    """
    bgzip -c dedup_raw_snp_vcf > dedup_raw_snp_bgzip.gz
    bcftools index dedup_raw_snp_bgzip.gz
    bgzip -c dedup_raw_indel_vcf > dedup_raw_indel_bgzip.gz
    bcftools index dedup_raw_indel_bgzip.gz
    bcftools concat -a -D dedup_raw_snp_bgzip.gz dedup_raw_indel_bgzip.gz -o ${dedup_raw_snp_indel_vcf}    
    
    # Filter >= 20.0% 
    java -jar /opt/VarScan.v2.3.10.jar filter ${dedup_raw_snp_indel_vcf} --min-var-freq 0.20 --output-file ${dedup_germline_snp_indel_vcf}

    """
}
      
process pindels_indels_standard {
    container '10.168.2.67:5000/pindel-yy:1.0'

    cpus params.threads/2
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
         pindel_process_mark == true

    input:
         path 'recal.bam'  from tumor_recal_bam_ch5
         path 'recal.bai' from tumor_recal_bai_ch5

    output:
         path "${pindel_out_vcf}" into pindel_vcf_ch

    script:
        pindel_tmp_bed = "${params.SampleId}_pindel.bed"
        pindel_out_vcf = "${params.SampleId}_pindel_raw.vcf"
        
    """
    cat ${indels_pindel_conf} | xargs \
    pindel -f ${ref_fa} -ii recal.bam -nn ${params.SampleId} -j ${indel_hotspots_bed} -o ${pindel_tmp_bed} -k -T ${task.cpus}
    pindel2vcf -r ${ref_fa} -R b37 -d 20101123 -P ${pindel_tmp_bed} -v ${pindel_out_vcf} -G
    """
}


process sv_svaba_paired {
     container '10.168.2.67:5000/svaba-yy:1.1.2'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus { params.threads/2 + 4 }
     memory { 32.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

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


     cpus { params.threads/2 + 4 }
     memory { 32.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3
     
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
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 2..140 ? 'retry' : 'ignore' }
     maxRetries 5

     when:
         params.ctl_fq1 != ""

     input:
         path "${params.SampleId}_tumor_recal.bam"  from tumor_recal_bam_ch3
         path "${params.SampleId}_tumor_recal.bai" from tumor_recal_bai_ch3
         path 'target.bed' from bed_target

     when:
         params.baseline =~ /^ZB*/

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


process statistic_sinotools_paired {
     container '10.168.2.67:5000/sinotools-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*_qc.txt'
     publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.cov'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 2

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
     maxRetries 2
     
     when:
         params.ctl_fq1 == ""

     input:
         path "${params.SampleId}_tumor_recal.bam"  from tumor_recal_bam_ch4
         path "${params.SampleId}_tumor_recal.bai" from tumor_recal_bai_ch4
         path 'target.bed' from st_bed_ch

     output:
         path "*_qc.txt" into statistic_tumor_ch
         path "*.cov" into coverage_tumor_ch

     script:
         tumor_statistics = "${params.SampleId}_tumor_qc.txt"
     """
     sinotools bam_qc -i ${params.SampleId}_tumor_recal.bam -bed target.bed -hotspot ${core_hotspots_bed} -o ${tumor_statistics}
     """
}


process anno_snpeff_paired {
     container '10.168.2.67:5000/snpeff-yy:4.4m'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus 4
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 != ""

     input:
         path 'raw_snp_indel.vcf' from raw_paired_vcf_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

     output:
         path '*.anno' into haplotype_anno_paired_ch
         path '*_snp_indel_hgvs.vcf' into all_anno_paired_ch

     script:
        prefix_anno_vcf = "${params.SampleId}_raw"

     """
     cat anno_snpeff.conf | xargs \
     java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}.haplotype.anno GRCh37.p13.RefSeq raw_snp_indel.vcf > ${prefix_anno_vcf}_snp_indel_hgvs.vcf

     """
}

process anno_snpeff_dedup_paired {
     container '10.168.2.67:5000/snpeff-yy:4.4m'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus 4
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 != "" && dedup_process_mark == true

     input:
         path 'raw_dedup_snp_indel.vcf' from dedup_raw_paired_vcf_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

     output:
         path '*_snp_indel_hgvs.vcf' into dedup_all_anno_paired_ch

     script:
        prefix_anno_vcf = "${params.SampleId}_raw_dedup"

     """
     cat anno_snpeff.conf | xargs \
     java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}_dedup.haplotype.anno GRCh37.p13.RefSeq raw_dedup_snp_indel.vcf > ${prefix_anno_vcf}_snp_indel_hgvs.vcf 
     """
}


process anno_snpeff_tumor_only {
     container '10.168.2.67:5000/snpeff-yy:4.4m'

     cpus 4
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 == ""

     input:
         path 'raw_snp_indel.vcf' from raw_tumor_vcf_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 
       
     output:
         path '*snp_indel_hgvs.vcf' into tumor_anno_vcf_ch

     script:
        prefix_anno_vcf = "${params.SampleId}_raw"


     """
     cat anno_snpeff.conf | xargs \
     java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}.haplotype.anno GRCh37.p13.RefSeq raw_snp_indel.vcf > ${prefix_anno_vcf}_snp_indel_hgvs.vcf

     """
}

process anno_tumor_dedup_snpeff {
     container '10.168.2.67:5000/snpeff-yy:4.4m'

     cpus 4
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 == "" && dedup_process_mark == true

     input:
         path 'raw_dedup_snp_indel.vcf' from dedup_raw_tumor_vcf_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 
       
     output:
         path '*snp_indel_hgvs.vcf' into dedup_tumor_anno_vcf_ch

     script:
        prefix_anno_vcf = "${params.SampleId}_raw_dedup"

     """
     cat anno_snpeff.conf | xargs \
     java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}_dedup.haplotype.anno GRCh37.p13.RefSeq raw_dedup_snp_indel.vcf > ${prefix_anno_vcf}_snp_indel_hgvs.vcf
     """
}


process anno_database_snpsift_paired {
    container '10.168.2.67:5000/snpsift-yy:1.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*cosmic_clinvar_gnomad_dbsnp_snv_indel_hgvs.vcf'

    cpus 1
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
         params.ctl_fq1 != ""

    input:
         path('hgvs_snp_indel.vcf') from all_anno_paired_ch

    output:
         path '*dbsnp_snv_indel_hgvs.vcf' into database_anno_vcf_paired_ch

    script:
        prefix_anno_vcf = "${params.SampleId}"
    """
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${cosmic_database} hgvs_snp_indel.vcf > ${prefix_anno_vcf}_cosmic_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${clinvar_database} ${prefix_anno_vcf}_cosmic_snv_indel_hgvs.vcf > ${prefix_anno_vcf}_cosmic_clinvar_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${gnomad_database} ${prefix_anno_vcf}_cosmic_clinvar_snv_indel_hgvs.vcf > ${prefix_anno_vcf}_cosmic_clinvar_gnomad_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${dbsnp_database} ${prefix_anno_vcf}_cosmic_clinvar_gnomad_snv_indel_hgvs.vcf > ${prefix_anno_vcf}_cosmic_clinvar_gnomad_dbsnp_snv_indel_hgvs.vcf

    """
}

process anno_dedup_database_snpsift_paired {
    container '10.168.2.67:5000/snpsift-yy:1.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*cosmic_clinvar_gnomad_dbsnp_snv_indel_hgvs.vcf'

    cpus 1
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
         params.ctl_fq1 != "" && dedup_process_mark == true

    input:
         path 'dedup_hgvs_snp_indel.vcf' from dedup_all_anno_paired_ch

    output:
         path '*dbsnp_snv_indel_hgvs.vcf' into dedup_database_anno_vcf_paired_ch

    script:
        prefix_anno_vcf = "${params.SampleId}"
        dedup_prefix_anno_vcf = "${params.SampleId}_dedup"
    """
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${cosmic_database} dedup_hgvs_snp_indel.vcf > ${dedup_prefix_anno_vcf}_cosmic_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${clinvar_database} ${dedup_prefix_anno_vcf}_cosmic_snv_indel_hgvs.vcf > ${dedup_prefix_anno_vcf}_cosmic_clinvar_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${gnomad_database} ${dedup_prefix_anno_vcf}_cosmic_clinvar_snv_indel_hgvs.vcf > ${dedup_prefix_anno_vcf}_cosmic_clinvar_gnomad_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${dbsnp_database} ${dedup_prefix_anno_vcf}_cosmic_clinvar_gnomad_snv_indel_hgvs.vcf > ${dedup_prefix_anno_vcf}_cosmic_clinvar_gnomad_dbsnp_snv_indel_hgvs.vcf       
    """
}
 

process anno_database_snpsift {
    container '10.168.2.67:5000/snpsift-yy:1.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*cosmic_clinvar_gnomad_dbsnp_snv_indel_hgvs.vcf'

    cpus 1
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
         params.ctl_fq1 == ""

    input:
         path('hgvs_snp_indel.vcf') from tumor_anno_vcf_ch

    output:
         path '*dbsnp_snv_indel_hgvs.vcf' into database_anno_vcf_ch

    script:
        prefix_anno_vcf = "${params.SampleId}"
    """
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${cosmic_database} hgvs_snp_indel.vcf > ${prefix_anno_vcf}_cosmic_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${clinvar_database} ${prefix_anno_vcf}_cosmic_snv_indel_hgvs.vcf > ${prefix_anno_vcf}_cosmic_clinvar_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${gnomad_database} ${prefix_anno_vcf}_cosmic_clinvar_snv_indel_hgvs.vcf > ${prefix_anno_vcf}_cosmic_clinvar_gnomad_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${dbsnp_database} ${prefix_anno_vcf}_cosmic_clinvar_gnomad_snv_indel_hgvs.vcf > ${prefix_anno_vcf}_cosmic_clinvar_gnomad_dbsnp_snv_indel_hgvs.vcf

    """
}

process anno_dedup_database_snpsift {
    container '10.168.2.67:5000/snpsift-yy:1.0'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*cosmic_clinvar_gnomad_dbsnp_snv_indel_hgvs.vcf'

    cpus 1
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 2

    when:
         params.ctl_fq1 == "" && dedup_process_mark == true

    input:
         path 'dedup_hgvs_snp_indel.vcf' from dedup_tumor_anno_vcf_ch

    output:
         path '*dbsnp_snv_indel_hgvs.vcf' into dedup_database_anno_vcf_ch

    script:
        prefix_anno_vcf = "${params.SampleId}"
        dedup_prefix_anno_vcf = "${params.SampleId}_dedup"
    """
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${cosmic_database} dedup_hgvs_snp_indel.vcf > ${dedup_prefix_anno_vcf}_cosmic_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${clinvar_database} ${dedup_prefix_anno_vcf}_cosmic_snv_indel_hgvs.vcf > ${dedup_prefix_anno_vcf}_cosmic_clinvar_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${gnomad_database} ${dedup_prefix_anno_vcf}_cosmic_clinvar_snv_indel_hgvs.vcf > ${dedup_prefix_anno_vcf}_cosmic_clinvar_gnomad_snv_indel_hgvs.vcf
    java -Xmx4g -jar /opt/SnpSift.jar annotate ${dbsnp_database} ${dedup_prefix_anno_vcf}_cosmic_clinvar_gnomad_snv_indel_hgvs.vcf > ${dedup_prefix_anno_vcf}_cosmic_clinvar_gnomad_dbsnp_snv_indel_hgvs.vcf       
    """
}
       
process anno_indels_snpeff {
     container '10.168.2.67:5000/snpeff-yy:4.4m'
     publishDir "${raw_unfilter_dir}", mode: 'copy'

     cpus 4
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 2

     input:
         path 'raw_indels.vcf' from pindel_vcf_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

     output:
         path '*_pindel_hgvs.vcf' into pindel_anno_vcf_ch

     script:
         prefix_anno_vcf = "${params.SampleId}"

     """
     cat anno_snpeff.conf | xargs \
     java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
             -haplotype-ouput ${prefix_anno_vcf}.haplotype.anno GRCh37.p13.RefSeq raw_indels.vcf > ${prefix_anno_vcf}_pindel_hgvs.vcf

     """
}


process filter_varscan_report_paired {
     container '10.168.2.67:5000/varscan-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 != ""
            
     input:
         path 'hgvs_snp_indel.vcf' from database_anno_vcf_paired_ch
         path 'variation_hotspots.conf' from variation_hotspots_conf

     output:
         path '*_review_for_report.xls' into small_var_report_paired

     script:
         snv_indel_report_output = "${params.SampleId}_review_for_report.xls"

     """
     python3 /opt/vcf2mut_match_normal.py hgvs_snp_indel.vcf ${variation_hotspots_conf} ${snv_indel_report_output} ${params.PanelName} 
     """
}
 
process filter_paired_dedup_varscan {
     container '10.168.2.67:5000/varscan-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 != "" && dedup_process_mark == true
            
     input:
         path 'dedup_database_hgvs_snp_indel.vcf' from dedup_database_anno_vcf_paired_ch
         path 'variation_hotspots.conf' from variation_hotspots_conf
      
     output:
         path '*_review_for_report.xls' into dedup_var_report_paired

     script:
         dedup_snv_indel_report_output = "${params.SampleId}_dedup_review_for_report.xls"

     """
     python3 /opt/vcf2mut_match_normal.py dedup_database_hgvs_snp_indel.vcf variation_hotspots.conf ${dedup_snv_indel_report_output} ${params.PanelName}
     """
}

process filter_varscan_report_tumor_only {
     container '10.168.2.67:5000/varscan-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 == ""
            
     input:
         path 'database_hgvs_snp_indel.vcf' from database_anno_vcf_ch
         path 'variation_hotspots.conf' from variation_hotspots_conf
      
     output:
         path '*_review_for_report.xls' into small_var_report_ch

     script:
         snv_indel_report_output = "${params.SampleId}_review_for_report.xls"

     """
     python3 /opt/vcf2mut_tumor_only.py database_hgvs_snp_indel.vcf variation_hotspots.conf ${snv_indel_report_output} ${params.PanelName}
     """
}

process filter_tumor_dedup_varscan {
     container '10.168.2.67:5000/varscan-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
     maxRetries 2

     when:
         params.ctl_fq1 == "" && dedup_process_mark == true
            
     input:
         path 'dedup_database_hgvs_snp_indel.vcf' from dedup_database_anno_vcf_ch
         path 'variation_hotspots.conf' from variation_hotspots_conf
      
     output:
         path '*_review_for_report.xls' into dedup_var_report_ch

     script:
         dedup_snv_indel_report_output = "${params.SampleId}_dedup_review_for_report.xls"

     """
     python3 /opt/vcf2mut_tumor_only.py dedup_database_hgvs_snp_indel.vcf variation_hotspots.conf ${dedup_snv_indel_report_output} ${params.PanelName}     
     """
}


process filter_pindel_to_report {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 2
     memory { 4.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 2

     input:
         path 'hgvs_indels.vcf' from pindel_anno_vcf_ch
         path 'indel_hotspots.bed' from indel_hotspots_bed

     output:
         path '*.xls' into pindel_report_ch

     script:
        prefix_anno_vcf = "${params.SampleId}"

     """
     python3 /opt/pindel2report.py hgvs_indels.vcf indel_hotspots.bed ${prefix_anno_vcf}_pindel_review.xls
     """
}


process filter_svaba_review_report {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
     maxRetries 2

     when:
         params.ctl_fq1 == ""

     input:
         path('unfiltered_sv.vcf') from sv_svaba_tumor_ch

     output:
         path '*' into sv_tumor_report_ch

     script:
         sv_report_output = "${params.SampleId}_sv_report.xls"
     """
     python3 /opt/sv2report_review.py unfiltered_sv.vcf ${sv_hotspots_conf} ${sv_partner_conf} ${sv_report_output}
     """
}


process filter_svaba_review_paired {
     container '10.168.2.67:5000/filter2report-yy:2.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'ignore' }
     maxRetries 2

     when:
         params.ctl_fq1 != ""

     input:
         path('somatic_sv.vcf') from sv_svaba_somatic_ch
         path('germline_sv.vcf') from sv_svaba_germline_ch

     output:
         path '*' into sv_paired_report_ch

     script:
         sv_report_output = "${params.SampleId}_sv_report.xls"
     """
     python3 /opt/sv2report_review.py somatic_sv.vcf germline_sv.vcf ${sv_hotspots_conf} ${sv_partner_conf} ${sv_report_output}
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
         path 'germline_control.vcf' from germline_paired_vcf_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

    output:
         path "${germline_anno_vcf}" into germline_anno_ch

    script:
        germline_anno_vcf = "${params.SampleId}_germline.vcf"
    """
    python3 /opt/vcf_multiple_split.py germline_control.vcf germline_control_splitted.vcf
    cat anno_snpeff.conf | xargs \
    java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
                     -haplotype-ouput haplotype.anno GRCh37.p13.RefSeq germline_control_splitted.vcf > ${germline_anno_vcf}
    """
}

process anno_germline_snpeff_tumor {
    container '10.168.2.67:5000/snpeff-yy:4.4m'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'
    publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*.xls'

    cpus 4
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 3

    when:
          params.ctl_fq1 == ""
    
    input:
         path 'germline_tumor.vcf' from germline_tumor_vcf_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

    output:
         path "${germline_tumor_anno_vcf}" into tumor_germline_anno_ch

    script:
        germline_tumor_anno_vcf = "${params.SampleId}_germline_tumor.vcf"
    """
    python3 /opt/vcf_multiple_split.py germline_tumor.vcf germline_tumor_splitted.vcf
    cat anno_snpeff.conf | xargs \
    java -Xmx4g -jar /opt/SnpEff-4.4-jar-with-dependencies.jar -c /opt/snpEff.config \
                     -haplotype-ouput haplotype.anno GRCh37.p13.RefSeq germline_tumor_splitted.vcf > ${germline_tumor_anno_vcf}
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
         path 'germline_hc.conf' from germline_calling_conf

     output:
         path '*.xls' into germline_clinvar_ch

     script:
         clinvar_anno_vcf = "${params.SampleId}_germline_clinvar.vcf"
         clinvar_anno_xls = "${params.SampleId}_germline_clinvar.xls"
     """
     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_database} germline.vcf > ${clinvar_anno_vcf}
     python3 /opt/germline_clinvar_filter2.py ${clinvar_anno_vcf} ${clinvar_anno_xls}
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
         path 'germline_hc.conf' from germline_calling_conf

     output:
         path '*.xls' into tumor_clinvar_ch

     script:
         clinvar_anno_vcf = "${params.SampleId}_tumor_clinvar.vcf"
         clinvar_anno_xls = "${params.SampleId}_tumor_clinvar.xls"
     """
     cat germline_hc.conf
     java -Xmx4g -jar /opt/SnpSift.jar annotate -info CLNSIG,CLNSIGCONF,CLNDN ${clinvar_database} pesudo_germline.vcf > ${clinvar_anno_vcf}
     python3 /opt/germline_clinvar_filter2.py ${clinvar_anno_vcf} ${clinvar_anno_xls}
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
         path('report.xls') from small_var_report_paired  

     output:
         path '*' into tmb_report_paired

     script:
         report_muts_output = "${params.SampleId}_review_muts_TMB_for_report.txt"
         report_tmb_output = "${params.SampleId}_review_TMB_for_report.txt"

     """
     python3 /opt/cal_ctDNA_TMB_match_normal.py report.xls ${params.codinglength} ${report_muts_output} ${report_tmb_output}
     """
}

process calc_paired_dedup_tmb {
     container '10.168.2.67:5000/calc2tmb-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
          params.codinglength != "" && params.ctl_fq1 != "" && dedup_process_mark == true 
       
     input:
         path 'dedup_report.xls' from dedup_var_report_paired

     output:
         path '*' into dedup_tmb_paired_ch

     script:
         dedup_vars_output = "${params.SampleId}_dedup_tmb_vars_review.txt"
         dedup_tmb_output = "${params.SampleId}_paired_dedup_tmb.txt"

     """
     python3 /opt/cal_ctDNA_TMB_match_normal.py dedup_report.xls ${params.codinglength} ${dedup_vars_output} ${dedup_tmb_output}
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
         path('report.xls') from small_var_report_ch

     output:
         path '*' into tmb_report_ch

     script:
         report_muts_output = "${params.SampleId}_review_muts_TMB_for_report.txt"
         report_tmb_output = "${params.SampleId}_review_TMB_for_report.txt"

     """
     python3 /opt/cal_ctDNA_tumor_only.py report.xls ${params.codinglength} ${report_muts_output} ${report_tmb_output}
     """
}

process calc_tumor_dedup_tmb {
     container '10.168.2.67:5000/calc2tmb-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 1
     memory { 2.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 3

     when:
         params.codinglength != "" && params.ctl_fq1 == "" && dedup_process_mark == true
       
     input:
         path('dedup_report.xls') from dedup_var_report_ch

     output:
         path '*' into dedup_tmb_ch

     script:
         dedup_vars_output = "${params.SampleId}_dedup_tmb_vars_review.txt"
         dedup_tmb_output = "${params.SampleId}_dedup_tmb.txt"
     """
     python3 /opt/cal_ctDNA_tumor_only.py dedup_report.xls ${params.codinglength} ${dedup_vars_output} ${dedup_tmb_output}
     """
}

process perbase_complex_standard_mgmt {
     container '10.168.2.67:5000/perbase-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 2
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
     maxRetries 2
     when:
         params.SampleId =~ /_M$/ 

     input:
         path 'realigned.bam' from tumor_recal_bam_ch9
         path 'realigned.bai' from tumor_recal_bai_ch9

     output:
         path "*" into perbase_mgmt_out_ch

     script:
         perbase_out_txt = "${params.SampleId}_perbase_complex_tips.txt"
         mgmt_sta_txt = "${params.SampleId}_mgmt_sta.xls"
         mgmt_methy_txt = "${params.SampleId}_mgmt_methy.xls"
     """
     cat ${perbase_conf} | xargs \
     perbase base-depth realigned.bam -r ${ref_fa} -t ${task.cpus} -o ${perbase_out_txt}
	 python3 /opt/mgmt.py ${ref_fa} ${perbase_out_txt} ${mgmt_sta_txt} ${mgmt_methy_txt}
     """
}

process perbase_complex_standard_1p19q {
     container '10.168.2.67:5000/perbase-yy:1.0'
     publishDir "${filtered_review_dir}", mode: 'copy'

     cpus 2
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
     maxRetries 2
     when:
         params.SampleId =~ /_N$/

     input:
         path 'realigned.bam' from tumor_recal_bam_ch9
         path 'realigned.bai' from tumor_recal_bai_ch9

     output:
         path "*" into perbase_1p19q_out_ch

     script:
         perbase_out_txt = "${params.SampleId}_perbase_complex_tips.txt"
         p19q_sta_txt = "${params.SampleId}_1p19q.xls"
     """
     cat ${perbase_conf} | xargs \
     perbase base-depth realigned.bam -r ${ref_fa} -t ${task.cpus} -o ${perbase_out_txt}
	 python3 /opt/1p19q.py ${rs_1p19q_bed} ${perbase_out_txt} ${p19q_sta_txt} 
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
          path 'recal.bam'  from tumor_recal_bam_ch7
          path 'recal.bai' from tumor_recal_bai_ch7

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

process germline_gatk4_for_tumor {
     container '10.168.2.67:5000/gatk-yy:4.2.0.0'

     cpus 2
     memory { 8.GB * task.attempt }

     errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
     maxRetries 5

     when:
          germline_process_mark == true

     input:
          path 'tumor_recal.bam'  from tumor_recal_bam_ch6
          path 'tumor_recal.bai' from tumor_recal_bai_ch6
          path 'chemo_radio_rs.bed' from chemo_radio_rs_bed

     output:
          path "${tumor_chemo_vcf}" into tumor_germline_ch

     script:
          tumor_chemo_vcf = "${params.SampleId}_germline_for_tumor.vcf"
     """
     gatk HaplotypeCaller -I tumor_recal.bam -R ${ref_fa} \
             -D ${gatk_db_dir}/dbsnp_138.b37.vcf --alleles ${chemo_radio_rs_vcf} \
             -L chemo_radio_rs.bed --max-reads-per-alignment-start 1000 \
             --dont-use-soft-clipped-bases true -stand-call-conf 10 \
             -A AlleleFraction -O ${tumor_chemo_vcf}
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
          path 'ctl_recal.bam' from ctl_recal_bam_ch5
          path 'ctl_recal.bai' from ctl_recal_bai_ch5
          path 'chemo_radio_rs.bed' from chemo_radio_rs_bed

     output:
          path "${germline_vcf}" into control_germline_ch

     script:
          germline_vcf = "${params.SampleId}_germline_for_control.vcf"
     """
     gatk HaplotypeCaller -I ctl_recal.bam -R ${ref_fa} \
             -D ${gatk_db_dir}/dbsnp_138.b37.vcf --alleles ${chemo_radio_rs_vcf} \
             -L chemo_radio_rs.bed --max-reads-per-alignment-start 1000 \
             --dont-use-soft-clipped-bases true -stand-call-conf 10 \
             -A AlleleFraction -O ${germline_vcf}
     """
}

process chemotherapy_anno_germline_snpeff_tumor {
    container '10.168.2.67:5000/snpeff-yy:4.4m'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'
    publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*.xls'

    cpus 4
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
         path 'germline_tumor.vcf' from tumor_germline_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

    output:
         path "${germline_tumor_anno_vcf}" into tumor_germline_rs_ch
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

process chemotherapy_anno_germline_snpeff_control {
    container '10.168.2.67:5000/snpeff-yy:4.4m'
    publishDir "${raw_unfilter_dir}", mode: 'copy', pattern: '*.vcf'
    publishDir "${filtered_review_dir}", mode: 'copy', pattern: '*.xls'

    cpus 4
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'terminate' }
    maxRetries 3

    when:
          params.ctl_fq1 != ""
    
    input:
         path 'germline_control.vcf' from control_germline_ch
         path 'anno_snpeff.conf' from anno_snpeff_conf 

    output:
         path "${germline_anno_vcf}" into germline_control_rs_ch
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


