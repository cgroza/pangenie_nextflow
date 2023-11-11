params.reads          = "reads.csv"
params.reference      = "reference.fa"
params.vcf            = "variants.vcf"
params.out            = "out"

process preprocess {
  cpus 40
  memory "180G"
  time "3h"

  input:
  path(vcf)
  path(fasta)

  output:
  path("index")

  script:
  """
  mkdir index
  PanGenie-index -v ${vcf} -r ${fasta} -t 10 -o index/processed
  """
}

process pangenie {
  cpus 20
  memory "90G"
  time "6h"
  publishDir "${params.out}/genotypes", mode: 'copy'

  input:
  tuple val(sample_name), path(sample_bam), path(ref), path("index/")

  output:
  path("${sample_name}_genotyping.vcf.gz")

  script:
  """
  PanGenie -t 20 -j 20 -s ${sample_name} -i <(samtools fastq ${sample_bam}) -f index/index/processed -o ${sample_name}
  bgzip ${sample_name}_genotyping.vcf
  """
}

workflow {
  Channel.fromPath(params.reference).set{ref_ch}
  Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path, checkIfExists:false)]}.set{reads_ch}

  Channel.fromPath(params.vcf).set{vcf_ch}

  preprocess(vcf_ch, ref_ch).set{index_ch}
  pangenie(reads_ch.combine(ref_ch).combine(index_ch))
}
