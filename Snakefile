builds = ['v2.1.1-b37', 'v3.1.1-b38']
ans = [1,10,20,30]

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from urllib.request import urlopen
from urllib.error import URLError

try:
  response = urlopen('https://www.google.com/', timeout=10)
  iconnect = True
except urllib.error.URLError as ex:
  iconnect = False

class dummyprovider:
  def remote(string_, allow_redirects = "foo"):
    return string_

HTTP = HTTPRemoteProvider() if iconnect else dummyprovider

localrules: download_gnomadv3, download_gnomadv2

rule all:
  input:
    expand('output/aj-enriched_{build}_an-ge{an}.png', build=builds, an=ans),
    expand('output/aj-enriched_{build}_an-ge{an}.stats.tsv',
           build=builds, an=ans)

rule download_gnomadv3:
  input:
    vcf = HTTP.remote('https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr{chrom}.vcf.bgz', allow_redirects=True),
    tbi = HTTP.remote('https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr{chrom}.vcf.bgz.tbi', allow_redirects=True)
  output:
    vcf = 'gnomad_AF/gnomad_v3.1.1-b38_AF.chr{chrom}.vcf.gz',
    tbi = 'gnomad_AF/gnomad_v3.1.1-b38_AF.chr{chrom}.vcf.gz.tbi'
  shell: '''
cp {input.vcf} {output.vcf}
cp {input.tbi} {output.tbi}
'''

rule download_gnomadv2:
  input:
    vcf = HTTP.remote('https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz', allow_redirects=True),
    tbi = HTTP.remote('https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz.tbi', allow_redirects=True)
  output:
    vcf = 'gnomad_AF/gnomad_v2.1.1-b37_AF.chr{chrom}.vcf.gz',
    tbi = 'gnomad_AF/gnomad_v2.1.1-b37_AF.chr{chrom}.vcf.gz.tbi'
  shell: '''
cp {input.vcf} {output.vcf}
cp {input.tbi} {output.tbi}
'''

rule get_enriched:
  input: "gnomad_AF/gnomad_{build}_AF.chr{chrom}.vcf.gz"
  output: "output/aj-enriched_{build}_chr{chrom}.vcf.gz"
  conda: 'vcfenv.yaml'
  shell: '''
bcftools view -v snps -f .,PASS \
 -i 'AF_asj >= 0.1 && AF_asj <= 0.9 && (AF_nfe < 0.01 || AF_nfe > 0.99)' \
 -Oz -o {output} {input}'''

rule tab_enriched:
  input: rules.get_enriched.output
  output: temp('output/aj-enriched_{build}_chr{chrom}.noheader.tsv')
  conda: 'vcfenv.yaml'
  shell: r'''
bcftools query -f '%CHROM:%POS:%REF:%ALT\t%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_asj\t%INFO/AF_nfe\t%INFO/AN_asj\t%INFO/popmax\n' \
  {input} > {output}
'''

rule merge_tabs:
  input: expand('output/aj-enriched_{{build}}_chr{chrom}.noheader.tsv', chrom=range(1,23))
  output: 'output/aj-enriched_{build}_chrall.tsv'
  shell: '''
awk 'BEGIN {{FS=OFS="\t"; print "CPRA","ID","CHROM","POS","REF","ALT","AF_asj","AF_nfe","AN_asj","popmax"}} 1 {{print}}' \
  {input} > {output}
'''

rule make_miami:
  input: rules.merge_tabs.output
  output: 'output/aj-enriched_{build}_an-ge{an}.png'
  conda: 'plotenv.yaml'
  script: 'miami.R'

rule make_stats:
  input: rules.merge_tabs.output
  output: 'output/aj-enriched_{build}_an-ge{an}.stats.tsv'
  conda: 'plotenv.yaml'
  script: 'stats.R'
