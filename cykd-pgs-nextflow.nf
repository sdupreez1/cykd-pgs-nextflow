include { getCohortData } from './modules/getCohortData.nf'
include { removeRelatedWithinCohorts; removeRelatedBetweenCohorts } from './modules/removeRelated.nf'
include { getGenomicPrincipalComponenets; ancestryMatching } from './modules/ancestryMatching.nf'
include { getLowMAFGenomicData } from './modules/getLowMAFGenomicData.nf'
include { finalCohortStats } from './modules/finalCohortStats.nf'
include { formatHarmonisedOPGS } from './modules/formatHarmonisedOPGS.nf'
include { calculatePGS; getRefPanelPlatekeys } from './modules/calculatePGS.nf'
include { chunkPvals } from './modules/chunkPvals.nf'
include { platformFDR } from './modules/platformFDR.nf'
include { filterSigPGSbyR2 } from './modules/filterSigPGSbyR2.nf'

workflow {
    def cohort_ch = Channel.of("control", "case")
    def filter_data_ch = Channel.fromPath("${projectDir}/misc/data_for_cohort_filters/*")
    getCohortData(cohort_ch, filter_data_ch.collect())

    removeRelatedWithinCohorts(getCohortData.out.related_cohort_platekeys)
    removeRelatedBetweenCohorts(removeRelatedWithinCohorts.out.collect())

    getGenomicPrincipalComponenets(removeRelatedBetweenCohorts.out)
    ancestryMatching(getGenomicPrincipalComponenets.out)

    finalCohortStats(ancestryMatching.out.platekeys)

    def extract_per_chr_ch = Channel.of(1..22, 'X').combine(ancestryMatching.out.platekeys)
    getLowMAFGenomicData(extract_per_chr_ch) 

    def platform_ch = Channel.of('Metabolon', 'Nightingale', 'Olink', 'RNAseq',
                                 'SomaScan', 'UKBEuropean', 'UKBMultiancestry')
    def unformatted_files_ch = platform_ch
               | flatMap { platform -> 
                  file("${params.OPGS_folders_dir}/${platform}/OPGS*")
                      .sort { a, b -> a.name <=> b.name }
                      .collate(params.chunk_size) 
                      .collect { files -> [platform, files] }
               }
               | map { platform, files ->
                  def filtered_files = files.findAll { file ->
                      def matches = file.name =~ /OPGS(\d{6})/
                      if (!matches) println "WARNING: Unexpected file name encountered and skipped ${file.name}"
                      matches    // This is the condition that findAll checks is truthy
                  }
                  def ids = filtered_files.collect { file ->
                      (file.name =~ /OPGS(\d{6})/)[0][1]
                  }
                  def from_file = ids.min()
                  def to_file   = ids.max()
                  def chunk_id  = "${platform}${from_file}to${to_file}"
                  [chunk_id, platform, filtered_files]
               }
    formatHarmonisedOPGS(unformatted_files_ch)

    formatHarmonisedOPGS.out.formatted_files 
        | combine(getLowMAFGenomicData.out.bfiles.collect() | map { genomic_data -> [genomic_data] }) 
        | calculatePGS

    def pgs = calculatePGS.out.first()
    getRefPanelPlatekeys(pgs, ancestryMatching.out.platekeys)

    def traits_with_chunk_id_ch = formatHarmonisedOPGS.out.traits
        | map {trait_file ->
            def chunk_id = (trait_file.name =~ /^([a-zA-Z]+\d+to\d+)_traits.tsv/)[0][1] 
            [chunk_id, trait_file]
        }

    calculatePGS.out
        | map { chunk_id ->  
            def matcher = (chunk_id =~ /^([A-Za-z]+)\d+to\d+$/)
            if (!matcher.matches()) {
                error "calculatePGS chunk_id '${chunk_id}' does not match expected pattern"
            }
            def platform = matcher[0][1]
            [chunk_id, platform]
        }
        | join(traits_with_chunk_id_ch)
        | combine(getRefPanelPlatekeys.out)
        | combine(ancestryMatching.out.platekeys)
        | combine(finalCohortStats.out.platekey_sex_age)
        | chunkPvals    // Outputs tuple val("${platform}"), path("${chunk_id}_pvals.tsv")
        | groupTuple    // Groups chunkPvals outputs by platform
        | platformFDR
   
    filterSigPGSbyR2(platformFDRandPlotting.out.sig_pgs.collect())
}
