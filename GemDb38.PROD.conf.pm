#!/usr/bin/perl -w
#/*************************************************************************
# *
# * TGen CONFIDENTIAL
# * __________________
# *
# *  [2010] - [2014] Translational Genomics Research Institute (TGen)
# *  All Rights Reserved.
# *
# * NOTICE:  All information contained herein is, and remains
# * the property of Translational Genomics Research Institute (TGen).
# * The intellectual and technical concepts contained herein are proprietary
# * to  TGen and may be covered by U.S. and Foreign Patents,
# * patents in process, and are protected by trade secret or copyright law.
# * Dissemination of this information, application or reproduction
# * is strictly forbidden unless prior written permission is obtained
# * from TGen.
# *
# * Author (s):
#    David Craig
#    Release 3/23/15
#/
######################################################################################
package Conf;
use Storable qw(dclone);

sub loadDefaults {
    my $function          = shift;
    my $RunParametersOrig = shift;
    my $bin               = $RunParametersOrig->{'binDir'};
    my $ensemblBuild = "GRCh38.86";
    my $startup           = $RunParametersOrig->{'StartupRunParameters'};
    print "\t+Loading Default Parameters\n";


    $RunParameters = {
###############################################################
        #  Mandatory Fields (Frequent Changes)
        #     UidName is array when project
###############################################################
        DatabaseName => "markers",    
        scrubEnabled         => 1,  ##  This will make sure that records aren't removed. 
        RulesFiles   => ['*'],           # All is all
        vcfPaths             => ["/labs/ngd-data/VCF38", "/labs/ngd-data/reports38"],
       #AddOnly=>"crdc", 
###############################################################        
    ###  HOST INFORMATION  #### Easy is to make all the same location
###############################################################
        MongodbHost    => "10.48.66.234:27223",
        AnnotateHost   => "10.48.66.234:27227",
        AnnotateDbName => "AnnotateBy",
        BufferHostList => [ { 'host'        => "diva.tgen.org:27056/buffer", 'connections' => 960 } ], #Computationally intensive
        #BufferHostList => [ { 'host'        => "localhost:27065/buffer", 'connections' => 960 } ], #Computationally intensive

###############################################################
###    LESS COMMON EDITS BELOW HERE                     #######
###############################################################

###############################################################
   #  reportableGenes
   #     Uses 'assay' to determine which genes can be reported on, defaults to first collectionName
###############################################################
        reportableGenes => {
            'A1STL'    => "$bin/resources/db/strexomeLite.v1.txt",
            'A1STX'    => "$bin/resources/db/tumor.txt",
            'K1STX'    => "$bin/resources/db/tumor.txt",
            'K1IDT'    => "$bin/resources/db/tumor.txt",
            'TSMRU'    => "$bin/resources/db/tumor.txt",
            'KAS5U'    => "$bin/resources/db/tumor.txt",
            'KHSCR'    => "$bin/resources/db/tumor.txt",
	    'KAWGL'    => "$bin/resources/db/tumor.txt",
	    'KHWGL'    => "$bin/resources/db/tumor.txt",
            'KHSTX'        => "$bin/resources/db/tumor.txt",
            'KBS5U'        => "$bin/resources/db/germline.txt",
            'KPSTX'        => "$bin/resources/db/tumor.txt",
            'tumor'    => "$bin/resources/db/tumor.txt",
            'germline' => "$bin/resources/db/germline.txt"
        },
        customAnnotation => { 'TERT'    => "$bin/resources/db/customTert.vcf" },
        'assay' => 'germline',
###############################################################
        #  filetypesFromFilename
        #     Requires 'regex' for matching to identify a file from a filename
        #     Requires 'valid' as 1/0 to determine if its a valid file for inserting
        #     Requires 'collectionDb' from Collections->CollectionName for specifying collection
###############################################################
        filetypesFromFilename => {
            'unknown' => {
                regex        => 'vcf$',
                valid        => 1,
                collectionDb => 'tumor'
            },
            'structuralVariantFile' => {
                regex        => '\.manta.somaticSV.vcf',
                valid        => 1,
                collectionDb => 'tumor'
            },
	    'copyNumberFile' => {
		regex	     => 'copy_number.vcf',
		valid 	     => 1,
		collectionDb => 'tumor'
            },
	    'cnvFile' => {
		regex	     => 'seg.vcf',
		valid 	     => 1,
		collectionDb => 'tumor'
            },
	    'freebayesFile' => {
		regex	     => 'tnFreebayes.filt.snpEff.vcf',
		valid 	     => 1,
		collectionDb => 'tumor'
            },
            'strelkaPointMutationFile' => {
                regex        => '(passed.somatic).*vcf',
                valid        => 1,
                collectionDb => 'tumor'
            },
            'mutectPointMutationFile' => {
                regex        => '(MuTect).*vcf$',
                valid        => 1,
                collectionDb => 'tumor'
            },
           'lumosVarPointMutationFile' => {
                regex        => '(lumosVarSNV).*vcf',
                valid        => 1,
                collectionDb => 'tumor'
            },
            'mergedVcfFile' => {
                regex        => 'merged.db.*vcf.gz$',
                valid        => 1,
                collectionDb => 'tumor'
            },
	    'germlineMergedFile' => {
	    	regex        => 'merged\.vcf',
		valid	     => 1,
		collectionDb => 'germline'
	    },
	    'deepVariantGermlineFile' => {
		regex	     => 'deepvariant.pass.*vcf',
		valid	     => 1,
		collectionDb => 'germline'
	    },
            'GATK_Haplocaller' => {
                regex        => '(HC)|(ermline).*vcf$',
                valid        => 1,
                collectionDb => 'germline'
            },
            'GATK_UGFile' => {
                regex        => '(UG)|(gatk)|(GATK).*vcf$',
                valid        => 1,
                collectionDb => 'germline'
            },
            'samtools' => {
                regex        => '\.(C|P|T)[1-9]\..*mpileup.*snp.*vcf$',
                valid        => 1,
                collectionDb => 'germline'
            },
            'thfFile' => {
                regex        => '.*thf\..*vcf$',
                valid        => 1,
                collectionDb => 'tumor'
            },
            'deseqFile' => {
                regex        => '.*deseq.*vcf$',
                valid        => 1,
                joinInName        => 'deseq',
                joinInfoFPKM1     => 2,
                joinGeneInMarkers => 1,
                joinInfoFPKM2     => 2,
                joinInfoLOG2FC     => 2,
                joinInfoPVALUE     => 2,
                joinInfoQVALUE     => 2,                
                collectionDb => 'tumor'
            },
	     'cuffdiffFile' => {
                regex        => '.*cuffdiff.*.vcf$',
                valid        => 1,
                joinInName        => 'cuffdiff',
                joinGeneInMarkers => 1,
                joinInfoLOG2FC     => 2,
                joinInfoPVALUE     => 2,
                joinInfoPADJ     => 2,                
                collectionDb => 'tumor'
            },
            'cnvFile' => {
                regex        => '.*seg.vcf$',
                valid        => 1,
                joinInName        => 'cnv',
                joinGeneInMarkers => 1,
                joinInfoLOG2FC     => 2,
                collectionDb => 'tumor'
            },
            'lumosVarCNVfile' => {
                regex        => '(lumosVarSeg).*vcf$',
                valid        => 1,
                collectionDb => 'tumor',
                joinInName   => 'lumosSeg',
                joinGeneInMarkers => 1,
                joinInfoLOG2FC    => 2
            }
        },
#######  DEFINITION OF REPORTABLE RANGE OCCURS HERE    ############
        
###############################################################
#  Filter Pipelines:
#     - FilterVariantPipeline (filters by variants array)
#     - FilterGeneInfoPipeline (filters by geneInfo sub-doc)
#     - FilterAnnotatePipeline (filters by annotate sub-doc)
#  Notes:
#     -Order should be least strict to most strict
#     -Within a 'type', if any are true the variant is marked
#
#     OPERATORS
#       >,<,<=,>=,== requires field,val1, and oper with numbers
#       NOT,EQUALS,MATCHES  requires field,val1, and oper
#       BETWEEN     requires field,val1, val2, and oper
#       DEF,UNDEF  requires field and oper
#
#     Secondary Operators (can combine Operators)
#       CONDITIONAL requires field,oper,conditOperator,conditVal1,
#                     fieldIfTrue,operIfTrue,val1IfTrue,*val1IfTrue,
#                           fieldIfFalse,val1IfFalse,*val2IfFalse
#       AND requires        field,operForField,val1forField,field2,
#                            operForField2,val1ForField2, *val2ForField2
#
###############################################################

###############################################################
#####  FilterVariantPipeline #####
###############################################################
        InitVariantStatus    => "PASS",
        FilterVariantPipeline => [
            {
                SetVariantToLOWQC => {
                    structuralVariantFile => [
                        {field => 'FILTER', oper => 'EQUALS', val1 =>'MinSomaticScore'},
                    ],
                    mergedVcfFile           => [
                        {field => 'CC', oper => '<=', val1 => 2}
                    ],
                    cnvFile => [
                        {field => 'LOG2FC',    oper => 'BETWEEN', val1 => -1.0, val2 => 1.0},
                    ],
                    lumosVarCNVfile => [
                        {field => 'LOG2FC',    oper => 'BETWEEN', val1 => -1.0, val2 => 1.0},
                        {field => 'SVLEN', oper => '>',       val1 => 25000000}
                    ],
                    GATK_Haplocaller    => [
                         {field => 'qual',oper => '<=', val1 => 500}
                    ],
 		    freebayesFile => [
                        {field => 'QUAL',oper => '<', val1 => 1}
                    ],     
 		    lumosVarPointMutationFile => [
                        {field => 'filter',oper => 'NOT', val1 => 'SomaticPASS'}
                    ],     
                    cuffdiffFile => [
                        {field => 'QVALUE',      oper => '>',   val1 => 0.05},
                        {field => 'SIGNIFICANT', oper => 'NOT', val1 => 'yes'},
                        {field => 'LOG2FC',  oper => 'BETWEEN', val1 => -2, val2 => 2}
                    ],
                    deseqFile => [
			{field => 'PADJ', oper => '>', val1 => 0.05},
                        {field => 'PADJ', oper => 'EQUALS',  val1 => 'NA'},
                        {field => 'LOG2FC',  oper => 'BETWEEN', val1 => -2, val2 => 2}
		    ],
                    thfFile   => [{field => 'qual',   oper => '<', val1 => 100}]
                }
            },
            {
                SetVariantToFAIL => {
                    structuralVariantFile => [
                        {field => 'qual',oper => '<=', val1 => 10},
                        {field => 'Gap', oper => '<=', val1 => 1500},
                        {field => 'NormalAlleleRatio', oper => '>=', val1 => 0.01}
                    ],
 		    freebayesFile => [
                        {field => 'QUAL',oper => '<', val1 => 0}
                    ],     
                    GATK_Haplocaller    => [
                         {field => 'qual',oper => '<=', val1 => 300}
                    ],
                    cnvFile => [
                        {field => 'LOG2FC',    oper => 'BETWEEN', val1 => -0.5, val2 => 0.5},
                        {field => 'cnvLength', oper => '>',       val1 => 25000000}
                    ],
		    lumosVarCNVfile => [
                        {field => 'SVTYPE',    oper => 'EQUALS', val1 => 'NONE'}
                    ],
		    mutectPointMutationFile => [{field => 'QUAL', oper => 'NOT', val1 => 'PASS'}],
                    lumosVarPointMutationFile => [
                        {field => 'filter',oper => 'EQUALS', val1 => 'REJECT'},
                        {field => 'filter',oper => 'EQUALS', val1 => 'GermlineHetPASS'},
                        {field => 'filter',oper => 'EQUALS', val1 => 'GermlineHetLowQC'},
                        {field => 'filter',oper => 'EQUALS', val1 => 'GermlineHomPASS'},
                        {field => 'filter',oper => 'EQUALS', val1 => 'GermlineHomLowQC'},
			{field => 'filter',oper => 'EQUALS', val1 => 'GermlineShiftPASS'},
			{field => 'filter',oper => 'EQUALS', val1 => 'GermlineShiftLowQC'}
		    ],
		    cuffdiffFile => [
                        {field => 'LOG2FC',      oper => 'BETWEEN', val1 => -1.9, val2 => 1.9},
                        {field => 'PVALUE',      oper => '>',       val1 => 0.01},
                        {field => 'LOG2FC', oper => 'EQUALS',  val1 => "inf"},
                        {field => 'LOG2FC', oper => 'EQUALS',  val1 => "-inf"}
                    ],
                    deseqFile => [
                        {field => 'PVALUE', oper => '>', val1 => 0.01},
                        {field => 'PVALUE', oper => 'EQUALS',  val1 => 'NA'},
			{field => 'LOG2FC', oper => 'BETWEEN', val1 => -1.9, val2 => 1.9},
                        {field => 'LOG2FC', oper => 'EQUALS',  val1 => "inf"},
                        {field => 'LOG2FC', oper => 'EQUALS',  val1 => "-inf"}
                    ],
		    thfFile   => [{field => 'qual',   oper => '<=', val1 => 10}]
                }
            }
        ],

        
###############################################################
#####  FilterAnnotatePipeline #####
###############################################################
        FilterAnnotatePipeline => [{
                SetBiomarkerToFAIL => {
                    DatabaseIfFilter      => "FAIL", #Mandatory
                    germline => [{field => 'maxPopFreq',oper => '>', val1 => 0.05}]
#                    #tumor    => [{field => 'maxPopFreq',oper => '>', val1 => 0.01}],
                }                           
       }],
       
###############################################################
#####  FilterGeneInfoPipeline #####
################################################################        FilterGeneInfoPipeline => [{
#                SetBiomarkerToFAIL => {
#                    DatabaseIfFilter      => "FAIL", #Mandatory
#                    germline => [{field => 'name',oper => 'DEF'}],
#                    tumor =>    [{field => 'name',oper => 'REGEX',val1=>'RAF'}],
#                }
#        }],

###############################################################
#####  snpeffValidAnnotations #####
#####  DEFINES REPORTABLE RANGE BY VARIANT TYPE (SNPEFF 4.1)   ######
###############################################################       
        snpeffValidAnnotations => {
            'ALL'                           => 0,
            'chromosome_number_variation'   => 1,
            'exon_loss_variant'             => 1,
            'frameshift_variant'            => 1,
            'stop_gained'                   => 1,
            'stop_lost'                     => 1,
            'start_lost'                    => 1,
            'splice_acceptor_variant'       => 1,
            'splice_donor_variant'          => 1,
            'rare_amino_acid_variant'       => 1,
            'missense_variant'              => 1,
            'inframe_insertion'             => 1,
            'disruptive_inframe_insertion'  => 1,
            'inframe_deletion'              => 1,
            'disruptive_inframe_deletion'   => 1,
            '5_prime_UTR_truncation+exon_loss_variant'=>1,
            '3_prime_UTR_truncation+exon_loss'=>1,
            'splice_branch_variant'         => 1,
            'splice_region_variant'         => 1,
            'splice_region_variant'         => 1,
            'stop_retained_variant'         => 1,
            'initiator_codon_variant'       => 1,
            'synonymous_variant'            => 1,
            'initiator_codon_variant+non_canonical_start_codon'=>1,
            'stop_retained_variant'         => 1,
            '5_prime_UTR_variant'           => 1,
            '3_prime_UTR_variant'           => 1,
            '5_prime_UTR_premature_start_codon_gain_variant'=>1,
            'upstream_gene_variant'         => 0,
            'downstream_gene_variant'       => 0,
            'TF_binding_site_variant'       => 1,
            'regulatory_region_variant'     => 1,
            'miRNA'                         => 1,
            'custom'                        => 0,
            'sequence_feature'              => 1,
            'conserved_intron_variant'      => 1,
            'intron_variant'                => 1,
            'intragenic_variant'            => 0,
            'conserved_intergenic_variant'  => 1,
            'intergenic_region'             => 1,
            'coding_sequence_variant'       => 1,
            'non_coding_exon_variant'       => 1,
            'nc_transcript_variant'         => 1,
            'gene_variant'                  => 1
        },
        reported_aberration_names => {
            'ALL'                           => 0,
            "Missense"                      => 1,
            "High Or Moderate Variant"      => 1,            
            "Initiator Codon Variant"       => 1,
            "Stop Retained Variant"         => 1,
            "Coding Sequence Variant"       => 1,
            "Loss Of Function"              => 1,
            "UTR Truncation"                => 1,
            'Low Impact Change'             => 1,
            'Frameshift'                    => 1,
            'Inframe Indel'                 => 1,
            "Synonymous"                    => 1,
            "Splice-site Loss"              => 1,          
            "Splicing Altered"              => 1,
            "Sequence Feature"              => 1,
            "Structural Variant Breakpoint" => 1,
            "Focal Copy Number Loss"        => 1,
            "Focal Copy Number Gain"        => 1,
            "Other Variant"                 => 1,
            "Large Insertion"               => 1,
            "Large Deletion"                => 1,
            "Over Expressed Gene"           => 1,
            "Under Expressed Gene"          => 1,
            "Fused Genes"                   => 1,
            "Transcript Variant"            => 1,
            'Splicing Altered'              => 1,
            'miRNA'                         => 1,
            "UTR"                           => 1,
            "Intron Variant"                => 1,
            "Upstream-Downstream"           => 1,
            "TF Binding Site"               => 1            
        },
###############################################################        
#######  DEFINITION OF INCLUDED INFO,ANNOTATE, GENEINFO    ####
# 0 - Do not include, 1 - String, 2- Double, 3 - Long
# ALL=1 includes unannotated variants
###############################################################
 
        IncludeInfoFields => {            
            'ALL'                           => 0,
            'DNA_ALT_ALLELE_TOTAL_FRACTION' => 2,
            'RNA_ALT_ALLELE_TOTAL_FRACTION' => 2,
            'LOF'                           => 1,
            'NMD'                           => 1,
            'TumorAlleleRatio'              => 2,
            'SVTYPE'                        => 1,
            'SvCoordinate'                  => 1,
            'TumorAlleleRatio'              => 2,
            'NormalAlleleRatio'             => 2,
            'AR1'                           => 2,
            'AR2'                           => 2,
            'LOG2FC'                        => 2,
            'Gap'                           => 3,
            'GAP'                           => 3,
            'gap'                           => 3,
            'FUSED_GENE'                    => 1,
            'LOF'                           => 1,
            'NMD'                           => 1,
            'SV'                            => 1,
            'AC'                            => 1,
            'AN'                            => 1,
            'exomes_Hom'		    => 2,
	    'exomes_AC'	  	    => 2,
	    'exomes_AF'		    => 2,
	    'exomes_AN'		    => 2,
	    'genomes_AC'	  	    => 2,
	    'genomes_AF'		    => 2,
	    'genomes_AN'		    => 2,
	    'genomes_Hom'		    => 2,
	    'gnomAD_total_AC'		    => 2,
	    'gnomAD_total_AN'		    => 2,
	    'gnomAD_total_HOM'		    => 2,
	    'gnomAD_total_AF'		    => 2,
	    'HOM'			    => 1,
	    'HET'			    => 1,
	    'AF'                            => 2,
            'DP'                            => 3,
            'DP2'                            => 3,
            'DP1'                            => 3,
            'MQ'                            => 1,
            "CNAME"                         => 1,
            'TRANSCRIPT'                    => 1,
            'FUSION'                        => 1,
            'PADJ'                          => 2,
            'FPKM1'                         => 2,
            'FPKM2'                         => 2,
            'PVALUE'                        => 2,
            'SIGNIFICANT'                   => 1,
            'QVALUE'                        => 2,
            'BASEMEANA'                     => 2,
            'BASEMEANB'                     => 2,
            'BPGENE2'                       => 1,
            'BPGENE1'                       => 1,
            'SEURAT_AR_NORMAL'		    => 2,
            'SEURAT_AR_TUMOR'		    => 2,
            'SEURAT_DNA_ALT_ALLELE_REVERSE_FRACTION'=>2,
            'SEURAT_RNA_ALT_ALLELE_TOTAL_FRACTION'=>2,
            'SEURAT_DP_NORMAL'		    =>3,
            'SEURAT_DP_TUMOR'		    =>3,
            'RNA_ALT_FREQ'	            =>2,
	    'RNA_ALT_COUNT'		    =>3,
	    'RNA_REF_COUNT'		    =>3,
	    'CALLERS_COUNT'		    =>3,
            'MUTECT'			    =>3,
            'SEURAT'			    =>3,
            'STRELKA'		 	    =>3,
            'NMD'			    =>1,
            'LOF'			    =>1,
	    'ARDIFF'			    =>2,
            'CloneId'                                   => 3,
            'JPS'                                       => 2,
            'CN'                                        => 2,
            'MACN'                                      => 2,
            'VSF'                                       => 2,
            'CNF'                                       => 2,
            'FT'                                        => 1,
            'ApopAF'                                    => 2,
            'BpopAF'                                    => 2,
            'PT'                                        => 2,        
            'HOMLEN'                           => 2,
            'BND_DEPTH'                           => 2,
            'MATE_BND_DEPTH'                           => 2,
            'CIPOS'                           => 1,
            'MATEID'                           => 1

},            
        
#######  AnnotateByGeneFields   ############                    
        AnnotateByCoordFields => {
        #1  Adds (not forcing of type - as type is defined in the database)        
            'ALL'                         => 1,
            '1000g'                       => 1
        },
        
#######  AnnotateByGeneFields   ############        
        #1  Adds (not forcing of type - as type is defined in the database)                
        AnnotateByGeneFields => {
            'ALL'                         => 0,
            'cgdIntervention'             => 1,
            MGI_mouse_phenotype=>1,
	    pLI=>1,            
	    Gene_full_name=>1,
	    geneFullName=>1,
	    PathwayKegg=>1,
            PathwayConsensus=>1,
            goCellularComponent=>1,
            goBiologicalProcess=>1,
            goMolecularFunction=>1,
	    goOntology=>1,
            refseq=>1,
            orphanet=>1,
            genecards=>1,
            cgdManifestation=>1,
            probabilityHaploinsufficiency=>1,
            probabilityRecessive=>1,
            PathwayConsensus=>1,
            BRCA_TCGA=>1,
            SKCM_TCGA=>1,
            LUAD_TCGA=>1,
            GBM_TCGA=>1,
            'diseaseDescription'          => 1,
            'functionalDescription'       => 1,
            'cgdRef'                      => 1,
            'isca'                        => 1,
            'cgdCondition'                => 1,
            'cgdInheritance'              => 1,
            'panels'                      => 1,
            'pathway'                     => 1,
            'pharmkgb'                    => 1,
            'name'                        => 1,
            'clinvar'                     => 1,
            'level2'                      => 1,
            'cgdAgeGroup'                 => 1,
            'gene'                        => 1
        },
        
###############################################################    
#######  MISC MANDATORY PARAMETERS   ############
###############################################################   
        #    writeVcfTo           => "$bin/GCF/",
        Collections  => [
            {CollectionName => 'tumor', GroupField1    => 'biomarker', GroupField2    => 'projectRun'},
            {CollectionName => 'germline', GroupField1    => 'biomarker', GroupField2    => 'vcfPath'}
        ],
        runDir               => "$bin/runDir",
        Snpeff4Path          => "$bin/resources/snpEff/",
        SnpeffParameter      => "-canon  $ensemblBuild -quiet -noStats -noLog",
        StartupRunParameters => $startup,
        defaultAssay         => 'tumor',
        geneSpecific         => 1,         # Skips locations for SNPs unlikely to be functional by snpeff
        rareFreqThresh       => 0.05,
        gemGeneLocPath       => "$bin/resources/db/GemDb.locations.txt",
        varGlobalMax         => int(50),
        maxInsertsPerRun     => int(1000000),
        largeRun=>1,
        InitVariantStatus    => "PASS",        
        uidName              => "projectRun",
        UidName              => "ProjectRun",    #UID for when 'multiple uids per document allowed        
        loadLocations        => int(1),
        AnnotateByGeneCodon  => 1,
        maxProcesses         => int(8),
        sync                 => 1,
        multipledVariant     => 1    # Allow multiple projectRuns per order
    };
    
    
    my $temp = dclone($RunParameters);
    foreach my $key (keys %{$RunParametersOrig}) {
        $RunParameters->{$key} = $RunParametersOrig->{$key};
    }
    
    $RunParameters = $temp;
    Maintain->reloadStartupRunParameters($RunParameters);
    Maintain->validStartCheck($RunParameters);
    return $RunParameters;
}
1;
