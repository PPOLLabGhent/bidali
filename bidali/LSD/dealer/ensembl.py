# -*- coding: utf-8 -*-
"""ensembl

Reference: https://www.ensembl.org/
"""
from bidali import LSD
from bidali.LSD import retrieveSources,cacheableTable,processedDataStorage,datadir
import os, gzip, pandas as pd
from io import TextIOWrapper, StringIO
from urllib.parse import parse_qsl

## Biomart
@cacheableTable
def get_biomart(atts=None,dataset='hsapiens_gene_ensembl'):
    """
    Get biomart id mappings

    Attributes to specify can be:

    name description
    1 ensembl_gene_id Ensembl Gene ID
    2 ensembl_transcript_id Ensembl Transcript ID
    3 ensembl_peptide_id Ensembl Protein ID
    4 canonical_transcript_stable_id Canonical transcript stable ID(s)
    5 description Description
    6 chromosome_name Chromosome Name
    7 start_position Gene Start (bp)
    8 end_position Gene End (bp)
    9 strand Strand
    10 band Band
    11 transcript_start Transcript Start (bp)
    12 transcript_end Transcript End (bp)
    13 external_gene_id Associated Gene Name
    14 external_transcript_id Associated Transcript Name
    15 external_gene_db Associated Gene DB
    16 transcript_db_name Associated Transcript DB
    17 transcript_count Transcript count
    18 percentage_gc_content % GC content
    19 gene_biotype Gene Biotype
    20 transcript_biotype Transcript Biotype
    21 source Source
    22 status Status (gene)
    23 transcript_status Status (transcript)
    24 go_biological_process_id GO Term Accession (bp)
    25 name_1006 GO Term Name (bp)
    26 definition_1006 GO Term Definition (bp)
    27 go_biological_process_linkage_type GO Term Evidence Code (bp)
    28 go_cellular_component_id GO Term Accession (cc)
    29 go_cellular_component__dm_name_1006 GO Term Name (cc)
    30 go_cellular_component__dm_definition_1006 GO Term Definition (cc)
    31 go_cellular_component_linkage_type GO Term Evidence Code (cc)
    32 go_molecular_function_id GO Term Accession
    33 go_molecular_function__dm_name_1006 GO Term Name (mf)
    34 go_molecular_function__dm_definition_1006 GO Term Definition (mf)
    35 go_molecular_function_linkage_type GO Term Evidence Code (mf)
    36 goslim_goa_accession GOSlim GOA Accession(s)
    37 goslim_goa_description GOSlim GOA Description
    38 ucsc UCSC ID
    39 pdb PDB ID
    40 clone_based_ensembl_gene_name Clone based Ensembl gene name
    41 clone_based_ensembl_transcript_name Clone based Ensembl transcript name
    42 clone_based_vega_gene_name Clone based VEGA gene name
    43 clone_based_vega_transcript_name Clone based VEGA transcript name
    44 ccds CCDS ID
    45 embl EMBL (Genbank) ID
    46 ox_ens_lrg_transcript__dm_dbprimary_acc_1074 Ensembl LRG transcript
    47 entrezgene EntrezGene ID
    48 ottt VEGA transcript ID(s) (OTTT)
    49 ottg VEGA gene ID(s) (OTTG)
    50 shares_cds_with_enst Ensembl transcript (where OTTT shares CDS with ENST)
    51 shares_cds_with_ottt HAVANA transcript (where ENST shares CDS with OTTT)
    52 shares_cds_and_utr_with_ottt HAVANA transcript (where ENST identical to OTTT)
    53 hgnc_id HGNC ID
    54 hgnc_symbol HGNC symbol
    55 hgnc_automatic_gene_name HGNC automatic gene name
    56 hgnc_curated_gene_name HGNC curated gene name
    57 hgnc_automatic_transcript_name HGNC automatic transcript name
    58 hgnc_curated_transcript_name HGNC curated transcript name
    59 hgnc_mb001 HGNC mb001 ID
    60 ipi IPI ID
    61 merops MEROPS ID
    62 mim_morbid_accession MIM Morbid Accession
    63 mim_morbid_description MIM Morbid Description
    64 mim_gene_accession MIM Gene Accession
    65 mim_gene_description MIM Gene Description
    66 mirbase_accession miRBase Accession(s)
    67 mirbase_id miRBase ID(s)
    68 protein_id Protein (Genbank) ID
    69 refseq_dna RefSeq DNA ID
    70 refseq_dna_predicted RefSeq Predicted DNA ID
    71 refseq_peptide RefSeq Protein ID
    72 refseq_peptide_predicted RefSeq Predicted Protein ID
    73 refseq_genomic RefSeq Genomic ID(s)
    74 rfam Rfam ID
    75 unigene Unigene ID
    76 uniprot_sptrembl UniProt/TrEMBL Accession
    77 uniprot_swissprot_accession UniProt/SwissProt Accession
    78 wikigene_name WikiGene name
    79 wikigene_description WikiGene description
    80 hpa Human Protein Atlas Antibody ID
    81 dbass3_id Database of Aberrant 3' Splice Sites (DBASS3) IDs
    82 dbass3_name DBASS3 Gene Name
    83 dbass5_id Database of Aberrant 5' Splice Sites (DBASS5) IDs
    84 dbass5_name DBASS5 Gene Name
    85 affy_hc_g110 Affy HC G110
    86 affy_hg_focus Affy HG FOCUS
    87 affy_hg_u133_plus_2 Affy HG U133-PLUS-2
    88 affy_hg_u133a_2 Affy HG U133A_2
    89 affy_hg_u133a Affy HG U133A
    90 affy_hg_u133b Affy HG U133B
    91 affy_hg_u95av2 Affy HG U95AV2
    92 affy_hg_u95b Affy HG U95B
    93 affy_hg_u95c Affy HG U95C
    94 affy_hg_u95d Affy HG U95D
    95 affy_hg_u95e Affy HG U95E
    96 affy_hg_u95a Affy HG U95A
    97 affy_hugenefl Affy HuGene FL
    98 affy_huex_1_0_st_v2 Affy HuEx 1_0 st v2
    99 affy_hugene_1_0_st_v1 Affy HuGene 1_0 st v1
    100 affy_u133_x3p Affy U133 X3P
    101 agilent_cgh_44b Agilent CGH 44b
    102 agilent_wholegenome Agilent WholeGenome
    103 codelink Codelink
    104 illumina_humanwg_6_v1 Illumina HumanWG 6 v1
    105 illumina_humanwg_6_v2 Illumina HumanWG 6 v2
    106 illumina_humanwg_6_v3 Illumina HumanWG 6 v3
    107 illumina_humanht_12 Illumina Human HT 12
    108 phalanx_onearray Phalanx OneArray
    109 anatomical_system Anatomical System (egenetics)
    110 development_stage Development Stage (egenetics)
    111 cell_type Cell Type (egenetics)
    112 pathology Pathology (egenetics)
    113 anatomical_system_gnf Anatomical System (gnf)
    114 development_stage_gnf Development Stage (gnf)
    115 cell_type_gnf Cell Type (gnf)
    116 pathology_gnf Pathology (gnf)
    117 family_description Ensembl Family Description
    118 family Ensembl Protein Family ID(s)
    119 pirsf PIRSF SuperFamily ID
    120 superfamily Superfamily ID
    121 smart SMART ID
    122 profile PROFILE ID
    123 prosite PROSITE ID
    124 prints PRINTS ID
    125 pfam PFAM ID
    126 tigrfam TIGRFam ID
    127 interpro Interpro ID
    128 interpro_short_description Interpro Short Description
    129 interpro_description Interpro Description
    130 transmembrane_domain Transmembrane domain
    131 signal_domain Signal domain
    132 ncoils
    
    Info:
        For further possible attributes, see attribute table in
        https://www.bioconductor.org/packages/2.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf


    >>> bmrt = get_biomart()
    """
    import biomart
    server = biomart.BiomartServer("http://www.ensembl.org/biomart")
    # default atts
    if not atts:
        atts = ['external_gene_name','external_gene_source','ensembl_gene_id',
                'ensembl_transcript_id','ensembl_peptide_id']
    seda = server.datasets[dataset]
    s = seda.search({'attributes': atts}, header=1)
    data = pd.read_table(StringIO(s.content.decode()))
    return data

## Human resources
def get_ensembl(onlyGeneLabeled=True,onlyInChromosomes=None):
    import warnings
    warnings.warn("deprecated, use get_biomart instead", DeprecationWarning)
    ensembl = get_biomart(
        atts=[
            'ensembl_gene_id','ensembl_transcript_id','start_position','end_position',
            'percentage_gene_gc_content','chromosome_name','strand','transcript_start',
            'gene_biotype','transcript_biotype','external_gene_name','entrezgene'
            ]
        )
    ensembl.columns = ['egid','etid','start','stop','gcC','chr','strand','TSS','typeg','typet','gene_label','entrez']
    ensembl.set_index('egid',inplace=True)
    ensembl = ensembl[~ensembl.index.duplicated()]
    if onlyGeneLabeled: ensembl = ensembl[ensembl.gene_label.isnull().apply(lambda x: not(x))]
    ensembl.chr = ensembl.chr.apply(lambda x: 'chr'+x.replace('MT','M'))
    if onlyInChromosomes: ensembl = ensembl[ensembl.chr.isin(onlyInChromosomes)]
    return ensembl

@retrieveSources
def get_ensemblGeneannot():
    """
    Info: http://www.ensembl.org/info/data/ftp/index.html
    Source: ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz
    """
    import gffutils
    try: db = gffutils.FeatureDB(processedDataStorage+'Homo_sapiens.GRCh38.88.sqlite3')
    except ValueError:
        if not os.path.exists(processedDataStorage+'Homo_sapiens.GRCh38.88.gtf.gz'):
            raise FileNotFoundError
        db = gffutils.create_db(processedDataStorage+'Homo_sapiens.GRCh38.88.gtf.gz',
                                processedDataStorage+'Homo_sapiens.GRCh38.88.sqlite3',
                                disable_infer_genes=True,disable_infer_transcripts=True)
    return db

## Cross species resources
### Mouse
@retrieveSources
def get_mouseEnsemblSet():
    """
    Source: ftp://ftp.ensembl.org/pub/release-90/gff3/mus_musculus/Mus_musculus.GRCm38.90.gff3.gz
    """
    db = pd.read_table(TextIOWrapper(gzip.open(processedDataStorage+'Mus_musculus.GRCm38.90.gff3.gz')),comment='#',
                       names='seqid,source,type,start,end,score,strand,phase,attribute'.split(','),low_memory=False)
    genes = db[db.type == 'gene'].copy()
    transcripts = db[db.type == 'mRNA'].copy()
    del db
    transcripts['TRANSCRIPT_ID'] = transcripts.attribute.apply(lambda x: parse_qsl(x)[0][1].split(':')[1])
    transcripts['GENE_ID'] = transcripts.attribute.apply(lambda x: parse_qsl(x)[1][1].split(':')[1])
    transcripts.set_index('TRANSCRIPT_ID',inplace=True)
    return (genes,transcripts)
