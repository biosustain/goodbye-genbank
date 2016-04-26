DEFAULT_SO_TERM = 'sequence_feature'

# The following feature keys have been discontinued:
GENBANK_DISCONTINUED_FEATURE_KEY_SO_TERMS = {
    # Release Notes For GenBank Release 205
    # http://www.ncbi.nlm.nih.gov/genbank/release/205/
    #     1.3.2 New 'regulatory' feature to replace eleven existing feature types.
    #   As of the December 2015 GenBank release, eleven different features that
    # describe various aspects of regulation have been brought under the umbrella
    # of a single new feature: regulatory . This new feature has a mandatory
    # qualifier (/regulatory_class) to indicate the nature of its regulatory
    # activity.
    #
    #   The existing features that have been replaced include: enhancer, promoter,
    # CAAT_signal, TATA_signal, -35_signal, -10_signal, RBS, GC_signal,
    # polyA_signal, attenuator, and terminator.
    'enhancer': 'enhancer',
    'promoter': 'promoter',
    'CAAT_signal': 'CAAT_signal',
    'TATA_signal': 'TATA_box',
    '-10_signal': 'minus_10_signal',
    '-35_signal': 'minus_35_signal',
    'RBS': 'ribosome_entry_site',
    'GC_signal': 'GC_rich_promoter_region',
    'polyA_signal': 'polyA_signal_sequence',
    'attenuator': 'attenuator',
    'terminator': 'terminator'
}

GENBANK_FEATURE_KEY_SO_TERMS = {
    "3'UTR": 'three_prime_UTR',
    "5'UTR": 'five_prime_UTR',
    'CDS': 'CDS',
    'D-loop': 'D_loop',
    'LTR': 'long_terminal_repeat',
    'STS': 'STS',

    # Immune system related:
    'N_region': 'N_region',
    'S_region': 'S_region',
    'C_region': 'C_region',
    'V_region': 'V_region',
    'J_segment': 'J_gene_segment',
    'D_segment': 'D_gene_segment',
    'V_segment': 'V_gene_segment',

    # XXX what is the difference?
    'assembly_gap': 'gap',

    'conflict': 'sequence_conflict',

    'centromere': 'centromere',

    'exon': 'exon',
    'gap': 'gap',
    'gene': 'gene',
    'iDNA': 'iDNA',
    'intron': 'intron',
    'mRNA': 'mRNA',
    'mat_peptide': 'mature_protein_region',
    'misc_RNA': 'RNA',
    'misc_binding': 'binding_site',
    'misc_difference': 'sequence_difference',
    'misc_feature': 'sequence_feature',
    'misc_recomb': 'recombination_feature',

    # misc_signal
    # GenBank definition:
    #  any region containing a signal controlling or altering gene
    #  function or expression that cannot be described by other
    #  signal keys (promoter, CAAT_signal, TATA_signal,
    #  -35_signal, -10_signal, GC_signal, RBS, polyA_signal,
    #  enhancer, attenuator, terminator, and rep_origin).
    #
    # Using regulatory_region, but unsure.
    # regulatory_region: A region of sequence that is involved in the control of a biological process.
    'misc_signal': 'regulatory_region',

    # misc_structure
    # GenBank definition:
    #  any secondary or tertiary nucleotide structure or
    #  conformation that cannot be described by other Structure
    #  keys (stem_loop and D-loop);

    # XXX Sequence Ontology term used here may be incorrect.
    #     It can be secondary structure of the sequence or tertiary structure of the protein.
    'misc_structure': 'sequence_secondary_structure',

    # mobile_element
    # GenBank definition:
    #  region of genome containing mobile elements;
    #
    # mobile_element could also be mobile_element_insertion.
    'mobile_element': 'mobile_genetic_element',
    'modified_base': 'modified_DNA_base',
    'ncRNA': 'ncRNA_gene',

    #
    # the presented sequence revises a previous version of the sequence at this location;
    # Undefined in Sequence Ontology.
    'old_sequence': 'polypeptide_sequencing_information',

    'operon': 'operon',
    'oriT': 'oriT',
    'polyA_site': 'polyA_site',
    'precursor_RNA': 'primary_transcript',
    'prim_transcript': 'primary_transcript',
    'primer_bind': 'primer_binding_site',
    'protein_bind': 'protein_binding_site',
    'rRNA': 'rRNA',
    'rep_origin': 'origin_of_replication',
    'repeat_region': 'repeat_region',
    'sig_peptide': 'signal_peptide',

    # Might want to exclude this one?
    'source': 'databank_entry',

    'stem_loop': 'stem_loop',
    'tRNA': 'tRNA',
    'telomere': 'telomere',
    'tmRNA': 'tmRNA',
    'transit_peptide': 'transit_peptide',

    # unsure
    # GenBank definition:
    #  a small region of sequenced bases, generally 10 or fewer in its length, which
    #  could not be confidently identified. Such a region might contain called bases
    #  (A, T, G, or C), or a mixture of called-bases and uncalled-bases ('N').
    #  The unsure feature should not be used when annotating gaps in genome assemblies.
    #  Please refer to assembly_gap feature for gaps within the sequence of an assembled
    #  genome. For annotation of gaps in other sequences than assembled genomes use the
    #  gap feature.
    #
    # sequence_uncertainty
    # Sequence Ontology definition:
    #  Describes the positions in a sequence where the authors are unsure about the
    #  sequence assignment.
    'unsure': 'sequence_uncertainty',

    # variation
    # GenBank definition:
    #  a related strain contains stable mutations from the same
    #  gene (e.g., RFLPs, polymorphisms, etc.) which differ
    #  from the presented sequence at this location (and
    #  possibly others);
    #
    # sequence_variant
    # Sequence Ontology definition:
    #  A sequence_variant is a non exact copy of a sequence_feature
    #  or genome exhibiting one or more sequence_alteration.
    'variation': 'sequence_variant',
}

SO_TERM_GENBANK_FEATURE_KEYS = {v: k for k, v in GENBANK_FEATURE_KEY_SO_TERMS.items()}

GENBANK_FEATURE_KEYS = tuple(GENBANK_FEATURE_KEY_SO_TERMS.keys())

UNAMBIGUOUS_INVALID_KEY_SO_TERMS = {
    'primer': 'primer_binding_site',
}

GENBANK_FEATURE_KEY_SO_TERMS.update(GENBANK_DISCONTINUED_FEATURE_KEY_SO_TERMS)

GENBANK_REGULATORY_DEFAULT_SO_TERM = 'regulatory_region'

GENBANK_REGULATORY_CLASS_SO_TERMS = {
    "attenuator": "attenuator",
    "enhancer": "enhancer",

    # GB: transcriptional cis regulatory region that prevents the enhancer from modulating the expression of the gene,
    #     but is not an insulator.
    "enhancer_blocking_element": None,  # TODO what about this one?

    "CAAT_signal": "CAAT_signal",
    'GC_signal': "GC_rich_promoter_region",

    "imprinting_control_region": None,  # TODO what about this one?

    "insulator": "insulator",
    "locus_control_region": "locus_control_region",

    "minus_35_signal": "minus_35_signal",
    "minus_10_signal": "minus_10_signal",

    "polyA_signal_sequence": "polyA_signal_sequence",
    "promoter": "promoter",

    "ribosome_binding_site": "ribosome_entry_site",

    "riboswitch": "riboswitch",
    "silencer": "silencer",

    "TATA_box": "TATA_box",
    "terminator": "terminator",
    "other": GENBANK_REGULATORY_DEFAULT_SO_TERM
}

SO_TERM_GENBANK_FEATURE_KEYS.update({v: GENBANK_REGULATORY_DEFAULT_SO_TERM
                                     for k, v in GENBANK_REGULATORY_CLASS_SO_TERMS.items()})

GENBANK_NC_RNA_DEFAULT_SO_TERM = "ncRNA_gene"

GENBANK_NC_RNA_CLASS_SO_TERMS = {
    "snoRNA": "snoRNA_gene",
    "gRNA": "gRNA_gene",
    "tmRNA": "tmRNA_gene",
    "lincRNA": "lincRNA_gene",
    "tRNA": "tRNA_gene",
    "rRNA": "rRNA_gene",
    "miRNA": "miRNA_gene",
    "piRNA": "piRNA_gene",
    "scRNA": "scRNA_gene",

    "SRP_RNA": "SRP_RNA_gene",
    "RNase_P_RNA": "RNase_P_RNA_gene",
    "RNase_MRP_RNA": "RNase_MRP_RNA_gene",
    "telomerase_RNAe": "telomerase_RNA_gene",

    # XXX Genes for these are missing from Sequence Ontology:
    "antisense_RNA": "antisense_RNA",
    "siRNA": "siRNA",

    "other": GENBANK_NC_RNA_DEFAULT_SO_TERM
}

SO_TERM_GENBANK_FEATURE_KEYS.update({v: GENBANK_NC_RNA_DEFAULT_SO_TERM
                                     for k, v in GENBANK_NC_RNA_CLASS_SO_TERMS.items()})

GENBANK_PSEUDOGENE_DEFAULT_SO_TERM = "pseudogene"

GENBANK_PSEUDOGENE_TYPE_SO_TERMS = {
    "processed": "processed_pseudogene",
    "unprocessed": "non_processed_pseudogene",
    "unitary": "unitary_pseudogene",
    "allelic": "allelically_excluded_gene",
    "unknown": GENBANK_PSEUDOGENE_DEFAULT_SO_TERM,
}

SO_TERM_GENBANK_FEATURE_KEYS.update({v: GENBANK_PSEUDOGENE_DEFAULT_SO_TERM
                                     for k, v in GENBANK_PSEUDOGENE_TYPE_SO_TERMS.items()})

REDUCED_SO_TERMS = {}

REDUCED_SO_TERMS.update({
                            term: GENBANK_REGULATORY_DEFAULT_SO_TERM
                            for term in GENBANK_REGULATORY_CLASS_SO_TERMS.values()
                            })

DEFAULT_GENBANK_FEATURE_KEY = 'misc_feature'
