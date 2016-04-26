import re

from gbgb.maps.features import GENBANK_FEATURE_KEYS
from gbgb.utils import single, as_set

RE_PROTEIN_ID = re.compile(r'^[A-Z]{3}[0-9]{5}\.[0-9]+$')


def remove_protein_id_and_add_to_xrefs(before, after):
    """
    Protein IDs are codes such as "AAF19666.1", which come from "International collaborators" and should all be
    on GenBank. This translation function removes the `/protein_id=""` qualifier and instead adds a "/db_xref=""`
    qualifier.

    Removes any malformed protein IDs as they are useless.

    :param before:
    :param after:
    :return:
    """
    if 'protein_id' in before:
        protein_id = single(before, 'protein_id')

        if RE_PROTEIN_ID.match(protein_id):
            protein_xref = 'GenBank:{}'.format(protein_id)
            after['db_xref'] = as_set(after.get('db_xref')) | {protein_xref}

        del after['protein_id']
    return after


def rename_label_to_note(before, after):
    """
    The `/label=""` qualifier was discontinued in 2010, but is still used frequently.

    See `GenBank Release 180 <http://www.ncbi.nlm.nih.gov/genbank/release/180/>`_.

    :param before:
    :param after:
    :return:
    """
    if 'label' in before:
        # TODO support cases where both a note and a label exist
        after['note'] = single(before, 'label')
    return after


GENBANK_SINGLE_QUALIFIERS = (
    'gene',
    'locus_tag',
    'codon_start',
    'protein_id',
    'note',
    'pseudogene',
    'product',
    'experiment',
    'transl_table',
    'rpt_family',
)


def format_single_qualifiers(before, after):
    """

    :param before:
    :param after:
    :return:
    """
    for name in before.keys():
        if name in GENBANK_SINGLE_QUALIFIERS:
            after[name] = single(before, name, on_multiple_ignore=True)
    return after


GENBANK_INTEGER_QUALIFIERS = (
    'transl_table',
    'codon_start'
)


def format_integer_qualifiers(before, after):
    """
    Formats qualifiers that must be integers as ``int`` objects.

    Removes any malformed qualifiers.

    :param before:
    :param after:
    :return:
    """
    for name in before.keys():
        if name in GENBANK_INTEGER_QUALIFIERS:
            try:
                after[name] = int(single(before, name, on_multiple_ignore=True))
            except ValueError:
                del after[name]
    return after


GENBANK_FLAG_QUALIFIERS = (
    'environmental_sample',
    'focus',
    'germline',
    'macronuclear',
    'partial',
    'proviral',
    'pseudo',
    'rearranged',
    'ribosomal_slippage',
    'transgenic',
    'trans_splicing',
)


def format_flag_qualifiers(before, after):
    for name, value in before.items():
        if name in GENBANK_FLAG_QUALIFIERS:
            after[name] = True
    return after


APE_QUALIFIER_NAMESPACE = 'ApEinfo_'


def remove_ape_a_plasmid_editor_qualifiers(before, after):
    for name, value in before.items():
        if name.startswith(APE_QUALIFIER_NAMESPACE):
            del after[name]
    return after


GENE_RELATED_FEATURES = [
    'C_region',
    'CDS',
    'D-loop',
    'D_segment',
    'exon',
    'gene',
    'iDNA',
    'intron',
    'J_segment',
    'LTR',
    'mat_peptide',
    'misc_binding',
    'misc_difference',
    'misc_feature',
    'misc_recomb',
    'misc_RNA',
    'misc_structure',
    'mobile_element',
    'modified_base',
    'mRNA',
    'ncRNA',
    'N_region',
    'old_sequence',
    'oriT',
    'polyA_site',
    'precursor_RNA',
    'prim_transcript',
    'primer_bind',
    'protein_bind',
    'regulatory',
    'repeat_region',
    'rep_origin',
    'rRNA',
    'S_region',
    'sig_peptide',
    'stem_loop',
    'STS',
    'tmRNA',
    'transit_peptide',
    'tRNA',
    'unsure',
    'V_region',
    'V_segment',
    'variation',
    "3'UTR",
    "5'UTR"
]

GENBANK_QUALIFIER_FEATURE_KEYS = {
    'allele': GENE_RELATED_FEATURES,
    'altitude': ['source'],
    'anticodon': ['tRNA'],
    'artificial_location': ['CDS', 'mRNA'],
    'bio_material': ['source'],
    'bound_moiety': ['misc_binding', 'protein_bind'],
    'cell_line': ['source'],
    'cell_type': ['source'],
    'chromosome': ['chromosome'],
    'citation': GENBANK_FEATURE_KEYS,
    'clone': ['source'],
    'clone_lib': ['source'],
    'codon_start': ['CDS'],
    'collected_by': ['source'],
    'collection_date': ['source'],
    'compare': ['misc_difference', 'unsure', 'old_sequence', 'variation'],
    'country': ['source'],
    'cultivar': ['source'],
    'culture_collection': ['source'],
    'db_xref': GENBANK_FEATURE_KEYS,
    'dev_stage': ['source'],
    'direction': ['oriT', 'rep_origin'],
    'EC_number': ['CDS'],
    'ecotype': ['source'],
    'environmental_sample': ['source'],
    'estimated_length': ['gap', 'assembly_gap'],
    'exception': ['CDS'],
    'experiment': GENBANK_FEATURE_KEYS,
    'focus': ['source'],
    'frequency': ['modified_base', 'variation'],
    'function': [
        'CDS',
        'exon',
        'gene',
        'iDNA',
        'intron',
        'LTR',
        'mat_peptide',
        'misc_binding',
        'misc_feature',
        'misc_RNA',
        'misc_structure',
        'mobile_element',
        'mRNA',
        'ncRNA',
        'operon',
        'precursor_RNA',
        'prim_transcript',
        'protein_bind',
        'regulatory',
        'repeat_region',
        'rRNA',
        'sig_peptide',
        'stem_loop',
        'tmRNA',
        'transit_peptide',
        'tRNA',
        "3'UTR",
        "5'UTR"
    ],
    'gap_type': ['assembly_gap'],
    'gene': GENE_RELATED_FEATURES,
    'gene_synonym': GENE_RELATED_FEATURES,
    'germline': ['source'],
    'haplogroup': ['source'],
    'haplotype': ['source'],
    'host': ['source'],
    'identified_by': ['source'],
    'inference': GENBANK_FEATURE_KEYS,
    'isolate': ['source'],
    'lab_host': ['source'],
    'lat_lon': ['source'],
    'linkage_evidence': ['assembly_gap'],
    'locus_tag': GENE_RELATED_FEATURES,
    'macronuclear': ['source'],
    'map': GENBANK_FEATURE_KEYS,
    'mating_type': ['source'],
    'mobile_element_type': ['mobile_element'],
    # TODO mobile element feature key
    #     Value format    "<mobile_element_type>[:<mobile_element_name>]" where
    #                 mobile_element_type is one of the following:
    #                 "transposon", "retrotransposon", "integron",
    #                 "insertion sequence", "non-LTR retrotransposon",
    #                 "SINE", "MITE", "LINE", "other".
    # Example         /mobile_element_type="transposon:Tnp9"
    #                     /mobile_element_type="insertion sequence:IS5I"
    'mod_base': ['modified_base'],
    'mol_type': ['source'],
    'ncRNA_class': ['ncRNA'],
    'note': GENBANK_FEATURE_KEYS,
    'number': ['CDS', 'exon', 'iDNA', 'intron', 'misc_feature'],
    'old_locus_tag': GENE_RELATED_FEATURES,
    'operon': [
        'CDS',
        'gene',
        'misc_RNA',
        'mRNA',
        'ncRNA',
        'operon',
        'precursor_RNA',
        'prim_transcript',
        'protein_bind',
        'regulatory',
        'rRNA',
        'stem_loop',
        'tRNA'],
    'organelle': ['source'],
    'organism': ['source'],
    'PCR_conditions': ['primer_bind'],
    'PCR_primers': ['primer_bind'],
    'phenotype': ['variation', 'regulatory'],
    'plasmid': ['source'],
    'pop_variant': ['source'],
    'product': ['CDS', 'mRNA', 'mat_peptide'],
    'protein_id': ['CDS'],
    'proviral': ['source'],
    'pseudo': GENE_RELATED_FEATURES,
    'pseudogene': GENE_RELATED_FEATURES,
    'rearranged': ['source'],
    'regulatory_class': ['regulatory'],
    'replace': ['misc_difference', 'old_sequence', 'unsure', 'variation'],
    'ribosomal_slippage': ['CDS', 'mRNA'],

    # NOTE: qualifiers /rpt_unit_range and /rpt_unit_seq replaced qualifier /rpt_unit in December 2005

    # TODO repeat_region!
    'rpt_family': ['repeat_region'],
    'rpt_type': ['repeat_region'],
    'rpt_unit_range': ['oriT', 'repeat_region'],
    'rpt_unit_seq': ['oriT', 'repeat_region'],

    'satellite': ['repeat_region'],
    'segment': ['source'],
    'serotype': ['source'],
    'serovar': ['source'],
    'sex': ['source'],
    'specimen_voucher': ['source'],
    'standard_name': GENBANK_FEATURE_KEYS,
    'strain': ['source'],
    'sub_clone': ['source'],
    'sub_species': ['source'],
    'sub_strain': ['source'],
    'tag_peptide': ['tmRNA'],
    'tissue_lib': ['source'],
    'tissue_type': ['source'],
    'transgenic': ['source'],
    'translation': ['CDS'],
    'transl_except': ['CDS'],
    'transl_table': ['CDS'],
    'trans_splicing': [
        'gene',
        'exon',
        'CDS',
        'mRNA',
        'precursor_RNA',
        'tRNA',
        'ncRNA',
        'misc_RNA',
        'intron'
    ],
    'type_material': ['source'],
    'variety': ['source']
}

GENBANK_QUALIFIERS = tuple(GENBANK_QUALIFIER_FEATURE_KEYS.keys())


def remove_unrecognized_qualifiers(before, after):
    for name, value in before.items():
        if name not in GENBANK_QUALIFIERS:
            del after[name]
    return after


def remove_qualifiers_inappropriate_for_feature(before, after, genbank_feature_key):
    for name, value in before.items():
        if name not in GENBANK_QUALIFIERS or genbank_feature_key not in GENBANK_QUALIFIER_FEATURE_KEYS[name]:
            del after[name]
    return after


DEFAULT_QUALIFIER_TRANSFORMERS = (
    format_single_qualifiers,
    format_integer_qualifiers,
    format_flag_qualifiers,
    remove_protein_id_and_add_to_xrefs,
    rename_label_to_note,
    remove_unrecognized_qualifiers,
)
