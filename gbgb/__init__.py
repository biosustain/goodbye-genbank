from Bio.SeqFeature import SeqFeature

from gbgb.maps.features import GENBANK_REGULATORY_DEFAULT_SO_TERM, GENBANK_REGULATORY_CLASS_SO_TERMS, \
    UNAMBIGUOUS_INVALID_KEY_SO_TERMS, GENBANK_PSEUDOGENE_TYPE_SO_TERMS, GENBANK_NC_RNA_DEFAULT_SO_TERM, \
    GENBANK_PSEUDOGENE_DEFAULT_SO_TERM, DEFAULT_SO_TERM, SO_TERM_GENBANK_FEATURE_KEYS, DEFAULT_GENBANK_FEATURE_KEY, \
    GENBANK_FEATURE_KEYS, GENBANK_FEATURE_KEY_SO_TERMS
from gbgb.maps.qualifiers import single, DEFAULT_QUALIFIER_TRANSFORMERS, remove_qualifiers_inappropriate_for_feature

# https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/VERSION-INFO
# https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/so-xp-simple.obo
__sequence_ontology_version__ = '2015-06-22'

# http://www.ncbi.nlm.nih.gov/genbank/release/
__genbank_release__ = 211
__genbank_release_date__ = '2015-12'


def convert_feature_type(feature):
    """
    Finds a Sequence Ontology term for a GenBank feature.

    This function requires a :class:`SeqFeature` as opposed to just a GenBank feature key, since the type of a GenBank
    feature is not always fully described by its feature key. For example a `regulatory` GenBank feature could have
    a `/regulatory_class="promoter"` qualifier.

    :param SeqFeature feature:
    :return: a Sequence Ontology term for the type of this feature
    """
    type_ = feature.type

    if type_ == 'regulatory':
        regulatory_class = single(feature.qualifiers, 'regulatory_class', on_multiple_ignore=True)

        if regulatory_class is None:
            return GENBANK_REGULATORY_DEFAULT_SO_TERM
        return GENBANK_REGULATORY_CLASS_SO_TERMS.get(regulatory_class, GENBANK_REGULATORY_DEFAULT_SO_TERM)

    elif type_ == 'ncRNA':
        nc_rna_class = single(feature.qualifiers, 'ncRNA_class')

        if nc_rna_class is None:
            return GENBANK_NC_RNA_DEFAULT_SO_TERM
        return GENBANK_REGULATORY_CLASS_SO_TERMS.get(nc_rna_class, GENBANK_NC_RNA_DEFAULT_SO_TERM)

    else:
        if 'pseudo' in feature.qualifiers:
            pass
            # "The qualifier /pseudo should be used to describe non-functional
            #  genes that are not formally described as pseudogenes, e.g. CDS
            #  has no translation due to other reasons than pseudogenisation events.
            #  Other reasons may include sequencing or assembly errors.
            #  In order to annotate pseudogenes the qualifier /pseudogene= must be
            #  used indicating the TYPE which can be taken from the INSDC controlled vocabulary
            #  for pseudogenes."

            # FIXME what to do about these?
            # if 'pseudogene' not in feature.qualifiers:
            #     raise NotImplementedError()

        #  /pseudo and /pseudogene="" used with gene
        if 'pseudogene' in feature.qualifiers:
            pseudogene = single(feature.qualifiers, 'pseudogene', on_multiple_ignore=True)

            if type_ == 'gene':
                return GENBANK_PSEUDOGENE_TYPE_SO_TERMS.get(pseudogene, GENBANK_PSEUDOGENE_DEFAULT_SO_TERM)

        # The /ribosomal_slippage qualifier is used with genes that have a translational frameshift.
        if 'ribosomal_slippage' in feature.qualifiers:
            if type_ == 'gene':
                return 'gene_with_mRNA_with_frameshift'
            elif type_ == 'mRNA':
                return 'mRNA_with_frameshift'
            # TODO CDS with frameshift, ..

        # /trans_splicing is used on features such as CDS, mRNA and other features that are produced as
        # a result of a trans-splicing event.
        if 'trans_splicing' in feature.qualifiers:
            if type_ == 'mRNA':
                return 'trans_spliced_mRNA'
            # TODO trans-spliced CDS, tRNA, ..

        # TODO repeat_region, mobile_element

        try:
            return GENBANK_FEATURE_KEY_SO_TERMS[type_]
        except KeyError:
            pass

        if type_ in SO_TERM_GENBANK_FEATURE_KEYS:
            return type_

        # TODO use a sequence ontology and allow any term that inherits from sequence_feature

        try:
            return UNAMBIGUOUS_INVALID_KEY_SO_TERMS[type_]
        except KeyError:
            return DEFAULT_SO_TERM


def convert_feature(feature, qualifier_transformers=DEFAULT_QUALIFIER_TRANSFORMERS):
    """

    :param SeqFeature feature:
    :param qualifier_transformers: a set of transformation functions used to clean up the qualifiers.
    :return: a :class:`SeqFeature` with valid GenBank qualifiers and a feature type that is a Sequence Ontology term.
    """
    type_ = convert_feature_type(feature)

    before, after = feature.qualifiers, feature.qualifiers

    for transformer in qualifier_transformers:
        before = after = transformer(before, dict(after))

    # finally, remove all qualifiers that do not belong to a certain feature key
    qualifiers = remove_qualifiers_inappropriate_for_feature(before, dict(after), genbank_feature_key(type_))

    return SeqFeature(location=feature.location,
                      location_operator=feature.location_operator,
                      type=type_,
                      id=feature.id,
                      qualifiers=qualifiers,
                      ref=feature.ref,
                      ref_db=feature.ref_db)


def _unconvert_qualifier(qualifier):
    if qualifier is True:
        return []
    if isinstance(qualifier, (list, tuple, set)):
        return list(qualifier)
    return [qualifier]


def unconvert_feature(feature):
    """
    Restores a feature to its sad former shape.

    :param SeqFeature feature:
    :return: a :class:`SeqFeature` with each qualifier being a list and having a GenBank feature key
    """
    return SeqFeature(location=feature.location,
                      location_operator=feature.location_operator,
                      type=genbank_feature_key(feature),
                      id=feature.id,
                      qualifiers={name: _unconvert_qualifier(value) for name, value in feature.qualifiers.items()},
                      ref=feature.ref,
                      ref_db=feature.ref_db)


def genbank_feature_key(sequence_feature_term, idempotent=True):
    """

    :param str sequence_feature_term:
    :param bool idempotent: whether to pass through GenBank feature keys
    :return: a GenBank feature key for this Sequence Ontology term
    """
    if idempotent and sequence_feature_term in GENBANK_FEATURE_KEYS:
        return sequence_feature_term

    if sequence_feature_term in SO_TERM_GENBANK_FEATURE_KEYS:
        return SO_TERM_GENBANK_FEATURE_KEYS[sequence_feature_term]

    # TODO use a sequence ontology to look up the nearest feature for any term that inherits from sequence_feature

    return DEFAULT_GENBANK_FEATURE_KEY
