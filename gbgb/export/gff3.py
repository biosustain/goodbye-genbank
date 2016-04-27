from gbgb import single


class GFF3Writer(object):

    def __init__(self):
        pass

    def write(self, record, fp, include_sequence=False):
        pass

    def write_feature(self, feature, fp, sequence=None):
        pass

    def _feature_id(self, feature):
        """

        Uses the following steps to find a unique feature ID:

        1. Looks for a unique identifier in `SeqFeature.id`
        2. Looks for a unique identifier in the qualifiers. How this is done depends
           on the feature type.
        3. Generates a unique identifier from the feature type and an incrementing counter.
           e.g. exon00001

        :param feature:
        :return:
        """
        pass

    # TODO handle trans-splicing:
    # http://sequenceontology.org/so_wiki/index.php/Discontinuous_features#Example_of_NCBI_Trans-splicing


    # TODO simplify feature types, but then use Ontology_term with a attribute (quality):
    # e.g. mRNA_with_frameshift becomes mRNA with quality frameshift
    def _feature_type_quality(self, feature):
        try:
            return {
                'mRNA_with_frameshift': ('mRNA', 'frameshift'),
                'gene_with_mRNA_with_frameshift': ('gene', 'frameshift'),
                # ...
            }[feature.type]
        except KeyError:
            return feature.type, None

    def _feature_parent(self, feature, record):
        """

        Attempts to find a parent or parents for a feature.

        Where multiple features share a `/locus_tag=` or `/gene=` qualifier, parents are determined by feature type:

        - A "gene" is a top level feature; all other features are children of it.
        - A "CDS" or "exon" is a child of "mRNA" at the same location, or of the "gene" if no "mRNA" at that position
          exists.

        This corresponds to "part_of" relationships of the Sequence Ontology.

        .. todo::

            Automatically

        :param feature:
        :param record:
        :return:
        """
        pass

    def _feature_display_name(self, feature):

        if 'standard_name' in feature.qualifiers:
            return single(feature.qualifiers, 'standard_name')

        # TODO should include pseudogene etc.
        if feature.type == 'gene' and 'gene' in feature.qualifiers:
            return single(feature.qualifiers, 'gene')

        # TODO if number could name exons/CDDs as gene.number.

        return None

    def _feature_note(self, feature):
        return feature.qualifiers.get('note', None)

    def _CDS_feature_phase(self, feature):
        return feature.qualifiers.get('codon_start', 1) - 1

    # TODO source feature Is_circular. If no source feature exists, generate a "region" feature with the same ID as
    #  sequence ID and Is_circular=true

    # TODO handling of gaps.

