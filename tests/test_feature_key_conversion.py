import unittest

from gbgb import genbank_feature_key


class FeatureKeyConversionTestCase(unittest.TestCase):

    def test_genbank_feature_type(self):

        # simple feature keys:
        self.assertEqual('gene', genbank_feature_key('gene'))
        self.assertEqual('source', genbank_feature_key('databank_entry'))
        self.assertEqual('mat_peptide', genbank_feature_key('mature_protein_region'))

        # from feature key with qualifiers:
        self.assertEqual('regulatory', genbank_feature_key('TATA_box'))
        self.assertEqual('regulatory', genbank_feature_key('minus_10_signal'))

        self.assertEqual('primer_bind', genbank_feature_key('primer_binding_site'))

        self.assertEqual('regulatory', genbank_feature_key('GC_rich_promoter_region'))

        self.assertEqual('ncRNA', genbank_feature_key('miRNA_gene'))
        self.assertEqual('ncRNA', genbank_feature_key('antisense_RNA'))

        self.assertEqual('mobile_element', genbank_feature_key('transposable_element'))
        self.assertEqual('repeat_region', genbank_feature_key('tandem_repeat'))

        # TODO tes
        self.assertEqual('pseudogene', genbank_feature_key('allelically_excluded_gene'))

        # unknown terms:
        self.assertEqual('misc_feature', genbank_feature_key('thingamajig'))

        # TODO SO terms inheriting from sequence_feature

        # GenBank feature keys:
        self.assertEqual('source', genbank_feature_key('source'))
        self.assertEqual('primer_bind', genbank_feature_key('primer_bind'))
        self.assertEqual('misc_feature', genbank_feature_key('primer_bind', idempotent=False))

        # TODO exhaustive testing
