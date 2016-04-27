import os
import unittest


from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from six import StringIO

from gbgb import convert_feature, unconvert_feature

FILES_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'files/')

class FeatureConversionTestCase(unittest.TestCase):

    def assertFeatureLocationEqual(self, l1, l2):
        self.assertEqual(l1.__class__, l2.__class__)
        if isinstance(l1, CompoundLocation):
            self.assertEqual(l1.operator, l2.operator)
            self.assertEqual(len(l1.parts), len(l2.parts))
            for p1, p2 in zip(l1.parts, l2.parts):
                self.assertFeatureLocationEqual(p1, p2)
        else:
            self.assertEqual(l1.start, l2.start)
            self.assertEqual(l1.end, l2.end)
            self.assertEqual(l1.strand, l2.strand)
            self.assertEqual(l1.ref, l1.ref)
            self.assertEqual(l1.ref_db, l1.ref_db)

    def assertFeatureEqual(self, f1, f2):
        self.assertEqual(f1.__class__, f2.__class__)
        self.assertEqual(f1.type, f2.type)
        self.assertFeatureLocationEqual(f1.location, f2.location)
        self.assertEqual(f1.qualifiers, f2.qualifiers)

    def test_genbank_record_conversion(self):
        with open(os.path.join(FILES_PATH, 'sample.gb')) as f:
            first_record = next(SeqIO.parse(f, 'genbank'))

        self.assertEqual(4, len(first_record.features))

        EXPECTED_FEATURES = (
            SeqFeature(type='CDS',
                       location=FeatureLocation(24, 43, strand=-1),
                       qualifiers={
                           'note': 'CDS',
                           'locus_tag': 'b0001'
                       }),
            SeqFeature(type='sequence_uncertainty',
                       location=FeatureLocation(102, 106, strand=1),
                       qualifiers={
                           'note': 'What is this?'
                       }),
            SeqFeature(type='gene',
                       location=FeatureLocation(19, 90, strand=-1),
                       qualifiers={
                           'note': 'genA',
                           'gene': 'genA',
                           'gene_synonym': ['A', 'alpha'],
                           'locus_tag': 'b0001'
                       }),
            SeqFeature(type='sequence_feature',
                       location=FeatureLocation(106, 110, strand=1),
                       qualifiers={
                           'note': 'This is something.'
                       })
        )

        for f1, f2 in zip(EXPECTED_FEATURES, first_record.features):
            self.assertFeatureEqual(f1, convert_feature(f2))

    def test_genbank_record_fixing(self):
        with open(os.path.join(FILES_PATH, 'sample.gb')) as f:
            first_record = next(SeqIO.parse(f, 'genbank'))

        first_record.features = [unconvert_feature(convert_feature(f)) for f in first_record.features]

        output = StringIO()
        SeqIO.write(first_record, output, "genbank")
        with open(os.path.join(FILES_PATH, 'sample.fixed.gb')) as f:
            self.assertEqual(output.getvalue(), f.read())
        output.close()
