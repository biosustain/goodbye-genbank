import unittest

from gbgb.maps.qualifiers import remove_protein_id_and_add_to_xrefs, rename_label_to_note, format_flag_qualifiers, \
    remove_ape_a_plasmid_editor_qualifiers, format_single_qualifiers, format_integer_qualifiers, \
    remove_qualifiers_inappropriate_for_feature


class TranslationTestCase(unittest.TestCase):

    def test_remove_protein_id_and_add_to_xrefs(self):
        self.assertEqual({
            'db_xref': {'GenBank:XYZ12345.1'}
        }, remove_protein_id_and_add_to_xrefs({
            'protein_id': ['XYZ12345.1']
        }, {
            'protein_id': ['XYZ12345.1']
        }))

        self.assertEqual({
            'db_xref': {'GI:123456', 'GenBank:XYZ12345.1'}
        }, remove_protein_id_and_add_to_xrefs({
            'protein_id': ['XYZ12345.1']
        }, {
            'protein_id': ['XYZ12345.1'],
            'db_xref': ['GI:123456']
        }))

        self.assertEqual({
            'db_xref': {'GenBank:XYZ12345.1'}
        }, remove_protein_id_and_add_to_xrefs({
            'protein_id': ['XYZ12345.1']
        }, {
            'protein_id': ['XYZ12345.1'],
            'db_xref': ['GenBank:XYZ12345.1']
        }))

        # ignore malformed protein_ids
        self.assertEqual({
        }, remove_protein_id_and_add_to_xrefs({
            'protein_id': ['MyLittleProtein']
        }, {
            'protein_id': ['MyLittleProtein']
        }))

        # noop:
        self.assertEqual({
            'db_xref': {'GenBank:XYZ12345.1'}
        }, remove_protein_id_and_add_to_xrefs({
            'db_xref': {'GenBank:XYZ12345.1'}
        }, {
            'db_xref': {'GenBank:XYZ12345.1'}
        }))

    def test_rename_label_to_note(self):
        self.assertEqual({
            'note': 'This is a feature.'
        }, rename_label_to_note({
            'label': ['This is a feature.']
        }, {}))

        # noop:
        self.assertEqual({
            'note': 'This is a feature.'
        }, rename_label_to_note({
            'note': 'This is a feature.'
        }, {
            'note': 'This is a feature.'
        }))

    def test_format_single_qualifiers(self):
        self.assertEqual({
            'locus_tag': 'b2892',
            'gene': 'recJ',
            'gene_synonym': ['ECK2887', 'JW2860']
        }, format_single_qualifiers({
            'locus_tag': ['b2892'],
            'gene': ['recJ'],
            'gene_synonym': ['ECK2887', 'JW2860']
        }, {
            'locus_tag': ['b2892'],
            'gene': ['recJ'],
            'gene_synonym': ['ECK2887', 'JW2860']
        }))

        # remove duplicates because they cannot be interpreted properly:
        # XXX should maybe remove instead of set to None
        self.assertEqual({
            'locus_tag': None,
        }, format_single_qualifiers({
            'locus_tag': ['b1', 'b2']
        }, {
            'locus_tag': ['b1', 'b2']
        }))

    def test_format_integer_qualifiers(self):
        self.assertEqual({
            'codon_start': 1
        }, format_integer_qualifiers({
            'codon_start': ['1']
        }, {
            'codon_start': ['1']
        }))

        self.assertEqual({}, format_integer_qualifiers({
            'codon_start': ['second-codon']
        }, {
            'codon_start': ['second-codon']
        }))

    def test_format_flag_qualifiers(self):
        self.assertEqual({
            'pseudo': True
        }, format_flag_qualifiers({
            'pseudo': ['']
        }, {
            'pseudo': ['']
        }))

    def test_remove_ape_a_plasmid_editor_qualifiers(self):
        self.assertEqual({
            'note': 'This feature, of course, is blue.',
        }, remove_ape_a_plasmid_editor_qualifiers({
            'note': 'This feature, of course, is blue.',
            'ApEinfo_some_color': ['#0000ff']
        }, {
            'note': 'This feature, of course, is blue.',
            'ApEinfo_some_color': ['#0000ff']
        }))

    def test_remove_qualifiers_inappropriate_for_feature(self):
        self.assertEqual({
            'note': 'This feature, of course, is blue.',
            'focus': True,
        }, remove_qualifiers_inappropriate_for_feature({
            'note': 'This feature, of course, is blue.',
            'focus': True,
            'ApEinfo_some_color': ['#0000ff']
        }, {
            'note': 'This feature, of course, is blue.',
            'focus': True,
            'ApEinfo_some_color': ['#0000ff']
        }, genbank_feature_key='source'))

        self.assertEqual({
            'gene': 'recJ',
            'gene_synonym': ['ECK2887', 'JW2860']
        }, remove_qualifiers_inappropriate_for_feature({
            'gene': 'recJ',
            'gene_synonym': ['ECK2887', 'JW2860'],
            'country': 'Denmark'
        }, {
            'gene': 'recJ',
            'gene_synonym': ['ECK2887', 'JW2860'],
            'country': 'Denmark'
        }, genbank_feature_key='gene'))

        self.assertEqual({
            'gene': 'genA',
            'regulatory_class': 'attenuator'
        }, remove_qualifiers_inappropriate_for_feature({
            'gene': 'genA',
            'label': 'This attenuator, attenuates genA.',
            'regulatory_class': 'attenuator'
        }, {
            'gene': 'genA',
            'label': 'This attenuator, attenuates genA.',
            'regulatory_class': 'attenuator'
        }, genbank_feature_key='regulatory'))

        self.assertEqual({
            'experiment': 'DESCRIPTION: sequenced by a monkey with a typewriter'
        }, remove_qualifiers_inappropriate_for_feature({
            'locus_key': 'b0001',
            'experiment': 'DESCRIPTION: sequenced by a monkey with a typewriter',
            'invalid': True
        }, {
            'locus_key': 'b0001',
            'experiment': 'DESCRIPTION: sequenced by a monkey with a typewriter',
            'invalid': True
        }, genbank_feature_key='misc_feature'))

        # TODO exhaustive testing

