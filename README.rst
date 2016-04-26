================
Goodbye, GenBank
================

This package converts* `SeqFeature <http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html>`_ sequence
annotations from NCBI GenBank records to a common and simplified format. GenBank feature annotations have a
type and reasonably well defined qualifiers, but non-standard and discontinued feature types and qualifiers are commonly
used, and often the feature type is

This package converts most of these feature types to appropriate `Sequence Ontology <http://www.sequenceontology.org/>`_ terms used by GFF3 and SBOL.
 Non-standard qualifiers are repaired or removed.

*Goodbye, GenBank* is intended for those who wish to clean-up their GenBank files and then transition to a different format.
The philosophy of this project is to salvage what is salvageable and to discard what is not. GenBank feature types are translated
to Sequence Ontology terms; qualifiers are converted into a reduced set that contains only the parts that are not broken.

However, different options are available to configure what is kept and what is thrown away.

Example
-------

::

    >>> feature
    SeqFeature(FeatureLocation(ExactPosition(2931), ExactPosition(2936), strand=1), type='-10_signal')
    >>> feature.qualifiers
    {'ApEinfo_fwdcolor': ['pink'],
     'ApEinfo_graphicformat': ['arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0'],
     'ApEinfo_revcolor': ['pink'],
     'label': ['RNAII Promoter (-10 signal)']}
    >>>
    >>> from gbgb import convert_feature
    >>> feature = convert_feature(feature)
    >>> feature
    SeqFeature(FeatureLocation(ExactPosition(2931), ExactPosition(2936), strand=1), type='minus_10_signal')
    >>> feature.qualifiers
    {'note': 'RNAII Promoter (-10 signal)'}
    >>>
    >>> from gbgb import genbank_feature_key
    >>> genbank_feature_key('minus_10_signal')
    'regulatory'


Design considerations
---------------------

For the most part, *Goodbye, GenBank* attempts to be idempotent, i.e. features and their types/keys and qualifiers can be safely
transformed any number times with the same settings. The apparent mismatch between the conversion to Sequence Ontology feature
terms and valid/fixed GenBank qualifiers is to simplify downstream processing. It is up to the users which qualifiers they wish
to keep, but at least the choices they are given are reasonable.

Contributing
------------

If you have any questions or suggestions or if you have found a unique new specimen of GenBank files that you would like
to convert, please open an issue.


Issues
------

- SO Term: "regulatory" feature type with /regulatory_class="enhancer_blocking_element"

  There is apparently no matching Sequence Ontology term. An enhancer blocking element behaves like an insulator, but
  is not an insulator. It is a transcriptional cis regulatory region, but that description is too broad.

- SO Term: "misc_structure" feature type

  GenBank uses this feature type for secondary and tertiary nucleotide structures. There appears to be
  no matching Sequence Ontology term.

- SO Term: "assembly_gap" feature type

  GenBank has both "gap" and "assembly_gap" feature types, which appear to have slightly different meanings. However,
  SO only has a "gap" term, which refers to assembly gaps.

- GFF3 export

  There is no good GFF3 exporter out there, so why not write one?

  Skeleton code in gbgb.export.gff3

- Reduction of SO terms

  Allow users to specify a set of Sequence Ontology terms (inheriting from "sequence_feature"). Feature types will be
  reduced to the nearest Ontology term. This is to simplify downstream analysis.

- /pseudo qualifier without /pseudogene=""

  There is no matching Sequence Ontology term for this. Several GenBank files contain /pseudo without /pseudogene=""
  to mean pseudogene.

- Mandatory qualifiers

  These should be filled in using a reasonable guess or errors should be thrown when trying to convert a feature without
  its mandatory qualifiers.


Materials
---------

- `GenBank Feature Table Definition <http://www.insdc.org/documents/feature-table>`_

- `GenBank Release Notes (since December 1992) <http://www.ncbi.nlm.nih.gov/genbank/release/>`_

- `NCBI Prokaryotic Genome Annotation Guide <http://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/>`_

- `Sequence Ontology Wiki -- Discontinuous Features <http://sequenceontology.org/so_wiki/index.php/Discontinuous_features>`_.

  On trans-splicing.

- `GFF3 Specification <http://www.sequenceontology.org/resources/gff3.html>`_


