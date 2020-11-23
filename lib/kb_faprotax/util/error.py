class NoTaxonomyException(Exception): pass

class NoWsReferenceException(Exception): pass

class NonZeroReturnException(Exception): pass

class ValidationException(Exception): pass

msg_overwriteAttribute = (
'Overwriting attribute `%s` with source `%s` in row AttributeMapping with name `%s`')

msg_missingTaxonomy = (
"Input AmpliconSet is missing taxonomy for amplicon with id `%s`. (Taxonomic assignments can be obtained by first running kb_RDP_Classifier.)") # using amplicon id for testing

msg_dupGenomes = (
'Duplicate referenced Genomes in GenomeSet')
