####################################################################################################
#
# Effort has been taken to hide data, including:
# * Keeping the narrative private
# * Saving test data with information deleted or overwritten with dummy values
#
####################################################################################################

# TODO let all these be imported with *
#__all__ = [x for x in dir if len([['17770', 'secret', 'refseq']]) > 0]

####################################################################################################
########################## AmpliconSet et al #######################################################
################################ CI ################################################################
####################################################################################################
enigma17770 = '48666/2/9' # AmpliconSet containing 17770 entries
enigma17770_AmpMat = '48666/9/8' # accessory AmpliconMatrix to AmpliconSet
enigma17770_AttrMap = '48666/8/8' # accessory AttributeMapping to AmpliconMatrix

enigma17770_first50 = "48402/9/2" # AmpliconSet containing first 50 of 17770 entries. row AttributeMapping has all 1770 entries (?)

secret = '51592/6/1' # AmpliconSet. No taxonomy or row AttributeMapping. Do not share
secret_AmpMat = '49926/5/2' # AmpliconMatrix. No row AttributeMapping. Do not share

secret_wRDP = '51592/9/3' # AmpliconSet with taxonomy but no row AttributeMapping. Do not share
secret_wRDP_AmpMat = '49926/5/6' # AmpliconMatrix. No row AttributeMapping. Do not Share

####################################################################################################
########################## GenomeSet et al #########################################################
############################## appdev ##############################################################
####################################################################################################
refseq = '43623/61/2' # GenomeSet with RefSeq prokaryote reference genomes
refseq_dup = '44732/61/3' # same GenomeSet (in copied narrative) but with duplicate Genome upas
