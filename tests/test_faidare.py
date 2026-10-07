'''
FAIDARE: the entry description counts the cluster's genes and species instead of listing them.
'''
import pss_export.entity_classes as ec
from pss_export.faidare.faidare_api import description, plural


def test_description_counts_genes_and_species():
    cluster = ec.Node('WRKY33[fc00166]', labels=['Node', 'FunctionalCluster', 'Plant', 'PlantCoding'],
                      short_name='WRKY33', description='WRKY transcription factor 33', pathway='P1',
                      homologues={'ath': ['AT2G38470'], 'stu': ['S1', 'S2'], 'mdo': []})
    text = description('AT2G38470', cluster, {'protein activation': {'MPK3,6'}})
    assert 'FunctionalCluster WRKY33[fc00166]' in text
    assert 'includes 3 genes across 2 species.' in text
    assert 'S1' not in text and 'protein activation with MPK3,6' in text


def test_plural():
    assert (plural(1, 'gene'), plural(2, 'gene'), plural(1, 'species', 'species')) == ('1 gene', '2 genes', '1 species')
