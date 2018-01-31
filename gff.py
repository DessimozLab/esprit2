import collections
import tempfile

import gffutils
import os
from gffutils.gffwriter import GFFWriter
from gffutils.attributes import dict_class as gff_attr_class
import logging

logger = logging.getLogger(__name__)


class MergeError(Exception):
    pass


def unique(seq):
    """Return the elements of a list uniquely while preserving the order
    :param list seq: a list of hashable elements
    :returns: new list with first occurence of elements of seq"""
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]

class GFFOperations(object):
    def __init__(self, file, in_memory=False, db_path=None):
        self.remove_db = False
        if in_memory:
            self.dbfn = ":memory:"
        elif db_path is not None:
            self.dbfn = db_path
        else:
            with tempfile.NamedTemporaryFile() as t:
                self.dbfn = t.name
            self.remove_db = True
        self.db = gffutils.create_db(file, self.dbfn, merge_strategy="error")

    def __del__(self):
        if self.remove_db:
            try:
                os.remove(self.dbfn)
            except IOError:
                pass

    def is_valid(self):
        # for now, once file could be loaded, we can assume it is valid as well.
        return True

    def combine_attributes(self, attrs):
        comb_attrs = collections.defaultdict(list)
        for attr in attrs:
            for k, v in attr.items():
                comb_attrs[k].extend(v)
        for k, v in comb_attrs.items():
            comb_attrs[k] = unique(v)
        if 'ID' in comb_attrs:
            comb_attrs['ID'] = "-".join(comb_attrs['ID'])
        return gff_attr_class(comb_attrs)

    def merge_genes(self, split_genes, allow_overlapping_merge=True, **kwargs):
        try:
            gene_features = [self.db[z] for z in split_genes]
        except gffutils.FeatureNotFoundError as e:
            raise KeyError('gene "{}" is not found in gff file'.format(e))
        merged_gene_feature = self.merge_features(gene_features, **kwargs)
        within_genes = list(self.db.region(merged_gene_feature, featuretype='gene', completely_within=True))
        if len(within_genes) > len(split_genes):
            cnt_overlaps = len(within_genes) - len(split_genes)
            if allow_overlapping_merge:
                logger.warning('merging genes that overlap with {} other genes'.format(cnt_overlaps))
            else:
                raise MergeError('merging would results in {} overlapped genes'.format(cnt_overlaps))

        self.db.update((merged_gene_feature,))
        for old_parent in gene_features:
            for child in self.db.children(old_parent, level=1):
                self.db.add_relation(merged_gene_feature, child, 1)
            self.db.delete(old_parent)

    def merge_features(self, features, ignore_strand=False, source=None, **kwargs):
        features = list(features)

        if len(features) == 0:
            return None

        # Sanity check to make sure all features are from the same chromosome.
        chroms = [i.chrom for i in features]
        if len(set(chroms)) > 1:
            raise NotImplementedError('Merging multiple chromosomes not '
                                      'implemented')

        # check all features are on the same strand
        strands = [i.strand for i in features]
        if len(set(strands)) > 1:
            if ignore_strand:
                strand = '.'
                logger.warning('merging features on opposite strands: '+str(strands))
            else:
                raise MergeError('not all features are on the same strand')
        else:
            strand = strands[0]

        merged_start = min(f.start for f in features)
        merged_stop = max(f.stop for f in features)
        attrs = self.combine_attributes([f.attributes for f in features])

        feature = features[0]
        merged_feature = dict(
            seqid=feature.chrom,
            source='.' if source is None else str(source),
            featuretype=feature.featuretype,
            start=merged_start,
            end=merged_stop,
            score=feature.score,
            strand=strand,
            frame=feature.frame,
            attributes=attrs)
        return self.db._feature_returner(**merged_feature)

    def write(self, out):
        writer = GFFWriter(out)
        writer.write_recs(self.db.all_features())
