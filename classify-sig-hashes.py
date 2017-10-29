#! /usr/bin/env python
"""

"""
import os
import argparse
import collections
import sys
import csv
from pickle import load, dump

import sourmash_lib.signature

import lca_json                      # from github.com/ctb/2017-sourmash-lca

LCA_DBs = ['db/genbank.lca.json']
SCALED=10000
THRESHOLD=5                               # how many counts of a taxid at min

want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--lca', nargs='+', default=LCA_DBs)
    p.add_argument('sig')
    p.add_argument('-o', '--output-csv')
    p.add_argument('--threshold', type=int, default=THRESHOLD,
                   help="minimum number of times a taxid must be present to count")
    p.add_argument('-X', '--output-unassigned', action='store_true')

    args = p.parse_args()

    if args.output_csv:
        output_filename = args.output_csv
    else:
        output_filename = os.path.basename(args.sig) + '.taxonomy.csv'

    # load the LCA databases from the JSON file(s)
    lca_db_list = []
    for lca_filename in args.lca:
        print('loading LCA database from {}'.format(lca_filename))
        lca_db = lca_json.LCA_Database(lca_filename)
        taxfoo, hashval_to_lca, _ = lca_db.get_database(args.ksize, SCALED)
        lca_db_list.append((taxfoo, hashval_to_lca))

    # load signature
    sig = sourmash_lib.signature.load_one_signature(args.sig,
                                                    ksize=args.ksize)
    hashes = sig.minhash.get_mins()

    # open output file.
    outfp = open(output_filename, 'wt')
    outw = csv.writer(outfp)
    outw.writerow(['hash', 'taxid', 'status', 'rank_info', 'lineage'])

    ###

    # track number of classifications at various rank levels
    classified_at = collections.defaultdict(int)

    # also track unassigned
    unassigned = set()

    # for each hash in the minhash signature, get its LCA.
    n_in_lca = 0
    taxid_counts = collections.defaultdict(int)
    lca_to_hashvals = collections.defaultdict(set)
    for hashval in hashes:
        # if a k-mer is present in multiple DBs, pull the
        # least-common-ancestor taxonomic node across all of the
        # DBs.

        this_hashval_taxids = set()
        for (_, hashval_to_lca) in lca_db_list:
            hashval_lca = hashval_to_lca.get(hashval)
            if hashval_lca is not None and hashval_lca != 1:
                this_hashval_taxids.add(hashval_lca)

        if this_hashval_taxids:
            this_hashval_lca = taxfoo.find_lca(this_hashval_taxids)
            if this_hashval_lca != None:
                taxid_counts[this_hashval_lca] += 1
                lca_to_hashvals[this_hashval_lca].add(hashval)
        else:
            unassigned.add(hashval)

    # filter on given threshold - only taxids that show up in this
    # signature more than THRESHOLD.
    abundant_taxids = set([k for (k, cnt) in taxid_counts.items() \
                           if cnt >= args.threshold])

    # remove root (taxid == 1) if it's in there:
    if 1 in abundant_taxids:
        abundant_taxids.remove(1)

    # now, output hashval classifications.
    n_classified = 0
    for lca_taxid in abundant_taxids:
        for hashval in lca_to_hashvals[lca_taxid]:
            status = 'match'
            status_rank = taxfoo.get_taxid_rank(lca_taxid)
            lineage = taxfoo.get_lineage(lca_taxid,
                                         want_taxonomy=want_taxonomy)
            lineage = ";".join(lineage)

            classified_at[status_rank] += 1
            n_classified += 1

            outw.writerow([str(hashval), str(lca_taxid),
                           status, status_rank, lineage])

    # output unassigned?
    if args.output_unassigned:
        for hashval in unassigned:
            status = 'nomatch'
            status_rank = 'unknown'
            lineage = 'UNKNOWN'

            outw.writerow([str(hashval), str(lca_taxid),
                           status, status_rank, lineage])

    print('')
    print('classified sourmash signature \'{}\''.format(args.sig))
    print('LCA databases used: {}'.format(', '.join(args.lca)))
    print('')

    print('total hash values: {}'.format(len(hashes)))
    print('num classified: {}'.format(n_classified))

    n_rare_taxids = sum([cnt for (k, cnt) in taxid_counts.items() \
                         if cnt < args.threshold ])
    print('n rare taxids not used: {}'.format(n_rare_taxids))
    print('unclassified: {}'.format(len(unassigned)))

    print('')
    print('number classified unambiguously, by lowest classification rank:')
    for rank in want_taxonomy:
        if classified_at.get(rank):
            print('\t{}: {}'.format(rank, classified_at.get(rank, 0)))

    print('')
    print('classification output as CSV, here: {}'.format(output_filename))


if __name__ == '__main__':
    main()
