#! /usr/bin/env python
"""
Classify a directory of sourmash signatures for individual genomes by
taxonomic rank, using a kraken-style least-common-ancestor analysis.

TODO: allow other LCA DBs. output in metacoder form.
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

# Bacteria, cellular organisms, Archaea
DEPRECATED_LCA_to_remove = set({ 2, 131567, 2157 })


def traverse_find_sigs(dirnames):
    """
    Find all the filenames ending with .sig under given directories.
    """
    for dirname in dirnames:
        for root, dirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.sig'):
                    fullname = os.path.join(root, name)
                    yield fullname


def load_all_signatures(dirname, ksize):
    """
    Load all signatures under given dirname with given ksize, return
    dictionary d[name] -> signature.

    Because this can be slow for many hundreds of signatures, cache dict
    using pickle.
    """
    sigd = {}
    filename_d = {}

    pickle_cache = dirname.rstrip('/') + '.pickle'
    if os.path.exists(pickle_cache):
        print('loading from cache:', pickle_cache)
        with open(pickle_cache, 'rb') as fp:
            filename_d, sigd = load(fp)

    loaded_new = False

    n = 0
    for filename in traverse_find_sigs([dirname]):
        if filename not in filename_d:
            loaded_new = True
            sig = sourmash_lib.signature.load_one_signature(filename,
                                                            select_ksize=ksize)
            filename_d[filename] = 1
            sigd[sig.name()] = sig

        n += 1

    if loaded_new:
        print('saving to cache:', pickle_cache)
        with open(pickle_cache, 'wb') as fp:
            dump((filename_d, sigd), fp)

    return sigd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('dir')
    p.add_argument('-o', '--output-csv')
    p.add_argument('--threshold', type=int, default=THRESHOLD,
                   help="minimum number of times a taxid must be present to count")
   
    args = p.parse_args()

    if args.output_csv:
        output_filename = args.output_csv
    else:
        output_filename = os.path.basename(args.dir) + '.taxonomy.csv'

    outfp = open(output_filename, 'wt')
    outw = csv.writer(outfp)
    outw.writerow(['name', 'taxid', 'status', 'rank_info', 'lineage'])

    # load the LCA databases from the JSON file(s)
    lca_db_list = []
    for lca_filename in LCA_DBs:
        print('loading LCA database from {}'.format(lca_filename))
        lca_db = lca_json.LCA_Database(lca_filename)
        taxfoo, hashval_to_lca, _ = lca_db.get_database(args.ksize, SCALED)
        lca_db_list.append((taxfoo, hashval_to_lca))
    
    print('loading all signatures in directory:', args.dir)
    sigdict = load_all_signatures(args.dir, args.ksize)
    print('...loaded {} signatures at k={}'.format(len(sigdict), args.ksize))

    ###

    # track number of disagreements at various rank levels
    disagree_at = collections.defaultdict(int)

    # for each minhash signature in the directory,
    n_in_lca = 0
    for name, sig in sigdict.items():

        # for each k-mer in each minhash signature, collect assigned taxids
        # across all databases (& count).
        taxid_set = collections.defaultdict(int)
        
        for hashval in sig.minhash.get_mins():

            # if a k-mer is present in multiple DBs, pull the
            # least-common-ancestor taxonomic node across all of the
            # DBs.
            
            this_hashval_taxids = set()
            for (_, hashval_to_lca) in lca_db_list:
                hashval_lca = hashval_to_lca.get(hashval)
                if hashval_lca is not None and hashval_lca != 1:
                    this_hashval_taxids.add(hashval_lca)

            this_hashval_lca = taxfoo.find_lca(this_hashval_taxids)
            taxid_set[this_hashval_lca] += 1

        # filter on given threshold - only taxids that show up in this
        # signature more than THRESHOLD.
        abundant_taxids = set([k for (k, cnt) in taxid_set.items() \
                               if cnt >= args.threshold])

        # remove root (taxid == 1) if it's in there:
        if 1 in abundant_taxids:
            abundant_taxids.remove(1)

        # ok - out of the loop, got our LCAs, ...are there any left?
        if abundant_taxids:
            # increment number that are classifiable at *some* rank.
            n_in_lca += 1

            disagree_rank, disagree_taxids = \
              taxfoo.get_lineage_first_disagreement(abundant_taxids,
                                                    want_taxonomy)

            # we found a disagreement - report the rank *at* the disagreement,
            # the lineage *above* the disagreement.
            if disagree_rank:
                # global record of disagreements
                disagree_at[disagree_rank] += 1

                list_at_rank = [taxfoo.get_taxid_name(r) for r in disagree_taxids]
                list_at_rank = ", ".join(list_at_rank)

                print('{} has multiple LCA at {}: \'{}\''.format(name,
                                                                 disagree_rank,
                                                                 list_at_rank))

                # set output
                status = 'disagree'
                status_rank = disagree_rank
                taxid = taxfoo.find_lca(abundant_taxids)
            else:
                # found unambiguous! yay.
                status = 'found'

                taxid = taxfoo.get_lowest_lineage(abundant_taxids,
                                                  want_taxonomy)
                status_rank = taxfoo.get_taxid_rank(taxid)
                status = 'found'
        else:
            # nothing found. boo.
            status = 'nomatch'
            status_rank = ''
            taxid = 0

        if taxid != 0:
            lineage_found = taxfoo.get_lineage(taxid,
                                               want_taxonomy=want_taxonomy)
            lineage_found = ";".join(lineage_found)
        else:
            lineage_found = ""
        
        outfp.write('{},{},{},{},{}\n'.format(name, taxid, status, status_rank, lineage_found))

    print('')
    print('classified sourmash signatures in directory: \'{}\''.format(args.dir))
    print('LCA databases used: {}'.format(', '.join(LCA_DBs)))
    print('')
              
    print('total signatures found: {}'.format(len(sigdict)))
    print('no classification information: {}'.format(len(sigdict) - n_in_lca))
    print('')
    print('could classify {}'.format(n_in_lca))
    print('of those, {} disagree at some rank.'.format(sum(disagree_at.values())))

    print('disagreements by rank:')
    for rank in want_taxonomy:
        if disagree_at.get(rank):
            print('\t{}: {}'.format(rank, disagree_at.get(rank, 0)))

    print('')
    print('classification output as CSV, here: {}'.format(output_filename))


if __name__ == '__main__':
    main()
