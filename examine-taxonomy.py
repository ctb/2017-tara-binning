#! /usr/bin/env python
import os
import argparse
import collections
import sys
from pickle import load, dump
import sourmash_lib.signature

sys.path.insert(0, '../sourmash/lca/')
import lca_json                      # from github.com/ctb/2017-sourmash-lca

LCA_DBs = ['db/delmont.lca.json', 'db/tully.lca.json',
           '../sourmash/genbank/genbank.lca.json']
#LCA_DBs = ['db/tully.lca.json']
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
    args = p.parse_args()

    # load all the LCA JSON files
    lca_db_list = []
    for lca_filename in LCA_DBs:
        lca_db = lca_json.LCA_Database(lca_filename)
        taxfoo, hashval_to_lca, _ = lca_db.get_database(args.ksize, SCALED)
        lca_db_list.append((taxfoo, hashval_to_lca))
    
    print('loading all signatures:', args.dir)
    sigdict = load_all_signatures(args.dir, args.ksize)
    print('...loaded {} signatures at k={}'.format(len(sigdict), args.ksize))

    ###

    disagree_at = collections.defaultdict(int)

    n = 0
    for name, sig in sigdict.items():
        taxid_set = collections.defaultdict(int)
        for hashval in sig.minhash.get_mins():

            this_hashval_taxids = set()
            for (_, hashval_to_lca) in lca_db_list:
                hashval_lca = hashval_to_lca.get(hashval)
                if hashval_lca is not None:
                    this_hashval_taxids.add(hashval_lca)

            this_hashval_lca = taxfoo.find_lca(this_hashval_taxids)
            taxid_set[this_hashval_lca] += 1

        abundant_taxids = [k for k in taxid_set if taxid_set[k] >= THRESHOLD]

        if not abundant_taxids:
            continue

        n += 1

        ranks_found = collections.defaultdict(set)
        for taxid in abundant_taxids:
            d = taxfoo.get_lineage_as_dict(taxid)
            for k, v in d.items():
                ranks_found[k].add(v)

        found_disagree = False
        for rank in reversed(want_taxonomy):
            if len(ranks_found[rank]) > 1:
                disagree_at[rank] += 1
                found_disagree = True
                break

        if found_disagree:
            print('{} has multiple LCA at rank \'{}\': {}'.format(name,
                                                                  rank,
                                                                  ", ".join(ranks_found[rank])))

    print('for', args.dir, 'found', len(sigdict), 'signatures;')
    print('out of {} that could be classified, {} disagree at some rank.'.format(n, sum(disagree_at.values())))

    for rank in want_taxonomy:
        if disagree_at.get(rank):
            print('\t{}: {}'.format(rank, disagree_at.get(rank, 0)))


if __name__ == '__main__':
    main()
