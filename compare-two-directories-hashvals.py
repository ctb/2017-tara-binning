#! /usr/bin/env python
"""
Compare two directories of signatures, for fun and profit.

This version just looks at distinct hash values.
"""
import argparse
import os
import sourmash_lib
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import (SigLeaf, search_minhashes,
                                SearchMinHashesFindBest)
from pickle import load, dump


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
    p.add_argument('dir1')
    p.add_argument('dir2')
    p.add_argument('-k', '--ksize', type=int, default=31)
    args = p.parse_args()

    print('loading all signatures:', args.dir1)
    sigdict1 = load_all_signatures(args.dir1, args.ksize)
    print('...loaded {} signatures at k={}'.format(len(sigdict1), args.ksize))

    print('loading all signatures:', args.dir2)
    sigdict2 = load_all_signatures(args.dir2, args.ksize)
    print('...loaded {} signatures at k={}'.format(len(sigdict2), args.ksize))

    # collect all hash values
    hashvals_1 = set()
    for v in sigdict1.values():
        hashvals_1.update(v.minhash.get_mins())

    hashvals_2 = set()
    for v in sigdict2.values():
        hashvals_2.update(v.minhash.get_mins())

    total = len(hashvals_1.union(hashvals_2))
    print('total hashvals:', total)
    print('both: {:.1f}%'.format(len(hashvals_1.intersection(hashvals_2)) / total * 100))
    print('unique to {}: {:.1f}%'.format(args.dir1, len(hashvals_1 - hashvals_2) / total * 100))
    print('unique to {}: {:.1f}%'.format(args.dir2, len(hashvals_2 - hashvals_1) / total * 100))

                    
if __name__ == '__main__':
    main()
