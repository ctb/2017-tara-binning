#! /usr/bin/env python
"""
Compare two directories of signatures, for fun and profit.
"""
import argparse
import os
import sourmash_lib
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import (SigLeaf, search_minhashes,
                                SearchMinHashesFindBest)


def traverse_find_sigs(dirnames):
    for dirname in dirnames:
        for root, dirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.sig'):
                    fullname = os.path.join(root, name)
                    yield fullname


def make_all_matches(sigdict, tree, threshold):
    match_d = {}
    search_fn = lambda: SearchMinHashesFindBest().search

    for query in sigdict.values():
        matching_sig = None
        for leaf in tree.find(search_fn(), query, threshold):

            # deal with bug? in this search_fn; thresholds not always met.
            similarity = leaf.data.similarity(query)
            if similarity >= threshold:
                matching_sig = leaf.data
                print('match:', query.name(), matching_sig.name(), similarity)
                match_d[query.name()] = (matching_sig.name(), similarity)

        if not matching_sig:
            print('no match found:', query.name())

    return match_d


def load_all_signatures(dirname, ksize):
    sigd = {}

    n = 0
    for filename in traverse_find_sigs([dirname]):
        sig = sourmash_lib.signature.load_one_signature(filename,
                                                        select_ksize=ksize)
        sigd[sig.name()] = sig

        n += 1
        if n > 20:
            break

    return sigd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('dir1')
    p.add_argument('sbt1')
    p.add_argument('dir2')
    p.add_argument('sbt2')
    p.add_argument('-k', '--ksize', type=int, default=31)
    args = p.parse_args()

    print('loading all signatures:', args.dir1)
    sigdict1 = load_all_signatures(args.dir1, args.ksize)
    tree1 = SBT.load(args.sbt1, leaf_loader=SigLeaf.load)
    print('...loaded {} signatures at k={}'.format(len(sigdict1), args.ksize))

    print('loading all signatures:', args.dir2)
    sigdict2 = load_all_signatures(args.dir2, args.ksize)
    tree2 = SBT.load(args.sbt2, leaf_loader=SigLeaf.load)
    print('...loaded {} signatures at k={}'.format(len(sigdict2), args.ksize))

    # first, find all matches in 2 for 1, and 1 for 2
    THRESHOLD=0.05
    matches_1_in_2 = make_all_matches(sigdict1, tree2, THRESHOLD)
    matches_2_in_1 = make_all_matches(sigdict2, tree1, THRESHOLD)

    # now, do containment, too.
    contained_1_in_2 = {}
    for query_name, (match_name, similarity) in matches_1_in_2.items():
        if query_name not in sigdict1 or match_name not in sigdict2:
            # during testing with subsets, this may be true
            continue

        query = sigdict1[query_name]
        match = sigdict2[match_name]

        print('got one!')
        

                    
if __name__ == '__main__':
    main()
