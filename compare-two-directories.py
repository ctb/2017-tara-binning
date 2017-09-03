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


def make_all_matches(sigdict, tree, threshold):
    """
    Find all the matches between a dictionary of signatures and an
    SBT (search tree), at or above given threshold.

    Return a dictionary of d[signame] -> (match_name, similarity)
    """
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


def containment(matches, query_sigdict, match_sigdict):
    """
    Compute containment in query for its match.

    Returns dict d[query_name] -> containment.
    """
    contained = dict()
    for query_name, (match_name, similarity) in matches.items():
        query = query_sigdict[query_name]
        match = match_sigdict[match_name]

        cont = query.contained_by(match)
        contained[query_name] = cont

    return contained


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

    # now, do containment
    contained_1_in_2 = containment(matches_1_in_2, sigdict1, sigdict2)
    contained_2_in_1 = containment(matches_2_in_1, sigdict2, sigdict1)

    # summary stats
    CONTAIN_THRESHOLD=0.95
    IDENT_THRESHOLD=0.80
    
    print('thresholds:')
    print('min Jaccard similarity for any match:', THRESHOLD)
    print('to score as identical, similarity must be >=', IDENT_THRESHOLD)
    print('to score as contained, containment must be >=', CONTAIN_THRESHOLD)

    # 1 in 2
    c_ident = 0
    c_match = 0
    c_contain = 0
    c_no_match = 0
    c_no_contain = 0
    for query_name in sigdict1:
        best_match = None
        similarity = 0.0
        cont = 0.0
        
        if query_name in matches_1_in_2:
            (best_match, similarity) = matches_1_in_2[query_name]
        if query_name in contained_1_in_2:
            cont = contained_1_in_2[query_name]

        if not best_match:
            c_no_match += 1
        else:
            c_match += 1
            
        if cont < CONTAIN_THRESHOLD:
            c_no_contain += 1
        else:
            c_contain += 1

        if similarity > IDENT_THRESHOLD:
            c_ident += 1

    print('----')
    print('{} vs {}: {} signatures'.format(args.dir1, args.dir2, len(sigdict1)))
    print('identical count:', c_ident)
    print('containment count:', c_contain)
    print('matches:', c_match)
    
    print('no match:', c_no_match)
    print('no contain:', c_no_contain)

    # 2 in 1
    c_ident = 0
    c_match = 0
    c_contain = 0
    c_no_match = 0
    c_no_contain = 0
    for query_name in sigdict2:
        best_match = None
        similarity = 0.0
        cont = 0.0
        
        if query_name in matches_2_in_1:
            (best_match, similarity) = matches_2_in_1[query_name]
        if query_name in contained_2_in_1:
            cont = contained_2_in_1[query_name]

        if not best_match:
            c_no_match += 1
        else:
            c_match += 1
            
        if cont < CONTAIN_THRESHOLD:
            c_no_contain += 1
        else:
            c_contain += 1

        if similarity > IDENT_THRESHOLD:
            c_ident += 1

    print('----')
    print('{} vs {}: {} signatures'.format(args.dir2, args.dir1, len(sigdict2)))
    print('identical count:', c_ident)
    print('containment count:', c_contain)

    print('matches:', c_match)
    print('no match:', c_no_match)
    print('no contain:', c_no_contain)

                    
if __name__ == '__main__':
    main()
