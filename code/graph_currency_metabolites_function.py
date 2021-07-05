# -*- coding: utf-8 -*-
import itertools
import networkx as nx

def get_metpair(rea, pi_pairs1, h_pairs1, pi_pairs2, h_pairs2, nh4_pairs, other_pairs, currency_mets):
    # get the metabolite links for a reaction, excluding links through currency metabolites
    # processing in the order of P & H transfer, N transfer and other transfers
    sub_pro = []
    mark = 0  # mark if there is currency metabolite pairs for P/H transfer in the reaction
    c2mark = 0  # mark if there is currency metabolite pairs for P/H transfer in the reaction, second batch
    nmark = 0  # mark if there is nh4 transfer currency metabolite pairs in the reaction, deal seperately
    omark = 0  # mark if there is other group transfer currency metabolite pairs in the reaction, deal seperately
    cmet = []  # temorary currency metabolite list
    c2met = []  # temorary currency metabolite list, second batch
    ncmet = []  # temorary currency metabolite list for N transfer
    ocmet = []  # temorary currency metabolite list for other pair transfer
    phmet = []  # metabolite for P/H transfer, to be excluded again in the last step if >1 pairs still remaining, especially for UDP ADPsugars
    ex_pairs1 = pi_pairs1+h_pairs1  # excluded currency metabolite pairs, first batch
    ex_pairs2 = pi_pairs2+h_pairs2  # excluded currency metabolite pairs, second batch
    for sp in ex_pairs1:
        phmet.append(sp[0])
        phmet.append(sp[1])
    for sp in ex_pairs2:
        phmet.append(sp[0])
        phmet.append(sp[1])
    phmet = list(set(phmet))  # remove repeats
    # also check if CoA is a currency metabolite in the last step
    phmet += ['coa_c', 'coa_p', 'coa_e']
    subs = [
        m.id for m in rea.reactants if 'C' in m.elements if m.id not in currency_mets]
    pros = [m.id for m in rea.products if 'C' in m.elements if m.id not in currency_mets]
    for s, p in itertools.product(subs, pros):
        if (s, p) in ex_pairs1 or (p, s) in ex_pairs1:  # need to consider direction
            mark = 1
            cmet.append(s)
            cmet.append(p)
        if (s, p) in ex_pairs2 or (p, s) in ex_pairs2:  # need to consider direction
            c2mark = 1
            c2met.append(s)
            c2met.append(p)
        if (s, p) in nh4_pairs or (p, s) in nh4_pairs:  # for nh4 transfer
            nmark = 1
            ncmet.append(s)
            ncmet.append(p)
        if (s, p) in other_pairs or (p, s) in other_pairs:  # for other pair transfer
            omark = 1
            ocmet.append(s)
            ocmet.append(p)
    # if len(sub_pro)>1 and mark==1:
    if mark == 1:  # process in order
        subs = [m for m in subs if m not in cmet]
        pros = [m for m in pros if m not in cmet]
    if c2mark == 1:  # process in order
        subsn = [m for m in subs if m not in c2met]
        prosn = [m for m in pros if m not in c2met]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if nmark == 1:
        subsn = [m for m in subs if m not in ncmet]
        prosn = [m for m in pros if m not in ncmet]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if omark == 1:
        subsn = [m for m in subs if m not in ocmet]
        prosn = [m for m in pros if m not in ocmet]
        if subsn and prosn:  # proceed if only there are still other metabolites in the reactant and product list
            subs = subsn
            pros = prosn
    if len(subs) > 1:  # to remove UTP, UDP etc.
        subsn = subs
        for m in subsn:
            if m in phmet:
                subs.remove(m)
                if len(subs) == 1:
                    break
    if len(pros) > 1:  # to remove UTP, UDP etc.
        prosn = pros
        for m in prosn:
            if m in phmet:
                pros.remove(m)
                if len(pros) == 1:
                    break
    for s, p in itertools.product(subs, pros):
        sub_pro.append((s, p))
    return sub_pro

def search_path(G, result_save_path, substrate, target,n):#substrate and target
    allrecID = [] 
    G1 = G
    tm = open(result_save_path + substrate+'+'+target+'_paths_shortest' + '.txt', 'w')  #pathway file
    if (substrate in G1.nodes()) and (target in G1.nodes()):
        try:
            i=0
            allpaths = nx.all_simple_paths(G1, substrate, target, n)
            for eachpath in allpaths:  #eachpath is a list
                i+=1
                allrecID.append(eachpath)
                tmpath = '\t'.join(eachpath)
                tm.write(tmpath + '\n')  # write all pathways to allpaths.txt file
                tm.flush()
        except nx.NetworkXNoPath:
            pass
    print('finished!')
    return allrecID

