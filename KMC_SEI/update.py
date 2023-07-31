from initialize import *
from support import *
from collection import *
###############################################################################
"""
update function
"""
def update(ev,t,lattice):
    # e1/e2 = event([["Surf"],["EC"]],[["Surf"],["P"])
    if (ev.reactant == e1.reactant) or (ev.reactant == e2.reactant):
        ev.new_coord[1].species = "P"
        return lattice
    
    # e3/e7 = event([["Surf"],["P"]],[["Surf"],["F"])
    elif (ev.reactant == e3.reactant and ev.product == e3.product) or (ev.reactant == e7.reactant and ev.product == e7.product):
        ev.new_coord[0].species = "F"
        return lattice

    # e4 = event([["P"],["P"]],[["O"],["S"]])
    elif ev.reactant == e4.reactant:
        # pick one randomly
        pick = [ev.new_coord[0], ev.new_coord[1]]
        p = choice(pick)
        p.species = "O"
        g = [i for i in pick if i != p]
        g[0].species = "S"
        g[0].bonds = []
        p.bonds = []
        p.status = [[p.coordinate, t], [],'', [0,'']]
        return lattice

    # e5 = event([["O"],["O"]],[["OO"]])
    elif ev.reactant == e5.reactant :
        ev.new_coord[0].bonds = []
        ev.new_coord[1].bonds = []
        ev.new_coord[0].bonds = [ev.new_coord[1]]
        ev.new_coord[1].bonds = [ev.new_coord[0]]
        return lattice
    
    # e6 = event([["O"],["OO"]],[["C"]])
    elif ev.reactant == e6.reactant:
        s1 = ev.new_coord[0]
        s2 = ev.new_coord[1]
        sb = s1.bonds[0]
        s1.bonds = []
        s2.bonds = []
        sb.bonds = []
        s2.status = [s2.status[0], [s2.coordinate, t], num, [0,'']]
        s1.status = [s1.status[0], [s1.coordinate, t], num, [0,'']]
        sb.status = [sb.status[0], [sb.coordinate, t], num, [0,'']]
        s1.species = "C"
        s2.species = "C"
        sb.species = "C"
        if s2 in s1.nbr:
            s1.bonds = [s2,sb]
            s2.bonds = [s1]
            sb.bonds = [s1]
        else:
            sb.bonds = [s1,s2]
            s1.bonds = [sb]
            s2.bonds = [sb]
        # SEI flag check
        new = [s1,s2,sb]
        all_nbr = [j for j in list(dict.fromkeys(list(itertools.chain(*([i.nbr for i in new]))))) if j not in new]
        named_list = [i.species for i in all_nbr]
        if "E" in named_list or "F" in named_list:
            for i in new:
                i.status[3][1] = "SEI"
        else:
            for i in new:
                i.status[3][1] = ''
        return lattice

    # e8= event([["O"],["C"]],[["C"]])
    elif ev.reactant == e8.reactant:
        c = ev.new_coord[0] #C
        o = ev.new_coord[1] # O

        cb = cluster_x(c)
        o.species = "C"
        c.bonds   = c.bonds + [o]
        cnb = [i for i in o.nbr if i in cb]
        o.bonds   = cnb
        o.status  = c.status
        for item in cnb:
            item.bonds = item.bonds + [o]
        o.status = [o.status[0], [o.coordinate, t], c.status[2], c.status[3]]

        new = cluster_x(c)
        all_nbr = [j for j in list(dict.fromkeys(list(itertools.chain(*([i.nbr for i in new]))))) if j not in new]
        named_list = [i.species for i in all_nbr]
        if "E" in named_list or "F" in named_list:
            for i in new:
                i.status[3][1] = "SEI"
        else:
            for i in new:
                i.status[3][1] = ''
        return lattice

    # e9= event([["OO"],["C"]],[["C"]])
    elif ev.reactant == e9.reactant:
        c = ev.new_coord[0] #C
        clust = cluster_x(c)
        o2 = ev.new_coord[1] #O2
        ob = o2.bonds[0]
        o2.species = "C"
        ob.species = "C"
        ob.status = c.status
        o2.status = c.status
        o2b =  [j for j in o2.nbr if j in clust]
        o2.bonds = o2.bonds + o2b
        for items in o2b:
                items.bonds = items.bonds + [o2]
        obb =  [j for j in ob.nbr if j in clust]
        ob.bonds = ob.bonds + obb
        for items in obb:
                items.bonds = items.bonds + [ob]

        o2.status = [o2.status[0], [o2.coordinate, t], c.status[2], c.status[3]]
        ob.status = [ob.status[0], [ob.coordinate, t], c.status[2], c.status[3]]

        new = cluster_x(c)
        all_nbr = [j for j in list(dict.fromkeys(list(itertools.chain(*([i.nbr for i in new]))))) if j not in new]
        named_list = [i.species for i in all_nbr]
        
        if "E" in named_list or "F" in named_list:
            for i in new:
                i.status[3][1] = "SEI"
        else:
            for i in new:
                i.status[3][1] = ''
        return lattice

    # e10= event([["OO"],["OO"]],[["C"]])
    elif ev.reactant == e10.reactant:
        o1  = ev.new_coord[0]
        o2  = ev.new_coord[1]
        o1b = o1.bonds[0]
        o2b = o2.bonds[0]
        for i in [o2,o2b,o1,o1b]:
            i.species = "C"
            i.bonds = [j for j in [o2,o2b,o1,o1b] if j in i.nbr and j!=i]

        o1.status = [o1.status[0], [o1.coordinate, t], num, [0,'']]
        o2.status = [o2.status[0], [o2.coordinate, t], num, [0,'']]
        o1b.status = [o1b.status[0], [o1b.coordinate, t], num, [0,'']]
        o2b.status = [o2b.status[0], [o2b.coordinate, t], num, [0,'']]

        new = cluster_x(o1)
        all_nbr = [j for j in list(dict.fromkeys(list(itertools.chain(*([i.nbr for i in new]))))) if j not in new]
        named_list = [i.species for i in all_nbr]
        if "E" in named_list or "F" in named_list:
            for i in new:
                i.status[3][1] = "SEI"
        else:
            for i in new:
                i.status[3][1] = ''
        return lattice

    # e11 C+C---> C 
    elif ev.reactant == e11.reactant:
        c1 = ev.new_coord[0]
        c2 = ev.new_coord[1]
        clust1 = cluster_x(c1)
        clust2 = cluster_x(c2)
        for items in clust1:
            bonded = [j for j in items.nbr if j in clust2]
            items.bonds = items.bonds + bonded
            for k in bonded:
                k.bonds = k.bonds + [items]
        cb = cluster_x(c1)
        if "SEI" in [i.status[3][1] for i in cb]:
            for j in cb:
                j.status[3][1]="SEI"
        return lattice

    #e12 = event([["OO"],["S"]],[["S"],["OO"]])      
    elif ev.reactant == e12.reactant:
        O1 = ev.new_coord[0]  # OO
        n1 = ev.new_coord[1]  # EC/IS/S
        O2 = O1.bonds[0]
        n2 = ev.new_coord[2]  # S/EC/[]
        O1.bonds = []
        O2.bonds = []
        n1.bonds = []
        n2.bonds = []
        O1.species = "S"
        O2.species = "S"
        n1.species = "O"
        n2.species = "O"
        n1.status = O1.status
        n2.status = O2.status
        O1.status = [[], [],'', [0,'']]
        O2.status = [[], [],'',[0,'']]
        n1.bonds = [n2]
        n2.bonds = [n1]
        return lattice

    # e13= event([["O"],["S"]],[["S"],["O"]]) 
    elif ev.reactant == e13.reactant:
        s1 = ev.new_coord[0]
        s2 = ev.new_coord[1]
        s1.species = "S"
        s2.species = "O"
        s2.status = s1.status
        s1.status = [[], [], '', [0,'']]
        s1.bonds = []
        s2.bonds = []
        return lattice


    #  e14= event([["P"],["S"]],[["S"],["P"]])
    elif ev.reactant == e14.reactant:
        s1 = ev.new_coord[0]
        s2 = ev.new_coord[1]
        s1.species = "S"
        s2.species = "P"
        s1.bonds = []
        s2.bonds = []
        return lattice
    

    # e15 C---> C
    elif ev.reactant == e15.reactant and ev.product == e15.product:
        
        s = ev.new_coord[0]
        # the direction vector
        d = ev.new_coord[1]
        sorted_cb = []
        sb = s.bonds
        cb = cluster_x(s)
        for i in cb:
            i.bonds = []
            i.species = "S"

        if d == [1, 0]:
            sort_cb = sorted([i.coordinate for i in cb], key=lambda k: [k[0], k[1]], reverse=True)
            for j in sort_cb:
                temp = [i for i in cb if i.coordinate == j]
                sorted_cb.append(temp[0])
            new = []
            for i in sorted_cb:
                # i.species = ["S"]
                jcoord = [i.coordinate[0] + d[0], i.coordinate[1] + d[1]]
                test = [k for k in i.nbr if k.coordinate == jcoord]
                test[0].species = "C"
                test[0].bonds = []
                test[0].status = [i.status[0], i.status[1], i.status[2], [i.status[3][0] + 1,'']]
                i.status = [[], [], '', [0,'']]
                new.append(test[0])
        elif d == [-1, 0]:
            sort_cb = sorted([i.coordinate for i in cb], key=lambda k: [k[0], k[1]], reverse=False)
            for j in sort_cb:
                temp = [i for i in cb if i.coordinate == j]
                sorted_cb.append(temp[0])
            new = []
            for i in sorted_cb:
                # i.species = ["S"]
                jcoord = [i.coordinate[0] + d[0], i.coordinate[1] + d[1]]
                test = [k for k in i.nbr if k.coordinate == jcoord]
                test[0].species = "C"
                test[0].bonds = []
                test[0].status = [i.status[0], i.status[1], i.status[2], [i.status[3][0] + 1,'']]
                i.status = [[], [], '', [0,'']]
                new.append(test[0])
        elif d == [0, 1]:
            sort_cb = sorted([i.coordinate for i in cb], key=lambda k: [k[1], k[0]], reverse=True)
            for j in sort_cb:
                temp = [i for i in cb if i.coordinate == j]
                sorted_cb.append(temp[0])
            new = []
            for i in sorted_cb:
                # i.species = ["S"]
                jcoord = [i.coordinate[0] + d[0], i.coordinate[1] + d[1]]
                test = [k for k in i.nbr if k.coordinate == jcoord]
                test[0].species = "C"
                test[0].bonds = []
                test[0].status = [i.status[0], i.status[1], i.status[2], [i.status[3][0] + 1,'']]
                i.status = [[], [], '', [0,'']]
                new.append(test[0])
        elif d == [0, -1]:
            sort_cb = sorted([i.coordinate for i in cb], key=lambda k: [k[1], k[0]], reverse=False)
            for j in sort_cb:
                temp = [i for i in cb if i.coordinate == j]
                sorted_cb.append(temp[0])
            new = []
            for i in sorted_cb:
                # i.species = ["S"]
                jcoord = [i.coordinate[0] + d[0], i.coordinate[1] + d[1]]
                test = [k for k in i.nbr if k.coordinate == jcoord]
                test[0].species = "C"
                test[0].status = [i.status[0], i.status[1], i.status[2], [i.status[3][0] + 1,'']]
                i.status = [[], [], '', [0,'']]
                test[0].bonds = []
                new.append(test[0])

        for cc in new:
            cc.bonds = [j for j in cc.nbr if j in new and j.species == "C"]
        cb = cluster_x(new[0])
        '''needed for crack procedure'''
        all_nbr = [j for j in list(dict.fromkeys(list(itertools.chain(*([i.nbr for i in cb]))))) if j not in cb]
        named_list = [i.species for i in all_nbr]
        if "E" in named_list or "F" in named_list:
            for i in new:
                i.status[3][1] = "SEI"
        else:
            for i in new:
                i.status[3][1] = ''

        return lattice


    # Escapes
    elif ev.reactant == e19.reactant or ev.reactant == e18.reactant or ev.reactant == e16.reactant \
            or ev.reactant == e17.reactant :
        '''if we need full space with solvent we have to change this part to '''
        s = ev.new_coord[0]
        if "P" in ev.reactant:
            s.species = "S"
            s.status =  [[], [], '', [0, '']]
            s.bonds = []
            # Counter[0] += 1
        elif "O" in ev.reactant:
            s.species = "S"
            s.status =  [[], [], '', [0, '']]
            s.bonds = []
            # Counter[1] += 1
        elif "OO" in ev.reactant:
            s.species = "S"
            s.bonds[0].species = "S"
            s.status =  [[], [], '', [0, '']]
            s.bonds[0].status  =  [[], [], '', [0, '']]
            s.bonds[0].bonds   = []
            s.bonds = []
            # Counter[2] += 2
        elif ev.reactant == e19.reactant:
            cb = ev.new_coord[2]
            bd = len(cb)
            for i in cb:
                i.bonds = []
                i.status =  [[], [], '', [0, '']]
                i.species = "S"
            # Counter[3] += bd

        return lattice
