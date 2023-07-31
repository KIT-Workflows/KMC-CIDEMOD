from initialize import *
from support import *

###############################################################################
# reading the input file from xml file as yaml
with open('rendered_wano.yml') as file:
    wano_file = yaml.full_load(file)

xdim     = wano_file["Lattice size"]['Xdim']
ydim     = wano_file["Lattice size"]['Ydim']
T        = wano_file["Temperature (K)"]['T']
SaveStep = wano_file["Save step"]['SaveStep']
 #barriers E1 to E17
E1      = wano_file["Activation energy barriers (eV)"]['E1']
E2      = wano_file["Activation energy barriers (eV)"]['E2']
E3      = wano_file["Activation energy barriers (eV)"]['E3']
E4      = wano_file["Activation energy barriers (eV)"]['E4']
E5      = wano_file["Activation energy barriers (eV)"]['E5']
E6      = wano_file["Activation energy barriers (eV)"]['E6']
E7      = wano_file["Activation energy barriers (eV)"]['E7']
E8      = wano_file["Activation energy barriers (eV)"]['E8']
E9      = wano_file["Activation energy barriers (eV)"]['E9']
E10     = wano_file["Activation energy barriers (eV)"]['E10']
E11     = wano_file["Activation energy barriers (eV)"]['E11']
E12     = wano_file["Activation energy barriers (eV)"]['E12']
E13     = wano_file["Activation energy barriers (eV)"]['E13']
E13     = wano_file["Activation energy barriers (eV)"]['E13']
E14     = wano_file["Activation energy barriers (eV)"]['E14']
E15     = wano_file["Activation energy barriers (eV)"]['E15']   
E16     = wano_file["Activation energy barriers (eV)"]['E16']
E17     = wano_file["Activation energy barriers (eV)"]['E17']
E18     = wano_file["Activation energy barriers (eV)"]['E18']
E19     = wano_file["Activation energy barriers (eV)"]['E19']


'''think about the holes in the structure ---> different colors!!!'''
# ###############################################################################
'''
combine the events in one for similar reactions and name them something that you named before <= 4nm

'''
'''Events'''
e1  = event(["E", "S"], ["E", "P"], 0.24,0) # First electron reduction next to Electrode
e2  = event(["F", "S"], ["F", "P"],0.36,0)  # First electron reduction next to SEI (<=4nm)
e3  = event(["E", "P"], ["E", "F"],0.36,0)  # Second electron reduction next to Electrode
e4  = event(["P", "P"], ["O","S"],0.372,0)
e5  = event(["O", "O"], ["OO"],0.43,0)
e6  = event(["O", "OO"], ["C"], 0.44,0)
e7  = event(["F", "P"], ["F", "F"],0.36,0)  # second electron reduction next to SEI (<=4nm)
e8  = event(["O", "C"], ["C"], 0.49,0)
e9  = event(["OO", "C"], ["C"], 0.51,0)
e10 = event(["OO", "OO"], ["C"], 0.47,0)
e11 = event(["C", "C"], ["C"], 0.49,0)
e12 = event(["OO", "S"], ["S", "OO"],0.3770,0)
e13 = event(["O", "S"], ["S", "O"], 0.3755,0)
e14 = event(["P", "S"], ["S", "P"],0.374 ,0)
e15 = event(["C","S"], ["S", "C"],0.38,0)
e16 = event(["P", "A"], ["A","S"], 0.01,0)
e17 = event(["O", "A"], ["A", "S"], 0.01,0)
e18 = event(["OO", "A"], ["A", "S"],0.01,0)
e19 = event(["C", "A"], ["A", "S"],0.01,0)
# list of events
Events = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19]
#diretion list which helps with cluster diffusion
direction = [[1, 0], [-1, 0], [0, 1], [0, -1]]

###############################################################################
E_barrier = [E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12,E13, E14, E15, E16, E17, E18, E19]
# load the barriers for sample
for i in Events:
    i.barrier = E_barrier[Events.index(i)]
#calculate rates
for i in Events:
    i.rate = bar_rate(i.barrier,T)
#=================================================================================
"""
collection function
"""
"""
collection function
"""
def pre_event(top_list):
    pre_list = []
    for s in top_list:
            # determine the type of sites
            ss = Full(s)
            if ss != "C":
                NBR = [i for i in s.nbr if i not in s.bonds]
                for nbr in NBR:
                        ns = Full(nbr)
                        # to divide the list into noncluster elements
                        '''to all diffusions'''
                        if ss == "E" and ns == "S":
                                # no need to check the distance
                                    pre_list.append(new_event(e1.reactant, e1.product, e1.barrier, e1.rate, [s, nbr, Events.index(e1)]))
                                # here we need to check the distance, close to the 4nm
                        if ss == "F" and  ns == "S" and nbr.coordinate[1] <= 5:
                                    pre_list.append(new_event(e2.reactant, e2.product, e2.barrier, e2.rate, [s, nbr, Events.index(e2)]))
                    

                        if ss == "P" and s.coordinate[1] <= 4 and nbr.coordinate[1] <= 4 and ns in ["E", "F", "C"]:
                                # for <=4 nm distance it could be F or C (should be deposited)
                                if nbr.species == "C" and nbr.status[3][1] == "SEI":
                                    C = cluster_x(nbr)
                                    coordc = [i for i in C if i.coordinate[1] <= 4]

                                    if coordc:
                                        ev = e7
                                        pre_list.append(new_event(ev.reactant, ev.product, ev.barrier,ev.rate,
                                                                [s, nbr, Events.index(ev)]))
                                # this is the case for F <=4
                                elif nbr.species == "F":
                                    ev = e7
                                    pre_list.append(new_event(ev.reactant, ev.product, ev.barrier,ev.rate,
                                                        [s, nbr, Events.index(ev)]))
                                
                                # this is the case next to the electrode
                                elif nbr.species == "E":
                                    ev = e3
                                    pre_list.append(new_event(ev.reactant, ev.product, ev.barrier,ev.rate,
                                                        [s, nbr, Events.index(ev)]))

                        if ss == "P" and ns == "P":
                            ev = e4
                            pre_list.append(
                                new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, Events.index(ev)]))

                        if ss == "O" and ns == "O":
                            ev = e5
                            pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, Events.index(ev)]))

                        if (ss == "O" and ns == "OO") or (ss == "OO" and ns == "O") :
                            ev = e6
                            o2 = [i for i in [s,nbr] if i.bonds]
                            o = [i for i in [s,nbr] if i.bonds == []]
                            pre_list.append( new_event(ev.reactant, ev.product,ev.barrier,ev.rate, [o2[0], o[0], Events.index(ev)]))

                            
                        if ss == "O"  and ns == "C":
                                ev = e8
                                pre_list.append(
                                    new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [nbr, s, Events.index(ev)]))
                                
                        if ss == "OO" and ns == "C":
                                ev = e9
                                pre_list.append(
                                    new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [nbr, s, Events.index(ev)]))
        
                        if ss == "OO" and ns == "OO":
                            ev = e10
                            pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, Events.index(ev)]))


                        if ss == "OO" and ns == "S":
                                ev = e12
                                new = [i for i in nbr.nbr if i.species == "S"]
                                if new:
                                    nw_nbr = choice(new)
                                    pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, ev.rate,
                                                      [s, nbr, nw_nbr, Events.index(ev)]))

                        if ss == "O" and ns == "S":
                            ev = e13
                            pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, Events.index(ev)]))

                        if (ss == "P" and ns == "S"):
                            ev = e14
                            pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, Events.index(ev)]))

                        if ss == "P" and ns == "A":
                            ev = e16
                            pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, Events.index(ev)]))

                        if ss == "O" and ns == "A":
                            ev = e17
                            pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, Events.index(ev)]))

                        if ss == "OO" and ns == "A":
                            ev = e18
                            pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, Events.index(ev)]))


            elif ss == "C":               
                                           
                # this is the case for C <=4 ,first electron reduction
                NBR = [i for i in s.nbr if i.species == "S"]
                for nbr in NBR:
                    if s.status[3][1] == "SEI" and s.coordinate[1] <= 4 and nbr.coordinate[1] <= 5:
                        C = cluster_x(nbr)
                        coordc = [i for i in C if i.coordinate[1] <= 4]
                        if coordc:
                            pre_list.append(new_event(e2.reactant, e2.product, e2.barrier, e2.rate,[s, nbr, Events.index(e2)]))
                # this is for the second electron reduction
                NBR = [i for i in s.nbr if i.species == "P"]
                for nbr in NBR:
                    if s.status[3][1] == "SEI" and s.coordinate[1] <= 4 and nbr.coordinate[1] <= 4:
                        C = cluster_x(nbr)
                        coordc = [i for i in C if i.coordinate[1] <= 4]
                        if coordc:
                            pre_list.append(new_event(e7.reactant, e7.product, e7.barrier, e7.rate,[nbr, s, Events.index(e7)]))
                             
                # When it is close to electrode or an inorganic layer
                cb = cluster_x(s)
                cb_cord = [i.coordinate for i in cb]
                all_nbr = [j for j in list(dict.fromkeys(list(itertools.chain(*([i.nbr for i in cb]))))) if j not in cb]
                named_list = [i.species for i in all_nbr]
                if "O" in named_list:
                    n_spec = [i for i in all_nbr if i.species == "O"]
                    if Full(n_spec[0]) == "O":
                        pre_list.append(
                            new_event(e8.reactant, e8.product, e8.barrier, e8.rate, [s, n_spec[0], Events.index(e8)]))
                    elif Full(n_spec[0]) == "OO":
                        pre_list.append(
                            new_event(e9.reactant, e9.product, e9.barrier, e9.rate, [s, n_spec[0], Events.index(e9)]))
                if "A" in named_list:
                    if "E" not in named_list and "F" not in named_list:
                        pre_list.append(
                            new_event(e19.reactant, e19.product, e19.barrier, e19.rate, [s, s, cb, Events.index(e19)]))

                # single diff of particles:
                #not for attached C
                if s.status[3][1] != "SEI":
                    # directions lists
                    if xdim not in [i.coordinate[0] for i in cb]:
                        rs = [i for i in all_nbr if [(i.coordinate[0] - 1), i.coordinate[1]] in cb_cord]
                    else:
                        rs = []

                    if 1 not in [i.coordinate[0] for i in cb]:
                        ls = [i for i in all_nbr if [(i.coordinate[0] + 1), i.coordinate[1]] in cb_cord]
                    else:
                        ls = []

                    if ydim not in [i.coordinate[1] for i in cb]:
                        us = [i for i in all_nbr if [i.coordinate[0], (i.coordinate[1] - 1)] in cb_cord]
                    else:
                        us = []

                    if 1 not in [i.coordinate[1] for i in cb]:
                        ds = [i for i in all_nbr if [i.coordinate[0], (i.coordinate[1] + 1)] in cb_cord]
                    else:
                        ds = []
                    dir_list = [rs, ls, us, ds]

                    for i in dir_list:  # 0 1 2 3
                        if i:
                            if len(i) == len([j for j in i if j.species == "S"]):
                                # not more move when reached to F
                                ev = e15
                                ind = dir_list.index(i)
                                '''diffusion of clusters depends on the size of them !!!!!!!!!!!!!!!!!!!!'''
                                # pass the direction of diffusion using the direction vector RATE DEPENDS ON SIZE
                                evra = ev.rate #* exp(-(len(cb)) / 30)
                                pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, evra,
                                                        [s, direction[ind], i, all_nbr, cb, Events.index(ev)]))

                            C_test = [j for j in i if j.species == "C"]
                            if C_test:
                                k = C_test[0]
                                ev = e11
                                # cb2 = sb + [s]
                                evra = ev.rate
                                '''do we need to add this weihght?????????'''
                                # evra = ev.rate * exp(-(len(cb)) / 30)
                                # all surounding sisepciestes
                                all_nbr2 = all_nbr2 = [j for j in list(dict.fromkeys(list(itertools.chain(*([i.nbr for i in cb]))))) if j not in cb]
                                pre_list.append(new_event(ev.reactant, ev.product, ev.barrier, evra,
                                                        [s, k, all_nbr, all_nbr2, Events.index(ev)]))
    return pre_list