# required packages
from support import *
from initialize import *
from collection import *
from update import *
"""
************************************************* MAIN BODY ************************************************************
"""
#record the inital time
start_time = time.time()
'''
this loop determines termination condition
 True  : up to end on no more move
 < step: at desired step
 <t    : at desired time
'''
while True:
    temp = []
    del temp[:]
    '''get single representative of each cluster'''
    spec_in_top = [i for i in lattice if i.species == "C" and i.status[3][1] == '']
    top_list = [i for i in top_list if i not in spec_in_top]
    catch = []
    for i in spec_in_top:
        if i not in catch:
            temp = cluster_x(i)
            catch = catch + temp
            top_list.append(i)
    top_list =  list(dict.fromkeys((top_list)))
    tt = time.time()
    temp = pre_event(top_list)
    pre = temp + pre
    if len(pre) == 0:
        spec_in_top = [i for i in lattice if i.species not in ["S","A"]]
        catch = []
        spec_top = []
        for i in spec_in_top:
            if i not in catch:
                temp = cluster_x(i)
                catch = catch + temp
                spec_top.append(i)
        spec_top = list(dict.fromkeys((spec_top +  [j for j in lattice if j.status == 'SEI' and "S" in [k.species for k in j.nbr]])))
        pre = pre_event(spec_top)
        del spec_top[:]
        del spec_in_top[:]
        if len(pre) == 0:
            stat = [time.time() - start_time, t, num]
            stored = thickness(xdim, lattice)
            tracking['thickness'] = float(stored[0])
            tracking['porosity'] = float(stored[1])
            try:
                with open('../kmc_data.yml', 'w') as out:
                    yaml.dump(tracking, out, default_flow_style=False)
            except IOError:
                print("I/O error")
            break

    else:
        '''selection scheme'''
        R_index = []
        del R_index[:]
        # t1 = time.time()
        Events_per_type = [[i for i in pre if i.new_coord[-1] == j] for j in range(19)]
        giant_pre = [list(accumulate([i.rate for i in Events_per_type[j]])) for j in range(len(Events_per_type))]
        # print(giant_pre)
        giant_pre_index = [giant_pre.index(i) for i in giant_pre if i != []]
        giant_pre = [i for i in giant_pre if i != []]
        R_index = [0] + [i[-1] for i in giant_pre]
        R_tot = list(accumulate(R_index))
        cum_rate = R_tot[-1]
        '''Event type selection among 19 event types'''
        while True:
            rho1 = random()
            r = rho1 * (R_tot[-1])
            if (R_tot[0] < r and r <= R_tot[-1]):
                break
        rate = []
        for i in range(1, len(R_tot)):
            temp = ratee(i, r, R_tot)
            rate.append(temp)
        which = [i for i in rate if i]
        # index
        type_index = which[0][1]
        # which one in main list
        type_index2 = giant_pre_index[type_index]
        '''the second selection same as the first one'''
        rho2 = random()
        event_at = int(len(Events_per_type[type_index2]) * rho2)
        event = Events_per_type[type_index2][event_at]
        """
        Time increment
        """
        dt = -(np.log(random()) / cum_rate)
        t = t + dt
        """
        UPDATING the list of sites using the function considering the type of selected event.
        """
        '''lattice update'''
        lattice = update(event, t, lattice)
        """
        removing old events related to selected sites and nbrs
        """
        # it would be usefull later to keep C around eliminated F
        if event.new_coord[-1] == 14:
            all_nbr = event.new_coord[3]
            main_before = [event.new_coord[0]] + event.new_coord[0].bonds + all_nbr
            '''AND or OR . taking care of attached SEI'''
            if "SEI" in [i.status[3][1]  for i in all_nbr]:
                pre = [i for i in pre if event.new_coord[0] not in all_nbr]
                pre = [i for i in pre if event.new_coord[1] not in all_nbr]
                '''put flag for fixed C'''
            pre = [i for i in pre if i.new_coord[0] not in main_before and i.new_coord[1] not in main_before]
        else:
            s1 = event.new_coord[0]
            s2 = event.new_coord[1]
            main_ = [s1] + s1.nbr + [s2] + s2.bonds + s2.nbr
            main_before = list(dict.fromkeys(main_))
            pre = [i for i in pre if i.new_coord[0] not in main_before and i.new_coord[1] not in main_before]
            test = [i.bonds for i in main_ if i.species == "C"]
            if test:
                all_nbr = list(itertools.chain(*(test)))
                main_before_ = [j for j in list(dict.fromkeys(list(itertools.chain(*([i.nbr for i in all_nbr]))))) if
                                j not in all_nbr]
                main_before = main_before_ + main_before
                pre = [i for i in pre if i.new_coord[0] not in main_before and i.new_coord[1] not in main_before]



        pre = [i for i in pre if i.new_coord[-1] not in [10,14,18]]
        main_before = list(dict.fromkeys(main_before))
        pre1 = [i for i in pre if i.new_coord[-1] == 11]
        pre2 = [i for i in pre if i not in pre1]
        pre1 = [i for i in pre1 if i.new_coord[2] not in main_before]
        pre = pre1 + pre2
        counter = event.new_coord[-1]
        num = num + 1
        """
        Top-list changes 
        """
        del top_list[:]
        '''just kick S and A out of the  list'''
        top_list = [i for i in main_before if i.species not in ["A", "S"]]
        if num % frac_save == 0:
            stat = [time.time() - start_time, t, num]
            # creating the dataframe from traj
            tracking["Status"] = round_float_list(stat, 30)
            stored = thickness(xdim, lattice)
            tracking['thickness'] = float(stored[0])
            tracking['porosity'] = float(stored[1])
            try:
                with open('../kmc_data.yml', 'w') as out:
                    yaml.dump(tracking, out, default_flow_style=False)
            except IOError:
                print("I/O error")
            keep = num


