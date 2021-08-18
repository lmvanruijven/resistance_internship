import numpy as np

def get_i_j(N, V, flux2):
    flux = np.array(flux2)
    flux = flux[flux != 0]


    ################Get (i,j)
    i_j = np.array(V)


    ################Get i
    XV_i = np.argwhere(N<0)
    XV_i = XV_i.astype(str)
    XV_i = np.transpose(
                    np.core.defchararray.add(["X", "V"], XV_i))
    #print( XV_i, "=X and V for i\n")
    #put np array in a dict
    dict = {}
    t=0
    for v in XV_i[1]:
        dict[v] = ""
    for v in XV_i[1]:
        current_dict_v = dict[v]
        dict[v] = current_dict_v + XV_i[0][t]
        t=t+1
    #print(dict)
    i_distinct = np.array(list(dict.values()))
    #print(i_distinct)
    i_total = np.repeat(i_distinct, flux)
    i_text = np.array(np.unique(i_total, return_counts=True))
    i = i_text[1].astype(int)
    i = i / np.sum(i)
    #print("i probabilities are: ", i, "\n with following i's:", i_text[0], "\n")


    ################Get j
    XV_j = np.argwhere(N>0)
    XV_j = XV_j.astype(str)
    XV_j = np.transpose(
                    np.core.defchararray.add(["X", "V"], XV_j))
    #print( XV_j, "=X and V for j\n")
    #put np array jn a djct
    dict = {}
    t=0
    for v in XV_j[1]:
        dict[v] = ""
    for v in XV_j[1]:
        current_dict_v = dict[v]
        dict[v] = current_dict_v + XV_j[0][t]
        t=t+1
    j_distinct = np.array(list(dict.values()))
    j_total = np.repeat(j_distinct, flux)
    j_text = np.array(np.unique(j_total, return_counts=True))
    j = j_text[1].astype(int)
    j = j / np.sum(j)
    #print("j probabilities are: ", j, "\n with following j's:", j_text[0])


    return i,j,i_j
