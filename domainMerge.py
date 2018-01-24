from itertools import product
import numpy as np
class cluster(object):
    def __init__(self, elems):
        self.elems = elems
    def merge(self, cluster):
        self.elems.extend(cluster.elems)
        return self
    def dist(self, cluster):
        dmin=1.0
        for i,j in product(self.elems, cluster.elems):
            d = domainCluster._similarity(i,j)
            if d < dmin:
                dmin = d
        return dmin
class domainCluster(object):
    def __init__(self, inlist):
        self.inlist= inlist
    @staticmethod
    def _similarity(list1, list2):
        #assert list1 and list2 are 2-tuple(or 2-list)
        def overlap(list1, list2):
            ret = min(list1[1], list2[1]) - max(list1[0], list2[0]) + 1
            if ret > 0:
                return ret
            else:
                return 0
        def dist(l):
            return l[1] - l[0] + 1
        o = overlap(list1, list2)
        d1 = float(dist(list1))
        d2 = float(dist(list2))
        return (o/d1) * (o/d2)
    def hierarchical(self, threshold=0.9):
        def maxRange(c):
            start = 20000
            end = 0
            for l in c.elems:
                if l[0] < start:
                    start = l[0]
                if l[1] > end:
                    end  = l[1]
            return (start, end)
        def make_arr(cc):
            arr = np.zeros((len(cc),len(cc)))
            for i in xrange(len(cc)):
                for j in xrange(len(cc)):
                    if j < i:
                        arr[i][j] = arr[j][i]
                        continue
                    arr[i][j] = cc[i].dist(cc[j])
            return arr
            
            
        # a list of cluster
        cc = [cluster([l]) for l in self.inlist]
        #generate pairwise similarity
        arr = make_arr(cc)
        while(len(cc) > 1): #
            #calculate maximum similarity among clusters
            it = np.nditer(arr, flags=['multi_index'])
            merge = (0,0)
            max_sim = 0.0
            while not it.finished:
                if it.multi_index[1] == it.multi_index[0]:
                    it.iternext()
                    continue
                if it[0] > max_sim:
                    max_sim = it[0]
                    merge = it.multi_index
                it.iternext()

            if(max_sim < threshold ):
                break
            else:
                #merge 2 cluster
                c1 = cc[merge[0]]
                c2 = cc[merge[1]]
                cc = [c for i,c in enumerate(cc) if i not in merge]
                cc.append(c1.merge(c2))
                #update pairwise similarity
                arr = make_arr(cc)
        #arrange result and return
        ret = [maxRange(c) for c in cc]
        return ret
        
