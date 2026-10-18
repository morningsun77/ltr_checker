from Bio import SeqIO
# define a dict for genome (keys: chromosome's id, values: chromosome's sequence)
def genome_dict(file):
    seqfile = [x for x in SeqIO.parse(file, 'fasta')]
    genome = dict()
    for i in seqfile:
        genome[str(i.id)] = str(i.seq)
    return genome

# Splitting consecutive integer subsequences
def split_num_l(num_lst):
    num_lst_tmp = [int(n) for n in num_lst]
    sort_lst = sorted(num_lst_tmp)  # ascending
    len_lst = len(sort_lst)
    i = 0
    split_lst = []
    
    tmp_lst = [sort_lst[i]]
    while True:
        if i + 1 == len_lst:
            break
        next_n = sort_lst[i+1]
        if sort_lst[i] + 1 == next_n:
            tmp_lst.append(next_n)
        else:
            split_lst.append(tmp_lst)
            tmp_lst = [next_n]
        i += 1
    split_lst.append(tmp_lst)
    return split_lst

def split_ltr_range(list):
    result = []
    for txt in list:
        txt = txt.split('\t')
        End = int(txt[2])
        Start = int(txt[1])
        if End - Start > 100000000:
            distance = End - Start
            num = distance // 10
            for Num in range(0, num*10, num):
                Start = Start + Num
                End = Num + num
                result.append(txt[0] + '\t' + str(Start) + '\t' + str(End))
                Start = int(txt[1])
                End = int(txt[2])
        elif End - Start < 100000000 and End - Start > 50000000:
            distance = End - Start
            num = distance // 2
            for Num in range(0, num*2, num):
                Start = Start + Num
                End = num + Num
                result.append(txt[0] + '\t' + str(Start) + '\t' + str(End))
                Start = int(txt[1])
                End = int(txt[2])
        else:
            result.append(txt[0] + '\t' + str(Start) + '\t' + str(End))
    return result

# transform predict's results to position information
def ltr_range(chr_id,list_result,length,slide):
    new_result = sum(list_result,[])
    num = 0
    temp = []
    for i in new_result:
        if i == 1:
            temp.append(num)
        num += 1
    res = []
    if len(temp) == 1:
        start = temp[0] * slide 
        end = temp[-1] * slide + 50000
        l = chr_id + "\t" + str(start) + "\t" + str(end)
        res.append(l)
    elif len(temp) == 0:
        
        l = chr_id + "has no LTR_RT."
    elif len(temp) > 1:
        consistent_list = split_num_l(temp)
        if len(consistent_list) >= 2:
            for i in range(1,len(consistent_list),1):
                if consistent_list[i][0] - consistent_list[i-1][-1] > int(50000 // slide):
                    start = consistent_list[i-1][0] * slide 
                    end = consistent_list[i-1][-1] * slide + 50000
                    l = chr_id + "\t" + str(start) + "\t" + str(end)
                    res.append(l)
                elif consistent_list[i][0] - consistent_list[i-1][-1] <= int(50000 // slide) and consistent_list[i][0] - consistent_list[i-1][-1] >=2:
                    start = consistent_list[i-1][0] * slide
                    # Bug fix: extend by the full window width (50000), not one stride
                    # (10000). The previous +10000 truncated 40 kb off the region end
                    # whenever two positive-window runs were 2-5 windows apart, clipping
                    # real intact LTR-RTs (~11% of elements) that the CNN had flagged.
                    end = consistent_list[i-1][-1] * slide + 50000
                    l = chr_id + "\t" + str(start) + "\t" + str(end)
                    res.append(l)
            start = consistent_list[-1][0] * slide
            end = consistent_list[-1][-1] * slide + 50000
            if i == len(consistent_list) - 1:
                end = length
            l = chr_id + "\t" + str(start) + "\t" + str(end)
            res.append(l)
        elif len(consistent_list) == 1:
            start = consistent_list[0][0] * 10000
            end = consistent_list[0][-1] * 10000 + 50000
            l = chr_id + "\t" + str(start) + "\t" + str(end)
            res.append(l)
    if res == []:
        result = []
    else:
        # LTR_FINDER_parallel does its own fixed 5Mb splitting on each candidate region,
        # so the legacy split_ltr_range() large-block chopping is redundant here. It is
        # also skipped because it produces malformed (start>end) coordinates on very
        # large regions (e.g. a high-recall gate flagging ~100% of a chromosome), which
        # would drop those regions. Hand the merged regions straight to the annotator.
        #
        # Merge overlapping/adjacent regions before returning: after each positive-window
        # run is extended by the full 50 kb window, neighbouring runs overlap, and an
        # unmerged mask would extract the same sub-sequence more than once -> duplicate
        # candidates. Merging keeps the union coverage identical while removing the
        # redundancy. Coordinates are also clamped to [0, length].
        ivs = []
        for l in res:
            _c, s, e = l.split("\t")
            s = max(0, int(s)); e = min(length, int(e))
            if e > s:
                ivs.append((s, e))
        ivs.sort()
        merged = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        result = [chr_id + "\t" + str(s) + "\t" + str(e) for s, e in merged]
    return result