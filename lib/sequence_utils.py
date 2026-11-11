from Bio import SeqIO
# define a dict for genome (keys: chromosome's id, values: chromosome's sequence)
#
# NOTE (known limitation): this loads the ENTIRE genome into memory as Python strings,
# and the CNN gate (cnn_gate.py) holds it for the whole run. Peak memory therefore
# scales with genome size (~1 byte/bp plus Python overhead), so a multi-Gb plant genome
# can exhaust RAM. For such genomes, run the gate per chromosome (split the FASTA and
# run ltr_checker on each) rather than on the whole assembly at once. See README.
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

# transform the CNN window predictions into genome position ranges
def ltr_range(chr_id,list_result,length,slide,win=50000,pad=0):
    # `win` is the CNN detection window width (bp). A positive window at index i
    # spans [i*slide, i*slide + win], so a run of positive windows extends by the
    # full window width past its last start. `win`/`slide` were previously hard-coded
    # to 50000/10000; they are now parameters so the 15/20/25/30 kb retrained gates
    # can be used. `pad` (bp) is a symmetric safety margin added to BOTH sides of
    # every block before clamping/merging: the CNN only guarantees an element midpoint
    # lies in a positive window, not that the whole element (+ TSD/TG..CA flank the
    # tools need) fits the block; pad = max_LTR_separation//2 fully encloses any such
    # element. default pad=0 = legacy pure-window-union behaviour.
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
        end = temp[-1] * slide + win
        l = chr_id + "\t" + str(start) + "\t" + str(end)
        res.append(l)
    elif len(temp) == 0:
        pass  # no positive window on this chromosome -> no regions (res stays [])
    elif len(temp) > 1:
        consistent_list = split_num_l(temp)
        if len(consistent_list) >= 2:
            for i in range(1,len(consistent_list),1):
                if consistent_list[i][0] - consistent_list[i-1][-1] > int(win // slide):
                    start = consistent_list[i-1][0] * slide
                    end = consistent_list[i-1][-1] * slide + win
                    l = chr_id + "\t" + str(start) + "\t" + str(end)
                    res.append(l)
                elif consistent_list[i][0] - consistent_list[i-1][-1] <= int(win // slide) and consistent_list[i][0] - consistent_list[i-1][-1] >=2:
                    start = consistent_list[i-1][0] * slide
                    # Extend by the full window width (win): the run's last positive
                    # window itself spans [last*slide, last*slide+win], so truncating
                    # to +slide would clip real intact LTR-RTs the CNN flagged.
                    end = consistent_list[i-1][-1] * slide + win
                    l = chr_id + "\t" + str(start) + "\t" + str(end)
                    res.append(l)
            # Final run: end at its last window's right edge like every other run.
            # (Previously forced to `length`, stretching the last block unbounded to
            #  the chromosome end -- removed; pad below gives the safety margin.)
            start = consistent_list[-1][0] * slide
            end = consistent_list[-1][-1] * slide + win
            l = chr_id + "\t" + str(start) + "\t" + str(end)
            res.append(l)
        elif len(consistent_list) == 1:
            start = consistent_list[0][0] * slide
            end = consistent_list[0][-1] * slide + win
            l = chr_id + "\t" + str(start) + "\t" + str(end)
            res.append(l)
    if res == []:
        result = []
    else:
        # Merge overlapping/adjacent regions before returning: after each positive-window
        # run is extended by the full 50 kb window, neighbouring runs overlap, and an
        # unmerged mask would extract the same sub-sequence more than once -> duplicate
        # candidates. Merging keeps the union coverage identical while removing the
        # redundancy. Coordinates are also clamped to [0, length].
        ivs = []
        for l in res:
            _c, s, e = l.split("\t")
            s = max(0, int(s) - pad); e = min(length, int(e) + pad)
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