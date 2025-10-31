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
                    end = consistent_list[i-1][-1] * slide + 10000
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
        result = split_ltr_range(res)
        length_list =[]
        for txt in result:
            txt = txt.split('\t')
            length_list.append(int(txt[2]) - int(txt[1]))
            max_length = max(length_list)
        while max_length > 100000000:
            print(max_length)
            result = split_ltr_range(result)
            length_list = []
            for txt in result:
                txt = txt.split('\t')
                length_list.append(int(txt[2]) - int(txt[1]))
            max_length = max(length_list)
    return result