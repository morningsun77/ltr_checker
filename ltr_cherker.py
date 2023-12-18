import torch 
import torch.utils.data as data_utils
import numpy as np
from random import seed,randint
import subprocess
from time import sleep
import argparse
from numpy import array
import os
import multiprocessing
from Bio import SeqIO
import time


# Model's structure
class DETECT(torch.nn.Module):
    def __init__(self):
        super(DETECT,self).__init__()
        self.conv_unit = torch.nn.Sequential(
            torch.nn.Conv2d(1,32,kernel_size=(5,20),stride=(1,1),bias=True),
            torch.nn.ReLU(),
            torch.nn.MaxPool2d(kernel_size=(1,10),stride=(1,10)),
            torch.nn.Conv2d(32,64,kernel_size=(1,20),stride=(1,1)),
            torch.nn.ReLU(),
            torch.nn.MaxPool2d(kernel_size=(1,15),stride=(1,15)),
            torch.nn.Conv2d(64,128,kernel_size=(1,35),stride=(1,1)),
            torch.nn.ReLU(),
            torch.nn.MaxPool2d(kernel_size=(1,15),stride=(1,15))
        )

        self.fc_unit = torch.nn.Sequential(
            torch.nn.Flatten(),
            torch.nn.Linear(2432,1000),
            torch.nn.ReLU(),
            torch.nn.Linear(1000,500),
            torch.nn.ReLU(),
            torch.nn.Linear(500,2),
            torch.nn.Sigmoid()
        )
        
    def forward(self,x):
        x = self.conv_unit(x)
        x = self.fc_unit(x)
        return x

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

# determine which sequence of 50000 bp has LTR-RTs
def predict(dataloader,model):
    result = []
    for i, (images,labels) in enumerate(dataloader):
        images = torch.unsqueeze(images,dim=1)
        outputs = model(images)
        _,pred = outputs.max(1)
        result.append(pred)

        new_result = []
    for tensor_reuult in result:
        new_result.append(tensor_reuult.cpu().numpy().tolist())
    new_result = sum(new_result,[])
    return new_result

# transform sequence to one-hot encoding
def fasta2one_hot(sequence,length):
    langu = ['A', 'C', 'G', 'T', 'N']
    posNucl = 0
    rep2d = np.zeros((1, 5, length), dtype=int)

    for nucl in sequence:
        posLang = langu.index(nucl.upper())
        rep2d[0][posLang][posNucl] = 1
        posNucl += 1
    return rep2d

# define a dict for genome (keys: chromosome's id, values: chromosome's sequence)
def genome_dict(file):
    seqfile = [x for x in SeqIO.parse(file, 'fasta')]
    genome = dict()
    for i in seqfile:
        genome[str(i.id)] = str(i.seq)
    return genome

# predict regions which have LTR-RTs
def fasta_to_onehot_to_predict(genome_seq,chr_id,threads,stride,model):
    start = 0
    end = 50000
    result = []
    num = len(genome_seq) // threads
    if num >= 50000:
        for i in range(num,len(genome_seq),num):
            list = []
            while end < i:
                seq = fasta2one_hot(genome_seq[start:end],50000)
                start += stride
                end += stride
                list.append(seq)
            seqence = array(list,dtype=int)
            seqence = np.squeeze(seqence,axis=1)
            seq_label = np.zeros(len(seqence),dtype=np.float64)
            X = torch.tensor(seqence,dtype=torch.float32)
            Y = torch.tensor(seq_label,dtype=torch.float32)
            dataset = data_utils.TensorDataset(X,Y)
            dataloader = data_utils.DataLoader(dataset,batch_size=64,shuffle=False)
            res = predict(dataloader,model)
            result.append(res)
            X = None
            Y = None
            seqence = None
            seq_label = None
            dataset = None
            dataloader = None
            list = None
            res = None
            chr_res = dict()
            chr_res[chr_id] = result 
    elif num < 50000:
        list = []
        while end < len(genome_seq):
            seq = fasta2one_hot(genome_seq[start:end],50000)
            start += stride
            end += stride
            list.append(seq)
        seqence = array(list,dtype=int)
        seqence = np.squeeze(seqence,axis=1)
        seq_label = np.zeros(len(seqence),dtype=np.float64)
        X = torch.tensor(seqence,dtype=torch.float32)
        Y = torch.tensor(seq_label,dtype=torch.float32)
        dataset = data_utils.TensorDataset(X,Y)
        dataloader = data_utils.DataLoader(dataset,batch_size=64,shuffle=False)
        res = predict(dataloader,model)
        result.append(res)
        X = None
        Y = None
        seqence = None
        seq_label = None
        dataset = None
        dataloader = None
        list = None
        res = None
        chr_res = dict()
        chr_res[chr_id] = result 
    return chr_res

# transform predict's results to position information
def ltr_range(chr_id,list_result,length,slide):
    new_result = sum(list_result,[])
    num = 0
    temp = []
    for i in new_result:
        if i == 1:
            temp.append(num)
        num += num
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
    result = []
    for txt in res:
        txt = txt.split('\t')
        End = int(txt[2])
        Start = int(txt[1])
        if End - Start > 100000000:
            distance = End - Start
            num = distance // 10
            for Num in range(0, num*10, num):
                Start = Start + Num
                End = End + Num
                result.append(txt[0] + '\t' + str(Start) + '\t' + str(End))
        elif End - Start < 100000000 and End - Start > 50000000:
            distance = End - Start
            num = distance // 2
            for Num in range(0, num*2, num):
                Start = Start + Num
                End = End + Num
                result.append(txt[0] + '\t' + str(Start) + '\t' + str(End))
        else:
            result.append(txt[0] + '\t' + str(Start) + '\t' + str(End))
    return result


## annote total LTR by LTR_FINDER
def ltr(max_len_threshold,min_len_threshold, outputDir,tg_ca, TSD,dict,seq_id_file):
    if tg_ca:
        finder_filter = '1111'
    else:
        finder_filter = '0000'
    if TSD:
        finder_filter += '1'
    else:
        finder_filter += '0'
    finder_filter += '000000'
    # finder_filter += '111001'
    path = outputDir + '/' + seq_id_file
    f = open(outputDir + '/' + seq_id_file + '.txt','r')
    fasta = open(outputDir + '/' + seq_id_file + '.fasta','w')
    for seq_id in f:
        seq_id = seq_id.split('\n')[0]
        fasta.write('>' + seq_id+'\n')
        fasta.write(dict[seq_id]+'\n')
    fasta.close()
    sleep(0.01)
    if os.path.exists(path) != True:
        os.mkdir(path)
    sleep(0.01)
    print(seq_id_file + " is start")
    print(path)
    num = int(seq_id_file.split('temp')[1])
    sleep(num)
    command = 'ltr_finder -F ' + finder_filter + ' -D ' + str(max_len_threshold) + ' -d ' + str(min_len_threshold) + ' -w2 -C -p 20 -M 0.85 -L 7000 -l 100 ' + outputDir + '/' + seq_id_file + '.fasta'
    output = subprocess.Popen(command,stdout=subprocess.PIPE,cwd=path,shell=True)
    sleep(0.01)
    res = output.stdout
    res = res.read().decode().split('\n')
    print(seq_id_file + " is finish!")
    # bestHits = []
    Res = []
    for i in res:
        if '[' in i:
            Res.append(i)
    if os.path.exists(outputDir + '/' + seq_id_file + '.fasta'):
        os.remove(outputDir + '/' + seq_id_file + '.fasta')
    return Res

parser = argparse.ArgumentParser()

parser.add_argument('--genome',type=str,required=True,help="genome file")
parser.add_argument('--threads',type=int,required=False,default=10,help='number of threads')
parser.add_argument('--output',type=str,required=True,help="output path")
parser.add_argument('--stride',type=int,required=False,default=10000,help="stride")
parser.add_argument('--max',type=int,required=False,default=15000,help="max separation distance of two LTR")
parser.add_argument('--min',type=int,required=False,default=1000,help="min separation distance of two LTR")
parser.add_argument('--tgca',type=str,required=False,default="no",help="whether need TGCA")
parser.add_argument('--tsd',type=str,required=False,default="no",help="whether need TSD")
parser.add_argument('--model',type=str,required=True,help="path of model")
parser.add_argument('--split',type=int,required=False,default=100,help="chromosome segmentation number")

args = parser.parse_args()
if __name__ == '__main__':
    print("#"*41)
    print("#" + " "*39 + "#")
    print("#" + " "*39 + "#")
    print("#" + " "*14 + "ltr_checker" + " "*14 + "#")
    print("#" + " "*39 + "#")
    print("#" + " "*39 + "#")
    print("#"*41)
    start_time = time.time()
    ## argparse
    file = args.genome
    file = os.path.abspath(file)
    threads = args.threads
    stride = args.stride
    max_len_threshold = args.max
    min_len_threshold = args.min
    tg_ca = args.tgca
    split = args.split
    TSD = args.tsd
    total_win_len = 50000
    outputDir = args.output
    outputDir = os.path.abspath(outputDir)
    model_path = args.model
    model_path = os.path.abspath(model_path)
    model = torch.load(model_path)
    model = model.cpu()
    genome = genome_dict(file)

    ## result of model's predction
    print("Begin to dectect LTR-RTs by model.")
    model_res = []
    for k,v in genome.items():
        chr = fasta_to_onehot_to_predict(v,k,split,10000,model)
        model_res.append(chr)
    model = None
    
    ## confirm LTR-RT's region

    region = []
    for Dict in model_res:
        for k,v in Dict.items():
            cc = ltr_range(k,v,len(genome[k]),10000)
            region.append(cc)
    Region = sum(region,[])

    ## define environment variable

    for i in range(1,threads + 1,1):
        locals()["temp" + str(i)] = []

    ## ditribute LTR-RT's region to environment variables 

    if len(Region) < threads:
        threads = len(Region)
        start = 0
        end = len(Region) // threads
        num = (len(Region) // threads)
        for i in range(1,threads + 1,1):
            if i < threads:
                for index in range(start,end,1):
                    locals()["temp" + str(i)].append(Region[index])
                start += num
                end += num
            elif i == threads:
                for index in range(start,len(Region),1):
                    locals()["temp" + str(i)].append(Region[index])

    else:
        start = 0
        end = len(Region) // threads
        for i in range(1,threads + 1,1):
            num = (len(Region) // threads)
            if i < threads:
                for index in range(start,end,1):
                    locals()["temp" + str(i)].append(Region[index])
                start += num
                end += num
            elif i == threads:
                for index in range(start,len(Region)-1,1):
                    locals()["temp" + str(i)].append(Region[index])


    seq_id_file = []
    for i in range(1,threads + 1,1):
        with open(outputDir + '/' +"temp" +str(i) +'.txt','w') as f:
            seq_id_file.append("temp" +str(i))
            for seq in locals()["temp" + str(i)]:
                seq = seq.split('\t')
                Seq = seq[0] + ":" + str(seq[1]) + "-" + str(seq[2])
                Seq = Seq + '\n'
                f.writelines(Seq)

    seq_id_file = []
    BED = []
    for i in range(1,threads + 1,1):
        seq = "temp" + str(i)
        seq_id_file.append(seq)
        f = open(outputDir + '/' +seq + '.txt','r')
        for txt in f:
            txt = txt.split('\n')[0]
            txt = txt.split(':')
            bed = txt[0] + '\t' + txt[1]
            bed = bed.split('-')
            final_txt = bed[0] + '\t' + bed[1]
            BED.append(final_txt)
    
    ltr_seq = dict()
    for i in BED:
        i = i.split("\t")
        ltr_seq[i[0] + ":" + str(i[1]) + "-" + str(i[2])] = genome[i[0]][int(i[1]):int(i[2])]
    
    ##################################################################################
    #  annotate LTR-RTs by LTR_FINDER
    print(" "*10)
    print("Begin to annotate LTR-TRs by LTR_FINDER.")
    pool = multiprocessing.Pool(processes=threads)

    localresults = [pool.apply_async(ltr,args=[max_len_threshold,min_len_threshold, outputDir,tg_ca, TSD,ltr_seq,seq_id]) for seq_id in seq_id_file]
    pool.close()
    pool.join()
    localTables = [p.get() for p in localresults]

    total_txt = sum(localTables,[])
    for i in total_txt:
        line = i.split('\t')
        line = line[1:]
        element_1 = line[0].replace(':','-').split('-')
        element_2 = line[1].split('-')
        element_2 = str(int(element_2[0]) + int(element_1[1])) + '-' + str(int(element_2[1]) + int(element_1[1]))
        element_1 = element_1[0]
        line[0] = element_1
        line[1] = element_2
        line = '\t'.join(str(i) for i in line)
        Index = total_txt.index(i)
        total_txt[Index] = line

    title = 'SeqID\tLocation\tLTR len\tInserted element len\tTSR\tPBS\tPPT\tRT\tIN (core)\tIN (c-term)\tRH\tStrand\tScore\tSharpness\tSimilarity'
    total_txt.insert(0,title)

    f = open(outputDir + '/' + 'ltr_check.out','w')
    for i in total_txt:
        i = i+ '\n'
        f.writelines(i)
    for i in range(1,threads + 1,1):
        if os.path.exists(outputDir + '/temp' + str(i) + '.txt'):
            os.remove(outputDir + '/temp' + str(i) + '.txt')
        if os.path.exists(outputDir + '/temp' + str(i)):
            os.removedirs(outputDir + '/temp' + str(i))
    end_time = time.time()

    print("ltr_checker cost {} s.".format(end_time - start_time))
