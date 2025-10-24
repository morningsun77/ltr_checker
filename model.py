import torch 
import torch.utils.data as data_utils
import numpy as np
from numpy import array
from tqdm import tqdm

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

# transform sequence to one-hot encoding
def fasta2one_hot(sequence,length):
    langu = ['A', 'C', 'G', 'T', 'N']
    posNucl = 0
    rep2d = np.zeros((1, 5, length), dtype=int)

    for nucl in sequence:
        nucl = nucl.upper()
        if nucl not in langu:
            nucl = 'N'  # Treat any other character as 'N'
        posLang = langu.index(nucl)
        rep2d[0][posLang][posNucl] = 1
        posNucl += 1
    return rep2d

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

# predict regions which have LTR-RTs
def fasta_to_onehot_to_predict(genome_seq,chr_id,split,stride,model):
    start = 0
    end = 50000
    result = []
    num = len(genome_seq) // split
    if num >= 50000:
        for i in tqdm(range(num,len(genome_seq),num)):
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
    elif num < 50000 and len(genome_seq) >= 50000:
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
    elif len(genome_seq) < 50000:
        list = []
        genome_seq = genome_seq + (50000 - len(genome_seq)) * "N"
        seq = fasta2one_hot(genome_seq,50000)
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