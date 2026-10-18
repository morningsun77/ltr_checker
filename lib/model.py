import torch
import torch.utils.data as data_utils
import torch.nn.functional as F
import numpy as np
from numpy import array
from tqdm import tqdm

# ASCII base -> code (A,C,G,T,N -> 0..4 ; any other char -> 4/N). Upper- and
# lower-case both map. Used by the vectorised one-hot encoder below.
_BASE2CODE = np.full(256, 4, dtype=np.int64)
for _ch, _code in zip("ACGTN", range(5)):
    _BASE2CODE[ord(_ch)] = _code
    _BASE2CODE[ord(_ch.lower())] = _code

def _encode(sequence):
    # non-ascii bytes -> '?' (63) -> code 4/N, matching "any other char -> N"
    return _BASE2CODE[np.frombuffer(sequence.encode("ascii", "replace"), dtype=np.uint8)]

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
    # Vectorised: (1, 5, length) one-hot with A,C,G,T,N->0..4 and any other base->N.
    # Semantically identical to the original per-character loop, ~100x faster.
    codes = _encode(sequence)
    n = min(len(codes), length)
    rep2d = np.zeros((1, 5, length), dtype=int)
    if n:
        rep2d[0, codes[:n], np.arange(n)] = 1
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
def fasta_to_onehot_to_predict(genome_seq, chr_id, split, stride, model, win=50000, batch=128):
    """Slide a `win`-bp window at `stride` across the sequence, run the CNN on every
    window, and return {chr_id: [preds]} where preds is a flat 0/1 list (1 = window
    predicted to contain an LTR-RT).

    The chromosome is one-hot encoded once (vectorised) and inference runs on the
    model's own device, so this uses the GPU whenever the model lives on cuda. This
    replaces the old per-character Python one-hot loop (~28 min on maize chr1 -> ~10s).

    `split` is kept for signature compatibility and is no longer used: windows are
    enumerated directly as range(0, L-win+1, stride), i.e. every full window including
    the tail one (matches cnn_predict.py; the old chunked loop dropped a tail window).
    """
    device = next(model.parameters()).device
    L = len(genome_seq)
    preds = []
    with torch.no_grad():
        if L >= win:
            codes = _encode(genome_seq)
            starts = list(range(0, L - win + 1, stride))  # full windows only (old code: end<len)
            for b0 in tqdm(range(0, len(starts), batch)):
                idx = starts[b0:b0 + batch]
                mat = np.stack([codes[s:s + win] for s in idx])            # (B, win) int
                t = torch.from_numpy(mat).to(device)
                oh = F.one_hot(t, num_classes=5).permute(0, 2, 1).float()  # (B, 5, win)
                oh = oh.unsqueeze(1)                                       # (B, 1, 5, win)
                out = model(oh)
                preds.extend(out.max(1)[1].cpu().numpy().tolist())
        else:
            # pad the short sequence up to one window, exactly like the original branch
            codes = _encode(genome_seq + "N" * (win - L))
            t = torch.from_numpy(codes[None, :win]).to(device)
            oh = F.one_hot(t, num_classes=5).permute(0, 2, 1).float().unsqueeze(1)
            out = model(oh)
            preds.extend(out.max(1)[1].cpu().numpy().tolist())
    return {chr_id: [preds]}