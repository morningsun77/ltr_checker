import torch


def _conv1d(w, k, s):
    """Width out of a 1-D conv/pool: floor((w-k)/s)+1 (PyTorch default, no padding)."""
    return (w - k) // s + 1

def conv_out_width(L):
    """Final conv-stack width for input length L, mirroring DETECT.conv_unit.
    Used to size the first Linear layer so a model trained at window W loads into
    an inference DETECT(win=W). W=50000->19 (flat 2432); 30000->10; 25000->8;
    20000->6; 15000->4."""
    w = _conv1d(L, 20, 1)    # Conv2d k(5,20) s1
    w = _conv1d(w, 10, 10)   # MaxPool (1,10) s10
    w = _conv1d(w, 20, 1)    # Conv2d k(1,20) s1
    w = _conv1d(w, 15, 15)   # MaxPool (1,15) s15
    w = _conv1d(w, 35, 1)    # Conv2d k(1,35) s1
    w = _conv1d(w, 15, 15)   # MaxPool (1,15) s15
    return w

# Model's structure
class DETECT(torch.nn.Module):
    def __init__(self, win=50000):
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

        # First Linear input dim scales with the window width (win=50000 -> 2432,
        # the original bundled-model value). Anything else loads a win-matched
        # retrained checkpoint (retrain/result_win{W}/detect_best_testf1.pth).
        flat = 128 * conv_out_width(win)
        self.fc_unit = torch.nn.Sequential(
            torch.nn.Flatten(),
            torch.nn.Linear(flat,1000),
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
