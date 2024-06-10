import torch 
import torchvision.datasets as dataset
import torchvision.transforms as transforms
import torch.utils.data as data_utils
import numpy as np
from random import seed,randint
import argparse

# Define parameter
parser = argparse.ArgumentParser()

parser.add_argument('--dataset_file',type=str,required=True,help="supply dataset file's path")
parser.add_argument('--num',type=int,required=True,default=100000,help='the number of samples to trian and test')
parser.add_argument('--output',type=str,required=True,help="output path")
parser.add_argument('--proportion',type=int,default=0.7,help="sequence of length")

arg = parser.parse_args()

# Define model's structural
class DETECT(torch.nn.Module):
    def __init__(self):
        super(DETECT,self).__init__()
        self.conv_unit = torch.nn.Sequential(
            torch.nn.Conv2d(1,32,kernel_size=(5,20),stride=(1,1),bias=True),
            torch.nn.ReLU(),
            torch.nn.AvgPool2d(kernel_size=(1,10),stride=(1,10)),
            torch.nn.Conv2d(32,64,kernel_size=(1,20),stride=(1,1)),
            torch.nn.ReLU(),
            torch.nn.AvgPool2d(kernel_size=(1,15),stride=(1,15)),
            torch.nn.Conv2d(64,128,kernel_size=(1,35),stride=(1,1)),
            torch.nn.ReLU(),
            torch.nn.AvgPool2d(kernel_size=(1,15),stride=(1,15))
        )

        self.fc_unit = torch.nn.Sequential(
            torch.nn.Flatten(),
            # torch.nn.Linear(5376,1000),
            torch.nn.Linear(2432,1000),
            # torch.nn.Linear(1152,1000),
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


if __name__ == '__main__':

    # read parametre
    dataset_file = arg.dataset_file
    num = arg.num
    proportion = arg.proportion
    output = arg.output
    if output[-1] == '/':
        output = output
    else:
        output = output + '/'
    # Assign a model variable
    dec = DETECT()
    dec = dec.cuda()

    # Definde optmizer
    opt = torch.optim.SGD(dec.parameters(),lr=0.005)

    # Define loss function
    loss_cn = torch.nn.CrossEntropyLoss()

    # Load dataset
    if dataset_file[-1] == '/':
        X = np.load(dataset_file + 'dataset.npy').astype(bool)
        Y = np.load(dataset_file + 'dataset_labels.npy').astype(np.int32)
    else:
        X = np.load(dataset_file + '/' + 'dataset.npy').astype(bool)
        Y = np.load(dataset_file + '/' + 'dataset_labels.npy').astype(np.int32)

    number_dataset = X.shape[0]
    if number_dataset <= num:
        Num = number_dataset
    elif number_dataset > num:
        Num = num
    # extract sequences randomly
    # already = []
    # X_new = np.zeros((Num, X.shape[1], X.shape[2]), dtype=np.int8)
    # Y_new = np.zeros((Num, Y.shape[1]), dtype=np.float64)
    # i = 0
    # seed(1)
    # while i < Num:
    #     value = randint(0, X.shape[0]-1)
    #     if value not in already:
    #         X_new[i] = X[value, :, :]
    #         Y_new[i, 0] = Y[value, 0]
    #         Y_new[i, 1] = Y[value, 1]
    #         i += 1
    #         already.append(value)
    # X = X_new
    # Y = Y_new[:, 0]
    # Y[Y > 0] = 1


    # 
    train_num = round(Num * proportion)
    X_train = X[0:train_num,:,:]
    Y_train = Y[0:train_num]
    X_test = X[train_num:,:,:]
    Y_test = Y[train_num:]
    print(X_train.shape)
    print(Y_train.shape)
    print(X_test.shape)
    print(Y_test.shape)
    X = None
    Y = None


    X_train = torch.tensor(X_train,dtype=torch.float32)
    Y_train = torch.tensor(Y_train,dtype=torch.float32)
    X_test = torch.tensor(X_test,dtype=torch.float32)
    Y_test = torch.tensor(Y_test,dtype=torch.float32)



    train_dataset = data_utils.TensorDataset(X_train,Y_train)
    test_dataset = data_utils.TensorDataset(X_test,Y_test)

    #batch_size = 64
    train_loader = data_utils.DataLoader(dataset=train_dataset,batch_size=64,shuffle=True)
    test_loader = data_utils.DataLoader(dataset=test_dataset,batch_size=64,shuffle=True)



    X_train = None
    Y_train = None
    X_test = None
    Y_test = None




    # Training model
    test_accuracy = []
    train_accuracy = []
    train_confusion_matrix = []
    test_confusion_matrix = [] 
    LOSS = []
    for epoch in range(100):
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        for i, (images,labels) in enumerate(train_loader):
            images = torch.unsqueeze(images,dim=1)
            images = images.cuda()
            labels = labels.cuda()
            outputs = dec(images)
            loss = loss_cn(outputs,labels.long())
            _,Pred = outputs.max(1)
            # Accuracy += (Pred == labels).sum().item()
            # loss = - loss
            opt.zero_grad()
            loss.backward()
            opt.step()
            a = Pred.cpu().numpy() == labels.cpu().numpy()
            for index in range(1,len(a),1):
                sample = str(a[index]) + str(labels.cpu().numpy()[index])
                # print(sample)
                
                if sample == 'True0.0':
                    TN += 1
                elif sample == 'True1.0':
                    TP += 1
                elif sample == 'False0.0':
                    FN += 1
                elif sample == 'False1.0':
                    FP += 1
        Accuracy = (TP + TN) / (TP + TN + FP + FN)
        Accuracy = Accuracy * 100
        train_accuracy.append(Accuracy)
        metrics = (TP, TN, FP, FN)
        train_confusion_matrix.append(metrics)
        LOSS.append(loss.item())
        if (epoch + 1) % 100 == 0:
            torch.save(dec,output + 'saved-model-{}-{}.pkl'.format(epoch + 1,Accuracy))
        print("Train: epoch is {},Accuracy is {},loss is {}.".format(epoch+1,Accuracy,loss.item()))
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        
        for i,(images,labels) in enumerate(test_loader):
            images = torch.unsqueeze(images,dim=1)
            images = images.cuda()
            labels = labels.cuda()

            outputs = dec(images)
            # loss_test += loss_cn(outputs,labels.long())

            _,pred = outputs.max(1)
            a = pred.cpu().numpy() == labels.cpu().numpy()
            for index in range(1,len(a),1):
                sample = str(a[index]) + str(labels.cpu().numpy()[index])
                # print(sample)
                if sample == 'True0.0':
                    TN += 1
                elif sample == 'True1.0':
                    TP += 1
                elif sample == 'False0.0':
                    FN += 1
                elif sample == 'False1.0':
                    FP += 1
        Accuracy = (TP + TN) / (TP + TN + FP + FN)
        Accuracy = Accuracy * 100
        metrics = (TP, TN, FP, FN)
        test_confusion_matrix.append(metrics)
        print("Train: epoch is {},Accuracy is {}.".format(epoch+1,Accuracy))
    train_confusion_matrix_array = np.array(train_confusion_matrix)
    np.save(output + '/Metrics_array.npy',train_confusion_matrix_array)
    test_confusion_matrix_array = np.array(test_confusion_matrix)
    np.save(output + '/metrics_array.npy',test_confusion_matrix_array)
    loss_array = np.array(LOSS)
    np.save(output + '/loss.npy',loss_array)