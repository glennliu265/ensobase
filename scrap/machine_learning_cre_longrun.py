#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Long Run Version for many Epochs, based on machine_learning_cre

Try FNN, use CCF as predictors. Target_global mean CRE

inputs        : CCFs (6x)
outputs       : global mean cre
total samples : 23 * 12 = 276 samples :(

or can increase the number of samples globally... (focus on maritime points)

Created on Mon Apr 20 10:58:39 2026

@author: gliu

"""

import sys
import numpy as np
import os
import time
from tqdm import tqdm

import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset,Dataset
import torch.optim as optim

import importlib
import xarray as xr
import copy

import matplotlib.pyplot as plt

#%% Import Modules


amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

pamvpath = "/home/niu4/gliu8/scripts/predict_nasst"
sys.path.append(pamvpath)
import amv_dataloader as dl
import amvmod as am

import predict_amv_params as pparams


#%%
def train_ResNet_regression(model,loss_fn,optimizer,dataloaders,
                 max_epochs,early_stop=False,verbose=True,
                 reduceLR=False,LRpatience=3,checkgpu=True,debug=True):
    """
    Copied above function, but modified for regression problem
    
    inputs:
        model       - Resnet model
        loss_fn     - (torch.nn) loss function
        opt         - tuple of [optimizer_name, learning_rate, weight_decay] for updating the weights
                      currently supports "Adadelta" and "SGD" optimizers
        dataloaders - List of (torch.utils.data.DataLoader) [train, test, val]
            - trainloader - (torch.utils.data.DataLoader) for training dataset
            - testloader  - (torch.utils.data.DataLoader) for testing dataset
            - valloader   - (torch.utils.data.DataLoader) for validation dataset
        max_epochs  - number of training epochs
        early_stop  - BOOL or INT, Stop training after N epochs of increasing validation error
                     (set to False to stop at max epoch, or INT for number of epochs)
        verbose     - set to True to display training messages
        reduceLR    - BOOL, set to true to use LR scheduler
        LRpatience  - INT, patience for LR scheduler
    
    output:
    
    dependencies:
        from torch import nn,optim

    """
    # Check if there is GPU
    if checkgpu:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device('cpu')
    model = model.to(device)
    
    # Get list of params to update
    params_to_update = []
    for name,param in model.named_parameters():
        if param.requires_grad == True:
            params_to_update.append(param)
            # if verbose:
            #     print("Params to learn:")
            #     print("\t",name)
    
    # Set optimizer
    if optimizer[0] == "Adadelta":
        opt = optim.Adadelta(model.parameters(),lr=optimizer[1],weight_decay=optimizer[2])
    elif optimizer[0] == "SGD":
        opt = optim.SGD(model.parameters(),lr=optimizer[1],weight_decay=optimizer[2])
    elif optimizer[0] == 'Adam':
        opt = optim.Adam(model.parameters(),lr=optimizer[1],weight_decay=optimizer[2])
    elif optimizer[0] == 'Adamax':
        opt = optim.Adamax(model.parameters(),lr=optimizer[1],weight_decay=optimizer[2])
    
    # Set up loaders
    mode_names = ["train","test","val"]
    mode_loop  = [(mode_names[i],dataloaders[i]) for i in range(len(dataloaders))]
    val_flag = False
    if len(dataloaders) > 2:
        val_flag = True
    
    # Add Scheduler
    if reduceLR:
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, patience=LRpatience)
    
    # Set early stopping threshold and counter
    if early_stop is False:
        i_thres = max_epochs
    else:
        i_thres = early_stop
    i_incr    = 0 # Number of epochs for which the validation loss increases
    bestloss  = np.inf
    
    # Preallocation (in the future can allocate this to [3 x max_epoch] array)
    losses = {'train': np.full((max_epochs),np.nan), 'test' : np.full((max_epochs),np.nan),'val' : np.full((max_epochs),np.nan)}
    #accs   = {'train': np.full((max_epochs),np.nan), 'test' : np.full((max_epochs),np.nan),'val' : np.full((max_epochs),np.nan)}
    
    # Main Loop
    for epoch in tqdm(range(max_epochs)): # loop by epoch
        for mode,data_loader in mode_loop: # train/test for each epoch
            if mode == 'train':  # Training, update weights
                model.train()
            else: # Testing/Validation, freeze weights
                model.eval()
            
            runningloss = 0
            #runningmean = 0
            #correct     = 0
            #total       = 0
            for i,data in enumerate(data_loader):
                # Get mini batch
                batch_x, batch_y = data
                batch_x = batch_x.to(device)
                batch_y = batch_y.to(device)
                # Set gradients to zero
                opt.zero_grad()
                
                # Forward pass
                pred_y = model(batch_x)
                #_,predicted = torch.max(pred_y.data,1) # Not classification, so dont track acc
                
                # Calculate loss
                loss = loss_fn(pred_y,batch_y)
                
                # Track accuracy
                #total   += batch_y.size(0)
                #correct += (predicted == batch_y[:,0]).sum().item()
                #print("Total is now %.2f, Correct is now %.2f" % (total,correct))
                
                # Update weights
                if mode == 'train':
                    loss.backward() # Backward pass to calculate gradients w.r.t. loss
                    opt.step()      # Update weights using optimizer
                elif (mode == 'val') or (val_flag is False and mode == "test"):  # update scheduler after 1st epoch for validation
                    if reduceLR:
                        scheduler.step(loss)
                
                runningloss += float(loss.item()) # Accumulate Loss
                #runningmean += correct/total
            
            # Compute the Mean Loss Across Mini-batches
            meanloss_batch = runningloss/len(data_loader) 
            #meanacc_batch  = runningmean/len(data_loader) # correct/total
            
            if verbose: # Print progress message
                print('{} Set: Epoch {:02d}. loss: {:3f}.'.format(mode, epoch+1, \
                                                meanloss_batch))#,meanacc_batch*100))
            
            # Save model if this is the best loss
            if (meanloss_batch < bestloss) and ((mode == 'val') or (val_flag is False and mode == "test")):
                bestloss  = meanloss_batch
                bestmodel = copy.deepcopy(model)
                if verbose:
                    print("Best Loss of %f at epoch %i"% (bestloss,epoch+1))
            
            # Save running loss values for the epoch
            losses[mode][epoch] = meanloss_batch
            #accs[mode][epoch]   = meanacc_batch
            
            # Evaluate if early stopping is needed
            if mode == 'val' or (val_flag is False and mode == "test"):
                if epoch == 0: # Save previous loss
                    lossprev = meanloss_batch
                else: # Add to counter if validation loss increases
                    if meanloss_batch > lossprev:
                        i_incr += 1 # Add to counter
                        if verbose:
                            print("Validation loss has increased at epoch %i, count=%i"%(epoch+1,i_incr))
                    else:
                        i_incr = 0 # Zero out counter
                    lossprev = meanloss_batch
                if (epoch != 0) and (i_incr >= i_thres): # Apply Early stopping and exit script
                    print("\tEarly stop at epoch %i "% (epoch+1))
                    # Decompress dicts (At some point, edit this so that you just return the dictionary...)
                    loss_arr = [losses[md] for md in mode_names]
                    train_loss,test_loss,val_loss = loss_arr
                    #acc_arr = [accs[md] for md in mode_names]
                    #train_acc,test_acc,val_acc = acc_arr
                    return bestmodel,train_loss,test_loss,val_loss
            # Clear some memory
            #print("Before clearing in epoch %i mode %s, memory is %i"%(epoch,mode,torch.cuda.memory_allocated(device)))
            del batch_x
            del batch_y
            torch.cuda.empty_cache()
            
            #print("After clearing in epoch %i mode %s, memory is %i"%(epoch,mode,torch.cuda.memory_allocated(device)))
            # <End Train/Test/Val Mode Loop>
        # <End Epoch Loop>
    
    #bestmodel.load_state_dict(best_model_wts)
    # Decompress dicts (At some point, edit this so that you just return the dictionary...)
    loss_arr = [losses[md] for md in mode_names]
    train_loss,test_loss,val_loss = loss_arr
    #acc_arr  = [accs[md] for md in mode_names]
    #train_acc,test_acc,val_acc = acc_arr
    return bestmodel,train_loss,test_loss,val_loss#train_acc,test_acc,val_acc

#%% Experiment and Training Settings

# Set Output Path
expname   = "FNN6_128_Pilot_LongRun"
outpath   = "/home/niu4/gliu8/projects/ccfs/machine_learning/%s/" % expname
proc.makedir(outpath)


eparams     = dict(
    
    # Cross-validation Parameters
    cv_loop=False, # Repeat for cross-validation
    cv_offset=0,   # Set cv option. Default is test size chunk
    
    # Network Hyperparameters
    netname="FNN4_128",       # Key for Architecture Hyperparameters
    opt=['Adamax',1e-1,0],      # [Optimizer Name, Learning Rate, Weight Decay]
    loss_fn=nn.MSELoss(),     # Loss Function
    use_softmax=False,        # Set to true to change final layer to softmax
    reduceLR=False,          # Set to true to use LR scheduler
    LRpatience=False,        # Set patience for LR scheduler
    
    # Test/Train/Validate and Sampling Splits
    nsamples        =300   , # Number of samples of each class to train with (for classification)
    percent_train   =0.80  , # Training Percentage
    percent_val     =0.00  , # Validation Percentage
    shuffle_class   =False , # Set to True to sample DIFFERENT subsets prior to class subsetting
    shuffle_trainsplit=False,  # Set to False to maintain same set for train/test/val split
    
    # Regularization and Training
    early_stop      =50,       # # of epochs were validation loss increases before stopping
    max_epochs      =1000,     # Maximnum number of Epochs to train for
    batch_size      =256,      # Number of batches for training
    unfreeze_all    =True      # Set to true to unfreeze all layers, false to only unfreeze last layer
    
    )

#% Model and Parameter Settings

verbose     = True
checkgpu    = True

FNN6_dict={
    "nlayers"     : 6,
    "nunits"      : [128,128,128,128,128,128],
    "activations" : [nn.ReLU(),nn.ReLU(),nn.ReLU(),nn.ReLU(),nn.ReLU(),nn.ReLU()],
    "dropout"     : 0.5}
netname="FNN6_128"
nn_params=FNN6_dict

#%%


# Set Output Path
expname   = "FNN2_test"
outpath   = "/home/niu4/gliu8/projects/ccfs/machine_learning/%s/" % expname
proc.makedir(outpath)
niter     = 10


eparams     = dict(
    
    # Cross-validation Parameters
    cv_loop=False, # Repeat for cross-validation
    cv_offset=0,   # Set cv option. Default is test size chunk
    
    # Network Hyperparameters
    netname="FNN4_128",       # Key for Architecture Hyperparameters
    opt=['Adamax',1e-1,0],      # [Optimizer Name, Learning Rate, Weight Decay]
    loss_fn=nn.MSELoss(),     # Loss Function
    use_softmax=False,        # Set to true to change final layer to softmax
    reduceLR=False,          # Set to true to use LR scheduler
    LRpatience=False,        # Set patience for LR scheduler
    
    # Test/Train/Validate and Sampling Splits
    nsamples        =300   , # Number of samples of each class to train with (for classification)
    percent_train   =0.80  , # Training Percentage
    percent_val     =0.00  , # Validation Percentage
    shuffle_class   =False , # Set to True to sample DIFFERENT subsets prior to class subsetting
    shuffle_trainsplit=False,  # Set to False to maintain same set for train/test/val split
    
    # Regularization and Training
    early_stop      =1,       # # of epochs were validation loss increases before stopping
    max_epochs      =1,     # Maximnum number of Epochs to train for
    batch_size      =256,      # Number of batches for training
    unfreeze_all    =True      # Set to true to unfreeze all layers, false to only unfreeze last layer
    
    )

#% Model and Parameter Settings

verbose     = True
checkgpu    = True

FNN6_dict={
    "nlayers"     : 6,
    "nunits"      : [128,128,128,128,128,128],
    "activations" : [nn.ReLU(),nn.ReLU(),nn.ReLU(),nn.ReLU(),nn.ReLU(),nn.ReLU()],
    "dropout"     : 0.5}
netname="FNN6_128"
nn_params=FNN6_dict


#%% Load Predictors and Target

# Before starting, limit years to withhold test set
tstart      = '2001-01-01'
tend        = '2015-01-01'#'2015-12-31'
region_crop = [0,360,-60,60]

# Load Predictors
predictor_path  = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5_2001_2024/anom_detrend1/"
predictor_names = ["sst","eis","Tadv","r700","w700","ws10"]
npredictors     = len(predictor_names)
predictors      = []
for cc in tqdm(range(npredictors)):
    predname    = predictor_names[cc]
    ncname      = "%s%s.nc" % (predictor_path,predname)
    ds          = xr.open_dataset(ncname)[predname].load().squeeze()
    ds          = ut.standardize_names(ds)
    predictors.append(ds)

# Load Target
target_path     = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5_2001_2024/anom_detrend1/"
target_name     = "cre"
nctarget        = "%s%s.nc" % (target_path,target_name)
target          = xr.open_dataset(nctarget)[target_name].load()
target          = ut.standardize_names(target)

# Make Land/Sea Mask and apply
li_mask         = proc.make_mask([p.mean('time') for p in predictors])
predictors      = [p*li_mask for p in predictors]
target          = target * li_mask

# Crop to region if needed
if region_crop is not None:
    predictors =[proc.sel_region_xr(ds,region_crop) for ds in predictors]
    target     = proc.sel_region_xr(target,region_crop)
    li_mask    = proc.sel_region_xr(li_mask,region_crop)


# Withold test set restrict to just training + validation
predictors      = [p.sel(time=slice(tstart,tend)) for p in predictors]
target          = target.sel(time=slice(tstart,tend))

# ======================================================
#%% Approach (1): Train a model using all Non-NaN points
# ======================================================

#  Part 1.1 Reshape Data

# Get Indices for non NaN points
ntime,nlat,nlon = predictors[0].shape
limask_reshape  = li_mask.data.reshape(nlat*nlon)[None,:]
nan_dict        = proc.find_nan(limask_reshape,0,return_dict=True)
okid            = nan_dict['ok_indices']
nsamples        = ntime * len(nan_dict['cleaned_data'].squeeze())

# Get non-NaN Points for Predictor
def reshape_get_non_nan(predictor,ok_indices):
    predictor       = predictor.transpose('time','lat','lon')
    ntime,nlat,nlon = predictor.shape
    predictor       = predictor.data.reshape(ntime,nlat*nlon)
    predictor_ok    = predictor[:,ok_indices]
    
    print("NaNs Found: " + str(np.any(np.isnan(predictor_ok))))
    return predictor_ok
predictors_reshape  = np.array([reshape_get_non_nan(p,okid) for p in predictors]) # [Predictor,Time,Space]
predictors_reshape  = predictors_reshape.reshape(npredictors,nsamples) # [Predictor x Sample]

# Get non-NaN Points for Target
target_ok           = reshape_get_non_nan(target,okid) # [Time x space]
target_reshape      = target_ok.reshape(nsamples) # [Sample]

# This leaves us about 7,364,340 samples... (for training and validation)

# Also Standardize predictors (along sample dimension)
predictors_reshape  = predictors_reshape / predictors_reshape.std(1)[:,None]

#%% Part 1.2: Train/Validate/Test Split and place into data loaders


inputsize   = npredictors
outsize     = 1 # Regression Problem

Xin = predictors_reshape.T # [samples x predictor]
yin = target_reshape[:,None] # [samples x 1]

# Apply Train Test Split (alreadt converts to Tensor?)
debug                    = True
X_subsets,y_subsets      = am.train_test_split(Xin,yin,eparams['percent_train'],
                                               percent_val=eparams['percent_val'],
                                               debug=debug,offset=eparams['cv_offset'])

# Convert to Tensors
X_subsets = [torch.from_numpy(X).type(dtype=torch.float32) for X in X_subsets]
y_subsets = [torch.from_numpy(y).type(dtype=torch.float32) for y in y_subsets]

# # Put into pytorch dataloaders
data_loaders = [DataLoader(TensorDataset(X_subsets[iset],y_subsets[iset]),shuffle=True,
                           batch_size=eparams['batch_size']) for iset in range(len(X_subsets))]
if len(data_loaders) == 2:
    if debug:
        print("There is no validation portion. Validation perc is set to %.2f" % (eparams['percent_val']))
    train_loader,test_loader = data_loaders
else:
    train_loader,test_loader,val_loader = data_loaders

for modeliter in range(niter):
    
    # Build the Model
    layers    = am.build_FNN_simple(inputsize,outsize,nn_params['nlayers'],nn_params['nunits'],nn_params['activations'],
                              dropout=nn_params['dropout'],use_softmax=eparams['use_softmax'])
    pmodel    = nn.Sequential(*layers)




    # Train/Validate Model
    model,trainloss,testloss,valloss = train_ResNet_regression(pmodel,eparams['loss_fn'],eparams['opt'],
                                                                               data_loaders,
                                                                               eparams['max_epochs'],early_stop=eparams['early_stop'],
                                                                               verbose=verbose,reduceLR=eparams['reduceLR'],
                                                                               LRpatience=eparams['LRpatience'],checkgpu=checkgpu,debug=debug)
    
    # ------------------------------
    # 13. Save the model and metrics
    # ------------------------------
    modout = "%s/model%04i.pt" %(outpath,modeliter)
    torch.save(model.state_dict(),modout)
    
    # Save Metrics
    savename = "%s/model%04i_metrics.npz" %(outpath,modeliter)
    np.savez(savename,**{
              'train_loss'     : trainloss,
              'test_loss'      : testloss,
              }
              )
    
    # Clear some memory
    del model
    torch.cuda.empty_cache()  # Save some memory
    
    print("Completed Training for model %04i" % modeliter)
    print("============================================")
    print("\n")








