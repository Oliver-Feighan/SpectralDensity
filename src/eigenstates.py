import numpy as np

def follow_previous_eigenstate(eigvec, start):
    previous = start

    args = []
    max_overlaps = []

    for f in range(len(eigvec)):
        overlaps = [np.abs(np.dot(eigvec[max(f-1, 0)][:, previous], eigvec[f][:,i])) for i in range(28)]     
        
        arg = np.argmax(overlaps)

        args.append(arg)
        max_overlaps.append(overlaps[arg])
        
        previous = arg
        
    return np.array(args), np.array(max_overlaps)
    

def follow_start_eigenstate(eigvec, start):

    args = []
    max_overlaps = []

    for f in range(len(eigvec)):
        overlaps = [np.abs(np.dot(eigvec[0][:, start], eigvec[f][:,i])) for i in range(28)]     

        arg = np.argmax(overlaps)

        args.append(arg)
        max_overlaps.append(overlaps[arg])
        
    return np.array(args), np.array(max_overlaps)