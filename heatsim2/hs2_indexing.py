# Note C-style indexing, order z,y,x 
# ** NOT CURRENTLY USED **

# This is just for debugging:
xy_scaling_log2=[0,0,1,1,1,1,2,2,2,2,2,2,2,2]; # indexed by z position

# index manipulation for compressed geometry

def compressed_index(uncompressed_index,xy_scaling_log2):
    assert(len(compressed_index)==3)
    
    res=[]
    for axis in range(len(uncompressed_index)):
        if axis > 0:
            res.append(uncompressed_index[axis]//(1 << xy_scaling_log2[uncompressed_index[0]]))
            pass
        else:
            res.append(uncompressed_index[axis])
        pass
    return tuple(res)

def compressed_single_index(compressed_index,xy_scaling_log2,uncompressed_shape):
    # origindex =uncompressed_index[2] + uncompressed_index[1]*uncompressed_shape[2] + uncompressed_index[0]*uncompressed_shape[2]*uncompressed_shape[1]
    # index = compressed_index[2] + compressed_index[1]*compressed_shape[2] + sum from 0 to res2-1 of compressed_shape[2]i*compressed_shape[1]i

    # Calculate sum of RHS above first
    index=0
    for zpos in range(compressed_index[0]):
        index+=(uncompressed_shape[2]//(1 << xy_scaling_log2[zpos]))*(uncompressed_shape[1]//(1<<xy_scaling_log2[zpos]))
        pass
     
    index +=compressed_index[2] + compressed_index[1]*(uncompressed_shape[2]//(1<<xy_scaling_log2[uncompressed_index[0]]))
    return index


def compressed_index_from_single(single_index,xy_scaling_log2,uncompressed_shape):
    # First determine zpos
    zposindex=0;
    for zpos in range(uncompressed_shape[0]):
        newzposindex=zposindex+(uncompressed_shape[2]//(1 << xy_scaling_log2[zpos]))*(uncompressed_shape[1]//(1<<xy_scaling_log2[zpos]))
        if newzposindex >= single_index:
            break
        zposindex=newzposindex
        pass

    indexremainder=single_index-zposindex
    rowsize=(uncompressed_shape[2]//(1<<xy_scaling_log2[uncompressed_index[0]]))
    ypos=indexremainder // rowsize
    xpos=indexremainder % rowsize
    
def uncompressed_index_range(compressed_index,xy_scaling_log2):
    res=[]
    
    for axis in range(len(compressed_index)):
        if axis > 0:
            scalefactor=(1 << xy_scaling_log2[uncompressed_index[0]])
            res.append((compressed_index[axis]*scalefactor,(compressed_index[axis]+1)*scalefactor))
            pass
        else:
            res.append((compressed_index[axis],compressed_index[axis]+1))
            pass
        pass
    pass

