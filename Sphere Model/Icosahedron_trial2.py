# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 17:03:51 2017

@author: tigan_5ytncvu
"""
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import numpy as np

t = (1.0 + (5.0)**(0.5) )/ 2.0

icosahedron_vertices = [[-1,  t,  0],
                            [ 1,  t,  0],
                            [-1, -t,  0],
                            [1, -t,  0],
                            [ 0, -1,  t],
                            [0,  1,  t],
                            [ 0, -1, -t],
                            [ 0,  1, -t],
                            [ t,  0, -1],
                            [t,  0,  1],
                            [-t,  0, -1],
                            [-t,  0,  1]]

faces = np.array(
            [[0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1]])



########################################################################

#Remove duplicate vertexes and faces. (I don't know if we really have to do this)

#Find a way to find the pentagons (use the simplices of the Convex Hull maybe?)

#Define Horizontal and vertical connections on the sphere using the Z coordinates
#I was thinking to use some sort of reference lines and if a reference line goes through 

#Define full graph

########################################################################

def normalize_v3(arr):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''

    lens = arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 
    lens = np.sqrt(lens)
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens        
    return arr
    
def normalize(arr):
    ''' Normalize components '''
 
    lens = arr[0]**2 + arr[1]**2 + arr[2]**2 
    lens = np.sqrt(lens)
    arr[0] /= lens
    arr[1] /= lens
    arr[2] /= lens          
    return arr


def splitter(vertices, faces):

    new_faces = []
    vert_list = []
    face_counter = 0
    
    for face in faces: 
        face_counter += 1
        #print("facecount =", face_counter)
        v0 = face[0]
        v1 = face[1]
        v2 = face[2]
        
        a = ( v0+v2 ) * 0.5
        b = ( v0+v1 ) * 0.5
        c = ( v1+v2 ) * 0.5  
        
        normalize( a )
        normalize( b )
        normalize( c )
        """
        if np.any(v0== i for i in vert_list):
            v0marker = 1    
            #print(i, v0)
        else:
            v0marker = 0
            vert_list.append(v0)
        
        if np.any(v1== i for i in vert_list):
            v1marker = 1    
        else:
            v1marker = 0
            vert_list.append(v1)
            
        if np.any(v2== i for i in vert_list):
            v2marker = 1    
        else:
            v2marker = 0
            vert_list.append(v2)
            
        if np.any(a== i for i in vert_list):
            amarker = 1    
        else:
            amarker = 0
            vert_list.append(a)
        
        if np.any(b== i for i in vert_list):
            bmarker = 1    
        else:
            bmarker = 0
            vert_list.append(b)
            
        if np.any(c== i for i in vert_list):
            cmarker = 1    
        else:
            cmarker = 0
            vert_list.append(c)
        """
        vert_list.append(v0)   
        vert_list.append(v1)
        vert_list.append(v2)
        vert_list.append(a)
        vert_list.append(b)
        vert_list.append(c)
        #if v0marker + bmarker + amarker == 0:
        new_faces.append([v0, b, a])
        #if bmarker + v1marker + cmarker == 0:
        new_faces.append([b,v1,c])
        #if amarker + bmarker + cmarker == 0
        new_faces.append([a,b,c])
        new_faces.append([a,c,v2])

    new_faces = np.asarray(new_faces)
    vert_list = np.asarray(vert_list)

    return  vert_list, new_faces


def create_unit_sphere_vert( recursion_level=2 ):
    vertex_array, index_array = icosahedron_vertices, faces
    for i in range( recursion_level - 1 ):
        vertex_array, index_array = splitter(vertex_array, index_array)
    return vertex_array, index_array


def sph_polar_convert(xyz):
    """converts xyz array into r,theta,phi array"""
    ptsnew = np.zeros(xyz.shape)
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,0] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    ptsnew[:,2] = np.arctan2(xyz[:,1], xyz[:,0])
    return ptsnew  


def find_adjacent_faces(face,faces):
    #check that you first erase face from faces
    
    adjacent_index=[]
    shape_faces=faces.shape
    nvertices=shape_faces[0]*shape_faces[1]
    verticesflat=faces.reshape((nvertices,3))  #flattens the faces array in an array of vertices
    vertices_index=range(nvertices)     #creates a list of indices
    
    
   
    bins=[ [vertices_index,verticesflat,face] ]  # a variable which is a list where each element is a list as follows [list of indices of vertices, flattened list of vertices, face to stuudy]]
    
    while len(bins)>0:   #repeats binary tree until bins in empty

        rnd=np.random.randint(0,len(bins)) 
        
        new_binary=binary_search(bins[rnd])   #uses binary search on a random element in bins
  
        if new_binary[0]!=True :  #when the binary creates two lists, splitting the input in two
            bins=bins+new_binary
        elif new_binary[0]==True:#when the binary search arrives at only one element and this element is a repeated vertex
            if np.array_equal(verticesflat[new_binary[1]],face[0]) :
                adjacent_index.append(int(new_binary[1]/3.)) #divides by three because it finds the index of faces not verticesflat
            elif np.array_equal(verticesflat[new_binary[1]],face[1]):
                adjacent_index.append(int(new_binary[1]/3.))
            elif np.array_equal(verticesflat[new_binary[1]],face[2]):
                adjacent_index.append(int(new_binary[1]/3.))
        
        del bins[rnd]  #erases the element studied from the list
        
    
    
    return adjacent_index  #gives back the index of the faces list whose faces have a common vertex
 
    
def binary_search(v):                               
    

    face=v[2]
    nvertices=len(v[0])
    
    if nvertices>=2: #if the indices are more than two then it splits the original list into two
        one_branch=v[1][0:int(nvertices/2.)]
        two_branch=v[1][int(nvertices/2.) : ]
        

        index_one_branch=v[0][0:int(nvertices/2.)]
        index_two_branch=v[0][int(nvertices/2.):]
        
      
        uniqueone_branch=(np.unique(one_branch.view([('',one_branch.dtype)]*3)))
        uniqueone_branch=uniqueone_branch.view(one_branch.dtype).reshape((uniqueone_branch.shape[0], 3))
        size_uniqueone_branch=int(np.size(uniqueone_branch)/3.)
        
        
        uniqueone_branch_face=np.vstack((face,one_branch))
        uniqueone_branch_face=(np.unique(one_branch.view([('',one_branch.dtype)]*3)))
        uniqueone_branch_face=uniqueone_branch_face.view(one_branch.dtype).reshape((uniqueone_branch_face.shape[0], 3))
        size_uniqueone_branch_face=int(np.size(uniqueone_branch_face)/3.)
        
        uniquetwo_branch=(np.unique(two_branch.view([('',two_branch.dtype)]*3)))
        uniquetwo_branch=uniquetwo_branch.view(two_branch.dtype).reshape((uniquetwo_branch.shape[0], 3))
        size_uniquetwo_branch=int(np.size(uniquetwo_branch)/3.)
        
        
        uniquetwo_branch_face=np.vstack((face,two_branch))
        uniquetwo_branch_face=(np.unique(two_branch.view([('',two_branch.dtype)]*3)))
        uniquetwo_branch_face=uniquetwo_branch_face.view(two_branch.dtype).reshape((uniquetwo_branch_face.shape[0], 3))
        size_uniquetwo_branch_face=int(np.size(uniquetwo_branch_face)/3.)
        
        
        if size_uniqueone_branch_face- size_uniqueone_branch <3  and size_uniquetwo_branch_face-size_uniquetwo_branch <3: #if both of them contain an element which is also in face
            
            return  [[index_one_branch, one_branch ,face],[index_two_branch, two_branch ,face] ]
        
        elif size_uniqueone_branch_face- size_uniqueone_branch <3: #if only branch one contains an element which is also in face
            
            return [ [index_one_branch, one_branch ,face]  ]
        
        elif size_uniquetwo_branch_face- size_uniquetwo_branch<3:
            
        
            return [ [index_two_branch, two_branch ,face] ]
    else: #if the input is only one vertex the binary search ends
        return [True,v[0][0]]
        
   

    
    
#def unique_rows(a):
#    a = np.ascontiguousarray(a)
#    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
#    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))    
#######################################
####MAKING THE ICOSOHEDRON############

icosahedron_vertices = normalize_v3(np.array(icosahedron_vertices))
icosahedron_vertices.tolist()
face_cache = []
for face in faces:
    new_face_vert = [icosahedron_vertices[face[0]], icosahedron_vertices[face[1]], icosahedron_vertices[face[2]]]
    face_cache.append(new_face_vert)
faces = np.asarray(face_cache )

#######################################
###### RECURSIVELY MAKING THE ICOSPHERE###

sph_vert, sph_tri = create_unit_sphere_vert(4)
x, y, z = sph_vert[:,0], sph_vert[:,1], sph_vert[:,2] #Arrays corresponting the the cartesian coords
ch = ConvexHull(sph_vert)


"""
hull_indices = np.unique(ch.simplices.flat)
hull_pts = sph_vert[hull_indices, :]
print("len hull points", len(hull_pts))
print(len(ch.simplices))

z_pts = hull_pts[:,2]

zmax = np.max(z_pts)
neighmax = np.max(ch.neighbors)
print("z max", zmax)
print("neigh max", neighmax)
"""
# ch.simplices contains the information of the vertices of each of the faces     
#      
#     3   
#    /\        
#   /  \       
#  /____\     
# 1      2       


vertex1, vertex2, vertex3 = ch.points[ch.simplices[:,0]],ch.points[ch.simplices[:,1]], ch.points[ch.simplices[:,2]]# hull_pts[ch.simplices[:,1]], hull_pts[ch.simplices[:,2]]

#CHANGE TO SPHERICAL POLARS
sph_v_1 = sph_polar_convert(vertex1)
sph_v_2 = sph_polar_convert(vertex2)
sph_v_3 = sph_polar_convert(vertex3)







colours = sph_v_3[:,2] #Making the colors of the faces correspond to the value of the angle

print("sph v 1")
print(sph_v_1)


fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_trisurf(x,y,z, triangles=ch.simplices, cmap=plt.cm.Greys_r)
surf.set_array(colours)


#ax.scatter(hull_pts[:,0], hull_pts[:,1], hull_pts[:,2])
"""
print("ch simplices")
print(ch.simplices.flatten())
print("ch neighbours")
print(ch.neighbors)
print("chvertices")
print(ch.vertices)
"""
plt.show()

