# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:36:07 2022

@author: Busra Bulut
"""

##Version numero 2
import tifffile
import skimage.measure
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import pandas as pd
from munkres import Munkres, print_matrix
import networkx as nx
from more_itertools import set_partitions
from itertools import combinations
import scipy.ndimage
from matplotlib.backends.backend_pdf import PdfPages
from random import *
#import essainumber2
import sys

plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [6, 4]
#Experiment parameters 
#longueur caractéristique d'une bactérie
#L_c=22 #totalement à l'oeil
#la longueur que gagne une bactérie en moyenne durant un frame
l=2 
#le seuil d'acceptation pour la zone d'int
r_vois=1.5

z_s=3

#seuil pour check si c'est trop petit pour être maman
seuil_too_small_to_be_mother=2.5

#pour afficher tous les élements 
#pd.set_option('display.max_columns',None)
randomCmapVals = np.random.rand(256, 3)
randomCmapVals[0, :] = np.array([1.0, 1.0, 1.0])
randomCmapVals[-1, :] = np.array([0.0, 0.0, 0.0])
randomCmap = matplotlib.colors.ListedColormap(randomCmapVals)






#Select experiment we want to analyse (type : .tif)


#Save all frames of the experiment in a list
def Save_frames(file_name,number , figure=False):
    
    """ ________________Description________________
    
        Take the file name.tif of  the experiment, figure a bool to know if print or not the experiment, and the number of the experiment
        Return a list in which Frames[k]=properties of cells in the frame k-1
    """
    lab = tifffile.imread(file_name)
        
    pp = PdfPages('Exp'+str(number)+'.pdf')   

    Frames=[]
    for t in range(lab.shape[0]): #on parcours tous les frames 
    
        #skimage.measure.regiopropos_table save properties in a table
        regionProps = skimage.measure.regionprops_table(lab[t], properties=("centroid", "orientation", "major_axis_length"))
        
        Frames.append(regionProps)
        
        if figure:
            
            fig=plt.figure()
            plt.imshow(lab[t],cmap='tab20b',interpolation='none')
            plt.title('frame : %s' %t)
            pp.savefig(fig)
            plt.show()
    pp.close() 
    return(Frames)

            
           

    
 





class Cell:
    
    def __init__(self,x_centre,y_centre,orientation,longueur_ellipse, indice ):
        self.x=x_centre
        self.y=y_centre
        self.orientation=orientation
        self.major_axis_length=longueur_ellipse
        #on va utiliser l'étiquetage global pour conserver l'information maman/fille, pour 
        self.label=indice 
        
        
    #Redefinition of print function
    
    def __repr__(self): 
        return "Cell Number :%s \n x:% s y:% s angle : % s, length : %s. " % (self.label,self.x, self.y,self.orientation,self.major_axis_length) 
        


class Scenario:
    def __init__(self,list_indice,list_cell_divises,cost,name,cost_graph,path,framek,framekk):
        
        #name of scenario framek-framek+_NumberofScenario
        self.name=name
        
        #list of cells matching between frame k and k+1
        self.matching=list_indice
        
        #list of cells which divide
        self.cells_div=list_cell_divises
        
        #cost of the matching proposed by this Scenario
        self.cost=cost
        
        #Cost added because of the path choosen, 
        #To apply dynamic programming in a graph of Scenarios
        self.cost_graph=cost_graph
        
        
        self.path=path
        
        #List of cells in frame k, to have acces to cell properties
        self.framek=framek
        
        #List of cells in frame k+1, to have acces to cell properties
        self.framekk=framekk
    def __repr__(self): 
        return "Scenario : %s, the matching is : %s, so divided cells are : %s, cost : %s, le noeud voisin : %s, " % (self.name,self.matching,self.cells_div,self.cost,self.path)




def cells_from_frame(frame):
    """________________Description______________
    
    Take frame=Save_frames[k] 
    Return a list iterated over cells 
    """

    Liste_cells=[]
    
    for i in range (np.size(frame['centroid-0'])):
        
        #- for the label by default
        Liste_cells.append(Cell(frame['centroid-1'][i],frame['centroid-0'][i],frame['orientation'][i],frame['major_axis_length'][i],-1))
        
    return(Liste_cells)




 
#cell is the cell in frame k that we want to compute its neighbors
#list_cell_frame is the list of cells in a frame, it can be the frame k or frame k+1
def Neighbors(mother_cell, list_cells_frame,L_c,r_vois):
    
    '''______________________________Description______________________________
    
        Takes an object Cell : mother_cell, and a list of Cell Objects (it can be in the same frame or not than the mother cell)
        
        Return a dictionnary of neighbors cell and the number of cells in the neighborhood of the mother_cell
    '''
    c=0 #Number of cells in the neighbors
    
    #dict to stock neighboors cells of the cell mother given 
    dict_neighbours={}
    
    for i in range(len(list_cells_frame)):
        
        #distance between barycenters of two cells
        r=math.sqrt((list_cells_frame[i].x - mother_cell.x)**2 + (list_cells_frame[i].y - mother_cell.y)**2)
        
        if( (r/L_c) >r_vois )or ( 0<=abs(math.cos(list_cells_frame[i].orientation - mother_cell.orientation))<0.3):
            c +=0
        else:
            #ça veut dire que c'est un potentiel voisin
            dict_neighbours[str(i)]=list_cells_frame[i]
            c+=1
    return((dict_neighbours,c))





def pointer(Voisins,Frame_cells,Frame_cells_ancien,t_image,lab,i_maman):
    
    """_____________Description_______________
    
    Take a list of Neighbors cells in frame_cells of cell number i_mom in frame_cells_ancien
    Return : we plot the position of neighbors in a frame[t_image] 
    
    """
    
    centroid_0=[]
    centroid_1=[]
    orientation=[]
    major_axis_length=[]
    for v in Voisins:
        if int(v) != i_maman:
            centroid_1.append(Frame_cells[int(v)].x)
            centroid_0.append(Frame_cells[int(v)].y)
            orientation.append(Frame_cells[int(v)].orientation)
            major_axis_length.append(Frame_cells[int(v)].major_axis_length)
    plt.imshow(lab[t_image])
    plt.plot(Frame_cells_ancien[i_maman].x, Frame_cells_ancien[i_maman].y, 'bx')
    xg,yg=barycentre(Voisins)
    plt.plot(xg, yg, 'kx')
    plt.plot(centroid_1, centroid_0, 'rx')
    plt.show()
    
    
    


def cost_without_division(cell_mere, cell_candidat,L_c,lambdaa):        
    """ ___________________________Description_____________________________
    Take two cell objects, the value Lc, average cell length, and lst of coefficients Lambda
    
    Return : The cost to evualate how cell_mere and cell candidat matches without division

    """
    
    terme1=abs(cell_mere.major_axis_length-cell_candidat.major_axis_length)/L_c
    terme2=abs(math.cos(cell_mere.orientation-cell_candidat.orientation))
    terme3=abs(cell_mere.x-cell_candidat.x)/cell_mere.x
    terme4=abs(cell_mere.y-cell_candidat.y)/cell_mere.y
    return(lambdaa[0]*(terme1**2)+lambdaa[1]*((1-terme2)**2)+lambdaa[2]*(terme3**2+terme4**2)) 




def proba_without_division(cell_mere, Dict_Voisins,cle_maman,L_c,lambdaa):
            
    """ ___________________________Description_____________________________
    Take an object cell, a dictionnary of cell objects, an index, Lc, and lambdaa
    
    Return : a line to the table final for the adjacency matrix for indivisible cells
    """
    Dict_proba=Dict_Voisins.copy()
    #ici on sait que ça s'est pas divisé du coup on va garder dans les voisins que ceux qui ont une longueur supérieur  la cellule mère
    
    seuil=1.5
    for v in Dict_Voisins:
        if Dict_Voisins[v].major_axis_length < cell_mere.major_axis_length -seuil:
            del Dict_proba[v]
           
        else:
            Dict_proba[v]=abs(cost_without_division(cell_mere, Dict_Voisins[v],L_c,lambdaa))
    Dict_proba['cle']=cle_maman
    return(Dict_proba)


def is_plus(chaine):
    """ ___________________________Description_____________________________
    Take a string
    Check there is a + in the string
    Return : a bool
    """
    for c in chaine:
        if c == '+':
            return True
    return False

#renvoie un string
def candidat1(chaine):
    """ ___________________________Description_____________________________
    Take a string
    a+b => a 
    Return : a string
    """
    candidat=''
    for c in chaine:
        if c=='+':
            return(candidat)
        else:
            candidat+=c
    return(candidat)

def candidat2(chaine):
    """ ___________________________Description_____________________________
    Take a string
    a+b => b
    Return : a string
    """
    candidat=''
    for i in range(len(chaine)):
        if chaine[i]=='+':
            for j in range(i+1,len(chaine)):
                candidat+=chaine[j]
    return(candidat)

def find_index_cle(liste_columns):
    
    for i in range(len(liste_columns)):
        if liste_columns[i]=='cle':
            return i
    print("pas de colonne clé")
    return -1


        ##generateur
def all_pairs_filter(lst,Frame_cells_1,L_c,r_vois):
    if len(lst) < 2:
        yield []
        return
    if len(lst) % 2 == 1:
        for i in range(len(lst)):
            for result in all_pairs_filter(lst[:i] + lst[i+1:],Frame_cells_1,L_c,r_vois):
                yield result
    else:
        a = lst[0]
        for i in range(1,len(lst)):
            #ici on enlève les partitions qui contient un couple non admissibles
            r=math.sqrt((Frame_cells_1[a].x - Frame_cells_1[lst[i]].x)**2 + (Frame_cells_1[a].y - Frame_cells_1[lst[i]].y)**2)
            if (r/L_c) <=r_vois:
                  pair = (a,lst[i])
                  for rest in all_pairs_filter(lst[1:i]+lst[i+1:],Frame_cells_1,L_c,r_vois):
                      yield [pair] + rest
                      
                      

             
                      
def all_pairs(lst,L_c):
    if len(lst) < 2:
        yield []
        return
    if len(lst) % 2 == 1:
        for i in range(len(lst)):
            for result in all_pairs(lst[:i] + lst[i+1:]):
                yield result
    else:
        a = lst[0]
        for i in range(1,len(lst)):
            #Ùr=math.sqrt((Frame_cells_1[a].x - Frame_cells_1[lst[i]].x)**2 + (Frame_cells_1[a].y - Frame_cells_1[lst[i]].y)**2)
            #if (r/L_c) <=1.5:
                  pair = (a,lst[i])
                  for rest in all_pairs(lst[1:i]+lst[i+1:]):
                      yield [pair] + rest
                             
      
#On travaille entre frame k et frame k+1
#renvoie liste_cell_divisé
#liste_2d_filles c'est une liste de 2d cellules filles sélectionnées, liste d'int
def trier_mother_cells(Frame_cells_0, Frame_cells_1,L_c,r_vois,lambdaa):
    #liste d'objets cellules qui peuvent se diviser
    liste_mother=[]
    #liste d'objets cellules qui peuvent sûrement pas se diviser
    liste_not_mother=[]
    #tableau_sans_div
    tableau_without_div=pd.DataFrame()
   # mm=0
    #i ici c'est la numérotation locale à une image dans l'odre où on a tracké
    if len(Frame_cells_0)==len(Frame_cells_1):
        liste_not_mother=[i for i in range(len(Frame_cells_0))]
        for i in range(len(Frame_cells_0)):
            cell_mere=Frame_cells_0[i]
            V_0,c_0=Neighbors(cell_mere,Frame_cells_0,L_c,r_vois)
            V_1,c_1=Neighbors(cell_mere,Frame_cells_1,L_c,r_vois)
            tableau_without_div=tableau_without_div.append(proba_without_division(cell_mere,V_1,i,L_c,lambdaa),ignore_index=True)
        
        return( liste_mother, liste_not_mother,tableau_without_div)
    for i in range(len(Frame_cells_0)):
        cell_mere=Frame_cells_0[i]
        V_0,c_0=Neighbors(cell_mere,Frame_cells_0,L_c,r_vois)
        V_1,c_1=Neighbors(cell_mere,Frame_cells_1,L_c,r_vois)
        #ici on enlève les 2d cellules sûres qu'on a choisit pour l'affectation
        #for k in liste_2d_filles:
          #  del V_1[str(i)]
          #ici pour l'instant on retire la condition de même nombre de voisins
        if  cell_mere.major_axis_length < (L_c/2)+seuil_too_small_to_be_mother :
        #ça signifie qu'il n'y a pas eu de divisions localement à côté de la bactérie
            liste_not_mother.append(i)
            tableau_without_div=tableau_without_div.append(proba_without_division(cell_mere,V_1,i,L_c,lambdaa),ignore_index=True)
            #mm+=1
        #elif c_0==c_1:
        if c_0==c_1:
                #on a le même nombre de voisins, soit il n'y a pas eu de division 
                #soit il y a eu une division, mais avec un décalage, 
                x_0,y_0=barycentre(V_0)
                x_1,y_1=barycentre(V_1) 
                if abs(x_0-x_1)>r_vois or  abs(y_0-y_1)>r_vois :
                    #un grand décalage, on pense qu'il y a eu substitution donc une possibilité de divison
                    #print("The cell %s interest zone barycenter moved" %i)
                    liste_mother.append(i)
                else:
                    liste_not_mother.append(i)
                    tableau_without_div=tableau_without_div.append(proba_without_division(cell_mere,V_1,i,L_c,lambdaa),ignore_index=True)
        else:
            liste_mother.append(i)
    return(liste_mother,liste_not_mother,tableau_without_div)



def barycentre(Dict_voisins):
    x_g=0
    y_g=0
    for v in Dict_voisins:
        x_g+=Dict_voisins[v].x
        y_g+=Dict_voisins[v].y
    return((x_g/len(Dict_voisins),y_g/len(Dict_voisins)))

def find_daughter_cells(Frame_cells_1,L_c):
    liste_daughter=[]
    for j in range(len(Frame_cells_1)):
        #small enough to be daughter
        if Frame_cells_1[j].major_axis_length < L_c - eps:
            liste_daughter.append(j)
    return(liste_daughter)

def len_2(p):
    for i in range(len(p)):
        if (len(p[i])!=2):
            return(False)
    return(True)

def len_2et_1(p):
    c1=0
    c2=0
    for i in range(len(p)):
        if (len(p[i])==2):
            c2+=1
        elif len(p[i])==1:
            c1+=1
    if c2==len(p)-1 and c1==1:
        return True
    return False
        
    


def get_partitions(lst):
    
    if len(lst) % 2 == 0:
        return list(
            map(tuple, [                
                p
                for p in set_partitions(lst, k=len(lst)//2)
                if len_2(p)==True]))
    else:
        return list(
             map(tuple, [
                 p
                 for p in set_partitions(lst, k=len(lst)//2+1)
                 if len_2et_1(p)==True]))
    
    
    

###On lit d'abord les cellules de la frame k+1

    
def filtrer_partitions(list_partitions, Frame_cells_1,L_c,r_vois):
    list_final=[]
    for i in range(len(list_partitions)):
        cpt=0
        p=list_partitions[i]
        for j in range(len(p)):
            if (len(p[j])==2):
                #print("p[j]" %p[j])
                #print(p[j])
                r=math.sqrt((Frame_cells_1[p[j][0]].x - Frame_cells_1[p[j][1]].x)**2 + (Frame_cells_1[p[j][0]].y - Frame_cells_1[p[j][1]].y)**2)
                if (r/L_c) <=r_vois:
                    #si il y a un couple qui pose problème, on enlève toute la partition
                   cpt+=1
            else:
                cpt+=1
        if cpt==len(p):
            list_final.append(p)
    return(list_final)



def colonne_cle(nom_data_frame):
    #print(nom_data_frame)
    cols=nom_data_frame.columns.tolist()
    #print(cols)
    index_cle=find_index_cle(cols)
    #print(cols)
    cols=[cols[index_cle]]+cols[:index_cle]+cols[index_cle+1:]
   # print(cols)
    nom_data_frame=nom_data_frame[cols]
    return(nom_data_frame)

def Nan_to_X(dataframe,x):
    for i in range(len(dataframe.index)):
        for j in range(len(dataframe.columns)):
            if pd.isna(dataframe.iloc[i,j]):
                dataframe.iloc[i,j]=x
    return(dataframe)


def Df_to_List(dataframe):
    Liste=[]
    M=dataframe.values
    for m in M:
        Liste.append(list(m))
    return(Liste)

#lst est une liste de liste, en autre une matrice sous forme d'une liste
def Become_Square(Lst):
    n=len(Lst)
    m=len(Lst[0])
    if n==m:
        #print("Square Matrix")
        return(Lst)
    elif n-m>0:
        #print("Row > columns")
        return(Lst)
    else:
        for k in range(m-n):
            Lst.append(m*[0])
    return(Lst)

def Hongrois(Liste_carrée,afficher,true_n):
    m = Munkres()
    indexes = m.compute(Liste_carrée)
    if afficher:
        #print("les indices de l'algorithme hongrois")
        #print(indexes)
        print ('cost=', sum([Liste_carrée[i[0]][i[1]] for i in indexes[:true_n]]))
    p=sum([Liste_carrée[i[0]][i[1]] for i in indexes[:true_n]])
    return(indexes,p)

def Print_Result_Hongrois(lst_ind,tableau_without_div,afficher):
    lst_matching=[]
    cell_prises=[]
    nmbr_ligne=len(tableau_without_div.index)
    for m in lst_ind:
        if(m[0]<nmbr_ligne): #si c'est pas une ligne artificielle
        #+1 car les indices snt calculés dans le tableau ans la colonne clé
            lst_matching.append([str(tableau_without_div['cle'][m[0]]),tableau_without_div.columns[m[1]+1]])
            cell_prises.append(tableau_without_div.columns[m[1]+1])
            if afficher:
                print("Cell Number :%s in frame %s is associated to Cell Number : %s in frame %s" %(tableau_without_div['cle'][m[0]],k,tableau_without_div.columns[m[1]+1],k+1) )
    return (lst_matching,cell_prises)

#liste ind est une liste de string
def elimine_elts_str(liste_Frame,liste_ind):
    new_liste=[]
    for i in range(len(liste_Frame)):
        if (str(i) not in liste_ind):
            new_liste.append(liste_Frame[i])
    return(new_liste)

#liste ind est une liste de int
def elimine_elts_int(liste_Frame,liste_ind):
    new_liste=[]
    for i in liste_Frame:
        if (i not in liste_ind):
            new_liste.append(i)
    return(new_liste)

#il faut que les d premières cases soient les d divisions
def proba_with_division(cell_mere, liste_cases,cle_mere,d,Dict_from_cells_1,L_c,lambdaa,gammaa):
    key_list = []
    for i in range(len(liste_cases)):
        if i <d:
        #if len(liste_cases[i])==2:
            #print(str(liste_cases[i][0])+'+'+str(liste_cases[i][1]))
            key_list.append(str(liste_cases[i][0])+'+'+str(liste_cases[i][1]))
        else:
            key_list.append(str(liste_cases[i]))
    #print("voici le key list ")
    #print(key_list)
    value_list = len(key_list)*[0]
    Dict_proba=dict(zip(key_list, value_list))
    for case in Dict_proba:
        #print(case)
        #Dict_proba[case]=0
        Dict_proba[case]=abs(cost_with_division(cell_mere, case,Dict_from_cells_1,L_c,lambdaa,gammaa))
        
    Dict_proba['cle']=cle_mere
    #print(Dict_proba)
    return(Dict_proba)


#ici le cell_candidat c'est une clé du dict de proba




def cost_with_division(cell_mere, case_candidate,Dict_from_cells_1,L_c,lambdaa,gammaa):
    #notre critère max va être la longueur de la bactérie
    if is_plus(case_candidate):
        #la case c'est une division
        cell1=Dict_from_cells_1[int(candidat1(case_candidate))]
        cell2=Dict_from_cells_1[int(candidat2(case_candidate))]
        t1=abs(cell1.major_axis_length + cell2.major_axis_length - cell_mere.major_axis_length)/L_c
        #t2=abs(math.cos(cell_mere.orientation-cell1.orientation))
        t2=1
        #•t3=abs(math.cos(cell_mere.orientation-cell2.orientation))
        t3=1
        t4=((cell_mere.x-cell1.x)**2+(cell_mere.y-cell1.y)**2)/(L_c)**2
        t5=((cell_mere.x-cell2.x)**2+(cell_mere.y-cell2.y)**2)/(L_c)**2
        t6=((cell2.x-cell1.x)**2+ (cell2.y-cell1.y)**2)/(L_c)**2
        t7=(((cell2.x+cell1.x)/2 - cell_mere.x)**2+ ((cell2.y+cell1.y)/2 - cell_mere.y)**2)/(L_c**2)
        return((1/cell_mere.major_axis_length)*(gammaa[0]*t1**2+gammaa[1]*(1-t2)**2+gammaa[2]*(1-t3)**2+gammaa[3]*t4+gammaa[4]*t5+gammaa[5]*t6 +gammaa[6]*t7))
    else:
        cell_candidat=Dict_from_cells_1[int(case_candidate)]
        return((1-1/cell_mere.major_axis_length)*cost_without_division(cell_mere, cell_candidat,L_c,lambdaa))
    
    
def smallest_cells(dict_from_frame_cells_1_trie,d,cells_prises,z_s):
    liste_indices_cell_filles=[]
    cpt=0
    
    #on augmente la amrge sur les 2d cellules prises
    if 2*d+z_s>len(dict_from_frame_cells_1_trie):
        for m in dict_from_frame_cells_1_trie:
            if cpt == len(dict_from_frame_cells_1_trie):
                return(liste_indices_cell_filles)
            elif str(m[0]) not in cells_prises:
                liste_indices_cell_filles.append(m[0]) #m[0] car on revoit une liste quand on trie
                cpt +=1
        #print("on a selectionné %s cellules filles" %(2*d+eps))
    else :
        for m in dict_from_frame_cells_1_trie:
            if cpt == 2*d+z_s:
                return(liste_indices_cell_filles)
            elif str(m[0]) not in cells_prises:
                liste_indices_cell_filles.append(m[0]) #m[0] car on revoit une liste quand on trie
                cpt +=1
    return(liste_indices_cell_filles)
    




def division(matching_with_div):
    liste_div=[]
    for p in matching_with_div:
        if is_plus(p[1]):
            liste_div.append(float(p[0]))
    return(liste_div)




#on teste si un générateur est vide ou non
def test_gen(gen):
    for i in gen:
        return False
    return True

 
def Create_scenarios(k,afficher,Frames,lab,z_s,r_vois,lambdaa,gammaa):
    Frame_cells_0=cells_from_frame(Frames[k])
    Frame_cells_1=cells_from_frame(Frames[k+1])
    L_c=(average_length(Frame_cells_0)+average_length(Frame_cells_1))/2
    if afficher :
        plt.imshow(lab[k])
        #plt.plot(Frame_cells_0[mere_a_voir].x, Frame_cells_0[mere_a_voir].y, 'bx')
        plt.show()
#on check l'image
        plt.imshow(lab[k+1])
        #plt.plot(Frame_cells_1[enfant].x, Frame_cells_1[enfant].y, 'bx')
        plt.show()
        
##Ici on transforme la liste Frame_cells_1 en dictinnaire car on veut garder le labelling
    key_list = [i for i in range(len(Frame_cells_1))]
    value_list = Frame_cells_1
    dict_from_frame_cells_1=dict(zip(key_list, value_list))
    dict_from_frame_cells_1_trie=sorted(dict_from_frame_cells_1.items(), key=lambda t:t[1].major_axis_length)

##Ici on transforme la liste Frame_cells_0 en dictinnaire car on veut garder le labelling
    key_list = [i for i in range(len(Frame_cells_0))]
    value_list = Frame_cells_0
    dict_from_frame_cells_0=dict(zip(key_list, value_list))
#le nombre de divisions
    d=len(Frame_cells_1)-len(Frame_cells_0)


#je charge la liste des cellules qui peuvent se diviser, qui peuvent forcément pas se diviser, et le tableau sans div
    liste_mother,liste_not_mother,tableau_without_div=trier_mother_cells(Frame_cells_0,Frame_cells_1,L_c,r_vois,lambdaa)

#Ici on met change l'emplacement de la colonne clé
    if not(tableau_without_div.empty):
      
        tableau_without_div=colonne_cle(tableau_without_div)
#on applique un algorithme pour affecter le meilleur matching des cellules mères que nous sommes sûres qui ne sont aps divisées. 

#ici on selectionne toutes les colonnes sauf la colonne clé
        tableau_without_div_reduit=tableau_without_div.iloc[:,1:]
#le dataframe initial est aussi modifié 

#je modifie les Nan en la valeur x que je veux ici 100
        tableau_without_div_reduit=Nan_to_X(tableau_without_div_reduit,100)
        
#on tarnsforme en liste pour appliquer facilement l'algorithme de python
        list_without_div=Df_to_List(tableau_without_div_reduit)
#on rend carrée la matrice, et on renomme le dataframe, car on veut garder l'initial
     
        without_div_carre=Become_Square(list_without_div)
        
        indices_without,cost_without=Hongrois(without_div_carre,False,len(tableau_without_div.index))
        
        matching_without_div,cells_prises=Print_Result_Hongrois(indices_without,tableau_without_div,False)
    else :
       # il y a aucune cellules dont on est sûres qu'elles ne sont pas divisées, de ce fait, tout est vide
        cells_prises=[]
        matching_without_div=[]
        cost_without=0

#cette fonction renvoie les 2d plus petits elements de dict privée de cells_prises
    liste_indices_cell_filles=smallest_cells(dict_from_frame_cells_1_trie,d,cells_prises,z_s)

    All_candidates=list(dict_from_frame_cells_1)
    resultat_hongrois=[]
    Liste_Scenarios=[]
    if d!=0:
        c_vide=0
        c_combi=0
        numero_scenario=1
        for combi in combinations(liste_indices_cell_filles,2*d):

            c_combi+=1
            liste_partitions_couples_filles_final=all_pairs_filter(combi,Frame_cells_1,L_c,r_vois)
            #attention ici il faut lancer un nvx générateur car sinon on reprend l'itératione n fonction du dernier
            if test_gen(all_pairs_filter(combi,Frame_cells_1,L_c,r_vois)):
              c_vide+=1
                #cela signifie que le generateur est vide

            for partition in liste_partitions_couples_filles_final:
      
                Candidates=elimine_elts_str(All_candidates,cells_prises)
                Candidates=list(partition)+elimine_elts_int(Candidates,combi)

                tableau_with_div=pd.DataFrame()
                for i in liste_mother:
                    cell_mere=dict_from_frame_cells_0[i]
                    tableau_with_div=tableau_with_div.append(proba_with_division(cell_mere, Candidates,i,d,dict_from_frame_cells_1,L_c,lambdaa,gammaa),ignore_index=True)

            #Ici on met change l'emplacement de la colonne clé
           
                tableau_with_div=colonne_cle(tableau_with_div)
    #on applique un algorithme pour affecter le meilleur matching des cellules mères que nous sommes sûres qui ne sont aps divisées. 
    #premiere etape, on remplace les Nan par des infinis
    
                tableau_with_div_reduit=tableau_with_div.iloc[:,1:]
    #le dataframe initial est aussi modifié 
                tableau_with_div_reduit=Nan_to_X(tableau_with_div_reduit,100)
    #on tarnsforme en liste
                list_with_div=Df_to_List(tableau_with_div_reduit)
    #on rend carré la matrice
                
                with_div_carre=Become_Square(list_with_div)
                indices,cost_with=Hongrois(with_div_carre,False,len(tableau_with_div.index))
       
                matching_with_div,n=Print_Result_Hongrois(indices,tableau_with_div,False)
                liste_div=division(matching_with_div)

                resultat_hongrois.append(matching_without_div+matching_with_div)
                
                
                Liste_Scenarios.append(Scenario(matching_without_div+matching_with_div,liste_div,cost_without+cost_with,str(k)+'_'+str(k+1)+'-'+str(numero_scenario),0,0,Frame_cells_0,Frame_cells_1))
                numero_scenario+=1
        if c_vide ==c_combi:
            #ça veut dire que c'était vide pour toutes les combi possibles
            print("Augmentez la valeur de z_small pour le frame %s "%k)
            
            
        
    else:
       
        matching_with_div=[]
        cost_with=0
        liste_div=[]
    
        Liste_Scenarios.append(Scenario(matching_without_div+matching_with_div,liste_div,cost_without+cost_with,str(k)+'_'+str(k+1)+'-'+str(1),0,0,Frame_cells_0,Frame_cells_1))
    return(Liste_Scenarios)

def average_length(Frame_cells_0):
    S=0
    for c in Frame_cells_0:
        S+=c.major_axis_length
    return(S/len(Frame_cells_0))




##il faut que cell_mere soit une str
def find_daug(cell_mere, liste_matching):
    ##we want to find the daughter matched l'indice de la mère
    for couple in liste_matching:
        if couple[0]==cell_mere:
            return(candidat1(couple[1]),candidat2(couple[1]))




def Print_results(lab,dict_layers,Frames,Name):
    pp = PdfPages(Name+'.pdf')
    for t in range(lab.shape[0]-1): #on parcours tous les frames 
        Frame_cells_0=cells_from_frame(Frames[t])
        Frame_cells_1=cells_from_frame(Frames[t+1])
        for sce in dict_layers[(t,t+1)]:
            if sce.cells_div:
                #il y a eu une division
                for cells in sce.cells_div:
                    enfant_1,enfant_2=find_daug(str(cells),sce.matching)
                    fig=plt.figure()
                    #f, axarr = plt.subplots(1,2) 
                    # use the created array to output your multiple images. In this case I have stacked 4 images vertically
                    #axarr[0].imshow(lab[t],cmap='tab20b',interpolation='none')
                    #plt.plot(Frame_cells_0[int(cells)].x, Frame_cells_0[int(cells)].y, 'bx')
                    #plt.title("frame %s => frame %s, Sce : %s "%(t,t+1, sce.name)) 
                    #axarr[1].imshow(lab[t+1],cmap='tab20b',interpolation='none')
                    plt.imshow(lab[t],cmap='tab20b',interpolation='none')
                    plt.plot(Frame_cells_0[int(cells)].x, Frame_cells_0[int(cells)].y, 'bx')
                    plt.title("frame %s, Sce : %s, Sce cost : %s, div number : %s"%(t, sce.name,round(sce.cost,5),len(sce.cells_div)))
                    pp.savefig(fig)
                    plt.show()
                    fig=plt.figure()
                    plt.imshow(lab[t+1],cmap='tab20b',interpolation='none')
                    plt.plot(Frame_cells_1[int(enfant_1)].x, Frame_cells_1[int(enfant_1)].y, 'rx')
                    plt.plot(Frame_cells_1[int(enfant_2)].x, Frame_cells_1[int(enfant_2)].y, 'rx')
                    #print("les couples filles choisies : %s, %s" %( enfant_1,enfant_2))
                    plt.title("frame %s, Sce : %s, cost : %s, div number : %s , index : %s, %s "%(t+1, sce.name,round(sce.cost,5),len(sce.cells_div), enfant_1,enfant_2))
                    pp.savefig(fig)
                    plt.show()
                
            else:
                    """
                    fig=plt.figure()
                    f, axarr = plt.subplots(1,2) 
                    # use the created array to output your multiple images. In this case I have stacked 4 images vertically
                    axarr[0].imshow(lab[t],cmap='tab20b',interpolation='none')
                    plt.title("frame %s => frame %s, Sce : %s "%(t,t+1, sce.name)) 
                 
                    axarr[1].imshow(lab[t+1],cmap='tab20b',interpolation='none')
                    pp.savefig(fig)
                    plt.show()
                    """
                    fig=plt.figure()
                    plt.imshow(lab[t],cmap='tab20b',interpolation='none')
                    plt.title("frame %s,  Sce : %s, cost:%s "%(t, sce.name,round(sce.cost,5)))
                    pp.savefig(fig)
                    plt.show()
                    fig=plt.figure()
                    plt.imshow(lab[t+1],cmap='tab20b',interpolation='none')
                    plt.title("frame %s, Sce : %s, cost : %s "%(t+1, sce.name,round(sce.cost,5)))
                    pp.savefig(fig)
                    plt.show()
    pp.close()
    

def Print_results_best(lab,dict_layers,c,name,Frames,file_name):
    pp = PdfPages('Result_best_exp30_z3'+name+'.pdf')
    for t in range(lab.shape[0]-1): #on parcours tous les frames 
        Frame_cells_0=cells_from_frame(Frames[t])
        Frame_cells_1=cells_from_frame(Frames[t+1])
        for sce in dict_layers[(t,t+1)]:
            if sce.name in c:
                if sce.cells_div:
                    #il y a eu une division
                    for cells in sce.cells_div:
                        enfant_1,enfant_2=find_daug(str(cells),sce.matching)
                        fig=plt.figure()
                        #f, axarr = plt.subplots(1,2) 
                        # use the created array to output your multiple images. In this case I have stacked 4 images vertically
                        #axarr[0].imshow(lab[t],cmap='tab20b',interpolation='none')
                        #plt.plot(Frame_cells_0[int(cells)].x, Frame_cells_0[int(cells)].y, 'bx')
                        #plt.title("frame %s => frame %s, Sce : %s "%(t,t+1, sce.name)) 
                        #axarr[1].imshow(lab[t+1],cmap='tab20b',interpolation='none')
                        plt.imshow(lab[t],cmap='tab20b',interpolation='none')
                        plt.plot(Frame_cells_0[int(cells)].x, Frame_cells_0[int(cells)].y, 'bx')
                        plt.title("frame %s, Sce : %s, Sce cost : %s, div number : %s"%(t, sce.name,round(sce.cost,5),len(sce.cells_div)))
                        pp.savefig(fig)
                        #plt.show()
                        fig=plt.figure()
                        plt.imshow(lab[t+1],cmap='tab20b',interpolation='none')
                        plt.plot(Frame_cells_1[int(enfant_1)].x, Frame_cells_1[int(enfant_1)].y, 'rx')
                        plt.plot(Frame_cells_1[int(enfant_2)].x, Frame_cells_1[int(enfant_2)].y, 'rx')
                        #print("les couples filles choisies : %s, %s" %( enfant_1,enfant_2))
                        plt.title("frame %s, Sce : %s, cost : %s, div number : %s , index : %s, %s "%(t+1, sce.name,round(sce.cost,5),len(sce.cells_div), enfant_1,enfant_2))
                        pp.savefig(fig)
                        plt.show()
                    
                else:
                        """
                        fig=plt.figure()
                        f, axarr = plt.subplots(1,2) 
                        # use the created array to output your multiple images. In this case I have stacked 4 images vertically
                        axarr[0].imshow(lab[t],cmap='tab20b',interpolation='none')
                        plt.title("frame %s => frame %s, Sce : %s "%(t,t+1, sce.name)) 
                     
                        axarr[1].imshow(lab[t+1],cmap='tab20b',interpolation='none')
                        pp.savefig(fig)
                        plt.show()
                        """
                        fig=plt.figure()
                        plt.imshow(lab[t],cmap='tab20b',interpolation='none')
                        plt.title("frame %s,  Sce : %s, cost:%s "%(t, sce.name,round(sce.cost,5)))
                        pp.savefig(fig)
                        #plt.show()
                        fig=plt.figure()
                        plt.imshow(lab[t+1],cmap='tab20b',interpolation='none')
                        plt.title("frame %s, Sce : %s, cost : %s "%(t+1, sce.name,round(sce.cost,5)))
                        pp.savefig(fig)
                        #plt.show()
    pp.close()




def find_ind_matching(sce,m):
    for k in range(len(sce.matching)):
        if int(float(sce.matching[k][0]))==m:
            return(k)
    else:
        print("pas de liens ligne 985")
        return -1
        
def cost_edges(sce0,sce1):
    S=0
    #variable qui indique si c'est possible de lier sce0 et sce1
    arr=True
    for k in range(len(sce0.matching)):
        #print("k : %s" %(k))
        if is_plus(sce0.matching[k][1]):
            m1=sce0.matching[k][0]
           # print("m1")
           # print(m1)
            x1=sce0.framek[int(float(m1))].x
           # print("x1 : %s"%x1)
            y1=sce0.framek[int(float(m1))].y
           # print("y1 : %s"%y1)
            #je trouve l'indice dans le frame k+1  de la cellule fille dans frame k associée à la maman avec le label k 
            c1=int(candidat1(sce0.matching[k][1]))
           # print("c1 : %s"%c1)
            c2=int(candidat2(sce0.matching[k][1]))
           # print("c2 : %s"%c2)
            indice1=find_ind_matching(sce1,c1)
           # print("indice1 : %s"%indice1)
            #print("indice1 de l'indice %s" %(k))
            indice2=find_ind_matching(sce1,c2)
           # print("indice2 : %s"%indice2)
            #print("indice2 de l'indice %s" %(k))
            mm2=sce1.matching[indice1][1]
           # print("mm2")
            #print(mm2)
            mmm2=sce1.matching[indice2][1]
           # print("mmm2")
           # print(mmm2)
            if is_plus(str(mm2)) or is_plus(str(mmm2)):
                print("il y a eu une succesion de division entre le %s et %s" %(sce0.name, sce1.name))
                return(math.sqrt(S)/len(sce0.matching),False)
            x2=(sce1.framekk[int(float(mmm2))].x+sce1.framekk[int(float(mm2))].x)/2
            y2=(sce1.framekk[int(float(mmm2))].y+sce1.framekk[int(float(mm2))].y)/2
            S+=(x1-x2)**2+(y1-y2)**2
            
        else:

            m1=sce0.matching[k][0]
            x1=sce0.framek[int(float(m1))].x
            y1=sce0.framek[int(float(m1))].y
            #je trouve l'indice dans le frame k+1  de la cellule fille dans frame k associée à la maman avec le label k 
            indice=find_ind_matching(sce1,int(sce0.matching[k][1]))
 
           
            m2=sce1.matching[indice][1]
            #on part du principe qu'on enchaine pas 2 div
            if is_plus(str(m2)):
                c1=int(candidat1(m2))
                c2=int(candidat2(m2))
                x2=(sce1.framekk[c1].x+sce1.framekk[c2].x)/2
                y2=(sce1.framekk[c1].y+sce1.framekk[c2].y)/2
                S+=(x1-x2)**2+(y1-y2)**2
                
            else:
                x2=sce1.framekk[int(m2)].x
                y2=sce1.framekk[int(m2)].y
                S+=(x1-x2)**2+(y1-y2)**2
            
    return(math.sqrt(S)/len(sce0.matching),arr)


def create_edge_list(L,name_final):
    edge_list=[]
    for i in range(len(L)-1):
        edge_list.append((L[i],L[i+1]))
    edge_list.append((L[0],name_final))
    return(edge_list)           


##dynamique programmation##

def find_number_sce(chaine):
    candidat=''
    for i in range(len(chaine)):
        if chaine[i]=='-':
            for j in range(i+1,len(chaine)):
                candidat+=chaine[j]
    return(candidat)

#G graph, dict_layers 
#là on modifie le dict_layers
def best_path(GG, dict_layers):
    G=GG.copy()
    min_name=''
    succesion_voisins=[]
    dict_succ_voisins={}
    dict_scenaris={}
    dict_best_voisin={}
    #la variable qui m'indique si un chemin de retour existe au moins
    if len(dict_layers[(len(dict_layers)-1,len(dict_layers))])==1:
         for temps in dict_layers:
             if temps[0]<len(dict_layers)-1:
                 #on saute le premier
                 for scenario in dict_layers[temps[0]+1,temps[1]+1]:
                     if scenario.name in G.nodes():
                         min_name=''
                         minn=math.inf
                         exist=False
                         #print("c'est l'etude du sce %s"%scenario.name)
                         cpt_neigh=0
                         for voisin in G.neighbors(scenario.name):
                             #parcequ'on minimise
                             
                            # print("les voisins sont : %s" %voisin)
                           
                             if int(find_layer(voisin)) == temps[0]:  
                              cpt_neigh+=1
                             
                             # print(G.get_edge_data(scenario.name,voisin)['weight']+dict_layers[(temps[0],temps[1])][int(find_number_sce(voisin))-1].cost_graph)
                              if G.get_edge_data(scenario.name,voisin)['weight']+dict_layers[(temps[0],temps[1])][int(find_number_sce(voisin))-1].cost_graph<minn: 
                               exist=True
                            
                               minn=G.get_edge_data(scenario.name,voisin)['weight']+dict_layers[(temps[0],temps[1])][int(find_number_sce(voisin))-1].cost_graph
    
                               min_name=voisin
        
                         if exist or cpt_neigh<2:
                            scenario.cost_graph=minn
                            #scenario.cost_graph=minn+dict_layers[(temps[0],temps[1])][int(find_number_sce(min_name))-1].cost_graph
                            scenario.path=min_name
                            #succesion_voisins.append(min_name)
                            dict_best_voisin[scenario.name]=min_name
                         else:
                             print("le noeud %s est supprimé" %(scenario.name))
                             G.remove_node(scenario.name)
                        
                
         ##backtracking
         #print(dict_best_voisin)
         temp=scenario.name
         while temp!=dict_layers[(0,1)][0].name:

             succesion_voisins.append(dict_best_voisin[temp])
             temp=dict_best_voisin[temp]
             
         return(G, dict_layers,succesion_voisins,scenario.cost_graph)

    else:

        for sce in dict_layers[(len(dict_layers)-1,len(dict_layers))]:
            if sce.name in G.nodes():
                 C=G.copy()
                 dictt=dict_layers.copy()
                 for scee in dict_layers[(len(dict_layers)-1,len(dict_layers))]:
                     if scee.name!=sce.name:
                         if scee.name in C.nodes():
                             C.remove_node(scee.name)
                 dictt[(len(dict_layers)-1,len(dict_layers))]=[]
                 dictt[(len(dict_layers)-1,len(dict_layers))].append(sce)
                 a,b,dict_succ_voisins[sce.name],dict_scenaris[sce.name]=best_path(C,dictt)
             
        return(G, dict_layers,dict_succ_voisins,dict_scenaris)
    

    
  

def find_layer(chaine):
    candidat=''
    for c in chaine:
        if c=='_':
            return(candidat)
        else:
            candidat+=c
    return(candidat)


def find_mother(matching,label_fille):
    for c in matching:
        if int(c[1])==label_fille:
            return(int(float(c[0])))
    print("une cellule fille n'a pas de mother")
    return -1

#ici c c'est déjà avec le noeud final choisi
def Extract_best_path(dict_layers,c,name):
    dict_final={}
    for k in dict_layers: #on parcours tous les frames 
        for sce in dict_layers[k]:
            if sce.name in c or sce.name==name:
                dict_final[k]=sce
    
    return(dict_final)

#entre
def define_color(dict_final, lab,name):
    ##initialization
    nCells_1=lab[0].max()
    ##on choisit les couleurs répartis sur le cercle pour initialiser 
    Hue = np.linspace(0, 1, nCells_1+1)[0:-1]
    key_list = [i for i in range(len(Hue))]
    dict_Hue=dict(zip(key_list, list(Hue)))
    newImageColoured = np.zeros((len(lab),
                                    lab.shape[1],
                                    lab.shape[2],
                                    3), dtype='u1')
    #ici on créé la première image juste en attribuant les couleurs comme ça
    for label in range(1, nCells_1+1):
        h = dict_Hue[label-1]
        rgbForThisLabel = matplotlib.colors.hsv_to_rgb([h, 1, 1])
        newImageColoured[0, lab[0]==label, :] = rgbForThisLabel*255
    ##maintenant on lance la machine pour avoir les couleurs matchant
    #on compte tout's les cellules en + qu'on a et qu'il faut caser sur le cercle
    compteur_d=0
    for k in dict_final:
        #on s'occupe de l'image k+1

        dict_Hue_old=dict_Hue.copy()

        for c in dict_final[k].matching:
            if is_plus(c[1]):
                compteur_d+=1
                dict_Hue[int(candidat1(c[1]))]=dict_Hue_old[int(float(c[0]))]+(1/nCells_1)/5**(compteur_d) #pas 255 parceque format
                dict_Hue[int(candidat2(c[1]))]=dict_Hue_old[int(float(c[0]))]-(1/nCells_1)/5**(compteur_d)
                #dict_Hue[int(candidat1(c[1]))]=dict_Hue_old[int(float(c[0]))] #pas 255 parceque format
                #dict_Hue[int(candidat2(c[1]))]=dict_Hue_old[int(float(c[0]))]
                #dict_Hue[int(candidat1(c[1]))]=dict_Hue_old[int(float(c[0]))]+(1/nCells_1)/3*(compteur_d) #pas 255 parceque format
                #dict_Hue[int(candidat2(c[1]))]=dict_Hue_old[int(float(c[0]))]-(1/nCells_1)/3*(compteur_d)
            else:
                dict_Hue[int(c[1])]=dict_Hue_old[int(float(c[0]))]
        #maintenant qu'on a rempli le tableau HUE, on colorie

        nn=len(dict_final[k].matching)+len(dict_final[k].cells_div)
        for label in range(1, nn+1):
            h = dict_Hue[label-1]
            rgbForThisLabel = matplotlib.colors.hsv_to_rgb([h, 1, 1])
            newImageColoured[k[1], lab[k[1]]==label, :] = rgbForThisLabel*255
    tifffile.imwrite(name+".tif", newImageColoured)
    return(newImageColoured)




def get_orientation(dict_final):
    dict_orientation={}
    for mere in range(len(dict_final[(0,1)].matching)):
        dict_orientation[mere]=[]
        mere_avant=mere
        for k in dict_final:
            #finir=False
            for c in dict_final[k].matching:
                if int(float(c[0]))==mere_avant:
                    if is_plus(c[1]):
                        o1=dict_final[k].framek[int(float(c[0]))].orientation-dict_final[k].framekk[int(float(candidat1(c[1])))].orientation
                        o2=dict_final[k].framek[int(float(c[0]))].orientation-dict_final[k].framekk[int(float(candidat2(c[1])))].orientation
                        if abs(o1)>abs(o2):
                            dict_orientation[mere].append(abs(o1))
                            mere_avant=int(float(candidat1(c[1])))
                        else:
                            dict_orientation[mere].append(abs(o2))
                            mere_avant=int(float(candidat2(c[1])))
                            
                    else:
                        o=dict_final[k].framek[int(float(c[0]))].orientation-dict_final[k].framekk[int(float(c[1]))].orientation
                        dict_orientation[mere].append(abs(o))
                        mere_avant=int(float(c[1]))
    return(dict_orientation)
                        

  
                    



def Create_Graph(Name_experiment,z_s,Save,Print_Graph,r_vois,lambdaa,gammaa):
    
    """______________________________Description______________________
    
    Takes the name of experiment to analyse, the value of z_s ie the 2*d+z_s smallest cells,
          the bool Save to save every scenarios as PDF, and bool Print_Graph to print graph of scenarios
          
    Return dict_layers a dictionnary of scenarios objet in every transition between frames,
            G the Graph, lab array of 0 and cells label to stock the image, Frames list of list of cell objects in each frame
    
    """

    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [6, 4]
    Frames=Save_frames(Name_experiment,False)
    lab = tifffile.imread(Name_experiment)
    liste_etage=[]
    for k in range(len(Frames)):
        if k <len(Frames)-1:
            if not(Create_scenarios(k,False,Frames,lab,z_s,r_vois,lambdaa,gammaa)):
                break
            liste_etage.append(Create_scenarios(k,False,Frames,lab,z_s,r_vois,lambdaa,gammaa))
       
                
    key_list = [(i,i+1) for i in range(len(Frames))]
    del key_list[len(Frames)-1]
    value_list = liste_etage
    dict_layers=dict(zip(key_list, value_list))
    
    if Save:
        Name='Result'+Name_experiment+'z_s'+str(z_s)
        Print_results(lab,dict_layers,cells_from_frame,Frames,Name)

    G = nx.Graph()
    

    ##construction du graph
    liste_colors=["gold","violet","limegreen","darkorange","cyan","green"]
    subset_colors=["gold"]
    for p in dict_layers:

        if  p[0]<len(dict_layers)-1:
            m=randint(0,len(liste_colors)-1)
            subset_colors.append(liste_colors[m])

            for sce_0 in dict_layers[p]:
                #subset_colors.append(liste_colors[m])
                for sce_1 in dict_layers[(p[0]+1,p[1]+1)]:
                    if cost_edges(sce_0,sce_1)[1]:
                        G.add_edge(sce_0.name,sce_1.name,weight=round(cost_edges(sce_0,sce_1)[0],3))
                        
                        
    for v, data in G.nodes(data=True):
        #print(v)
        #print(data)
        data['layer']=int(find_layer(v))
        

    pos = nx.multipartite_layout(G, subset_key="layer")
    color = [subset_colors[data["layer"]] for v, data in G.nodes(data=True)]

    if Print_Graph:
        plt.figure(figsize=(8, 8))
        nx.draw_networkx_nodes(G, pos, node_size=50,node_color=color)
        # edges
        nx.draw_networkx_edges(G, pos, width=0.5)
        nx.draw_networkx_edges(
            G, pos, width=0.2, alpha=0.5, edge_color="b", style=":")
    
        # node labels
        #nx.draw_networkx_labels(G, pos, font_size=3, font_family="sans-serif")
        edge_labels = nx.get_edge_attributes(G, "weight")
        nx.draw_networkx_edge_labels(G, pos, edge_labels,font_size=3)
        #plt.axis("equal")
        plt.show()
    
    return(dict_layers,G,lab,Frames)


def Create_lineage(dict_final,lab,Print_Graph):
    G=nx.Graph()
    liste_colors=["gold","violet","limegreen","darkorange","cyan","green","gold","violet","limegreen","darkorange","cyan","green","gold","violet","limegreen","darkorange","cyan","green","gold","violet","limegreen","darkorange","cyan","green","gold","violet","limegreen","darkorange","cyan","green","gold","violet","limegreen","darkorange","cyan","green"]
    subset_colors=["gold"]
    nCells_1=lab[0].max()
    first=0
    for k in dict_final:
        #on s'occupe de l'image k+1
        #print("k ")
        #print(k)

        for c in dict_final[k].matching:
            if is_plus(c[1]):
                G.add_edge(str(k[0])+'_'+str(int(float(c[0]))),str(k[1])+'_'+candidat1(c[1]),weight=0)
                G.add_edge(str(k[0])+'_'+str(int(float(c[0]))),str(k[1])+'_'+candidat2(c[1]),weight=0)
            else:
                G.add_edge(str(k[0])+'_'+str(int(float(c[0]))),str(k[1])+'_'+c[1],weight=0)
    k=k[0]+1

    for v, data in G.nodes(data=True):
        print(v)
        print(data)
    
        data['layer']=int(find_layer(v))
    pos = nx.multipartite_layout(G, subset_key="layer")
    color = [liste_colors[data["layer"]] for v, data in G.nodes(data=True)]

    if Print_Graph:
        plt.figure(figsize=(8, 8))
        nx.draw_networkx_nodes(G, pos, node_size=50,node_color=color)
        # edges
        nx.draw_networkx_edges(G, pos, width=0.5)
        nx.draw_networkx_edges(
            G, pos, width=0.2, alpha=0.5, edge_color="b", style=":")
    
        # node labels
        nx.draw_networkx_labels(G, pos, font_size=3, font_family="sans-serif")
    
        #plt.axis("equal")
        plt.show()
    
    return(G)

#CHECK The result
def Plot_Sce_Div(list_experiment,z_s,r_vois,lambdaa,gammaa,Plot_courbe=False,legend=str(z_s)):
    
    """______________________Description______________________
    
        Take list of experiment to analyse with the z_s critaria, and a bool to plot, legend :label for the plot
        Create a scenarios graph for a z_s for every experiment, and collect the number of scenarios for each division
        return X, Y : Number of division, number of scenarios, and save X Y in Sce_Div_z_s excel file
    """
    X=[]
    Y=[]
    for index in list_experiment :  
        print("Dealing with experiment %s" %index)
        if index//10>0:
            name="2021_08_29_SD_TS66_SHU1_FOS800_01_"+str(index)+"_R3D_labelled.tif"
        else:
            name="2021_08_29_SD_TS66_SHU1_FOS800_01_0"+str(index)+"_R3D_labelled.tif"
            
        dict_layers,G,lab,Frames=Create_Graph(name,z_s,False,False,r_vois,lambdaa,gammaa)
        for frame in dict_layers:
            #division number
           # XX.append(len(dict_layers[frame][0].cells_div))
            X.append(len(dict_layers[frame][0].cells_div))
            #secnarios number
            #YY.append(len(dict_layers[frame]))
            Y.append(len(dict_layers[frame]))
    if Plot_courbe:
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['figure.figsize'] = [6, 4]
        fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
        ax.plot(X,Y,'or', label=legend)  # Plot some data on the axes.
        ax.set_xlabel('Number of div')  # Add an x-label to the axes.
        ax.set_ylabel('Number of Sce')  # Add a y-label to the axes.
        ax.legend();  # Add a legend.
    col1 = "X"
    col2 = "Y"
    data = pd.DataFrame({col1:X,col2:Y})
    data.to_excel('Sce_Div'+legend+'.xlsx', sheet_name='sheet1', index=False)
    
    return(X,Y)








"""

a,b,c,cout_chem=best_path(G, dict_layers)

plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [6, 4]
#Print_results_best(lab,dict_layers,cells_from_frame,c,'')  
            
#afficher le meilleur chemin
##on itère sur les différents sommets

pos = nx.multipartite_layout(G, subset_key="layer")
color = [subset_colors[data["layer"]] for v, data in G.nodes(data=True)]

#plt.figure(figsize=(8, 8))
#nx.draw_networkx_nodes(G, pos, node_size=50,node_color=color)

# edges
#nx.draw_networkx_edges(G, pos, width=0.5)
#si plusieurs noeuds finaux

for sce in dict_layers[(len(dict_layers)-1,len(dict_layers))]:
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [6, 4]
   # Print_results_best(lab,dict_layers,cells_from_frame,c,sce.name)  
                
    #afficher le meilleur chemin
    ##on itère sur les différents sommets

    pos = nx.multipartite_layout(G, subset_key="layer")
    color = [subset_colors[data["layer"]] for v, data in G.nodes(data=True)]


    plt.figure(figsize=(8, 8))
    nx.draw_networkx_nodes(G, pos, node_size=50,node_color=color)
    nx.draw_networkx_edges(
    G, pos, width=0.5, alpha=0.5,edgelist=create_edge_list(c[sce.name],sce.name), edge_color="r")
    nx.draw_networkx_labels(G, pos, font_size=3, font_family="sans-serif")
    plt.title("Cost best path : %s, final node : %s"%(cout_chem[sce.name],sce.name))
    #plt.axis("equal")
    plt.show()
"""
#si un seul noeud
"""
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [6, 4]
sce=dict_layers[(len(dict_layers)-1,len(dict_layers))][0]
#Print_results_best(lab,dict_layers,cells_from_frame,c,sce.name)  
 
 #afficher le meilleur chemin
 ##on itère sur les différents sommets

pos = nx.multipartite_layout(G, subset_key="layer")
color = [subset_colors[data["layer"]] for v, data in G.nodes(data=True)]


plt.figure(figsize=(8, 8))
nx.draw_networkx_nodes(G, pos, node_size=50,node_color=color)
nx.draw_networkx_edges(G, pos, width=0.5, alpha=0.5,edgelist=create_edge_list(c,sce.name), edge_color="r")
nx.draw_networkx_labels(G, pos, font_size=3, font_family="sans-serif")
plt.title("Cost best path : %s, final node : %s"%(cout_chem,sce.name))
 #plt.axis("equal")
plt.show()
"""
         


#ici on extrait le meilleur path
"""
dictfinal=Extract_best_path(dict_layers,c['28_29-1'],'28_29-1')
new_image=define_color(dictfinal,lab)
                
dict_orientation=get_orientation(dictfinal)

for t in range(new_image.shape[0]):
        #print(new_image[t].shape)

        plt.imshow(new_image[t])
        plt.title('frame : %s' %t)
        plt.show()

"""

"""
plt.figure()
for mere in dict_orientation:
    mere=1
    x=[i for i in range(1,len(dict_orientation[mere])+1)]
    plt.plot(x,dict_orientation[mere],"o--",label="cell1")
    
plt.figure()
mere=0
x=[i for i in range(1,len(dict_orientation[mere])+1)]
plt.plot(x,dict_orientation[mere],"^",label="cell0")
mere=1
x=[i for i in range(1,len(dict_orientation[mere])+1)]
plt.plot(x,dict_orientation[mere],"o--",label="cell1")
plt.legend()
plt.show()
x=[i for i in range(1,len(dict_orientation[2])+1)]

plt.plot(x,dict_orientation[2],"o--",label="cell2")
x=[i for i in range(1,len(dict_orientation[3])+1)]
plt.plot(x,dict_orientation[3],"o:",label="cell3")
x=[i for i in range(1,len(dict_orientation[4])+1)]
plt.plot(x,dict_orientation[4],"o-.",label="cell4")
x=[i for i in range(1,len(dict_orientation[5])+1)]
plt.plot(x,dict_orientation[5],"o",label="cell5")
plt.legend()
plt.show()
"""

"""
#import sys
#t=6
#np.set_printoptions(threshold=sys.maxsize)
#print(new_image[t])
"""










  


"""
figure = plt.figure()
axes = figure.add_subplot(2, 1, 1)
axes.scatter(range(5), [x ** 2 for x in range(5)], s = 50, color = 'blue')
axes.margins(1, 0.5)
axes = figure.add_subplot(2, 1, 2)
axes.scatter(range(5), [x ** 2 for x in range(5)], s = 50, color = 'red')
axes.margins(0, 0)
"""
"""
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [6, 4]
#Print_results_best(lab,dict_layers,cells_from_frame,c,'')  
            
#afficher le meilleur chemin
##on itère sur les différents sommets
liste_colors=["gold","violet","limegreen","darkorange","cyan","green"]
subset_colors=["gold"]


#plt.figure(figsize=(8, 8))
#nx.draw_networkx_nodes(G, pos, node_size=50,node_color=color)

# edges
#nx.draw_networkx_edges(G, pos, width=0.5)
#si plusieurs noeuds finaux

for sce in dict_layers[(len(dict_layers)-1,len(dict_layers))]:
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [6, 4]
   # Print_results_best(lab,dict_layers,cells_from_frame,c,sce.name)  
                
    #afficher le meilleur chemin
    ##on itère sur les différents sommets

    pos = nx.multipartite_layout(G, subset_key="layer")
    color = [subset_colors[data["layer"]] for v, data in G.nodes(data=True)]


    plt.figure(figsize=(8, 8))
    nx.draw_networkx_nodes(G, pos, node_size=50,node_color=color)
    nx.draw_networkx_edges(
    G, pos, width=0.5, alpha=0.5,edgelist=create_edge_list(c[sce.name],sce.name), edge_color="r")
    nx.draw_networkx_labels(G, pos, font_size=3, font_family="sans-serif")
    plt.title("Cost best path : %s, final node : %s"%(cout_chem[sce.name],sce.name))
    #plt.axis("equal")
    plt.show()
"""