# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 17:32:08 2022

@author: Busra Bulut
"""
import Functions
import matplotlib.pyplot as plt
import pandas as pd
##MAIN##

""" _______________________Description__________________________________

    To plot with options, run this part before
"""
option_plot=['bo','r.','g*','m.','yo','ko','b*','r*','g*','m*','y*','k*']
option_plot_2=['o--','o:','o-','o','yo','ko','b*','r*','g*','m*','y*','k*']


""" ___________________________Description_____________________________
    To plot the number of scenarios over the number of divisions for different values of z_smallest and 
    to save it as Sce(Div)  
"""
r_vois=1.5
list_experiment=[1,2,5,8,10,11,14,15,17,18,19,21,22,24,27,28,29,30,31,32,33,34,35,36,37,39]
values_z_smallest=[0,1,2,3]
#Importance for the length, orientation, proximity for the cost function of without division
lambdaa=[0.05,0.03,0.92]
#coeff of the cost function with division
gammaa=[0.3,0.1,0.1,0.2,0.2,0.1,0.3]
fig=plt.figure()

for z in range(len(values_z_smallest)):
    
    ##we already have X and Y as excel 
    df=pd.read_excel('C:/Users/Busra Bulut/Documents/EPFL/Cell_lineage/05_labelled_images-20220530T153642Z-001/update190722/05_labelled_images/Sce_Div'+str(values_z_smallest[z])+'.xlsx')
    #plt.plot(df['X'],df['Y'],option_plot[z],label='z_smallest :'+str(z))
    df_zoom=df.loc[df['Y']<5000]
    plt.plot(df_zoom['X'],df_zoom['Y'],option_plot[z],label='z_smallest :'+str(z))
 

    ##we don't have 

    X,Y= Functions.Plot_Sce_Div(list_experiment,values_z_smallest[z],r_vois,lambdaa,gammaa,False,legend=str(values_z_smallest[z]))
    plt.plot(X,Y,option_plot[z],label='z_smallest :'+str(z))

plt.legend()
plt.xlabel('Number of Divisions')
plt.ylabel('Number of Scenarios')
#fig.savefig('Sce(Div)_5.jpg')
plt.show()



""" ___________________________Description_____________________________
    To plot the number of scenarios over the number of divisions for different values of r_vois and 
    to save it as Sce(Div)  
"""

X,Y=Functions.Plot_Sce_Div(list_experiment,2,3,False,legend='r_vois_3')

df1_2=pd.read_excel('C:/Users/Busra Bulut/Documents/EPFL/Cell_lineage/05_labelled_images-20220530T153642Z-001/update190722/05_labelled_images/Sce_Div_r_1_2.xlsx')
df1_5=pd.read_excel('C:/Users/Busra Bulut/Documents/EPFL/Cell_lineage/05_labelled_images-20220530T153642Z-001/update190722/05_labelled_images/Sce_Divr_vois_1_5.xlsx')
df0_9=pd.read_excel('C:/Users/Busra Bulut/Documents/EPFL/Cell_lineage/05_labelled_images-20220530T153642Z-001/update190722/05_labelled_images/Sce_Divr_vois_0_9.xlsx')


plt.figure()
plt.plot(df0_9['X'],df0_9['Y'],option_plot[0],label='r_vois : 0.9')
plt.plot(df1_2['X'],df1_2['Y'],option_plot[2],label='r_vois : 1.2')
plt.plot(df1_5['X'],df1_5['Y'],option_plot[1],label='r_vois : 1.5')
plt.legend()
plt.xlabel('Number of Divisions')
plt.ylabel('Number of Scenarios')
#fig.savefig('Sce(Div)_5.jpg')
plt.show()







""" ___________________________Description_____________________________
    To plot the min of the difference of orientation between cells picked by the best path
"""

##select the experiment ##
Name_experiment="2021_08_29_SD_TS66_SHU1_FOS800_01_10_R3D_labelled.tif"

#To save the experiment and plot the frames 
#Functions.Save_frames(Name_experiment, 10, True)
##create the graph of scenarios

#Parameters
z_smallest=3
r_vois=1.5

#Importance for the length, orientation, proximity for the cost function of without division
lambdaa=[0.05,0.03,0.92]

#coeff of the cost function with division
gammaa=[0.3,0.1,0.1,0.2,0.2,0.1,0.3]

#To construct the Graph
dict_layers, G,lab,Frames=Functions.Create_Graph(Name_experiment,z_smallest,False,True,r_vois,lambdaa,gammaa)



##To find the best path 
a,b,c,cout_chem=Functions.best_path(G, dict_layers)
#Here we have to test if c is a dict or a list, ie if several nodes for the last layer
#if several nodes, then we print cout_chem and pick the best one

#c.insert(0,'28_29-1')

#Let's choose the best final node
best_final_node = '28_29-9'
Functions.Print_results_best(lab,dict_layers,c[best_final_node],best_final_node,Frames,Name_experiment)

#Functions.Print_results_best(lab,dict_layers,cells_from_frame,c,sce.name)  
 
 #afficher le meilleur chemin

#ici c doit etre une liste 
dictfinal=Functions.Extract_best_path(dict_layers,c,'28_29-1')

G_lineage=Functions.Create_lineage(dictfinal, lab, True) 

dict_orientation=Functions.get_orientation(dictfinal)

#To plot the evolution of the diff of orientation over time (Frame)
plt.figure()
for mere in range(len(dict_orientation)):
    x=[i for i in range(1,len(dict_orientation[mere])+1)]
    plt.plot(x,dict_orientation[mere],option_plot_2[mere],label="cell"+str(mere))
plt.xlabel('Frames')
plt.ylabel('|Difference of orientation|')
plt.legend()
plt.show()


        
""" ___________________________Description_____________________________
    Create an illustration to show the cell lineage selected by the algorithm
"""
new_image=Functions.define_color(dictfinal,lab,'Image_exp_10_z_3')
for t in range(new_image.shape[0]):
        plt.imshow(new_image[t])
        plt.title('frame : %s' %t)
        plt.show()      
                   
import tifffile
import numpy as np

labb = tifffile.imread("Jet11.tif")

Fr=Functions.Save_frames("2021_08_29_SD_TS66_SHU1_FOS800_01_11_R3D_labelled.tif",False)
t=27
i_mam=4
frame_6=Functions.cells_from_frame(Fr[t])
Ll=Functions.average_length(frame_6)

plt.imshow(labb[t],interpolation='none')

 
theta = np.linspace( 0 , 2 * np.pi , 20)
 
radius = 2*Ll
 
a =frame_6[i_mam].x+ radius * np.cos( theta )
b = frame_6[i_mam].y + radius * np.sin( theta )
plt.plot(a,b,"r-")
plt.plot(frame_6[i_mam].x,frame_6[i_mam].y ,'bx')
plt.title('frame : %s' %t)
plt.show()


for t in range(labb.shape[0]): #on parcours tous les frames
    plt.imshow(labb[t],interpolation='none')
    plt.title('frame : %s' %t)
    plt.show()
       


from itertools import combinations
A = [3,6,8,9,11]
temp = combinations(A, 4)
for i in list(temp):
	print (i)
        
        