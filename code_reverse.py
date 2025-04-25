# -*- coding: utf-8 -*-
"""

@author: sanmuhe
"""

import itertools
import json
import statistics
from inspect import getsourcefile
import os
import pandas as pd
import numpy as np
import results as rs
import networkx as nx
import xlrd
import matplotlib.pyplot as plt
path = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))

data='data'       # data is the data set file name without extension
          # graph is the variable to store the network usin dictionary

output=path+'\\Output'
input=path+'\\Input'

#该函数是返回这两个列表的交集；
def Intersection(lst1, lst2): 
    return set(lst1).intersection(lst2)  

#返回两个元素的并集。
def Union(lst1, lst2): 
    final_list = lst1 + lst2 
    t=set(final_list)
    return t

def graphBuild(resultPath,dataPath):
    #其中的resultPath,和dataPath都是一个str类型的字符串；
    edgeF=open(input+'\\'+dataPath+'.txt')
    #就是打开input文件里面的data.txt文件；
    lvl=[]
    flag=0
    i=0
    graph={}
    #这里的graph是一个字典；
    while True:
        a=edgeF.readline()
        if a == '':
            break
        b=a.split()
        i+=1
        ##print (b)
        if flag == 0:
            lvl.append(b[1])
            graph[b[0]]=lvl.copy()
            #像这个字典里面添加键值对，键是b[0],值是后面的；
            lvl.clear()
            lvl.append(b[0])
            graph[b[1]]=lvl.copy()
            flag=1
            #添加一个互相的键值对；
        else:
            if b[0] in graph:
                graph[b[0]].append(b[1])
            #之所以这样，是因为一个键往往有很多个值，不是只有一个对；
            else:
                lvl.append(b[1])
                graph[b[0]]=lvl.copy()
                lvl.clear()
            if b[1] in graph:
                graph[b[1]].append(b[0])
            else:
                lvl.append(b[0])
                graph[b[1]]=lvl.copy()
            #这里只是给后面的元素添加键值对，每个元素都有，是互相的；
            #这里的形式和第一个if其实是一样的；只不过不是第一行了；
        lvl.clear()
    edgeF.close()

    file3=open(output+'\\details.txt','w+')
    file3.write('At begin :\n'+'node: '+ str(len(graph))+'  Edge: '+ str(i)+'\n')
    file3.close()

    file2=open(output+'\\'+resultPath+'.txt','w+')
    #这里的rusultPath里面的mainGraph,是在一开始传递时的那个形参；mainGraph.txt;
    file2.write(json.dumps(graph,indent=4, separators=(",",":")))
    #这里的意思是将这个graph对象转化为json文件，后面indent表示开始要缩进四格，后面是说的分隔符号；
    #最后写入file2文件中；
    file2.close()

    return graph
#本段代码主要就是将文件里面的每个元素建立其一个互相的图（双向的），最后返回的是一个字典，
#这个字典里面有的东西是每个数据的键值对；

def graphBuildList(resultPath,dataPath):
    #dataPath就是reduceEdge;,
    edgeF=open(output+'\\'+dataPath+'.txt')
    edgeList=json.load(edgeF)    
    lvl=[]
    lvl.clear()
    graph={}
    flag=0
    ##print('Edge in reduced graph: ',edgeList)
    for b in edgeList:
        ###print (b)
        if flag == 0:
            lvl.append(b[1])
            graph[b[0]]=lvl.copy()
            lvl.clear()
            lvl.append(b[0])  
            graph[b[1]]=lvl.copy()
            flag=1
        #这里的步骤和上面的一样都是建立一个双向的字典。
        else:
            if b[0] in graph:
                graph[b[0]].append(b[1])
            else:
                lvl.append(b[1])
                graph[b[0]]=lvl.copy()
                lvl.clear()
            if b[1] in graph:
                graph[b[1]].append(b[0])
            else:
                lvl.append(b[0])
                graph[b[1]]=lvl.copy()
        lvl.clear()
        #后面的这个也是相关的删除后面的要删除的边的过程，然后进行保留。
    file2=open(output+'\\'+resultPath+'.txt','w+')
    file2.write(json.dumps(graph, indent=4, separators=(",", ":")))
    file2.close()
    #最后就是生成一个删除的边的图；
    return graph
    
#
def uniqueProtein(filePath):
    uniqueProteinF=open(input+'\\'+filePath+'.txt')
    #这里同样是data;
    s1=set()
    while True:
        a=uniqueProteinF.readline()
        if a == '':
            break
        b=a.split()
    #这里的就是使用默认的分隔符号，也就是“”空格和换行符号等。将其分割成一个列表；
        s1.add(b[0])
        s1.add(b[1])
    ##print(s1)
    #由于s1是一个集合，所以最后返回的是所有的关键蛋白的种类，没有重复的，集合要求是没有重复的;
    uniqueProteinList=list(s1)
    uniqueProteinList.sort()
    #这个排序就是按照字典序来进行排列，字符串就是按照第一个字母进行比较，然后根据第二个来进行;
    #排序,最后来这样进行；
    #print('Unique protein list: ',uniqueProteinList);
    uniqueProteinF.close()
    return uniqueProteinList
#这里的data.txt里面所有的都是关键蛋白，不分b[0]和b[1];

def uniqueProteinL(filePath):
    uniqueProteinF=open(output+'\\'+filePath+'.txt')
    uniqueProteinEdgeList=json.load(uniqueProteinF)
    s1=set()
    for b in uniqueProteinEdgeList:
        s1.add(b[0])
        s1.add(b[1])
    ##print(s1)
    uniqueProteinList=list(s1)
    uniqueProteinList.sort()
    #print('Unique protein list: ',uniqueProteinList)
    uniqueProteinF.close()
    return uniqueProteinList

#   Non Essential Protein
def essentialNonEssential(k):
    reducedProteinF=open(output+'\\orthologous_IDC'+str(k)+'.txt')
    reducedProteinList=json.load(reducedProteinF)
    #print(reducedProteinList)
    totalProtein=len(reducedProteinList)
    print(totalProtein)
    essentialProteinNo=int(totalProtein*0.20)
    #取所有蛋白值的百分之20进行查准率的比较，是750*0.2左右的蛋白质；最后结果就是这个；
    essentialProteinList=[]
    nonEssentialProteinList=[]
    count=0
    for i in reducedProteinList:
        if count<=essentialProteinNo :
            essentialProteinList.append(i[0])
        else:
            nonEssentialProteinList.append(i[0])
        count+=1
    #print(s1)
    #这里的意思就是根据LIDC里面的计算方法，前百分之20是关键蛋白，后面的不是，因为前面已经根据LIDC的值，进行过排序，所以这样符合要求。
    essentialProteinF=open(output+'\\essentialProtein'+str(k)+'.txt','w+')
    essentialProteinF.write(json.dumps(essentialProteinList,indent=4 ,separators=(",",":")))
    essentialProteinF.close()
    
    nonEssentialProteinF=open(output+'\\nonEssentialProtein'+str(k)+'.txt','w+')
    nonEssentialProteinF.write(json.dumps(nonEssentialProteinList,indent=4 ,separators=(",",":")))
    nonEssentialProteinF.close()
    
#
#
def nodeWeight(data):
    graph=graphBuild('mainGraph',data)
    #这里时建立这个图的一个字典；graph是一个字典；
    ###print(graph)
    proteinList=uniqueProtein(data)
    #这个就是关键蛋白的列表；
    totalProtein=len(proteinList)
    ##print("totalProtein: ",totalProtein)
    ###print(proteinList)
    proteinDegList=[]
    #蛋白质度的列表；
    for p in proteinList:
        levList=graph[p]
        deg=0
        for child in levList:
            childList=graph[child]
            deg+=len(childList)
        proteinDegList.append(deg/totalProtein)
    #这里可以看成是每个节点的邻居节点/所有蛋白的比例；
    ##print(proteinDegList)
    degF=open(output+'\\nodeDeg'+'.txt','w+')
    degF.write(json.dumps(proteinDegList,indent=4, separators=(",",":")))
    degF.close()    
    alpha = statistics.mean(proteinDegList)
    #alpha是算数平均值，也就是所有节点平均度的算数平均值；
    SD = statistics.stdev(proteinDegList)
    #这里是算出这个list的标准偏差，这里也就是反映数据的分散程度，表示数据值与平均值的偏离程度；
    for k in range(1,4):
    #这里是从1到3,就是1，2，3不包括4，定义为低、介质、和高的否决，alpha被认为是边缘节点重点/体重的标准偏差值。
        print('started k = ',k)
        nodeDelete(alpha,SD,proteinList.copy(),proteinDegList.copy(),graph,k)
        print('finished k = ',k)
# SD = statistics.stdev(proteinDegList) 是一行 Python 代码，它使用 statistics 模块中的 stdev 函数来计算 proteinDegList 列表中所有元素的样本标准偏差，
# 并将结果赋值给变量 SD。
# statistics.stdev(data, xbar=None) 函数返回 data 的样本标准偏差，其中 data 是一个序列或迭代器。样本标准偏差是一种衡量数据集离散程度的统计量，
# 表示数据值与平均值的偏离程度。它是方差的平方根，用于描述数据的分散程度。与 pstdev() 方法不同，stdev() 方法计算的是样本标准偏差而不是总体标准偏差，
# 因此使用的是无偏估计。如果你有整个总体的数据，可以使用 pstdev() 方法计算总体标准偏差1。
# 例如，如果 proteinDegList = [1, 2, 3, 4, 5]，则执行 SD = statistics.stdev(proteinDegList) 后，变量 SD 的值将为 1.5811388300841898。

def nodeDelete(alpha,SD,proteinList,proteinDegList,graph,k):

    THk=alpha+k*SD*(1-(1/(1+(SD**2))))
    #公式；θk = α + k×σ×(1 − 1/( 1 + σ2))
    print('For k =',k,'Th = ',THk)
    #print(uniqueProteinList)
    #print(proteinDegList)
    graphTemp=graph.copy()
    #graph是一个双向图；
    for node,val in list(zip(proteinList,proteinDegList)):
#zip函数将两个列表proteinList和proteinDegList中的元素逐一配对，然后使用for循环遍历这些配对。
#然后返回一个迭代器，该迭代器生成由可地带对象中的元素组成的元组。而，node就是每个proteinList
#对应的元素而val就是proteinDegList中对应的值；
#print('Node :',node,'val: ',va0.2l)
        if val < 0*THk:
            graphTemp.pop(node)
            #这里会删除键值对，就是都删除的意思；
            #print('R Node :',node,'val: ',val)
            for parent in graphTemp:
                child=set(graphTemp[parent])
            #child = set(graphTemp[parent])是parent这个键对应的值。
            #它从字典graphTemp中获取键为parent的值然后使用set函数将其转化为一个。
            #集合，并将结果赋值给child;注意集合里面没有一样的值。
            #可以理解为每一个节点连接的其他节点，就是孩子节点；
                if node in child:
            #此时child就是图的度表的集合。
                    graphTemp[parent].remove(node)
            #如果这个节点在孩子节点里面，要删除所有节点的这个node点，因为其他所有的节点都有可能与这个节点由联系。
            #也就是这个节点对应字典中键值对的每个这个节点的值。
            proteinList.remove(node)
            #最后在彻底删除这个节点。
    proteinListSet=set(proteinList)
    #proteinList是关键蛋白的列表。注意是剪掉节点后的这个列表，上面是根据节点的截止值来删减。
    edge=open(input+'\\data.txt')
    reducedEdgeList=[]
    while True:
        a=edge.readline()
        if a == '':
            break
        i=a.split()
        if i[0] in proteinListSet and i[1] in proteinListSet:
            reducedEdgeList.append(i)
        #这里的i[0]和i[1]都在这个关键蛋白里面，就开始把这个边加到要剪的边里面。
    edge.close()
    df=pd.DataFrame(reducedEdgeList)
    df.to_csv(output+'\\edgesAfterNodeReduction'+str(k)+'.csv',sep=',',header=False,index=False)
    reducedEdgeF=open(output+'\\edgesAfterNodeReduction'+str(k)+'.txt','w+')
    reducedEdgeF.write(json.dumps(reducedEdgeList,indent=4,separators=(",",":")))
    reducedEdgeF.close()
    reducedGraphF=open(output+'\\reducedNodeGraph'+str(k)+'.txt','w+')
    reducedGraphF.write(json.dumps(graphTemp,indent=4,separators=(",",":")))
    reducedGraphF.close()
    file3=open(output+'\\details.txt','a+')
    file3.write('For k = '+ str(k) +'\n After node reduce :\n'+'node: '+ str(len(proteinList))+'  Edge: '+ str(len(reducedEdgeList))+'\n')
    file3.close()
    edgeWeight(k,graphTemp)

#   feature_pcc(k,graphTemp)
    #graphTemp是剪枝节点后的图。因为会进行三次，所以这里也是进行三次的。

def WeightgraphBulid(data):
    #第一步计算图的PCC；
    #第二部，得到每个图的边的权值；
    #最后根据每个节点的边的所有的权值进行剪枝，最后得到一个图；进行操作。
    graph = graphBuild('mainGraph',data)
    #这里的graph其实也是一个字典；
    protein_sum = open(input+'\\'+'data'+'.txt')
    gene_expressdata = open(input+'\\'+ 'yeastgenedata' +'.txt')
    #建立一个基因表达矩阵,
    s1 = set()
    gene_expressdatamatrix = {}
    Act_th = {}
    lvl = []
    lvl.clear()
    while True:
            a = gene_expressdata.readline()
            if a == '':
                break
            b = a.split()
            for i in range(2,38):
                lvl.append(b[i])
            gene_expressdatamatrix[b[1]] = lvl.copy()
            Act_th[b[1]] = 0
            lvl.clear()
    gene_expressdata.close()
    df = pd.DataFrame(gene_expressdatamatrix)
    df = df.T
    gene_df = df
    Act_th = pd.DataFrame([Act_th])
    Act_th = Act_th.T
    for gene,values in gene_expressdatamatrix.items():
        values = np.array(values)
        values = values.astype(float)
        mean = np.mean(values)
        std = np.std(values)
        Act_th.loc[gene] = mean + 3*std*(1-(1/(1+std*std)))

    for index,row in gene_df.iterrows():
        for col in gene_df.columns:
            val = np.array(row[col])
            val = val.astype(float)
            if val>=Act_th.loc[index,0]:
                gene_df.loc[index,col]=1
            else:
                gene_df.loc[index,col]=0
    AdjMatrix = pd.DataFrame(index=graph.keys(), columns=graph.keys())
    AdjMatrix = AdjMatrix.fillna(0).astype(float)
    # result = pd.DataFrame(0, index=df.index, columns=df.index)
    # for label_i in gene_df.index:
    #     for label_j in gene_df.index:
    #         # 计算每一位同时为1的数量
    #         if label_i!=label_j:
    #             count = np.sum((gene_df.loc[label_i,:] + gene_df.loc[label_j,:] ==2))
    #             if count>0:
    #         # 将结果存储在新的DataFrame中
    #                 result.loc[label_i, label_j] = count
    #                 result.loc[label_j, label_i] = count
    # print(result)
    # n = len(df.index)
    # for i in range(n):
    #     for j in range(i + 1, n):
    #         label_i = df.index[i]
    #         label_j = df.index[j]
    #         # 计算每一位同时为1的数量
    #         count = np.sum((gene_df.loc[label_i, :] + gene_df.loc[label_j, :] == 2))
    #         if label_i in AdjMatrix.index and label_j in AdjMatrix.index and count > 0:
    #             # 将结果存储在新的DataFrame中
    #                 AdjMatrix.loc[label_i, label_j] = count
    #                 AdjMatrix.loc[label_j, label_i] = count
    # print(AdjMatrix)
    gene_array = gene_df.values

    # 计算所有行对的和
    sums = gene_array[:, None, :] + gene_array[None, :, :]

    # 找出和为2的元素
    mask = sums == 2

    # 计算每一对行中和为2的元素的数量
    counts = np.sum(mask, axis=2)

    # 创建一个新的DataFrame来存储结果
    AdjMatrix = pd.DataFrame(counts, index=gene_df.index, columns=gene_df.index)

    # 将对角线上的元素设置为0，因为我们不关心行与其自身的比较
    np.fill_diagonal(AdjMatrix.values, 0)


    while True:
        a = protein_sum.readline()
        if a == '':
            break
        b = a.split()
        s1.add(b[0])
        s1.add(b[1])

    protein_sumList = list(s1)
    pro_total_length = len(protein_sumList)
    index = list(df.index)
    matrix = np.zeros((pro_total_length,36))
    for i in range(0,pro_total_length):
         if index[i] in protein_sumList:
             for j in range(0,35):
                 matrix[i,j] = df.iloc[i,j]
    #基因表达矩阵建立完成；
    standarddeviation = np.zeros((pro_total_length,1))
    meanvalue = np.zeros((pro_total_length,1))
    Act_temp = np.zeros((pro_total_length, 1))
    #计算方差和均值；
    for i in range(0,pro_total_length):
        standarddeviation[i] = np.std(matrix[i])

    for i in range(0, pro_total_length):
        meanvalue[i] = np.mean(matrix[i])
    standarddeviationhl = {}

    meanvaluehl = {}

    for index,value in zip(df.index,standarddeviation):
        standarddeviationhl[index] = value
    #给方差打个label;标签就是索引。
    #计算ECC
    for index,value in zip(df.index,meanvalue):
        meanvaluehl[index] = value
    #给均值打个label;标签就有了；
    ECC = np.zeros((pro_total_length,pro_total_length))
    PCC = {}
    #PCC数组；
    PCCMatrix = pd.DataFrame(index =graph.keys(),columns=graph.keys())
    PCCMatrix = PCCMatrix.fillna(0).astype(float)
    #averagePCC = 0
    # 计算PCC矩阵,这样的PCC矩阵可以使用字符串索引
    df = df.astype(np.float64)
    for key,values in graph.items():
            if key in standarddeviationhl:
                for value in values:
                    if value in standarddeviationhl:
                        if standarddeviationhl[key] != 0 and standarddeviationhl[value] != 0:
                            PCC[key] = 0
                            for i in range(0,36):
                                PCC[key] =  PCC[key]+(((df.loc[key,i])-meanvaluehl[key])/standarddeviationhl[key]*((df.loc[key,i])-meanvaluehl[value])/standarddeviationhl[value])
                            PCCMatrix.loc[key,value] = PCC[key]/35
                            PCCMatrix.loc[value,key] = PCC[key]/35
    array = PCCMatrix.to_numpy().flatten()

    #根据之前建立一个基因表达数据的dataframe有标签；的再建立一个基本的数据集合就是关于阈值计算的；


    #PCCstandarddeviation = np.std(array)
    #averagePCC = np.mean(array)
    #[0.21086401]
                            #要赋值两条相互的边；
                            #问题出在df.iloc里面的键值索引问题,
    # for i in range(1,pro_total_length):
    #            for j in range(1,36):
    #                #这里就是计算每条邻接矩阵的每一条边的方差都不为0才成立；
    #                if standarddeviation[graph[i]]!=0 and standarddeviation[i]!=0:

    #创建基因活性矩阵

    return PCCMatrix,array,AdjMatrix

def edgeWeight(k,graph):
    #print(data)    
    # print(data)
    PCCMatrix, array,AdjMatrix = WeightgraphBulid(data)
    j = 0
    # print(data)
    edgeF = open(output + '\\edgesAfterNodeReduction' + str(k) + '.txt')
    edgeList = json.load(edgeF)
    # print(edgeList)
    edgeWeightVal = []
    #print(edgeList[1][0])
    DegMatrix = pd.DataFrame(index=graph.keys(), columns=graph.keys())
    DegMatrix = PCCMatrix.fillna(0).astype(float)
    for j in range(len(edgeList)):
        lst1 = graph.get(edgeList[j][0])
        ##print(lst1)
        lst2 = graph.get(edgeList[j][1])
        ##print(lst2)
        w = len(Intersection(lst1, lst2)) / (len(lst1)*len(lst2))**0.5
        ##print(w)
        DegMatrix.loc[edgeList[j][0], edgeList[j][1]] = w
        edgeWeightVal.append(w)
        # lst1=lst2
    j = j + 1
    edgeF.close()
    # print('Edge weight: ',edgeWeightVal)
    # print('Edge: ',edgeList)
    alpha = statistics.mean(edgeWeightVal)
    SD = statistics.stdev(edgeWeightVal)
    ##print(alpha)
    ##print(SD)
    # temp = list()
    TH1 = edgeWeightVal
    T1 = edgeList
    THk = alpha +   k * SD * (1 - (1 / (1 + (SD ** 2))))
    AveragePCC = np.mean(array)
    PCCstandarddeviation = np.std(array)
    THk1 = AveragePCC +  k * PCCstandarddeviation * (1 - (1 / (1 + (PCCstandarddeviation ** 2))))

    # print('Threshhold for k = ',k,': ',THk)

    # new_T1 = []
    # for i in range(0, len(edgeList)):
    #    if edgeList[i][0] in PCCMatrix.index or edgeList[i][1] in PCCMatrix.index:
    #          if PCCMatrix.loc[edgeList[i][0], edgeList[i][1]] >=THk1 or DegMatrix.loc[edgeList[i][0], edgeList[i][1]] >=THk:
    #             new_T1.append(edgeList[i])
    # T1 = new_T1

    new_T1 = []
    for i in range(0, len(edgeList)):
       if edgeList[i][0] in AdjMatrix.index and edgeList[i][1] in AdjMatrix.index :
           if AdjMatrix.loc[edgeList[i][0],edgeList[i][1]]>0 and DegMatrix.loc[edgeList[i][0],edgeList[i][1]] >=THk :
                new_T1.append(edgeList[i])

    T1 = new_T1
    print('Edges after reduction: ', len(T1))

    file2 = open(output + '\\reducedEdge' + str(k) + '.txt', 'w+')
    file2.write(json.dumps(T1, indent=4, separators=(",", ":")))
    file2.close()

    proteinList = uniqueProteinL('reducedEdge' + str(k))

    reducedProteinF = open(output + '\\reducedNodeAfterEdgeReduce' + str(k) + '.txt', 'w+')
    reducedProteinF.write(json.dumps(proteinList, indent=4, separators=(",", ":")))
    reducedProteinF.close()

    file3 = open(output + '\\details.txt', 'a+')
    file3.write('For k = ' + str(k) + '\n After edge reduce :\n' + 'node: ' + str(len(proteinList)) + '  Edge: ' + str(len(T1)) + '\n')

    file3.close()
    df = pd.DataFrame(T1)
    df.to_csv(output + '\\edgesAfterEdgeReduction' + str(k) + '.csv', sep=',', header=False, index=False)
    LIDCDriver(k, graphBuildList('reducedEdgeGraph' + str(k), 'reducedEdge' + str(k)))
    essentialNonEssential(k)
#%%

def neighbour_affinity(k):
    edgeF = open(output + '\\reducedEdgeGraph'+str(k)+'.txt')
    edgeList = json.load(edgeF)
    #测量两个基因之间的重叠程度；
    df = pd.DataFrame(columns=edgeList.keys(), index=edgeList.keys())
    df = df.fillna(0)
    lst_1 = []
    lst_2 = []
    nodes = list(edgeList.keys())
    for node1, node2 in itertools.combinations(nodes,2):
            lst_1 = edgeList[node1]
            lst_2 = edgeList[node2]
            if  len(lst_1)!=0 and len(lst_2)!=0:
                lst_1 = set(lst_1)
                lst_2 = set(lst_2)
                ##print(lst1)
                ##print(lst2)
                w = (len(Intersection(lst_1, lst_2))**2) / len(lst_1)*len(lst_2)
                #计算一个lst_1和lst_2的聚类系数；
                ##print(w)
                #lst1_copy = list(edgeList[node1])
                #lst2_copy = list(edgeList[node1])
                df.loc[node1,node2] = w
                df.loc[node2,node1] = w
            else:
                continue
    return df
    #这个图就计算了
    # for j in range(len(edgeList)):
    #     lst1 = edgeList.get(edgeList[j][0])
    #     ##print(lst1)
    #     # 这里的get是接受一个键作为参数，并且返回改键对应的值；
    #     # 这里指的就是第一行所对应的数；
    #     lst2 = edgeList.get(edgeList[j][1])
    #     ##print(lst2)
    #     w = len(Intersection(lst1, lst2)) / len(Union(lst1, lst2))
    #     ##print(w)
    #     .append(w)
    #遍历字典类型，然后进行操作；
#     edgeList = json.load(edgeF)
#     edgeWeightVal = []
#     # print(edgeList)
#     ##print(edgeList[1][0])
#     # 这里就是对应所有的行数所做出的遍历。
#     # 这里的就是对每一对节点算的权值w就是交集/并集的总数量。
#     #for j in range(len(edgeList)):
#        # lst1 = graph.get(edgeList[j][0])
#         ##print(lst1)
#         # 这里的get是接受一个键作为参数，并且返回改键对应的值。
#         # 这里指的就是第一行所对应的数；
#       #  lst2 = graph.get(edgeList[j][1])
#         ##print(lst2)
#        # w = (len(Intersection(lst1, lst2)))**2/ (len(lst1)*len(lst2))
#         ##print(w)
#       #  edgeWeightVal.append(w)
#         # 这一句的j其实完全没有必要不知道为什么要这样设置，这里的j是另一个变量，不会影响循环里面的j;
#   # edgeF.close()
#   #  print(edgeWeightVal)
#     return 0

def orthologous():
    orthologous_data = open(input + '\\' + 'Nucleus_protein_allorthology' + '.txt')
    orthologous_datamatrix = {}
    lvl = []
    lvl.clear()
    while True:
        a = orthologous_data.readline()
        if a == '':
            break
        b = a.split()
        lvl.append(b[2])
        lvl.append(b[3])
        orthologous_datamatrix[b[1]] = lvl.copy()
        lvl.clear()

    orthologous_data.close()
    df = pd.DataFrame(orthologous_datamatrix)
    list = [0]
    orthologous_matrix = pd.DataFrame(index=df.keys(),columns=None)
    #print(df)
    #计算直系同源的值；
    #for i in range(0,5093):
    #print(orthologous_matrix.loc[])
    for key in df.keys():
        orthologous_matrix.loc[key,0]=int(df.loc[1,key])/99
    #计算出每个
       # print(df.loc[1,key])
    #print(orthologous_matrix)
    #计算的每个基因同源性评分；
    return orthologous_matrix
# def LID(k,graph):
#     node=set()
#     LIDList={}
#     for key in graph:
#         v=graph[key]
#         #这里的其实是这个里面的一个键对应的一个值，
#         s1=set(v)
#         edge=0
#         for i in v:
#             v1=graph[i]
#             #这里的v1也是一个键的集合；
#             edge+=len(set(v1).intersection(s1))
#             #计算v1(graph对应的字典中对应的许多值)，每次都会遍历它的邻居的邻居和这个节点的交集，如果有交集，则证明有条边，后面/2的操作是因为，
#             #这种交集的模式每条边都是会重复一次，所以要除以2；
#             node.update(set(v1).intersection(s1))
#             #这一句和这个不一样的就是node里面的不是数字，而是每个蛋白质的名称，就是节点，而上面一句其实加入的是一个键的的值的数量，也就是一个数字。
#             #j loop end
#         #i loop end
#         ##print(node)
#         edge=edge/2
#         if(len(node)>0):
#             w=edge/len(node)
#         else:
#             w=0
#         node.clear()
#         LIDList[key]=w
#     LIDListSorted=sorted(LIDList.items(), key=lambda t: t[1],reverse=True)
#     #这个代表是将LIDList（字典类型）里面的数值按照键值对里面的值进行降序排列，reverse=True代表降序，为false则为升序排列。
#     key = lambda t: t[1]
#     # 是一个可选参数，它指定了用于排序的键函数。在这种情况下，键函数是一个匿名函数（也称为
#     # lambda 函数），它接受一个元组 t
#     # 作为输入，并返回元组中的第二个元素
#     # t[1]（即值）。这意味着排序将根据字典中每个项的值进行。
#     LID=open(output+'\\LID'+str(k)+'.txt','w+')
#     LID.write(json.dumps(LIDList,indent=4 ,separators=(",",":")))
#     LID.close()
#     #print(LIDList)
#     return(LIDListSorted)
# #%%
def NC(num,graph):
    edges = pd.read_csv(output + '\\' + 'edgesAfterEdgeReduction' + str(num) + '.csv', header=None)
    G = nx.Graph()
    G.add_edges_from(edges.values)
    neigh_centrality = dict(list(nx.common_neighbor_centrality(G)))
    neigh_centrality_List = sorted(neigh_centrality.items(), key=lambda t: t[1], reverse=True)
    NC_List = open(output + '\\NC_List' + str(num) + '.txt', 'w+')
    NC_List.write(json.dumps(neigh_centrality_List, indent=4, separators=(",", ":")))
    NC_List.close()
    return 0
def SC(num,graph):
    edges = pd.read_csv(output + '\\' + 'edgesAfterEdgeReduction' +str(num) + '.csv', header=None)
    G = nx.Graph()
    G.add_edges_from(edges.values)
    subgraph_centrality = dict(nx.subgraph_centrality(G))
    subgraph_centrality_List = sorted(subgraph_centrality.items(), key=lambda t: t[1], reverse=True)
    SC_List = open(output + '\\SC_List' + str(num) + '.txt', 'w+')
    SC_List.write(json.dumps(subgraph_centrality_List, indent=4, separators=(",", ":")))
    SC_List.close()
    return 0
def EC(num,graph):
    edges = pd.read_csv(output + '\\' + 'edgesAfterEdgeReduction' +str(num) + '.csv', header=None)
    G = nx.Graph()
    G.add_edges_from(edges.values)
    eigenvector_centrality = dict(nx.eigenvector_centrality(G))
    eigenvector_centrality_List = sorted(eigenvector_centrality.items(), key=lambda t: t[1], reverse=True)
    EC_List = open(output + '\\EC_List' + str(num) + '.txt', 'w+')
    EC_List.write(json.dumps(eigenvector_centrality_List, indent=4, separators=(",", ":")))
    EC_List.close()
    return 0
def CC(num,graph):
    #贴进度中心性;
    edges = pd.read_csv(output + '\\' + 'edgesAfterEdgeReduction' + str(num) + '.csv', header=None)
    G = nx.Graph()
    G.add_edges_from(edges.values)
    closeness_centrality = dict(nx.closeness_centrality(G))
    closeness_centrality_List = sorted(closeness_centrality.items(), key=lambda t: t[1], reverse=True)
    CC_List = open(output + '\\CC_List' + str(num) + '.txt', 'w+')
    CC_List.write(json.dumps(closeness_centrality_List, indent=4, separators=(",", ":")))
    CC_List.close()
    return 0
def DC(num,graph):
    #使用networkx来写DC
    # ReducedEdgeGraph = open(output + '\\reducedEdgeGraph' + str(k) + '.txt')
    # ReducedEdgeGraph = json.load(ReducedEdgeGraph)
    edges = pd.read_csv(output + '\\' + 'edgesAfterEdgeReduction'+str(num) + '.csv', header=None)
    G = nx.Graph()
    G.add_edges_from(edges.values)
    degree_centrality = dict(nx.degree_centrality(G))
    degree_centrality_List = sorted(degree_centrality.items(), key=lambda t: t[1], reverse=True)
    DC_List = open(output + '\\DC_List' + str(num) + '.txt', 'w+')
    DC_List.write(json.dumps(degree_centrality_List, indent=4, separators=(",", ":")))
    DC_List.close()

    return 0
def BC(k,graph):
    edges = pd.read_csv(output + '\\' + 'data1' + '.csv', header=None)
    G = nx.Graph()
    G.add_edges_from(edges.values)
    betweenness_centrality = dict(nx.betweenness_centrality(G))
    betweenness_centrality_List = sorted(betweenness_centrality.items(), key=lambda t: t[1], reverse=True)
    BC_List = open(output + '\\BC_List' + str(k) + '.txt', 'w+')
    BC_List.write(json.dumps(betweenness_centrality_List, indent=4, separators=(",", ":")))
    BC_List.close()
    return 0
def IC(k,graph):
    edges = pd.read_csv(output + '\\' + 'edgesAfterEdgeReduction' +str(k)+ '.csv', header=None)
    G = nx.Graph()
    G.add_edges_from(edges.values)
    if not nx.is_connected(G):
        # 获取最大连通子图
        largest_connected_component = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_connected_component).copy()
    information_centrality = dict(nx.information_centrality(G))
    information_centrality_List = sorted(information_centrality.items(), key=lambda t: t[1], reverse=True)
    IC_List = open(output + '\\IC_List' + str(k) + '.txt', 'w+')
    IC_List.write(json.dumps(information_centrality_List, indent=4, separators=(",", ":")))
    IC_List.close()
    return 0
def NSC(num,graph):
    edges = pd.read_csv(output + '\\' + 'edgesAfterEdgeReduction'+ str(num)  +'.csv', header=None)
    G = nx.Graph()
    G.add_edges_from(edges.values)
    AdjMatrix = nx.adjacency_matrix(G).toarray()
    # 检查图是否是连通的
    # if not nx.is_connected(G):
    #     # 获取最大连通子图
    #     largest_connected_component = max(nx.connected_components(G), key=len)
    #     G = G.subgraph(largest_connected_component).copy()
    # subgraph_centrality = dict(nx.information_centrality(G))
    # max_centrality = max(subgraph_centrality.values())
    #sc = {node: centrality / max_centrality for node, centrality in subgraph_centrality.items()}
    degree = [deg for node, deg in G.degree()]
    number = len(G.nodes)
    nodes = list(G.nodes)  # 创建节点列表
    NSC = pd.DataFrame(np.zeros((number, number)), columns=nodes, index=nodes)
    for i in range(number):
        for j in range(number):
            neighbornumber = 0
            if AdjMatrix[i, j] != 0:
                for k in range(number):
                    if k != j and k != i and AdjMatrix[i, k] != 0 and AdjMatrix[j, k] != 0:
                        # 使用节点列表的索引
                        # if nodes[i] in sc and nodes[j] in sc and nodes[k] in sc:
                            # neighbornumber += sc[nodes[i]] + sc[nodes[j]] + sc[nodes[k]]  # 计算邻居数量
                              neighbornumber += 1
                if degree[k]>1:
                        neighbornumber = neighbornumber/degree[k]
                if degree[i] > 0 and degree[j] > 0:
                    NSC.iloc[i, j] = neighbornumber /min(degree[i],degree[j])
    SONSC = NSC.sum(axis=1)  # 对ECC的每行求和
    NSC_List = SONSC / np.max(SONSC)
    NSC_List = NSC_List.to_dict()
    NSC_F = open(output + '\\NSC_List' + str(num) + '.txt', 'w+')
    NSC_F.write(json.dumps(NSC_List, indent=4, separators=(",", ":")))
    NSC_F.close()
    NSC_Listsorted = sorted(NSC_List.items(), key=lambda t: t[1], reverse=True)
    NSC_R = open(output + '\\NSC_Listsorted' + str(num) + '.txt', 'w+')
    NSC_R.write(json.dumps(NSC_Listsorted, indent=4, separators=(",", ":")))
    NSC_R.close()
    return NSC_List
def NSO(num,graph):
    edges = pd.read_csv(output + '\\' + 'edgesAfterEdgeReduction' + str(num) + '.csv', header=None)
    G = nx.Graph()
    # 你可以添加边到G中
    G.add_edges_from(edges.values)
    # 获取邻接矩阵
    SparseAdjMatrix = nx.adjacency_matrix(G).toarray()

    degree = [deg for node, deg in G.degree()]
    max_degree = max(degree)
    min_degree = min(degree)

    # 归一化处理
    nor_degree = [(deg - min_degree) / (max_degree - min_degree) for deg in degree]

    number = len(G.nodes)
    nodes = list(G.nodes)  # 创建节点列表
    ECC = pd.DataFrame(np.zeros((number, number)), columns=nodes, index=nodes)

    for i in range(number):
        for j in range(number):
            neighbornumber = 0
            if SparseAdjMatrix[i, j] != 0:
                for k in range(number):
                    neighborK = 0
                    if k != j and k != i and SparseAdjMatrix[i, k] != 0 and SparseAdjMatrix[j, k] != 0:
                        y = np.nonzero(SparseAdjMatrix[k, :])[0]
                        for q in range(sum(SparseAdjMatrix[k, :])):
                            for p in range(sum(SparseAdjMatrix[k, :])):
                                if y[q] != y[p] and SparseAdjMatrix[y[q], y[p]] != 0:
                                    neighborK += 1  # 统计共邻蛋白质存在的三角形数
                        neighbornumber += (neighborK)/((degree[k]+1)**0.5)
                if degree[i] > 1 and degree[j] > 1:
                    ECC.iloc[i, j] = neighbornumber / min(sum(SparseAdjMatrix[i, :]),sum(SparseAdjMatrix[j, :]))
    SOECC = ECC.sum(axis=1)  # 对ECC的每行求和
    ND = SOECC / max(SOECC)
    NSO_List = ND.to_dict()
    NSO_F = open(output + '\\NSO_List' + str(num) + '.txt', 'w+')
    NSO_F.write(json.dumps(NSO_List, indent=4, separators=(",", ":")))
    NSO_F.close()
    NSO_Listsorted = sorted(NSO_List.items(), key=lambda t: t[1], reverse=True)
    NSO_R = open(output + '\\NSO_Listsorted' + str(num) + '.txt', 'w+')
    NSO_R.write(json.dumps(NSO_Listsorted, indent=4, separators=(",", ":")))
    NSO_R.close()
    return NSO_List
def ECC(num,graph):
    #这里首先进行了一个特征提取；
    edges = pd.read_csv(output + '\\' + 'data1'  + '.csv',header=None)
    # with open(output + '\\' + 'reducedNodeAfterEdgeReduce' + str(k) + '.txt') as f:
    #     edges = f.readlines()
    # 生成一个邻接矩阵
    # G = nx.Graph()
    # for edge in edges:
    #     node1, node2 = edge.strip().split(',')  # 假设边是由空格分隔的节点对
    #     G.add_edge(node1, node2)
    # 生成一个邻接矩阵
    G = nx.Graph()
    G.add_edges_from(edges.values)
    AdjMatrix = nx.adjacency_matrix(G).toarray()
    Dist = dict(nx.all_pairs_shortest_path_length(G))

    betweenness = nx.betweenness_centrality(G)
    #计算介数中心性；
    number = G.number_of_nodes()
    # 将无穷大的距离替换为number
    for k1 in Dist.keys():
        for k2 in Dist[k1].keys():
            if Dist[k1][k2] == float(10):
                Dist[k1][k2] = number
    # 获取节点的度
    degree = [deg for node, deg in G.degree()]
    cc = [(number - 1) / sum(Dist[i].values()) for i in Dist.keys()]
    # 初始化ECC为边聚集系数矩阵
    number = len(G.nodes)
    ECC = pd.DataFrame(np.zeros((number, number)), columns=G.nodes, index=G.nodes)
    for i in range(number):
        for j in range(number):
            neighbornumber = 0
            if AdjMatrix[i, j] != 0:
                for k in range(number):
                    if k != j and k != i and AdjMatrix[i, k] != 0 and AdjMatrix[j, k] != 0:
                        neighbornumber += cc[i]+cc[j]+cc[k]   # 计算邻居数量
                if degree[i] > 1 and degree[j] > 1:
                    ECC.iloc[i, j] = neighbornumber / min(np.sum(AdjMatrix[i, :]) - 1, np.sum(AdjMatrix[j, :]) - 1)

    SOECC = ECC.sum(axis=1)  # 对ECC的每行求和
    ECC_List = SOECC / np.max(SOECC)
    ECC_List = ECC_List.to_dict()
    ECC_F = open(output + '\\ECC_List' + str(num) + '.txt', 'w+')
    ECC_F.write(json.dumps(ECC_List, indent=4, separators=(",", ":")))
    ECC_F.close()
    ECC_Listsorted = sorted(ECC_List.items(), key=lambda t: t[1], reverse=True)
    ECC_R = open(output + '\\ECC_Listsorted' + str(num) + '.txt', 'w+')
    ECC_R.write(json.dumps(ECC_Listsorted, indent=4, separators=(",", ":")))
    ECC_R.close()
    return ECC_List

def subcell():
    sub_dict = {}
    with open('input\\subcellular_score.txt', 'r') as file:
        for line in file:
            key, value = line.strip().split('\t')
            sub_dict[key] = float(value)
    return sub_dict
def IDC(k,graph):
    #剪枝后的图
    path2 = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))+'\\Input\\pComplex.xlsx'
    #pComplex.xlsx是一个什么文件？根据论文是蛋白质复合物度中心性；
    reducedProteinF=open(output+'\\reducedNodeAfterEdgeReduce'+str(k)+'.txt')
    uniqueProteinSet=json.load(reducedProteinF)
    xl = pd.ExcelFile(path2)
    #这一句主要是通过这个方法将excel转化为padas里面的pd类型的文件
#    #print(xl.sheet_names)
    s=[]
    IDCList={}
    #temp
    s1=[]
    # Load a sheet into a DataFrame by name: df1
    df1 = xl.parse('proteinComplex',header=None)
#    #print(df1)
    #print(type(uniqueProteinSet))
    for row in df1.itertuples(index=False):
    #这个里面的数据是itertuples是pandas.DataFrame类型的一个方法
    #它用于按行迭代DataFrame对象。最后返回一个迭代器，可以用来逐行访问DataFrame对象。
    #每次迭代，都会返回一个命名的元组，其中包含当前行中所有列的数据。
    #index=False 是一个可选参数，它指定了是否在返回的命名元组中包含索引列。
    #当 index=False 时，返回的命名元组中不包含索引列。
            s1.clear()
            for x in row:        
                if x==x and x in uniqueProteinSet:
                #uniqueProteinSet就是reducedNodeAfterEdgeReduce+str(k);剪枝过边和节点的图；
                    s1.append(x)
            s.append(s1)
            for node in s1:
    #            #print(node)
                s2=set(graph[node])
                IDC=len(s2.intersection(set(s1)))
    #这里就是计算s2和s1的交集的大小，并将其存储再变量IDC。
    #检查字典IDList中是否存在键为node的项。
    #如果存在，则将该项的值增加IDC。
    #如果不存在，则在字典中创建一个新项，其键为node,值为IDC。
                if node in IDCList:
                    IDCList[node]+=IDC
                else:
                    IDCList[node]=IDC
    #print(IDCList)
    IDC=open(output+'\\IDC'+str(k)+'.txt','w+')
    IDC.write(json.dumps(IDCList,indent=4 ,separators=(",",":")))
    IDC.close()
    IDCdict_Sorted = sorted(IDCList.items(), key=lambda t: t[1], reverse=True)
    IDC_dic = open(output + '\\IDCdict_Sorted' + str(k) + '.txt', 'w+')
    IDC_dic.write(json.dumps(IDCdict_Sorted, indent=4, separators=(",", ":")))
    IDC_dic.close()
    return IDCList
#%%
#LIDC
def nonlinear_func(x):
    return np.exp(-x)
def rank_essentialprotein(IDCList,k):
    #根据直系同源和拓扑特征以及蛋白质复合物来识别关键蛋白；
    ReducedEdgeGraph = open(output+'\\reducedEdgeGraph'+str(k)+'.txt')
    ReducedEdgeGraph = json.load(ReducedEdgeGraph)

    # mainGraph = open(output + '\\mainGraph' + '.txt')
    # mainGraph = json.load(mainGraph)

    orthologous_matrix=orthologous()
    #feature_pccmatrix=neighbour_affinity(k)
    #最后一个边的系数暂时不用，先用这个IDC和直系同源；
    #首先将其中的三个元素进行标准化，比如说IDC;
    # ECC_dit=ECC(k,ReducedEdgeGraph)
    NSO_dit=NSO(k,ReducedEdgeGraph)
    Subcell_dit= subcell()
    # # #print(Subcell_dit)
    Subcell_ditListvalues = list(Subcell_dit.values())
    Subcell_ditListvalues_ave = statistics.mean(Subcell_ditListvalues)
    Subcell_ditList_std = statistics.stdev(Subcell_ditListvalues)
    for key, values in Subcell_dit.items():
        Subcell_dit[key] = (Subcell_dit[key] - Subcell_ditListvalues_ave) / Subcell_ditList_std

    IDC_values = list(IDCList.values())
    IDC_ave=statistics.mean(IDC_values)
    IDC_std=statistics.stdev(IDC_values)
    for key,values in IDCList.items():
         IDCList[key] = (IDCList[key] - IDC_ave)/IDC_std
    #接下来就是计算每个节点的同源蛋白的均值常数；
    #orthologous_ave = orthologous_matrix[0].mean()
    #orthologous_std = orthologous_matrix[0].var()
    orthologous_list=orthologous_matrix[0].to_dict()
    # for key,values in orthologous_list.items():
    #      orthologous_list[key] = (orthologous_list[key] -orthologous_ave)/orthologous_std
    Subcell_IDC_list = {}
    orthologous_ECC_list = {}
    for node in ReducedEdgeGraph:
         Subcell_IDC_list[node] = 0
    for node in ReducedEdgeGraph:
        orthologous_ECC_list[node] = 0
    # 接下来的一个思路就是压制模型结合非线性模型进行一个综合的评分从而得到一个良好的打分模型。
    # 这个要根据reducedEdgeGraph1里面的每个蛋白质的边来进行打分。
    # 是否可以图迭代；根据求IDC-IDC or orthlolgous-orthlolgous的最小值，
    orthologous_IDC_list1 = {}

    score_list = {}

    for node in ReducedEdgeGraph:
        orthologous_IDC_list1[node] = 0

    for node in ReducedEdgeGraph:
        score_list[node] = 0

    for node in ReducedEdgeGraph:
             flag = 0
             flag1 = 0
             alpha = 0
             beta = 0
             for neighbour in ReducedEdgeGraph[node]:
                     if node not in IDCList.keys():
                         IDCList[node] = 0.005
                     if neighbour not in IDCList.keys():
                         IDCList[neighbour] = 0.005
                     if(Subcell_dit[node]>=Subcell_dit[neighbour]):
                       alpha += Subcell_dit[neighbour]
                       flag = flag+1
                     if(Subcell_dit[node]<Subcell_dit[neighbour] ):
                       alpha += Subcell_dit[node]
                       flag = flag+1
                     if(NSO_dit[node]>=NSO_dit[neighbour]):
                       beta +=NSO_dit[neighbour]
                       flag1 = flag1+1
                     if(NSO_dit[node]<NSO_dit[neighbour]):
                       beta += NSO_dit[node]
                       flag1 =flag1+1
             if flag!=0:
                    alpha = alpha/flag
                    beta = beta/flag1
                    w1 = nonlinear_func(alpha)
                    w2 = nonlinear_func(beta)
                    Subcell_IDC_list[node] = (w1 * (NSO_dit[node]) + w2 * (Subcell_dit[node]))

    for node in ReducedEdgeGraph:
        if node not in ReducedEdgeGraph:
            Subcell_IDC_list[node] = 0

    orth = list(Subcell_IDC_list.values())
    orth_ave = statistics.mean(orth)
    orth_std = statistics.stdev(orth)
    for key, values in Subcell_IDC_list.items():
        Subcell_IDC_list[key] = (Subcell_IDC_list[key] - orth_ave) / orth_std


    # for node in mainGraph:
    #     for neig in mainGraph:
    #                 if(ECC_dit[node]>ECC_dit[neig] and orthologous_list[node]>orthologous_list[neig] ):
    #                         orthologous_ECC_list[node] = orthologous_ECC_list[node]+(ECC_dit[node]-ECC_dit[neig])+(orthologous_list[node]-orthologous_list[neig])
    #                 if (ECC_dit[node] < ECC_dit[neig] and orthologous_list[node] < orthologous_list[neig]):
    #                         orthologous_ECC_list[neig] = orthologous_ECC_list[neig]+(-ECC_dit[node]+ECC_dit[neig])+(-orthologous_list[node]+orthologous_list[neig])
    for node in ReducedEdgeGraph:
        for neig in ReducedEdgeGraph:
            # if node in NSC_dit and neig in NSC_dit and node in orthologous_list and neig in orthologous_list:
                if(orthologous_list[node]+NSO_dit[node]+Subcell_dit[node]> orthologous_list[neig]+NSO_dit[neig]+Subcell_dit[neig] ):
                    if( orthologous_list[node]>orthologous_list[neig] and NSO_dit[node]>NSO_dit[neig] and Subcell_dit[node]>Subcell_dit[neig]):
                        orthologous_ECC_list[node] = orthologous_ECC_list[node]+(orthologous_list[node]-orthologous_list[neig])*(NSO_dit[node]-NSO_dit[neig]+Subcell_dit[node]-Subcell_dit[neig])
                        score_list[node] =  (orthologous_list[node]-orthologous_list[neig]+NSO_dit[node]-NSO_dit[neig]+Subcell_dit[node]-Subcell_dit[neig])
                    elif( orthologous_list[node]>orthologous_list[neig] and NSO_dit[node]>NSO_dit[neig]):
                        orthologous_ECC_list[node] = orthologous_ECC_list[node]+
                        score_list[node] =  (orthologous_list[node]-orthologous_list[neig] + NSO_dit[node] - NSO_dit[neig])
                    elif(orthologous_list[node]>orthologous_list[neig] and Subcell_dit[node]>Subcell_dit[neig]):
                        orthologous_ECC_list[node] = (orthologous_ECC_list[node]+2)
                        score_list[node] = (orthologous_list[node] - orthologous_list[neig] + Subcell_dit[node]-Subcell_dit[neig])
                    elif(Subcell_dit[node]>Subcell_dit[neig] and NSO_dit[node]>NSO_dit[neig]):
                        orthologous_ECC_list[node] = (orthologous_ECC_list[node]+2)
                        score_list[node] = (Subcell_dit[node] - Subcell_dit[neig] + NSO_dit[node] -NSO_dit[neig])
                    elif(orthologous_list[node] > orthologous_list[neig]):
                        orthologous_ECC_list[node] = (orthologous_ECC_list[node] + 1)
                        score_list[node] = (orthologous_list[node] - orthologous_list[neig])
                    elif(NSO_dit[node]>NSO_dit[neig]):
                        orthologous_ECC_list[node] = (orthologous_ECC_list[node] + 1)
                        score_list[node] = (NSO_dit[node] - NSO_dit[neig])
                    elif(Subcell_dit[node]>Subcell_dit[neig]):
                        orthologous_ECC_list[node] = (orthologous_ECC_list[node] + 1)
                        score_list[node] = (Subcell_dit[node] - Subcell_dit[neig])
    # # OR_IDC_Listvalues = list(orthologous_IDC_list.values())
    # OR_IDC_Listvalues_ave = statistics.mean(OR_IDC_Listvalues)
    # OR_IDC_Liststd = statistics.stdev(OR_IDC_Listvalues)
    # for key, values in orthologous_IDC_list.items():
    #     orthologous_IDC_list[key] = (orthologous_IDC_list[key] - OR_IDC_Listvalues_ave) / OR_IDC_Liststd
    # orthologous_ECC_list_max = max(orthologous_ECC_list.values())
    # orthologous_ECC_list_min = min(orthologous_ECC_list.values())
    # for node in ReducedEdgeGraph:
    #     orthologous_ECC_list[node] = (orthologous_ECC_list[node] - orthologous_ECC_list_min)/(orthologous_ECC_list_max-orthologous_ECC_list_min)
    for node in ReducedEdgeGraph:
        orthologous_IDC_list1[node] = orthologous_ECC_list[node]

             # (orthologous_ECC_list[node])
                                      # *Subcell_IDC_list[node]

    # for node in ReducedEdgeGraph:
    #     # if node in orthologous_ECC_list:
    #         orthologous_IDC_list1[node] = orthologous_ECC_list[node]*Subcell_IDC_list[node]

                # (ECC_dit[node]-ECC_dit[neig])+(orthologous_list[node]-orthologous_list[neig])
                                #(-ECC_dit[node]+ECC_dit[neig])+-(orthologous_list[node]+orthologous_list[neig])
    # PF = {key:0 for key in ECC_dit}
    # for node in ReducedEdgeGraph:
    #         for neighbour in ReducedEdgeGraph[node]:
    #             if(node in ECC_dit.keys() and neighbour in ECC_dit.keys() and node in Subcell_dit.keys() and neighbour in Subcell_dit.keys()):
    #                 if(ECC_dit[node]>=ECC_dit[neighbour] or Subcell_dit[node]>=Subcell_dit[neighbour]):
    #                     PF[node] += (ECC_dit[node]-ECC_dit[neighbour]) + (Subcell_dit[node]-Subcell_dit[neighbour])
    #                     # orthologous_IDC_list[node] += (ECC_dit[node] + Subcell_dit[node]) * orthologous_list[node]
    #                 else:
    #                     PF[neighbour] += (ECC_dit[neighbour]-ECC_dit[node]) + (Subcell_dit[neighbour]-Subcell_dit[node])
                        # orthologous_IDC_list[neighbour] += (ECC_dit[neighbour] + Subcell_dit[neighbour])*orthologous_list[neighbour]
    # for node in ReducedEdgeGraph:
    #         flag = 0
    #         flag1 = 0
    #         alpha = 0
    #         beta = 0
    #         for neighbour in ReducedEdgeGraph[node]:
    #             if(node in orthologous_list.keys() and neighbour in orthologous_list.keys() and node in IDCList.keys() and neighbour in IDCList.keys()):
    #                 if(orthologous_list[node]>orthologous_list[neighbour] and IDCList[node]>IDCList[neighbour]):
    #                    alpha += (orthologous_list[node]-orthologous_list[neighbour])+(IDCList[node]-IDCList[neighbour])
    #                    flag = flag+1
    #                 if IDCList[node] < IDCList[neighbour] and orthologous_list[node] < orthologous_list[neighbour]:
    #                     beta += (orthologous_list[neighbour] - orthologous_list[node]) + (IDCList[neighbour] - IDCList[node])
    #                     flag = flag + 1
    #             if flag!=0:
    #                 alpha = alpha/flag
    #                 beta = beta/flag1
    #                 w1 = nonlinear_func(alpha)
    #                 w2 = nonlinear_func(beta)
    #                 orthologous_IDC_list[node] = w1 * (IDCList[node]) + w2 * (orthologous_list[node])
    # number = len(IDCList)
    # PF = {key: 0 for key in IDCList.keys()}
    # keys = list(IDCList.keys())
    # Subcell_dit = subcell()
    # #print(Subcell_dit)
    # Subcell_ditListvalues = list(Subcell_dit.values())
    # Subcell_ditListvalues_ave = statistics.mean(Subcell_ditListvalues)
    # Subcell_ditList_std = statistics.stdev(Subcell_ditListvalues)
    # for key, values in Subcell_dit.items():
    #     Subcell_dit[key] = (Subcell_dit[key] - Subcell_ditListvalues_ave) / Subcell_ditList_std

    # for i in range(number):
    #     for j in range(i + 1, number):
    #             if IDCList[keys[i]] > IDCList[keys[j]] and orthologous_list[keys[i]] > orthologous_list[keys[j]]:
    #                 PF[keys[i]] += (IDCList[keys[i]]-IDCList[keys[j]]) + (orthologous_list[keys[i]]-orthologous_IDC_list[keys[j]])
    #             if IDCList[keys[i]] < IDCList[keys[j]] and orthologous_list[keys[i]] < orthologous_list[keys[j]]:
    #                 PF[keys[j]] += (IDCList[keys[j]]-IDCList[keys[i]]) + (orthologous_list[keys[j]]-orthologous_IDC_list[keys[i]])
    #print(orthologous_IDC_list1)
    #ReducedEdgeGraph.close()
    # for key,values in IDCList.items():
    #      if (key in orthologous_list.keys()):
    #          alpha = 0.3
    #          w1 = nonlinear_func(alpha)
    #          w2 = nonlinear_func(1 - alpha)
    #          orthologous_IDC_list[key] = w1 * (IDCList[key]) + w2 * (orthologous_list[key])

    orthologous_IDC_ListSorted = sorted(orthologous_IDC_list1.items(), key=lambda t: t[1], reverse=True)
    orthologous_IDC = open(output + '\\orthologous_IDC' + str(k) + '.txt', 'w+')
    orthologous_IDC.write(json.dumps(orthologous_IDC_ListSorted, indent=4, separators=(",", ":")))
    orthologous_IDC.close()
    print('k= ', k, ' finished')
    return orthologous_IDC_ListSorted

# def LIDC(LIDList,IDCList,k):
#     rank=1
#     N=len(LIDList)
#     LIDCList={}
#     for tuple in LIDList:
#         #print(type(tuple))
#         LID=tuple[1]
#         #print(tuple[0])
#         if(tuple[0] in IDCList):
#             IDC=IDCList[tuple[0]]
#         else:
#             IDC=0
#         #print('LID: ',LID,'IDC: ',IDC,'Rank: ',rank,'N: ',rank/N)
#         LIDC=(LID*(1-(rank/N)))+(IDC*(rank/N))
#         #print('LIDC: ',LIDC)
#         LIDCList[tuple[0]]=LIDC
#         rank+=1
#     #print(LIDCList)
#     LIDCListSorted=sorted(LIDCList.items(), key=lambda t: t[1],reverse=True)
#     LIDC=open(output+'\\LIDC'+str(k)+'.txt','w+')
#     LIDC.write(json.dumps(LIDCListSorted,indent=4 ,separators=(",",":")))
#     LIDC.close()
#     print('k= ',k,' finished')

def optimality_method(ECCDit,OR_IDCList,k):
#这里的ECCList是字典，而OR_IDCList确实List:
    OR_IDCListvalues = list(OR_IDCList.values())
    OR_IDCListvalues_ave=statistics.mean(OR_IDCListvalues)
    OR_IDC_std=statistics.stdev(OR_IDCListvalues)
    for key,values in OR_IDCList.items():
         OR_IDCList[key] = (OR_IDCList[key] - OR_IDCListvalues_ave)/OR_IDC_std
    #print(OR_IDCList)
    ECCList = list(ECCDit.items())
    ECCList_maxtuple = max(ECCDit.values())
    ECCList_mintuple = min(ECCDit.values())
    OR_IDCList_maxtuple = max(OR_IDCListvalues)
    OR_IDCList_mintuple = min(OR_IDCListvalues)
    # print(OR_IDCList_mintuple)
    # print(OR_IDCList_maxtuple)
    # print(ECCList_mintuple)
    # print(ECCList_maxtuple)
    PF = {key:0 for key,_ in ECCList}
    ECC_dict = dict(ECCList)
    OR_IDCList_dict = dict(OR_IDCList)
# 遍历数据
#     for i in range(len(ECCList)):
#          for j in range(i+1, len(ECCList)):
#              if NCC_dict[ECCList[i][0]] >= NCC_dict[ECCList[j][0]] and OS_dict[OR_IDCList[i][0]] >= OS_dict[OR_IDCList[j][0]]:
#                  PF[ECCList[i][0]] + = (NCC_dict[ECCList[i][0]]-NCC_dict[ECCList[j][0]])+(OS_dict[OR_IDCList[i][0]]-OS_dict[OR_IDCList[j][0]])
#              if NCC_dict[ECCList[i][0]] < NCC_dict[ECCList[j][0]] and OS_dict[OR_IDCList[i][0]] < OS_dict[OR_IDCList[j][0]]:
#                  PF[ECCList[i][0]] + =
#     for i in range(len(ECCList)):
#         PF[ECCList[i][0]]=ECCDit[ECCList[i][0]]*OR_IDCList_dict[OR_IDCList[i][0]]
    keys = list(ECCDit.keys())
    number=len(ECCDit)
    for i in range(number):
        for j in range(i+1, number):
            if keys[i] in ECCDit and keys[i] in OR_IDCList_dict and keys[j] in ECCDit and keys[j] in OR_IDCList_dict:
                if ECCDit[keys[i]] > ECCDit[keys[j]] and OR_IDCList_dict[keys[i]] > OR_IDCList_dict[keys[j]]:
                    PF[keys[i]] += ((ECCDit[keys[i]]))/(ECCList_maxtuple-ECCList_mintuple)+((OR_IDCList_dict[keys[i]]))/(OR_IDCList_maxtuple-OR_IDCList_mintuple)
                if ECCDit[keys[i]] < ECCDit[keys[j]] and OR_IDCList_dict[keys[i]] < OR_IDCList_dict[keys[j]]:
                    PF[keys[j]] += ((ECCDit[keys[j]]))/(ECCList_maxtuple-ECCList_mintuple)+((OR_IDCList_dict[keys[j]]))/(OR_IDCList_maxtuple-OR_IDCList_mintuple)
    PF_ListSorted = sorted(PF.items(), key=lambda t: t[1], reverse=True)
    #rint(PF_ListSorted)
    PF_ListSorted_F = open(output + '\\PF_ListSorted' + str(k) + '.txt', 'w+')
    PF_ListSorted_F.write(json.dumps(PF_ListSorted, indent=4, separators=(",", ":")))
    PF_ListSorted_F.close()
    return PF_ListSorted
    #print(OR_IDCList)
def jackknife(k):
    TP = FP = TN = FN = 0
    essentialProteinData = []
    essentialProtein_Score = []
    essentialProtein_pr = []
    path2 = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0))) + '\\Input\\essential_non_essential_data.xls'
    xl = pd.ExcelFile(path2)
    # Load a sheet into a DataFrame by name: df1
    df1 = xl.parse('essential_proteins', header=None)
    # print(df1)
    for i in df1.itertuples(index=False):
        essentialProteinData.append(i[0])
    x = np.linspace(0, 600, 600)
    #这里的LIDC的值已经经过排序了；
    with open(output + '\\orthologous_IDC'  +str(k)+ '.txt', 'r') as file:
        data = json.load(file)
    sum = 0
        #这个是折刀法来进行均值的折刀分析方法来进行比较；
    flag = -1
    for i,row in enumerate(data):
        flag +=1
        if row[0] in essentialProteinData:
            sum+=1
        if flag<5094:
            essentialProtein_Score.append(sum)
            pr = sum/(flag+1)
            essentialProtein_pr.append(pr)
#    plt.plot(x,essentialProtein_Score)
#    plt.xlabel('Protein_nums')
#    plt.ylabel('Essential_protein')
    #plt.show()
    #画一个ROC曲线；
    Jaknife = open(output + '\\jaknife' + str(k) + '.txt', 'w+')
    Jaknife.write(json.dumps(essentialProtein_Score, indent=4, separators=(",", ":")))
    Jaknife.close()

    Pr = open(output + '\\pr_curve' + str(k) + '.txt', 'w+')
    Pr.write(json.dumps(essentialProtein_pr, indent=4, separators=(",", ":")))
    Pr.close()
def ml():

    return 0
def LIDCDriver(k,graph):
    IDCList=IDC(k,graph)
    #LIDList=LID(k,graph)
    #LIDC(LIDList,IDCList,k)
    # ECCDit = ECC(k,graph)
    # NSC(k,graph)
    #print(graph)

    # DC(k,graph)
    # BC(k,graph)
    # CC(k,graph)
    # EC(k,graph)
    # SC(k,graph)
    # NC(k,graph)
    # OR_IDCList = dict(OR_IDCList)
    #IC(k,graph)
    #subcell()
    NSO(k,graph) #
    rank_essentialprotein(IDCList,k)
    #optimality_method(ECCDit,OR_IDCList,k)
    jackknife(k)



def mainDriver():
    nodeWeight(data)
mainDriver()

rs.accuracy()

