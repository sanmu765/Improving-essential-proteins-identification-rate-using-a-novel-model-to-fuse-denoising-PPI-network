
import json

from inspect import getsourcefile
import os
import pandas as pd
path = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))
#获取当前文件所在的绝对路径；getsourcefile(lambda:0)会返回一个字符串，表示当前执行文件的路径，
#os.path.abspath会将这个路径转化为绝对路径。最后，os.dirname会返回这个绝对路径的目录名，也就是
#只返回这个绝对路径中的目录名，去除路径中的文件名，只返回文件所在的目录路径。
data='data'       # data is the data set file name without extension
          # graph is the variable to store the network usin dictionary

output=path+'\\Output'
input=path+'\\Input'
def uniqueProtein(filePath):
    print('in')
    uniqueProteinF=open(input+'\\'+filePath+'.txt')
    s1=set()
    #set是一种空的集合，也就是创建一个空的集合并将他赋值给s1;
    #set一种内置的数据类型，用于存储不重复的元素，并且可以进行并集，交集运算；
    while True:
        a=uniqueProteinF.readline()
        #读取一行的内容；
        if a == '':
            break
        b=a.split()
        #可以指定分割，里面什么都没有，就是按照换行符号进行分割。
        s1.add(b[0])
        s1.add(b[1])
        #add就是向集合里面添加一个元素；
    ##print(s1)
    uniqueProteinList=list(s1)
    #uniqueProteinList.sort()
    print('Unique protein list: ',uniqueProteinList)
    uniqueProteinF.close()
    return uniqueProteinList
L=uniqueProtein('data')
df=pd.DataFrame(L)
#DataFrame是一种表格型的数据结构。
df.to_csv(output+'\\uniqueProtein'+'.csv',sep=',',header=None,index=None)
file2=open(output+'\\'+'uniqueProtein'+'.txt','w+')
file2.write(json.dumps(L,indent=4, separators=("\n",":")))
file2.close()
