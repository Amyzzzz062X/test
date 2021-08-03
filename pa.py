# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 18:38:51 2021

@author: Amy
"""
import re
import requests
from bs4 import BeautifulSoup
def getHTML(url):
    try:
        r = requests.get(url,timeout=30)
        r.raise_for_status()
        r.encoding = 'utf-8'
        return r.text
    except:
        return ""
 
def getContent(url):
    html = getHTML(url)
    soup = BeautifulSoup(html,'html.parser')
    paras_tmp = soup.select('p')
    paras = paras_tmp[3:]
    return paras
 
def saveFile(text):
    f=open('novel.txt','w+',encoding="utf-8")
    for t in text:
        if len(t) > 0:
            f.writelines(t.get_text() + "\n")
    f.close()
    

url_1 = 'https://hznews.hangzhou.com.cn/kejiao/content/2021-07/13/content_8006604.htm'
text_1 = getContent(url_1)
print(text_1)
saveFile(text_1)

#url_2 = 'http://www.xinhuanet.com/politics/leaders/2021-07/13/c_1127651317.htm'
#text_2 = getContent(url_2)
#saveFile(text_2)