# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 09:07:29 2021

@author: Amy
"""

import requests
import bs4
import urllib.request

url_3 = 'https://www.163.com/news/article/GERTV2AC00019K82.html'

def get_img_text(url):
    response = requests.get(url)
    if response.status_code != 200 :
        print (response.status_code)
    response.encoding = 'utf-8'
    soup = bs4.BeautifulSoup(response.text, "html.parser")
    
    # get images
    imgs = soup.find_all('img')
    x = 1
    for i in imgs:
        # 获取src路径
        imgsrc = i.get('src')
        if imgsrc.startswith('https://nimg.ws.126.net/'):
            # 本地路径
            filename = 'C:/Users/Amy/Desktop/news/%s.jpg'%x
            # 将URL表示的网络对象复制到本地文件
            urllib.request.urlretrieve(imgsrc, filename)
            print('下载第%d张' % x)
            x += 1
    print('**下载完成!**')
    
    # get text
    paras_tmp = soup.select('p')
    paras = paras_tmp[3:]
    f=open('news.txt','w+',encoding="utf-8")
    for t in paras:
        if len(t) > 0:
            f.writelines(t.get_text() + "\n\n")
    f.close()

get_img_text(url_3)


