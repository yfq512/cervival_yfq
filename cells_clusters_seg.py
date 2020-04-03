import os,time
import cv2
import numpy as np
import matplotlib.pyplot as plt
import copy
import basefun as bf
from tqdm import tqdm
 
# 获得细胞坐标
def get_cells(img, imgpath, value, fit, limit_up, limit_down, side=5):
    ret_, thresh = cv2.threshold(img, value, 255, cv2.THRESH_BINARY_INV)
   # thresh = get_img_open(thresh,fit*0.1) #过滤杂质
    #cv2.imwrite(imgpath+'abc2.png', thresh)
    image_, contours, hierarchy_ = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
    cells_xoy = []
    for cnt in range(0,len(contours,)):
        area = cv2.contourArea(contours[cnt])
        x,y,w,h=cv2.boundingRect(contours[cnt])
        [x_,y_] = img.shape
        x1 = max(x-side,0) #向外扩展5个像素
        x2 = min(x+w+side,y_)
        y1 = max(y-side,0)
        y2 = min(y+h+side,x_)
        cell_xoy = [x1,y1,x2,y2]
        if area > limit_down and area < limit_up and w/h < 2 and w/h >1/2: #默认细胞长宽比在0.5~2之间浮动
            cells_xoy.append(cell_xoy)
    return cells_xoy
 
# 获得细胞团坐标
def get_clusters(img, value, fit, limit_down, side=5):
    ret_, thresh = cv2.threshold(img, value, 255, cv2.THRESH_BINARY_INV)
    thresh = bf.get_img_close(thresh,fit) #将密集的点连接成片
    fit2 = fit*0.1
    thresh = bf.get_img_open(thresh,fit2) #过滤孤立的点
    image_, contours, hierarchy_ = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
    clusters_xoy = []
    for cnt in range(0,len(contours,)):
        area = cv2.contourArea(contours[cnt])
        x,y,w,h=cv2.boundingRect(contours[cnt])
        [x_,y_] = img.shape
        x1 = max(x-side,0)
        x2 = min(x+w+side,x_)
        y1 = max(y-side,0)
        y2 = min(y+h+side,y_)
        cluster_xoy = [x1,y1,x2,y2]
        if area > limit_down:
            clusters_xoy.append(cluster_xoy)
    return clusters_xoy
 
def cells_save_sign(cell): #判断是否为细胞
    pass
 
def crop_from_fov(img, imgname, xoy, root_ = 'crop'): #扣细胞
    for n in xoy:
        [x1,y1,x2,y2] = n
        cropped = img[y1:y2,x1:x2,:]
        croppath = os.path.join(root_, imgname[:-4]+'_'+str(x1)+'_'+str(y1)+'_'+str(x2)+'_'+str(y2)+'.png')
        cv2.imwrite(croppath, cropped)
 
def plot_on_fov(img, imgpath, xoy):
    for n in xoy:
        [x1,y1,x2,y2] = n
        img =  cv2.rectangle(img, (x1,y1), (x2,y2), (0,255,0),2)
    #cv2.imwrite(imgpath+'abc.png', img)
  
def main2():
    dstroot = 'fovs'
    list_dst = os.listdir(dstroot)
    cnt = 0
    for n,i in zip(list_dst,tqdm(range(len(list_dst)))):
        imgpath = os.path.join(dstroot, n)
        #print(imgpath)
        img_org = cv2.imread(imgpath)
        img_gray = cv2.imread(imgpath,0)
        img_gray_fit = bf.get_fit_img(img_gray)
        value_1,value_2 = bf.get_2value(img_gray_fit) # 获取细胞核、背景阈值
        kernel_1 = np.ones((50,50),np.uint8)
        kernel_2 = np.ones((101,101),np.uint8)
        cells_xoy = get_cells(img_gray_fit, imgpath, value_2+10, kernel_1, 50000, 1600, side=5)
        clusters_xoy = get_clusters(img_gray_fit, value_1, kernel_2, 50000, side=5)
        crop_from_fov(img_org, n, cells_xoy)
        plot_on_fov(img_org, imgpath, cells_xoy)
        cnt = cnt + 1

#if __name__ == "__main__":
#    dstroot = 'fovs'
#    list_dst = os.listdir(dstroot)
#    for n in list_dst:
#        imgpath = os.path.join(dstroot, n)
#        print(imgpath)
#        img_org = cv2.imread(imgpath)
#        img_gray = cv2.imread(imgpath,0)
#        img_gray_fit = bf.get_fit_img(img_gray)
#        value_1,value_2 = bf.get_2value(img_gray_fit) # 获取细胞核、背景阈值
#        kernel_1 = np.ones((50,50),np.uint8)
#        kernel_2 = np.ones((101,101),np.uint8)
#        cells_xoy = get_cells(img_gray_fit, imgpath, value_2+10, kernel_1, 50000, 1600, side=5)
#        clusters_xoy = get_clusters(img_gray_fit, value_1, kernel_2, 50000, side=5)
#        crop_from_fov(img_org, n, cells_xoy)
#        plot_on_fov(img_org, imgpath, cells_xoy)
