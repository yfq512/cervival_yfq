import os,time,shutil
import cv2
import numpy as np
import matplotlib.pyplot as plt
import copy
import basefun as bf
from tqdm import tqdm
 
def get_cell_nuclei_mask(img, img_org, value):
    ret_, thresh = cv2.threshold(img, value, 255, cv2.THRESH_BINARY_INV)
    image_, contours, hierarchy_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if len(contours) == 0:
        return 0,0,0
    area_ = 0
    contours_max = []
    for cnt in range(0, len(contours)):
        area = cv2.contourArea(contours[cnt])
        if area > area_:
            contours_max = contours[cnt]
            area_ = area
    if not type(contours_max) == np.ndarray:
        return 0,0,0
    thresh_ = cv2.fillPoly(thresh, [contours_max], 255)
    img_temp = img * (thresh_/255)
    value_temp = img_temp.mean().mean()
    if value_temp < 11:
        return 0,0,0
    return 1, thresh_, area_
    
    
def get_cell_cytoplasm_mask(img, nuclei_mask, img_org, value):
    ret_, thresh = cv2.threshold(img, value, 255, cv2.THRESH_BINARY_INV)
    image_, contours, hierarchy_ = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)
    nuclei_mask_off = np.ones(nuclei_mask.shape)*255 - nuclei_mask
    thresh___ = thresh*(nuclei_mask_off/255)
    return thresh___

def main2():
    dstroot = 'clusters'
    listcells = os.listdir(dstroot)
    cellsinfo = []
    for n,i in zip(listcells,tqdm(range(len(listcells)))):
        cellinfo = {}
        cellpath = os.path.join(dstroot, n)
        img = cv2.imread(cellpath)
        img_gray = cv2.imread(cellpath, 0)
        img_gray = bf.get_fit_img(img_gray)
     #   cv2.imwrite(cellpath+'abc0.png',img_gray)
        value1_,value_2 = bf.get_2value(img_gray)
        value_1,_ = bf.get_2value(img_gray, grien = 3)
        sign_, cell_nuclei_mask1,area1 = get_cell_nuclei_mask(img_gray, img, value_1) #获取细胞核掩码、个数、面积、凸面积、周长、核形规则度、获取细胞>核深染程度
        sign_,cell_nuclei_mask2,area2 = get_cell_nuclei_mask(img_gray, img, value1_) #获取细胞核掩码、个数、面积、凸面积、周长、核形规则度、获取细胞核深染
        if sign_ == 0:
            continue
    #    cv2.imwrite(cellpath+'abc1.png',cell_nuclei_mask2)
        temp = area1/area2
        temp1 = str(temp)[0:5]
        if area1/area2 > 0.85:
            newpath = os.path.join('ab_clusters', n+'_'+temp1+'.png')
            shutil.copy(cellpath, newpath)


if __name__ == "__main__":
    dstroot = 'clusters'
    listcells = os.listdir(dstroot)
    cellsinfo = []
    for n,i in zip(listcells,tqdm(range(len(listcells)))):
        cellinfo = {}
        cellpath = os.path.join(dstroot, n)
        img = cv2.imread(cellpath)
        img_gray = cv2.imread(cellpath, 0)
        img_gray = bf.get_fit_img(img_gray)
     #   cv2.imwrite(cellpath+'abc0.png',img_gray)
        value1_,value_2 = bf.get_2value(img_gray)
        value_1,_ = bf.get_2value(img_gray, grien = 3)
        sign_, cell_nuclei_mask1,area1 = get_cell_nuclei_mask(img_gray, img, value_1) #获取细胞核掩码、个数、面积、凸面积、周长、核形规则度、获取细胞>核深染程度
        sign_,cell_nuclei_mask2,area2 = get_cell_nuclei_mask(img_gray, img, value1_) #获取细胞核掩码、个数、面积、凸面积、周长、核形规则度、获取细胞核深染
        if sign_ == 0:
            continue
    #    cv2.imwrite(cellpath+'abc1.png',cell_nuclei_mask2)
        temp = area1/area2
        temp1 = str(temp)[0:5]
        if area1/area2 > 0.85:
            newpath = os.path.join('ab_clusters', n+'_'+temp1+'.png')
            shutil.copy(cellpath, newpath)
