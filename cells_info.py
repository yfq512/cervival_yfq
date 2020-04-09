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
    [_h, _w] = img.shape
    xo = _w/2
    yo = _h/2
    long_ = 200
    contours_best = []
    area_max = 0
    _long = 0
    if len(contours) == 0:
        return 0, 0, 0, 0, 0, 0, 0, 0
    for cnt in range(0, len(contours)):
        x,y,w,h=cv2.boundingRect(contours[cnt])
        xo2 = x+w/2
        yo2 = y+h/2
        _long = ((xo-xo2)**2+(yo-yo2)**2)**(1/2)
        if _long < long_:
            contours_best = contours[cnt]
            long_ = _long
    [w_img,h_img] = thresh.shape
    if (0 in contours_best) or ((w_img-1) in contours_best) or ((h_img-1) in contours_best):
        return 0, 0, 0, 0, 0, 0, 0, 0
    if not type(contours_best) == np.ndarray:
        return 0, 0, 0, 0, 0, 0, 0, 0
    area = cv2.contourArea(contours_best)
    perimeter = cv2.arcLength(contours_best,True)
    hull = cv2.convexHull(contours_best)
    area2 = cv2.contourArea(hull)
    
    if area == 0:
        return 0, 0, 0, 0, 0, 0, 0, 0
    rule = area/area2
    thresh_ = thresh*0
    thresh__ = cv2.fillPoly(thresh_, [contours_best], 255)
    thresh__ = cv2.polylines(thresh__, [hull], True, 255, 1)
    img_b = img_org[:,:,2]*(thresh__/255)
    value_1 = bf.get_img_mean_value(img_b)
    return 1, thresh__, cnt+1, area, area2, perimeter, rule, value_1
    
    
def get_cell_cytoplasm_mask(img, nuclei_mask, nuclei_area, img_org, value):
    ret_, thresh = cv2.threshold(img, value, 255, cv2.THRESH_BINARY_INV)
    image_, contours, hierarchy_ = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)
    if len(contours) == 0:
        return 0,0,0,0,0,0,0,0
    area_ = 0
    contours_max = []
    for cnt in range(0, len(contours)):
        area = cv2.contourArea(contours[cnt])
        if area > area_:
            contours_max = contours[cnt]
            area_ = area
    if not type(contours_max) == np.ndarray:
        return 0,0,0,0,0,0,0,0
    hull = cv2.convexHull(contours_max)
    area2 = cv2.contourArea(hull)
    perimeter = cv2.arcLength(contours_max,True)
    if area_ == 0:
        return 0,0,0,0,0,0,0,0
    rule = area_/area2
    thresh_ = thresh*0
    thresh__ = cv2.fillPoly(thresh_, [contours_max], 255)
    thresh__1 = thresh__ * (nuclei_mask/255)
    _area = sum(sum(thresh__1))
    if _area == 0:
        return 0,0,0,0,0,0,0,0
    [w_img,h_img] = thresh__.shape
    if thresh__[0,0] == 255 or thresh__[w_img-1,0] == 255 or thresh__[0,h_img-1] == 255 or thresh__[w_img-1,h_img-1] == 255:
        return 0,0,0,0,0,0,0,0
    thresh__ = cv2.polylines(thresh__, [hull], True, 255, 1)
    nuclei_mask_off = np.ones(nuclei_mask.shape)*255 - nuclei_mask
    thresh___ = thresh__*(nuclei_mask_off/255)
    img_r = img_org[:,:,0]*(thresh___/255)
    value_1 = bf.get_img_mean_value(img_r)
    img_cytoplasm = img*(thresh___/255)
    cytoplasm_var = img_cytoplasm.var()
    return 1, thresh___, area_, area2, perimeter, rule, value_1, cytoplasm_var

def get_cell_save_sign(cellinfo):
    nuclei_area = []
    nuclei_hull_area = []
    cytoplasm_area = []
    cytoplasm_hull_area = []
    nuclei_area = cellinfo['nuclei_area']
    nuclei_hull_area = cellinfo['nuclei_hull_area']
    cytoplasm_area = cellinfo['cytoplasm_area']
    cytoplasm_hull_area = cellinfo['cytoplasm_hull_area']
    nuclei_va = nuclei_area/nuclei_hull_area
    cytoplasm_va = cytoplasm_area/cytoplasm_hull_area
    area_diff = cytoplasm_area - nuclei_area
    if nuclei_va < 0.95 or cytoplasm_va < 0.85 or area_diff <150:
        sign = 0
    else:
        sign = 1
    return sign

def main2():
    dstroot = 'crop'
    listcells = os.listdir(dstroot)
    cellsinfo = []
    for n,i in zip(listcells,tqdm(range(len(listcells)))):
        cellinfo = {}
        cellpath = os.path.join(dstroot, n)
        img = cv2.imread(cellpath)
        img_gray = cv2.imread(cellpath, 0)
        img_gray = bf.get_fit_img(img_gray)
     #   cv2.imwrite(cellpath+'abc0.png',img_gray)
        value_1,value_2 = bf.get_2value(img_gray)
        sign_nuclei, cell_nuclei_mask, nuclei_cnt, nuclei_area, nuclei_hull_area, nuclei_circ, nuclei_rule, cell_nuclei_value = get_cell_nuclei_mask(img_gray, img, value_1) #获取细胞核掩码、个数、面积、凸面积、周长、核形规则度、获取细胞核深染程度
     #   cv2.imwrite(cellpath+'abc1.png',cell_nuclei_mask)
        if sign_nuclei == 0:
            continue
        sign_cytoplasm, cell_cytoplasm_mask, cytoplasm_area, cytoplasm_hull_area, cytoplasm_circ, cytoplasm_rule, cell_cytoplasm_value, cytoplasm_var = get_cell_cytoplasm_mask(img_gray, cell_nuclei_mask, nuclei_area,img, value_2) #获取细胞质掩码、面积、凸面积、周长、细胞规则度、获取细胞质情况
     #   cv2.imwrite(cellpath+'abc2.png',cell_cytoplasm_mask)
        if sign_cytoplasm == 0 or (cytoplasm_area-nuclei_area) == 0:
            continue
        cell_N_C = nuclei_area/(cytoplasm_area-nuclei_area) #计算核质比
        cellinfo_keys = ['cellpath','nuclei_cnt','nuclei_area','nuclei_hull_area','nuclei_circ','nuclei_rule','cell_nuclei_value','cytoplasm_area','cytoplasm_hull_area','cytoplasm_circ','cytoplasm_rule','cell_cytoplasm_value','cytoplasm_var','cell_N_C']
        cellinfo_values = [cellpath,nuclei_cnt,nuclei_area,nuclei_hull_area,nuclei_circ,nuclei_rule,cell_nuclei_value,cytoplasm_area,cytoplasm_hull_area,cytoplasm_circ,cytoplasm_rule,cell_cytoplasm_value,cytoplasm_var,cell_N_C]
        cellinfo = dict(zip(cellinfo_keys, cellinfo_values))
        cell_save_sign = get_cell_save_sign(cellinfo)
        if cell_save_sign == 0:
            continue
        cellsinfo.append(cellinfo)
        newpath = os.path.join('valid', n)
        shutil.copy(cellpath, newpath)
    np.save("./cells_info/cells_info.npy", cellsinfo)

if __name__ == "__main__":
    dstroot = 'ab_cells'
    listcells = os.listdir(dstroot)
    cellsinfo = []
    for n,i in zip(listcells,tqdm(range(len(listcells)))):
        cellinfo = {}
        cellpath = os.path.join(dstroot, n)
#         print(cellpath)
        img = cv2.imread(cellpath)
        img_gray = cv2.imread(cellpath, 0)
        img_gray = bf.get_fit_img(img_gray)
        #cv2.imwrite(cellpath+'abc0.png',img_gray)
        value_1,value_2 = bf.get_2value(img_gray)
        sign_nuclei, cell_nuclei_mask, nuclei_cnt, nuclei_area, nuclei_hull_area, nuclei_circ, nuclei_rule, cell_nuclei_value = get_cell_nuclei_mask(img_gray, img, value_1) #获取细胞核掩码、个数、面积、凸面积、周长、核形规则度、获取细胞核深染程度
        #cv2.imwrite(cellpath+str(cell_nuclei_value)[0:5]+'abc1.png',cell_nuclei_mask)
        if sign_nuclei == 0:
            continue
        sign_cytoplasm, cell_cytoplasm_mask, cytoplasm_area, cytoplasm_hull_area, cytoplasm_circ, cytoplasm_rule, cell_cytoplasm_value,cytoplasm_var = get_cell_cytoplasm_mask(img_gray, cell_nuclei_mask, nuclei_area,img, value_2) #获取细胞质掩码、面积、凸面积、周长、细胞规则度、获取细胞质情况
        #cv2.imwrite(cellpath+'abc2.png',cell_cytoplasm_mask)
        if sign_cytoplasm == 0 or (cytoplasm_area-nuclei_area) == 0:
            continue
        cell_N_C = nuclei_area/(cytoplasm_area-nuclei_area) #计算核质比
        cellinfo_keys = ['cellpath','nuclei_cnt','nuclei_area','nuclei_hull_area','nuclei_circ','nuclei_rule','cell_nuclei_value','cytoplasm_area','cytoplasm_hull_area','cytoplasm_circ','cytoplasm_rule','cell_cytoplasm_value','cytoplasm_var','cell_N_C']
        cellinfo_values = [cellpath,nuclei_cnt,nuclei_area,nuclei_hull_area,nuclei_circ,nuclei_rule,cell_nuclei_value,cytoplasm_area,cytoplasm_hull_area,cytoplasm_circ,cytoplasm_rule,cell_cytoplasm_value,cytoplasm_var,cell_N_C]
        cellinfo = dict(zip(cellinfo_keys, cellinfo_values))
        cell_save_sign = get_cell_save_sign(cellinfo)
        if cell_save_sign == 0:
            continue
        cellsinfo.append(cellinfo)
        newpath = os.path.join('valid', n)
        shutil.copy(cellpath, newpath)
    np.save("./cells_info/cells_info.npy", cellsinfo)
