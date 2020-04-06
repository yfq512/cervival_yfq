import os,shutil
import cv2
import numpy as np
import basefun as bf
import cells_clusters_seg
import cells_info
import find_normal_cells
import clusters_info

def clean_file(dst):
    print('================clean dirs================')
    for n in dst:
        if os.path.exists(n):
            shutil.rmtree(n)
    print('================make dirs================')
    for m in dst:
        os.makedirs(m)

def seg_cells(fovs_root):
    print('================seg cells================')
    cells_clusters_seg.main2()

def get_valid_cells_info():
    print('================get valid cells info================')
    cells_info.main2()

def get_ab_clusters():
    print('================find abnormal clusters================')
    clusters_info.main2()

def get_abnormal_cells():
    print('================find abnormal cells================')
    find_normal_cells.main2()

if __name__ == "__main__":
    dst = ['crop','clusters','ab_clusters','valid','invalid','cells_info','ab_cells']
    clean_file(dst)
    fovs_root = 'fovs'
    list_fovs = os.listdir(fovs_root)
    if len(list_fovs) == 0:
        print('================can not find files in fovs================')
        exit()
    seg_cells(fovs_root)

    crop_root = 'crop'
    list_crop = os.listdir(crop_root)
    if len(list_crop) == 0:
        print('================can not find files in crop================')
        exit()
    get_valid_cells_info()

    clusters_info__root = 'clusters'
    list_clusters_info = os.listdir(clusters_info__root)
    if len(list_clusters_info) == 0:
        print('================can not find files in clusters================')
        exit()
    get_ab_clusters()

    cells_info__root = 'cells_info'
    list_cells_info = os.listdir(cells_info__root)
    if len(list_cells_info) == 0:
        print('================can not find files in cells_info================')
        exit()
    get_abnormal_cells()

