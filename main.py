import os,shutil
import cv2
import numpy as np
import basefun as bf
import cells_clusters_seg
import cells_info
import find_normal_cells

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
    activate_this = 'cells_clusters_seg.py'
    #execfile(activate_this, dict(__file__=activate_this))
    #exec(open(activate_this).read())
    cells_clusters_seg.main2()

def get_valid_cells_info():
    print('================get valid cells info================')
    activate_this = 'cells_info.py'
    #execfile(activate_this, dict(__file__=activate_this))
    #exec(open(activate_this).read())
    cells_info.main2()

def get_abnormal_cells():
    print('================find abnormal cells================')
    activate_this = 'find_normal_cells.py'
    #execfile(activate_this, dict(__file__=activate_this))
    #exec(open(activate_this).read())
    find_normal_cells.main2()

if __name__ == "__main__":
    dst = ['crop','valid','invalid','cells_info','cells_abnormal']
    clean_file(dst)
    fovs_root = 'fovs'
    list_fovs = os.listdir(fovs_root)
    if len(list_fovs) == 0:
        print('================can not find files in fovs================')
        exit()
    
    seg_cells(fovs_root)
    
    get_valid_cells_info()
    
    get_abnormal_cells()
