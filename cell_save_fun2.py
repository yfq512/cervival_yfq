import os, copy
import numpy as np
import shutil

def get_save_sign(cells_info_file):
    cells_info = np.load(cells_info_file)
    nuclei_area = []
    nuclei_hull_area = []
    cytoplasm_area = []
    cytoplasm_hull_area = []
    for n in cells_info:
        nuclei_area = n['nuclei_area']
        nuclei_hull_area = n['nuclei_hull_area']
        cytoplasm_area = n['cytoplasm_area']
        cytoplasm_hull_area = n['cytoplasm_hull_area']
        nuclei_va = nuclei_area/nuclei_hull_area
        cytoplasm_va = cytoplasm_area/cytoplasm_hull_area
 #       print(nuclei_va, cytoplasm_va)
        area_diff = cytoplasm_area - nuclei_area
        if nuclei_va < 0.9 or cytoplasm_va < 0.85 or area_diff <150:
            sign = 0
        else:
            sign = 1
        cellpath = n['cellpath']
        end_name = '_'+str(nuclei_va)[0:5]+'_'+str(cytoplasm_va)[0:5]+'_'+str(area_diff)[0:5]+'.png'
        if sign == 1:
            newpath = os.path.join('valid', cellpath.split('/')[1] + end_name)
            shutil.copy(cellpath, newpath)
        else:
            newpath = os.path.join('invalid', cellpath.split('/')[1] + end_name)
            shutil.copy(cellpath, newpath)
    return 0

if __name__ == "__main__":
    cells_info_file = 'cells_info/cells_info.npy'
    get_save_sign(cells_info_file)
#细胞核靠近图像边缘就pass（靠近用信心点坐标，或细胞核像素坐标贴在边缘）
#细胞质不包含细胞核pass
#细胞过滤就大功告成
