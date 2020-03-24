import os, copy
import numpy as np
import shutil 

def find_abnormal_data(temp,limit): #寻找不正常数据
    vaule_mean = temp.mean()
    cnt = 0
    temp_out = []
    for n in temp:
        if limit > 0:
            if (n-vaule_mean)/n > limit:
                temp_out.append(1)
            else:
                temp_out.append(0)
        else:
            if (n-vaule_mean)/n < limit:
                temp_out.append(1)
            else:
                temp_out.append(0)
        cnt = cnt + 1
    temp_out = np.array(temp_out)
    return temp_out

def get_sign(cells_info_file): #解析细胞信息
    cells_info = np.load(cells_info_file)
    nuclei_area_ = []
    nuclei_rule_ = []
    nuclei_value_ = []
    cytoplasm_area_ = []
    cytoplasm_rule_ = []
    cytoplasm_value_ = []
    cell_N_C_ = []
    for n in cells_info:
        nuclei_area = n['nuclei_area']
        nuclei_area_.append(nuclei_area)
        nuclei_rule = n['nuclei_rule']
        nuclei_rule_.append(nuclei_rule)
        nuclei_value = n['cell_nuclei_value']
        nuclei_value_.append(nuclei_value)
        cytoplasm_area = n['cytoplasm_area']
        cytoplasm_area_.append(cytoplasm_area)
        cytoplasm_rule = n['cytoplasm_rule']
        cytoplasm_rule_.append(cytoplasm_rule)
        cytoplasm_value = n['cell_cytoplasm_value']
        cytoplasm_value_.append(cytoplasm_value)
        cell_N_C = n['cell_N_C']
        cell_N_C_.append(cell_N_C)
    nuclei_area_ = np.array(nuclei_area_)
    nuclei_rule_ = np.array(nuclei_rule_)
    nuclei_value_ = np.array(nuclei_value_)
    cytoplasm_area_ = np.array(cytoplasm_area_)
    cytoplasm_rule_ = np.array(cytoplasm_rule_)
    cytoplasm_value_ = np.array(cytoplasm_value_)
    cell_N_C_ = np.array(cell_N_C_)
    
    nuclei_area_abnormal = find_abnormal_data(nuclei_area_, 2) #单项阈值
    nuclei_rule_abnormal = find_abnormal_data(nuclei_rule_, 0.5)
    nuclei_value_abnormal = find_abnormal_data(nuclei_value_, -0.2)
    cytoplasm_area_abnormal = find_abnormal_data(cytoplasm_area_, 0.1)
    cytoplasm_rule_abnormal = find_abnormal_data(cytoplasm_rule_, 0.1)
    cytoplasm_value_abnormal = find_abnormal_data(cytoplasm_value_, -0.2)
    cell_N_C_abnormal = find_abnormal_data(cell_N_C_, 4)

    cells_infos = []
    cells_infos = [nuclei_area_abnormal, nuclei_rule_abnormal, nuclei_value_abnormal, cytoplasm_area_abnormal, cytoplasm_rule_abnormal, cytoplasm_value_abnormal, cell_N_C_abnormal]
    cells_infos = np.array(cells_infos)
    cells_infos_2 = cells_infos.copy()
    cells_infos_2[0,:] = cells_infos[0,:] * 20 #综合权重
    cells_infos_2[1,:] = cells_infos[1,:] * 5
    cells_infos_2[2,:] = cells_infos[2,:] * 20
    cells_infos_2[3,:] = cells_infos[3,:] * 5
    cells_infos_2[4,:] = cells_infos[4,:] * 5
    cells_infos_2[5,:] = cells_infos[5,:] * 10
    cells_infos_2[6,:] = cells_infos[6,:] * 35
    cells_infos_2 = np.array(cells_infos_2)
    cnt_2 = 0
    for n in range(0, cells_infos_2.shape[1]):
        sum_ = sum(cells_infos_2[:,n])
        if sum_ >= 40:
            cellpath = cells_info[n]['cellpath']
            newpath = os.path.join('cells_abnormal', cellpath.split('/')[1])
            shutil.copy(cellpath, newpath)
        cnt_2 = cnt_2 + 1
    sign = 0
    if cnt_2 > 10: #不正常细胞个数大于１０认为为阳性病例
        sign = 1
    return sign

if __name__ == "__main__":
    cells_info_file = 'cells_info/cells_info.npy'
    get_sign(cells_info_file)
