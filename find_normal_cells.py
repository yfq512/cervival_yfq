import numpy as np
import matplotlib.pyplot as plt
import basefun as bf
import os,shutil

def find_normal_cells(cells_info_file):
    cells_info = np.load(cells_info_file)
    temp = []
    for n in cells_info:
        nuclei_area = n['nuclei_area']
        temp.append(nuclei_area)
    temp1 = np.array(temp)
    temp1 = (temp1/20).astype(int)*20
    temp2 = list(temp1)
    dict_ = {}
    for key in temp2:
        dict_[key] = dict_.get(key, 0) + 1
 #   print(dict_)
    value = dict_.values()
    key = dict_.keys()
    value = np.array(list(value))
 #   print('平滑矩阵',len(value), len(temp3))
    temp4 = np.argmax(value)
 #   print('最大索引',temp4, len(key))
    temp5 = list(key)[temp4]
 #   print('众数',temp5)
    limit = temp5*0.3
    cnt_ = 0
    cells_i = []
    for n in cells_info:
        nuclei_area = n['nuclei_area']
        nuclei_area = np.array(nuclei_area)
        if nuclei_area>(temp5-limit) and nuclei_area<(temp5+limit):
            cells_i.append(cnt_)
        cnt_ = cnt_ + 1
    return cells_i

def get_normal_cell_info(cells_info_file ,normal_cells_i):
    cells_info = np.load(cells_info_file)
    mar_nuclei_area = []
    mar_nuclei_rule = []
    mar_cell_nuclei_value = []
    mar_cytoplasm_area = []
    mar_cytoplasm_rule = []
    mar_cell_cytoplasm_value = []
    mar_cytoplasm_var = []
    mar_cell_N_C = []
    cnt = 0
    for n in cells_info:
        if cnt in normal_cells_i:
            mar_nuclei_area.append(n['nuclei_area'])
            mar_nuclei_rule.append(n['nuclei_rule'])
            mar_cell_nuclei_value.append(n['cell_nuclei_value'])
            mar_cytoplasm_area.append(n['cytoplasm_area'])
            mar_cytoplasm_rule.append(n['cytoplasm_rule'])
            mar_cell_cytoplasm_value.append(n['cell_cytoplasm_value'])
            mar_cytoplasm_var.append(n['cytoplasm_var'])
            mar_cell_N_C.append(n['cell_N_C'])
        cnt = cnt + 1
    mean_nuclei_area = np.array(mar_nuclei_area).mean()
    mean_nuclei_rule = np.array(mar_nuclei_rule).mean()
    mean_cell_nuclei_value = np.array(mar_cell_nuclei_value).mean()
    mean_cytoplasm_area = np.array(mar_cytoplasm_area).mean()
    mean_cytoplasm_rule = np.array(mar_cytoplasm_rule).mean()
    mean_cell_cytoplasm_value = np.array(mar_cell_cytoplasm_value).mean()
    mean_cytoplasm_var = np.array(mar_cytoplasm_var).mean()
    mean_cell_N_C = np.array(mar_cell_N_C).mean()
    cell_keys = ['mean_nuclei_area','mean_nuclei_rule','mean_cell_nuclei_value','mean_cytoplasm_area','mean_cytoplasm_rule','mean_cell_cytoplasm_value','mean_cytoplasm_var','mean_cell_N_C']
    cell_values = [mean_nuclei_area,mean_nuclei_rule,mean_cell_nuclei_value,mean_cytoplasm_area,mean_cytoplasm_rule,mean_cell_cytoplasm_value,mean_cytoplasm_var,mean_cell_N_C]
    cell_dict = dict(zip(cell_keys, cell_values))
    return cell_dict

def find_ab_fun(temp,vaule_mean,limit):
    cnt = 0
    temp_out = []
    for n in temp:
        if limit > 0:
            #print(n-vaule_mean, vaule_mean, limit)
            #exit()
            if (n-vaule_mean)/vaule_mean > limit:
                temp_out.append(1)
            else:
                temp_out.append(0)
        else:
            if (n-vaule_mean)/vaule_mean < limit:
                temp_out.append(1)
            else:
                temp_out.append(0)
        cnt = cnt + 1
    temp_out = np.array(temp_out)
    return temp_out

def find_abnormal_cells(cells_info_file, normal_cell_info):
    cells_info = np.load(cells_info_file)
    mar_nuclei_area = []
    mar_nuclei_rule = []
    mar_cell_nuclei_value = []
    mar_cytoplasm_area = []
    mar_cytoplasm_rule = []
    mar_cell_cytoplasm_value = []
    mar_cytoplasm_var = []
    mar_cell_N_C = []
    cnt = 0
    for n in cells_info:
        mar_nuclei_area.append(n['nuclei_area'])
        mar_nuclei_rule.append(n['nuclei_rule'])
        mar_cell_nuclei_value.append(n['cell_nuclei_value'])
        mar_cytoplasm_area.append(n['cytoplasm_area'])
        mar_cytoplasm_rule.append(n['cytoplasm_rule'])
        mar_cell_cytoplasm_value.append(n['cell_cytoplasm_value'])
        mar_cytoplasm_var.append(n['cytoplasm_var'])
        mar_cell_N_C.append(n['cell_N_C'])
        cnt = cnt + 1
    mar_nuclei_area = np.array(mar_nuclei_area)
    mar_nuclei_rule = np.array(mar_nuclei_rule)
    mar_cell_nuclei_value = np.array(mar_cell_nuclei_value)
    mar_cytoplasm_area = np.array(mar_cytoplasm_area)
    mar_cytoplasm_rule = np.array(mar_cytoplasm_rule)
    mar_cell_cytoplasm_value = np.array(mar_cell_cytoplasm_value)
    mar_cytoplasm_var = np.array(mar_cytoplasm_var)
    mar_cell_N_C = np.array(mar_cell_N_C)
    abnormal_nuclei_area = find_ab_fun(mar_nuclei_area,normal_cell_info['mean_nuclei_area'],2)
    abnormal_nuclei_rule = find_ab_fun(mar_nuclei_rule,normal_cell_info['mean_nuclei_rule'],-0.01)
    abnormal_cell_nuclei_value = find_ab_fun(mar_cell_nuclei_value,normal_cell_info['mean_cell_nuclei_value'],0.1)
    abnormal_cytoplasm_area = find_ab_fun(mar_cytoplasm_area,normal_cell_info['mean_cytoplasm_area'],-0.1)
    abnormal_cytoplasm_rule = find_ab_fun(mar_cytoplasm_rule,normal_cell_info['mean_cytoplasm_rule'],-0.02)
    abnormal_cell_cytoplasm_value = find_ab_fun(mar_cell_cytoplasm_value,normal_cell_info['mean_cell_cytoplasm_value'],0.1)
    abnormal_cytoplasm_var = find_ab_fun(mar_cytoplasm_rule,normal_cell_info['mean_cytoplasm_var'],0.1)
    abnormal_cell_N_C = find_ab_fun(mar_cell_N_C,normal_cell_info['mean_cell_N_C'],4)
    cells_infos = []
    cells_infos = [abnormal_nuclei_area,abnormal_nuclei_rule,abnormal_cell_nuclei_value,abnormal_cytoplasm_area,abnormal_cytoplasm_rule,abnormal_cell_cytoplasm_value,abnormal_cytoplasm_var,abnormal_cell_N_C]
    cells_infos = np.array(cells_infos)
    cells_infos_2 = cells_infos.copy()
    cells_infos_2[0,:] = cells_infos[0,:] * 10 #综合权重
    cells_infos_2[1,:] = cells_infos[1,:] * 0
    cells_infos_2[2,:] = cells_infos[2,:] * 0
    cells_infos_2[3,:] = cells_infos[3,:] * 0
    cells_infos_2[4,:] = cells_infos[4,:] * 0
    cells_infos_2[5,:] = cells_infos[5,:] * 0
    cells_infos_2[5,:] = cells_infos[6,:] * 0
    cells_infos_2[6,:] = cells_infos[7,:] * 10
    cells_infos_2 = np.array(cells_infos_2)
    cnt_2 = 0
    for n in range(0, cells_infos_2.shape[1]):
        sum_ = sum(cells_infos_2[:,n])
        #if sum_ >= 20:
        if True:
            cellpath = cells_info[n]['cellpath']
            nuclei_area_name = cells_info[n]['nuclei_area']
            nuclei_rule_name = cells_info[n]['nuclei_rule']
            cell_nuclei_value_name = cells_info[n]['cell_nuclei_value']
            cytoplasm_area_name = cells_info[n]['cytoplasm_area']
            cytoplasm_rule_name = cells_info[n]['cytoplasm_rule']
            cell_cytoplasm_value_name = cells_info[n]['cell_cytoplasm_value']
            cytoplasm_var_name = cells_info[n]['cytoplasm_var']
            cell_N_C_name = cells_info[n]['cell_N_C']
            newpath = os.path.join('ab_cells', cellpath.split('/')[1]+'_'+'nuarea'+str(nuclei_area_name)[0:5]+'_'+'nurule'+str(nuclei_rule_name)[0:5]+'_'+'nuvalue'+str(cell_nuclei_value_name)[0:5]+'_'+'cyarea'+str(cytoplasm_area_name)[0:5]+'_'+'cyrule'+str(cytoplasm_rule_name)[0:5]+'_'+'cyvalue'+str(cell_cytoplasm_value_name)[0:5]+'_'+'cyvar'+str(cytoplasm_var_name)[0:5]+'_'+'NC'+str(cell_N_C_name)[0:5]+'.png')
            shutil.copy(cellpath, newpath)
            cnt_2 = cnt_2 + 1
    sign = 0
    if cnt_2 > 20: #dd不正常细胞个数大于１０认为为阳性病例
        sign = 1
    return sign

def main2():
    cells_info_file = 'cells_info/cells_info.npy'
    normal_cells_i = find_normal_cells(cells_info_file)
    normal_cell_info = get_normal_cell_info(cells_info_file, normal_cells_i)
    print(normal_cell_info)
    sign = find_abnormal_cells(cells_info_file, normal_cell_info)
if __name__ == "__main__":
    cells_info_file = 'cells_info/cells_info.npy'
    normal_cells_i = find_normal_cells(cells_info_file)
    normal_cell_info = get_normal_cell_info(cells_info_file, normal_cells_i)
    print(normal_cell_info)
    sign = find_abnormal_cells(cells_info_file, normal_cell_info)
