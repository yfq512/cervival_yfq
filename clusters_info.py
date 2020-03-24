if __name__ == "__main__":
dstroot = 'clusters'
listcluster = os.listdir(dstroot)
clustersinfo = []
for n in listcluster:
clusterinfo = {}
clusterpath = os.path.join(dstroot, n)
img = cv2.imread(clusterpath)
img_gray = cv2.imread(clusterpath, 0)
#获取阈值
cells_mask, cells_cnt, cells_area, cells_circ, cells_rule = get_cluster_info_1(img_gray) #获取细胞团细胞核掩码、细胞核数、细胞面积，细胞周长、细胞规则度
cytoplasm_mask = get_cluster_mask(img_gray) #获取细胞质掩码
cluster_nuclei_value = get_cluster_nuclei_value(img, cells_mask) #获取细胞核深染程度
cluster_cytoplasm_value = get_cluster_cytoplasm_value(img, cytoplasm_mask) #获取细胞质信息
