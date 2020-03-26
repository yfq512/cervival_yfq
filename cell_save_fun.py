import cv2
import numpy as np
import os,shutil,math
import basefun as bf

def get_fit_img(img):
    dst = cv2.fastNlMeansDenoising(img,None,15,7,21)
    return dst
 
def get_mean(temp, long_=10):
    temp_1 = []
    for m in range(0+long_, 255-long_):
        temp_1.append(temp[m-long_:m+long_].mean())
    temp_2 = np.array(temp_1)
    temp_out = np.zeros(255)
    temp_out[0+long_: 255-long_] = temp_2
    return temp_out
 
def grien_value(temp):
    temp_1 = np.zeros(255)
    temp_1[1:255] = temp[0:254]
    temp_out = temp - temp_1
    temp_out = get_mean(temp_out, 10)
    return temp_out
 
def get_last2value(temp):
    sign1 = 0
    sign2 = 0
    cnt = 0
    value = 0
    temp_4 = temp[::-1]
    for k in range(0, len(temp_4)-1):
        if temp_4[k+1] > temp_4[k]:
            sign2 = sign1
            sign1 = 1
            if sign1 != sign2:
                cnt = cnt + 1
        else:
            sign2 = sign1
            sign1 = 0
            if sign1 != sign2:
                cnt = cnt + 1
        if cnt == 3:
            value = k
            break #!!!!!!!!!
    return 255-value
 
def get_2value(img, long_ = 10, chan = 0, mask = None):
    temp_2 = cv2.calcHist([img],[chan],mask,[256],[0,255])
    temp_4 = get_mean(temp_2, long_)
    temp_4 = get_mean(temp_4, 5)
    temp_4 = grien_value(temp_4) #使用二阶导数更容易检测到细胞核像素阈值
#    plt.figure("Image")
#    plt.imshow(img)
#    plt.figure("Image_value")
#    plt.plot(temp_4)
#    plt.figure("Image_value_gre")
#    plt.plot(grien_value(temp_4))
#    plt.figure("Image_value_gre2")
#    plt.plot(grien_value(grien_value(temp_4)))
#    plt.show()
    sign1 = 0
    sign2 = 0
    cnt = 0
    value_1 = 0
    value_2 = 0
    for k in range(0, len(temp_4)-1):
        if temp_4[k+1] > temp_4[k]:
            sign2 = sign1
            sign1 = 1
            if sign1 != sign2:
                cnt = cnt + 1
        else:
            sign2 = sign1
            sign1 = 0
            if sign1 != sign2:
                cnt = cnt + 1
        if cnt == 3:
            value_1 = k
            #print(temp_4[k])
    value_2 = get_last2value(temp_4)
    return cnt

def find_invalid_cells(dstroot,limit):
    list_ = os.listdir(dstroot)
    for n in list_:
        cellpath = os.path.join(dstroot, n)
        cell = cv2.imread(cellpath, 1)
        width,height = cell.shape[:2][::-1]
        img_resize = cv2.resize(cell,(int(width*1.0),int(height*1.0)),interpolation=cv2.INTER_CUBIC)
        img_gray = cv2.cvtColor(img_resize,cv2.COLOR_RGB2GRAY)
        imageVar = cv2.Laplacian(img_gray, cv2.CV_64F).var()
        if imageVar<limit:
            newpath = os.path.join('invalid',n)
            shutil.copy(cellpath,newpath)

def get_entropy(img_):
    x, y = img_.shape[0:2]
    img_ = cv2.resize(img_, (100, 100)) # 缩小的目的是加快计算速度
    tmp = []
    for i in range(256):
        tmp.append(0)
    val = 0
    k = 0
    res = 0
    img = np.array(img_)
    for i in range(len(img)):
        for j in range(len(img[i])):
            val = img[i][j]
            tmp[val] = float(tmp[val] + 1)
            k =  float(k + 1)
    for i in range(len(tmp)):
        tmp[i] = float(tmp[i] / k)
    for i in range(len(tmp)):
        if(tmp[i] == 0):
            res = res
        else:
            res = float(res - tmp[i] * (math.log(tmp[i]) / math.log(2.0)))
    return res

if __name__ == "__main__":
    dstroot = 'crop'
 #   find_invalid_cells(dstroot,50)
    list__ = os.listdir(dstroot)
    for n in list__:
        imgpath = os.path.join(dstroot, n)
        img = cv2.imread(imgpath, 0)
        img_fit = get_fit_img(img)
        cnt = get_2value(img_fit)
        res = get_entropy(img_fit)
        _, value2 = bf.get_2value(img_fit)
        ret_, thresh = cv2.threshold(img, value2, 255, cv2.THRESH_BINARY_INV)
        kernel_1 = np.ones((5,5),np.uint8)
        thresh = bf.get_img_open(thresh,kernel_1)
       # thresh = bf.get_img_close(thresh,kernel_1)
        image_, contours, hierarchy_ = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
        rule = 100
        if len(contours) == 0:
            rule = 100
        else:
            area_ = 0
            contours_max = []
            for cnt in range(0, len(contours)):
                area = cv2.contourArea(contours[cnt])
                if area > area_:
                    contours_max = contours[cnt]
                    area_ = area
            perimeter = cv2.arcLength(contours_max,True)
            if area_ == 0:
                rule = 100
            else:
                rule = perimeter/area_
        thresh2 = thresh*0
        thresh2 = cv2.fillConvexPoly(thresh2, contours_max, 255)
        print('规则度',rule)
        if rule < 0.14:
        #if cnt < 12 and cnt >=7 and rule <0.14:
 #       if res < 6:
            print(res)
            #print(imgpath)
            #print('abnormal',cnt)
            newpath = os.path.join('valid',n+'_'+str(cnt)+'_'+str(rule)+'.png')
            #cv2.imwrite(newpath+'_'+str(rule)+'_'+'abc.png', thresh2)
            #cv2.imwrite(newpath+'_'+str(rule)+'_'+'abc2.png', thresh)
            shutil.copy(imgpath,newpath)
        else:
            newpath = os.path.join('invalid',n+'_'+str(cnt)+'_'+str(rule)+'.png')
            shutil.copy(imgpath,newpath)
