import os,shutil,time
import main

def run_32(blpath, outpath):
    os.makedirs(outpath)
    fov = 'fovs'
    if os.path.exists(fov):
        shutil.rmtree(fov)
    shutil.copytree(blpath, fov)
    main.main2()
    outpath_ab_cells = os.path.join(outpath,'ab_cells')
    outpath_valid = os.path.join(outpath,'valid')
    outpath_cells_info = os.path.join(outpath,'cells_info')
    shutil.copytree('ab_cells', outpath_ab_cells)
    shutil.copytree('valid', outpath_valid)
    shutil.copytree('cells_info', outpath_cells_info)

    

if __name__ == "__main__":
    dstroot = '20200205'
    outroot = 'temp_ab_cells'
    dst_list = os.listdir(dstroot)
    for n in dst_list:
        blpath = os.path.join(dstroot, n, 'Images')
        outpath = os.path.join(outroot, n)
        run_32(blpath, outpath)

