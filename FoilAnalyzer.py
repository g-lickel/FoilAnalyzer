#-----------------------------------------------------------------------------------------------------------
# TH-GEM Foil Analyzer
# G. Lickel - v 1.0
# Last update: Mar. 05 2025
#-----------------------------------------------------------------------------------------------------------
# Packages
import matplotlib.image as img
import matplotlib.pyplot as plt
import numpy as np
import time as tm
import pandas as pd
#-----------------------------------------------------------------------------------------------------------
def main():

# Initial parameters:
    file = 'CMGL006_10X10_3600_Cor_24Bits'   # File name

    inDiam = 0.3    # expected internal diameter (mm)
    outDiam = 0.5    # expected external diameter (mm)
    pitch = 1.55  # expected pitch (mm)
    dpi = 3600  # image definition (dots per inch)

    firstCenter = 2.75   # Approximated first center of first column (mm)
#-----------------------------------------------------------------------------------------------------------
    inp = int(input('1 - Analyze image file ('+file+'.jpg).\n2 - Measure analyzed file ('+file+'_GS.png).\n3 - Plot histogram from '+file+'_hist.csv (WARNING: PyROOT required).\nWhat do you wish? '))

    if inp == 1:
        analysis(file)
    elif inp == 2:
        measure(inDiam, outDiam, pitch, dpi, file, firstCenter)
    elif inp == 3:
        plot(file, dpi)

    return
#-----------------------------------------------------------------------------------------------------------
def analysis(file):

    (copper, outhole, inhole) = data()  # Load colors data for reference

    print('Initializing Analysis...')
    t0 = tm.time()  # Takes initial time

    imag  = img.imread(file+'.jpg') # Load .jpg image

    image = imag.copy()

    (r, c, _) = image.shape # Define image resolution (rolls and columns)

    total = r*c
    j = 0

    for x in range(0, c):
        for y in range(0, r):
            dist = []
            dist.append(distPoint(image[y, x], copper))
            dist.append(distPoint(image[y, x], outhole))
            dist.append(distPoint(image[y, x], inhole))

            if min(dist) == dist[0]:
                image[y, x] = [0, 0, 0]
            elif min(dist) == dist[1]:
                image[y, x] = [189, 189, 189]
            else:
                image[y, x] = [255, 255, 255]

            j += 1
            p = j/total
            print_bar(t0, p) # Print progress bar

    t = tm.time() - t0
    if t < 60:
        print('Finished in {0:2.2f}s                                                                     '.format(t))
    elif t < 3600:
        print('Finished in {0:2.2f}min                                                                   '.format(t/60))
    else:
        print('Finished in {0:2.2f}h                                                                     '.format(t/3600))
    
    img.imsave(file+'_GS.png', image)
    plt.imshow(image)
    plt.show()

    return
#-----------------------------------------------------------------------------------------------------------
def measure(inDiam, outDiam, pitch, dpi, file, firstCenter):

    print('Initializing Measurements...')
    t0 = tm.time()  # Takes initial time

    imag  = img.imread(file+'_GS.png') # Load .jpg image

    image = imag.copy()

    (r, c, _) = image.shape # Define image resolution (rolls and columns)

    centers = [firstCenter*dpi/25.4]
    while True:
        ct = centers[-1]+(pitch*dpi/25.4)/2
        if ct+((pitch/2)*dpi/25.4)/2 <= c:
            centers.append(ct)
        else:
            break

    total = int(((pitch/2)*dpi/25.4)*r*len(centers))
    j = 0

    d = []
    m = 0

    for col in range(len(centers)):
        d.append([])
        for x in range(int(centers[col]-((pitch/2)*dpi/25.4)/2), int(centers[col]+((pitch/2)*dpi/25.4)/2)):
            n = -1
            for y in range(0, r):

                if image[y, x][0] == 1. and image[y, x][1] == 1. and image[y, x][2] == 1.:
                    m += 1
                    if m == int((inDiam*dpi/25.4)/2):
                        n += 1
                        if len(d[col]) == n:
                            d[col].append(0)
                else:
                    m = 0

                if n >= 0:
                    if m > d[col][n]:
                            d[col][n] = m

                j += 1
                p = j/total
                print_bar(t0, p) # Print progress bar

    t = tm.time() - t0
    if t < 60:
        print('Finished in {0:2.2f}s                                                                     '.format(t))
    elif t < 3600:
        print('Finished in {0:2.2f}min                                                                   '.format(t/60))
    else:
        print('Finished in {0:2.2f}h                                                                     '.format(t/3600))

    for lin in range(len(d)):
        for num in range(len(d[lin])):
            d[lin][num] = d[lin][num]*25.4/dpi
    
    with open(file+'_hist.csv', 'w') as file:
        file.write("\n")
        for lin in d:
            for num in lin:
                file.write(f"{num}\n")

    return
#-----------------------------------------------------------------------------------------------------------
def plot(file, dpi):

    from ROOT import gStyle, TCanvas, TH1D, TF1

    gStyle.SetOptFit(111) # superimpose fit results
    c1 = TCanvas(file,'Data' ,200 ,10 ,700 ,500) # make nice

    table = pd.read_csv(file+'_hist.csv',names=['x'])
    dados = table['x'].values

    N = len(dados)
    min_res = 25.4/dpi  # Min. resolution in mm
    b = int(min(np.sqrt(N)+1, (max(dados)-min(dados))/min_res))

    h1 = TH1D('Data', ';Diameter (mm);Entries', b, min(dados), max(dados))

    for i in dados:
        h1.Fill(i)

    # Define Gaussian function for fitting
    gauss_fit = TF1("gauss_fit", "gaus", min(dados), max(dados))
    h1.Fit(gauss_fit, "R")  # Perform fit in range

    h1.Draw()
    c1.Update()

    # request user action before ending (and deleting graphics window)
    input()

    return
#-----------------------------------------------------------------------------------------------------------
# Functions
def data():

    copper = [
        [0, 0, 0,],
        [99, 52, 8],
        [64, 32, 2],
        [63, 21, 0],
        [180, 119, 60],
        [119, 76, 15],
        [135, 90, 44],
        [255, 254, 187],
        [110, 71, 30],
        [50, 10, 0],
        [184, 126, 75],
        [172, 127, 86]
    ]

    outhole = [
        [120, 115, 72],
        [115, 147, 124],
        [126, 123, 78],
        [158, 156, 124],
        [116, 113, 64],
        [119, 112, 66],
        [156, 146, 111],
        [129, 119, 81],
        [144, 144, 111],
        [121, 117, 80],
        [99, 96, 60],
        [144, 136, 103]
    ]

    inhole = [
        [116, 110, 54],
        [126, 118, 65],
        [98, 89, 32],
        [149, 138, 98],
        [150, 137, 91],
        [118, 111, 63],
        [107, 103, 52],
        [115, 98, 48],
        [123, 99, 43],
        [118, 100, 47],
        [43, 34, 2],
        [87, 80, 36],
        [58, 50, 1],
        [89, 86, 33],
        [175, 163, 120]
    ]

    return (copper, outhole, inhole)

def print_bar(t0, p):

    print('Progress: [{0:4.0%}] ['.format(p), end='')

    carac = 58
    for i in range(0, int(p*carac)):
        print('#', end='')
    for i in range(0, carac-int(p*carac)):
        print('.', end='')
    t = tm.time() - t0
    prev = t/p*(1-p)

    h = prev//3600
    m = (prev%3600)//60
    s = (prev%3600)%60

    print('] {0:2.0f}:{1:2.0f}:{2:2.0f}  '.format(h, m, s),end='\r')

    return

def distPoint(point, list):

    dist = []

    for coord in list:
        d = np.sqrt((point[0]-coord[0])**2 + (point[1]-coord[1])**2 + (point[2]-coord[2])**2)
        dist.append(d)

    return min(dist)

#-----------------------------------------------------------------------------------------------------------
main()
