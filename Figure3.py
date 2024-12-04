"""
Produce the 'random landscape simulation' charts for Figure 5
"""
from math import pi
from pandas import DataFrame
from skimage.draw import disk
from numpy.random import uniform
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import subplots, savefig, show
from numpy import arange, hypot as np_hypot, meshgrid, nonzero, ones, \
    min as np_min, sum as np_sum, zeros, unique, nan, where


'''
*** WEIGHTING KEY ***
0 = Unweighted
1 = Away from focal location
2 = Towards focal location
'''
RADIUS = 100
ITERATIONS = 10000
Y_LIM = ITERATIONS * 1.2
PRINT_PROPORTIONS = False
''' ****** '''


def fab(d, r):
    """ Calculate FAB correction for a given distance """
    return r / d if d > 0 else 0

def idw(d, r):
    """ Calculate IDW for a given distance """
    return 1 - (d / r) if d > 0 else 0

def fab_idw(d, r):
    """ Calculate FAB-Corrected IDW for a given distance """
    return fab(d, r) * idw(d, r)

def euclidean_distance(arr):
    """
    Calculate Euclidean distance surface (in image space)
    """
    # create meshgrid
    y, x = arr.shape
    xx, yy = meshgrid(arange(x), arange(y))

    # create indices where arr is different from NoData and reshape them
    ind = nonzero(arr)
    ix = ind[1].reshape((-1, 1, 1))
    iy = ind[0].reshape((-1, 1, 1))

    # compute legs
    dx = abs(iy - yy)
    dy = abs(ix - xx)

    # return distance
    return np_min(np_hypot(dx, dy), axis=0)



#  set up circle mask
diameter = 2 * RADIUS
area = pi * RADIUS**2
mask = zeros((diameter, diameter))
rr, cc = disk((RADIUS, RADIUS), RADIUS)
mask[rr, cc] = 1

# get cell index for all zeros - these are removed before plotting
clipper = where(mask == 0)

# set up distance surface
distance = zeros((diameter, diameter))
distance[RADIUS, RADIUS] = 1
distance = euclidean_distance(distance) * mask


# plot
my_fig, my_ax = subplots(3, 3, figsize=(12, 18))

# do all weights and wighted / unweighted
srow = -1 
for WEIGHT in [2, 0, 1]:
    srow +=1
    for IDW in [False, True]:

        # get un-corrected and corrected surface for IDW
        if IDW:
            orig_layer = zeros(distance.shape)
            fab_layer = zeros(distance.shape)
            for d in unique(distance):
                orig_layer[distance == d] = idw(d, RADIUS)
                fab_layer[distance == d] = fab_idw(d, RADIUS)

        # ...or regular
        else:
            orig_layer = mask
            fab_layer = zeros(distance.shape)
            for d in unique(distance):
                fab_layer[distance == d] = fab(d, RADIUS)

        # get denominators
        fab_layer[RADIUS, RADIUS] = 1
        denom = orig_layer.sum()
        fab_denom = fab_layer.sum()

        # unweighted (0.5)
        if WEIGHT == 0:
            # just a surface of 0.5's
            p = ones(distance.shape, dtype=float) / 2

        # weight towards the edge (0 - MAX_P)
        elif WEIGHT == 1:
            p = 1 - (1 - (distance / RADIUS))*1.58

        # weight towards the focal point (0 - MAX_P)
        elif WEIGHT == 2:
            p = (1 - (distance / RADIUS))*1.58

        else:
            print("invalid weight setting")
            exit()

        # run n iterations
        n = 0
        outputs = []
        DONE = False
        while n < ITERATIONS:
            
            # get weighted noise (weighted based on p)
            noise = zeros((diameter, diameter))
            for r in range(noise.shape[0]):
                for c in range(noise.shape[1]):
                    noise[(r,c)] = 1 if uniform() < p[(r,c)] else 0
            noise *= mask

            # do things that should only happen once
            if not DONE:
                DONE = True

                # verify that the proportions of 1's and 0's in the noise surface are correct (runs once only)
                if PRINT_PROPORTIONS:
                    noise[clipper[0], clipper[1]] = nan
                    _, counts = unique(noise, return_counts=True)
                    print(f"0: {counts[0]/denom:.1f}\t1:{counts[1]/denom:.1f}")

                # remove areas outside of the radius
                noise[clipper[0], clipper[1]] = nan

                # add an example landcover map
                my_ax[srow][0].imshow(noise, cmap = LinearSegmentedColormap.from_list('binary', [(178/255, 223/255, 138/255, 1), (0.7, 0.7, 0.7, 1)], N=2))

                # add legend to map
                my_ax[srow][0].legend(
                    handles=[
                        Patch(facecolor=(0.7, 0.7, 0.7, 1), edgecolor=None, label=f'Urban'),
                        Patch(facecolor=(178/255, 223/255, 138/255, 1), edgecolor=None, label=f'Green'),
                    ], 
                    loc='lower left', 
                    bbox_to_anchor=(0.7, -0.2)
                    )
                my_ax[srow][0].set_title('Urban Scenario' if srow == 0 else 'Green Scenario' if srow == 2 else "Mixed Scenario")

                # disable axis
                my_ax[srow][0].set_axis_off()

            # store the proportions of natural land cover before & after and increment counter
            outputs.append({
                'Uncorrected': np_sum(noise * orig_layer) / denom,
                'Corrected': np_sum(noise * fab_layer) / fab_denom
            })
            n += 1

        # store dataframe
        df = DataFrame(outputs)
        print(f"mean uncorrected: {df.Uncorrected.mean():.2f}, mean corrected:{df.Corrected.mean():.2f}")

        # set range for histogram axes
        xlim = [0, 1]
        ylim = [0, Y_LIM]

        # plot historams
        splt = 1 if not IDW else 2
        bins = [n / 100 for n in list(range(0, 100, 5))]
        my_ax[srow][splt].hist([df.Uncorrected, df.Corrected], bins=bins, histtype='bar', stacked=True, alpha=0.67) #, color=['#1b9e77','#d95f02'])
        my_ax[srow][splt].grid(True)
        my_ax[srow][splt].legend(['Uncorrected', 'Corrected'], loc='upper right')
        my_ax[srow][splt].set_title("Unweighted"if not IDW else "IDW Weighted")
        my_ax[srow][splt].set_xlim(xlim)
        my_ax[srow][splt].set_ylim(ylim)

# write to file
savefig("Figure3.png", bbox_inches='tight', dpi=300)