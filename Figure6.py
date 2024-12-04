"""
Produce the 'approximate correction' charts for Figure 9
"""
from math import hypot
from skimage.draw import disk
from numpy import zeros, column_stack, array
from matplotlib.pyplot import subplots, savefig, rc, subplots_adjust


def fab(d, r):
    """ Calculate FAB Correction """
    return r / d if d > 0 else 0

def idw(d, r):
    """ Calculate IDW Weighting """
    return 1 - (d / r) if d > 0 else 0

# plot
fig, my_ax = subplots(1, 2, figsize=(17, 8))
FONT = 14
rc('font', size=FONT)          # controls default text sizes
rc('axes', labelsize=FONT)    # fontsize of the x and y labels
rc('legend', fontsize=FONT)    # legend fontsize

# set the radii
for i, radius in enumerate([50, 250]):

    # setup
    dim = radius * 2
    centre = (radius, radius)

    # calculate outer ring
    outer = disk(centre, radius)

    # initialise the landscape and outer ring
    landscape_bk = zeros((dim, dim))
    landscape_bk[outer] = 1

    # init
    x = [0]
    uncorrected = [0]
    corrected = [0]

    # increment 0.1-1 in 0.1 increments
    for inner_proportion in array(range(1, 51, 1)) / 50:

        # record x value
        x.append(inner_proportion)

        # refresh landscape for uncorrected and corrected versions
        uc_surf = landscape_bk.copy()
        c_surf = landscape_bk.copy()

        # calculate inner ring
        inner = disk(centre, radius * inner_proportion)

        # now apply weighting to entire disk
        for r, c in column_stack(outer):
            dist = hypot(c - centre[1], r - centre[0])
            c_surf[(r,c)] = fab(dist, radius)

        # mean value
        uncorrected.append(sum(uc_surf[inner]) / sum(uc_surf[outer]))
        corrected.append(sum(c_surf[inner]) / sum(c_surf[outer]))
            
    # add results from each function
    my_ax[i].plot(x, uncorrected, linewidth=3, color='#1f78b4')
    my_ax[i].plot(x, corrected, linewidth=3, color='#33a02c')
    my_ax[i].plot(x, x, linewidth=3, color='#000', linestyle=(5, (8, 4)))

# add labels
my_ax[1].legend([ 'Uncorrected', 'Approx. FAB Correction', 'FAB Correction'], loc='lower right')

# format plots
for n in [0, 1]:
    my_ax[n].grid(True)
    my_ax[n].set_xlabel('Radius')
    my_ax[n].set_ylabel('Effective Area')

# add labels
my_ax[0].annotate('A. 50 cell radius',
        xy=(0, 1), xycoords='axes fraction',
        xytext=(+0.5, -0.5), textcoords='offset fontsize',
        fontsize='large', verticalalignment='top', fontfamily='sans serif')
my_ax[1].annotate('B. 250 cell radius',
        xy=(0, 1), xycoords='axes fraction',
        xytext=(+0.5, -0.5), textcoords='offset fontsize',
        fontsize='large', verticalalignment='top', fontfamily='sans serif')

# reduce gaps
subplots_adjust(wspace=0.15)

# output
savefig(f"Figure7.png", bbox_inches='tight', dpi=300)