"""
Produce the assumption and 'geographical weighting' charts for Figure 2
"""
from math import hypot
from scipy.stats import norm
from skimage.draw import disk
from numpy import zeros, column_stack, array
from matplotlib.pyplot import subplots, savefig, subplots_adjust, rcParams


def fab(d, r):
    """ Calculate FAB Correction """
    return r / d if d > 0 else 0

def idw(d, r):
    """ Calculate IDW Weighting """
    return 1 - (d / r) if d > 0 else 0

def fab_idw(d, r):
    """ Calculate fab * IDW Weighting """
    return fab(d, r) * idw(d, r) if d > 0 else 0 

def idw2(d, r):
    """ Calculate IDW-squared Weighting """
    return (1 - (d / r))**2 if d > 0 else 0

def fab_idw2(d, r):
    """ Calculate fab * IDW-squared Weighting """
    return fab(d, r) * idw2(d, r) if d > 0 else 0 

def gauss(d, r, sigma=1):
    """ Calculate Gaussian Weighting """
    return norm.pdf(norm.ppf(d * 0.5 / r + 0.5, loc=0, scale=sigma), loc=0, scale=sigma) if d > 0 else 0 

def fab_gauss(d, r, sigma=1):
    """ Calculate fab * Gaussian Weighting """
    return fab(d, r) * gauss(d, r, sigma) if d > 0 else 0 


# setup
radius = 100
dim = radius * 2
centre = (radius, radius)

# set font size
rcParams.update({'font.size': 18})

# calculate outer ring
outer = disk(centre, radius)

# initialise the landscape and outer ring
landscape_bk = zeros((dim, dim))
landscape_bk[outer] = 1

# plot
fig, my_ax = subplots(1, 3, figsize=(25, 8))


''' FIGURE A '''

# loop through each pair of functions
for f1, f2 in [(None, fab)]:

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

            # apply uncorrected weighting
            if f1 is not None:
                uc_surf[(r,c)] = f1(dist, radius)
                
            # apply corrected weighting
            c_surf[(r,c)] = f2(dist, radius)

        # mean value
        uncorrected.append(sum(uc_surf[inner]) / sum(uc_surf[outer]))
        corrected.append(sum(c_surf[inner]) / sum(c_surf[outer]))
        
    # add results from each function
    my_ax[0].plot(x, uncorrected, linewidth=3, color='#1f78b4')
    my_ax[0].plot(x, x, linewidth=3, color='#33a02c')

# add labels
my_ax[0].legend(['Zonal Assumption', 'Focal Assumption'], loc='lower right')


''' FIGURE B & C '''

# add 1:1 lines
my_ax[1].plot(x, x, color='#000', linewidth=3, alpha=0.67, linestyle=(5, (10, 3)))
my_ax[2].plot(x, x, color='#000', linewidth=3, alpha=0.67, linestyle=(5, (10, 3)))

# loop through each pair of functions
for f1, f2, colour in [
    (idw, fab_idw, '#1b9e77'), 
    (idw2, fab_idw2, '#d95f02'), 
    (gauss, fab_gauss, '#7570b3')
    ]:

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

            # apply uncorrected weighting
            if f1 is not None:
                uc_surf[(r,c)] = f1(dist, radius)
                
            # apply corrected weighting
            c_surf[(r,c)] = f2(dist, radius)

        # mean value
        uncorrected.append(sum(uc_surf[inner]) / sum(uc_surf[outer]))
        corrected.append(sum(c_surf[inner]) / sum(c_surf[outer]))
        
    # add results from each function
    my_ax[1].plot(x, uncorrected, linewidth=3, color=colour)    
    my_ax[2].plot(x, corrected, linewidth=3, color=colour)

# add labels
my_ax[1].legend(['$IDW$', '$IDW^2$', '$Gaussian$ $(σ=1)$','$Unweighted$'], loc='lower right')
my_ax[2].legend(['$IDW$', '$IDW^2$', '$Gaussian$ $(σ=1)$','$Unweighted$'], loc='lower right')

# format plots
for n in range(3):
    my_ax[n].grid(True)
    my_ax[n].set_xlabel('Radius')
    my_ax[n].set_ylabel('Effective Area')

# add labels
# my_ax[0].annotate('A. Focal and Zonal Assumption)',
#         xy=(0, 1), xycoords='axes fraction',
#         xytext=(+0.5, -0.5), textcoords='offset fontsize',
#         fontsize='large', verticalalignment='top', fontfamily='sans serif')
# my_ax[1].annotate('B. Geographical Weighting (Zonal Assumption)',
#         xy=(0, 1), xycoords='axes fraction',
#         xytext=(+0.5, -0.5), textcoords='offset fontsize',
#         fontsize='large', verticalalignment='top', fontfamily='sans serif')
# my_ax[2].annotate('C. Geographical Weighting (Focal Assumption)',
#         xy=(0, 1), xycoords='axes fraction',
#         xytext=(+0.5, -0.5), textcoords='offset fontsize',
#         fontsize='large', verticalalignment='top', fontfamily='sans serif')
my_ax[0].set_title('A. Focal and Zonal Assumption')
my_ax[1].set_title('B. GW with Zonal Assumption')
my_ax[2].set_title('C. GW with Focal Assumption')

# reduce gaps
subplots_adjust(wspace=0.2)

# output plots
savefig("Figure2.png", bbox_inches='tight', dpi=300)