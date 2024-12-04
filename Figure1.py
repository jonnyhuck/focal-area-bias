"""
Produce the 'Barnoldswick' maps for Figure 1
"""
from affine import Affine
from geopandas import GeoSeries
from shapely.geometry import Point
from rasterio.windows import Window
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from rasterio import open as rio_open
from fab import FAB, euclidean_distance
from rasterio.features import rasterize
from rasterio.plot import show as rio_show
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.colors import LinearSegmentedColormap
from numpy import sum as np_sum, zeros, count_nonzero, isnan
from matplotlib.pyplot import subplots, savefig, subplots_adjust, text

# buffer size
MAX_BUFFER_D = 2000

# output image
fig, my_ax = subplots(1, 2, figsize=(12, 8))

# access land cover dataset
fab = None
outputs = []
with rio_open('gblcm10m2021.tif') as lcm:

    # for each town
    for i, coords in enumerate([(387854, 446942), (389270, 446930)]):
        
        # create point & buffer
        x, y = coords
        point = Point(x, y)
        poly = point.buffer(MAX_BUFFER_D)

        # get the top left & bottom right corner of ther AOI in image space
        bounds = point.buffer(MAX_BUFFER_D + (lcm.res[0]*5)).bounds
        tl_img = lcm.index(bounds[0], bounds[3])
        br_img = lcm.index(bounds[2], bounds[1])
        w, h = br_img[1]-tl_img[1], br_img[0]-tl_img[0]

        # read using window for speed / RAM and create affine transform
        lcm_band = lcm.read(1, window=Window(tl_img[1], tl_img[0], w, h))
        affine = Affine(lcm.res[0], 0, bounds[0], 0, -lcm.res[1], bounds[3])

        # reclassify the band
        lcm_band[isnan(lcm_band)] = 13  # set sea to saltwater
        lcm_band[lcm_band < 20] = 1     # natural land covers
        lcm_band[lcm_band >= 20] = 0    # suburban, urban

        # we only need to do the FAB setup once...
        if fab is None:

            # rasterize the buffer (we don't want to mask as it will involve operations on the lcm dataset directly)
            outer = rasterize([poly], (h, w), fill=0, transform=affine, default_value=1, all_touched=True) 
            
            # euclidean distance (re-scale distance to coordinate space and clip to outer buffer)
            point_raster = zeros((h,w))
            r, c = ~affine * (point.x, point.y)
            point_raster[(int(r), int(c))] = 1
            distance = euclidean_distance(point_raster, lcm.res[0], outer)

            # init FAB object
            fab = FAB(distance)

        # apply weighting to surface
        weighted = fab.get_fab_correction(lcm_band)

        # store the urbanness values
        print(f'before: {1 - (np_sum(lcm_band * outer) / count_nonzero(outer))}, after: {1-(np_sum(weighted) / fab.get_denominator())}')

        # add layer to map
        rio_show(
            lcm_band,
            ax=my_ax[i],
            transform=affine,
            cmap = LinearSegmentedColormap.from_list('binary', [(0.6, 0.6, 0.6, 1), (178/255, 223/255, 138/255, 1)], N=2)
            )
        
        GeoSeries(poly).plot(
            ax=my_ax[i],
            facecolor = "None",
            linewidth = 2,
            edgecolor = 'black',
            )
        
        GeoSeries(point).plot(
            ax=my_ax[i],
            facecolor = "black",
            edgecolor = 'None',
            )

# decorate maps
x, y, arrow_length = 0.97, 0.99, 0.1

# remove ticks
for the_ax in my_ax:
    the_ax.xaxis.set_tick_params(labelbottom=False)
    the_ax.yaxis.set_tick_params(labelleft=False)
    the_ax.set_xticks([])
    the_ax.set_yticks([])

# add north arrow
my_ax[1].annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
    arrowprops=dict(facecolor='black', width=5, headwidth=15),
    ha='center', va='center', fontsize=20, xycoords=my_ax[1].transAxes)

# add scalebar
my_ax[1].add_artist(ScaleBar(dx=1, units="m", location="lower right"))

my_ax[1].legend(
    handles=[
        Line2D([0], [0], marker='o', color=(1,1,1,0), label='Focal Location', markerfacecolor='black', markersize=8),
        Patch(facecolor=(1, 1, 1, 1), edgecolor=(0, 0, 0, 1), label=f'Focal Zone'),
        Patch(facecolor=(0.7, 0.7, 0.7, 1), edgecolor=None, label=f'Urban'),
        Patch(facecolor=(178/255, 223/255, 138/255, 1), edgecolor=None, label=f'Green'),
    ], 
    loc='lower left', 
    bbox_to_anchor=(1, 0)
    )

# reduce gaps
subplots_adjust(wspace=0.03)

# add labels
my_ax[0].annotate('A',
        xy=(0, 1), xycoords='axes fraction',
        xytext=(+0.5, -0.5), textcoords='offset fontsize',
        fontsize='large', verticalalignment='top', fontfamily='sans serif')
my_ax[1].annotate('B',
        xy=(0, 1), xycoords='axes fraction',
        xytext=(+0.5, -0.5), textcoords='offset fontsize',
        fontsize='large', verticalalignment='top', fontfamily='sans serif')

# save the result
savefig('./Figure1.png', bbox_inches='tight', dpi=300)
print("done!")