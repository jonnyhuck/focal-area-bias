"""
Produces the dataset (which is drawn in QGIS) - not the figure itself
"""
from affine import Affine
from pandas import DataFrame
from numpy.random import uniform
from shapely.geometry import Point
from rasterio.windows import Window
from matplotlib.lines import Line2D
from rasterio import open as rio_open
from fab import FAB, euclidean_distance
from rasterio.features import rasterize
from matplotlib.pyplot import subplots, savefig, imshow, show, rc
from numpy import sum as np_sum, zeros, count_nonzero, isnan, histogram


##### SETTINGS #####

ITERATIONS = 10000      # number of points to test
MAX_BUFFER_D = 2000     # maximum buffer distance
DEMO = False            # whether to show intermediate plots (kills it after 1)
OUTPUT_CSV = False      # should export results to CSV
OUTPUT_FIG = False       # should export results to chart
MORAN = True            # whether to do LISA analysis

####################

# p value threshold for LISA outputs
significance = 0.05

# convert LISA output to useful labels
quadList = ["NA", "HH", "LH", "LL", "HL"]

def getQuadrants(qs, sigs, acceptable_sig):
    """
    * Return list of quadrant codes depending upon specified significance level
    """
    # return quad code rather than number
    out = []
    for q in range(len(qs)):
        # overrride non-significant values as N/A
        if sigs[q] < acceptable_sig:
            out.append(quadList[qs[q]])
        else:
            out.append(quadList[0])
    return out


# access land cover dataset
fab = None
outputs = []
with rio_open('/Users/jonnyhuck/Documents/UK Land Cover Map/data/a22baa7c-5809-4a02-87e0-3cf87d4e223a/gblcm10m2021.tif') as lcm:
    
    # set AOI
    aoi_bounds = (341660, 371164, 416090, 431034)

    # start loop
    flag = False
    counter = 0
    while counter < ITERATIONS:

        # create a random point
        x = uniform(low=aoi_bounds[0], high=aoi_bounds[2])
        y = uniform(low=aoi_bounds[1], high=aoi_bounds[3])

        # create point & buffer
        point = Point(x, y)
        poly = point.buffer(MAX_BUFFER_D)

        # get the top left & bottom right corner of ther AOI in image space
        bounds = point.buffer(MAX_BUFFER_D + lcm.res[0]).bounds
        tl_img = lcm.index(bounds[0], bounds[3])
        br_img = lcm.index(bounds[2], bounds[1])
        w, h = br_img[1]-tl_img[1], br_img[0]-tl_img[0]

        # read using window for speed / RAM and create affine transform
        lcm_band = lcm.read(1, window=Window(tl_img[1], tl_img[0], w, h))
        affine = Affine(lcm.res[0], 0, bounds[0], 0, -lcm.res[1], bounds[3])

        # reclassify the band (be careful not to conflict!)
        lcm_band[isnan(lcm_band)] = 13  # set sea to saltwater
        lcm_band[lcm_band < 20] = 0     # natural land covers
        lcm_band[lcm_band >= 20] = 1    # suburban, urban
        # lcm_band[lcm_band > 1] = 0     # natural land covers
        if DEMO:
            imshow(lcm_band);show()

        # we only need to do the FAB setup once...
        if fab is None:

            # rasterize the buffer (we don't want to mask as it will involve operations on the lcm dataset directly)
            outer = rasterize([poly], (h, w), fill=0, transform=affine, default_value=1, all_touched=True) 
            denom = count_nonzero(outer)
            if DEMO:
                imshow(outer);show()
            
            # euclidean distance (re-scale distance to coordinate space and clip to outer buffer)
            point_raster = zeros((h,w))
            r, c = ~affine * (point.x, point.y)
            point_raster[(int(r), int(c))] = 1
            distance = euclidean_distance(point_raster, lcm.res[0], outer)
            if DEMO:
                imshow(distance);show()

            # init FAB object
            fab = FAB(distance)

        # don't use 1's or 0's
        before_val = np_sum(lcm_band * outer) / denom
        if before_val in [0,1]:
            continue
        counter +=1

        # apply correction to surface
        corrected = fab.get_fab_correction(lcm_band)

        # store the proportions of natural land cover before & after and increment counter
        outputs.append({
            'points': point,
            'before': before_val,
            'after': np_sum(corrected) / fab.get_denominator()
        })

# store dataframe
df = DataFrame(outputs)
df['difference'] = df.after - df.before             # difference between after & before
df['abs_diff'] = abs(df.difference)                 # absolute difference between after & before
df['magnitude'] = abs(df.difference) / df.before    # difference as a proportio of the original value


# report magnitudes
print("\nMagnitude:")
counts, bin_edges = histogram(df.magnitude, bins=[0, 0.1, 0.2, 0.3 ,0.4 ,0.5 ,0.6 ,0.7 ,0.8 ,0.9 ,1, 100])
to_sum = []
for c, b in zip(counts, bin_edges):
    to_sum.append(c/10000)
    s = sum(to_sum)
    print(f"{b:.6f}\t{s:.6f}\t{1-s:.6f}")
print()

# output data is required
if OUTPUT_CSV:
    df.to_csv(f"./outputs/lcsim_{ITERATIONS}.csv")

# output chart if required
if OUTPUT_FIG:

    # plot
    FONT = 14
    rc('font', size=FONT)          # controls default text sizes
    rc('axes', labelsize=FONT)    # fontsize of the x and y labels
    rc('legend', fontsize=FONT)    # legend fontsize
    my_fig, my_ax = subplots(1, 1, figsize=(8, 8))

    # set axis limits
    xlim = [0, 1]
    ylim = [0, int(ITERATIONS/2)]

    # difference between before & after
    df.difference.hist(ax=my_ax, facecolor="#ccebc5")
    my_ax.set_ylim(ylim)
    my_ax.axvline(df.abs_diff.median(), color='k', linestyle='dashed')
    my_ax.axvline(df.abs_diff.min(), color='#e41a1c', linestyle='dashed')
    my_ax.axvline(df.abs_diff.max(), color='#e41a1c', linestyle='dashed')
    my_ax.set_xlabel('Difference in urban-ness value')

    # legend
    my_ax.legend(handles=[
        Line2D([0], [1], color='k', linestyle='dashed', label='Median Value'),
        Line2D([0], [1], color='#e41a1c', linestyle='dashed', label='Min/Max Value')
        ], loc='upper right')

    # output figure
    savefig(f"Figure6-tmp.png", bbox_inches='tight')
    # show()

print(f"{df.abs_diff.min():.2f}, {df.abs_diff.median():.2f}, {df.abs_diff.max():.2f}")
print(f"{df.before.min():.2f}, {df.before.median():.2f}, {df.before.max():.2f}")
print(f"{df.after.min():.2f}, {df.after.median():.2f}, {df.after.max():.2f}")
print()

if MORAN:
    from geopandas import GeoDataFrame
    from esda.moran import Moran, Moran_Local
    from splot.esda import plot_moran, lisa_cluster
    from libpysal.weights import DistanceBand, min_threshold_distance

    # create GeoDataFrame
    gdf = GeoDataFrame(df, geometry=df['points'], crs=lcm.crs)

    # calculate min distance
    min_threshold = min_threshold_distance(list(zip(gdf.geometry.x, gdf.geometry.y))) * 1.1
    print(f'min threshold distance: {min_threshold}')

    # calculate and row standardise weights matrix
    W = DistanceBand.from_dataframe(gdf, min_threshold)
    W.transform = 'r'
    print(f"max neighbours: {W.max_neighbors}, min neighbours: {W.min_neighbors}")

    # calculate and report global Moran's I (Null Hypothesis: CSR)
    mi = Moran(
        gdf['abs_diff'], 
        W, 
        permutations=9999
        )
    print(f"\nGlobal Moran's I Results")
    print("I:\t\t\t", mi.I)					   # value of Moran's I
    print("Expected I:\t\t", mi.EI)			   # expected Moran's I
    print("Simulated p:\t\t", mi.p_sim, "\n")  # simulated p

    # plot output
    # plot_moran(mi, zstandard=True, figsize=(10,4))
    # savefig(f"./moran.png")

    # only do LISA if significant
    if mi.p_sim < 0.005:

        # calculate local I (Rate Adjusted for Population)
        lisa = Moran_Local(
            gdf['abs_diff'], 
            W, 
            transformation='R', 
            permutations=9999)

        # update GeoDataFrame
        gdf['Morans_I'] = lisa.Is                                        # value of Moran's I
        gdf['sig'] = lisa.p_sim                                          # simulated p
        gdf['quadrant'] = getQuadrants(lisa.q, lisa.p_sim, significance) # quadrant (HH, HL, LH, LL)

        # output to shapefile
        gdf.drop('points', axis=1).to_file('sim_lisa.shp')

        # plot local moran figure
        # lisa_cluster(lisa, gdf, significance)
        # savefig(f"./lisa.png", bbox_inches='tight')
    
    else:
        print("No LISA as I not significant")

    # report quick histograms of before, after and difference
    # g = gdf[gdf['quadrant'] == 'HH'][['before', 'after', 'difference', 'abs_diff']]
    # g.hist()
    # show()