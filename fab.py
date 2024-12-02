from scipy.ndimage import distance_transform_edt
from numpy import nonzero, sum as np_sum, max as np_max, zeros, unique, where, logical_and, logical_not


class FAB():
    """
    Class for generic Focal Area Bias correction.
    """

    def __init__(self, distance_raster, int_rounding=True, circle=True, resolution=None):
        """
        Constructor - pre-calculate all of the unique buffer areas required at runtime
        """
        # create blank array for output
        self.weighted = zeros(distance_raster.shape)

        # convert distances to integers if required, this makes it much faster with only a minor deccrease in precision
        if int_rounding:
            distance_raster = distance_raster.astype(int)  
        
        # get the max distance in the masked euclidean distance raster
        self.max_dist = np_max(distance_raster)
        
        # for the circle version
        if circle:

            # loop through all unique distances, calculate fab and load into output raster
            for d in unique(distance_raster):
                self.weighted[distance_raster == d] = self.fab_circle(d, self.max_dist)
        
        # for the generalised (polygon) version
        else:

            # get the area of one cell
            try:
                cell_area = resolution **2
            except:
                print("resolution is a required argument for polygon FAB calculation")
                exit()
            
            # get the max area of the buffer (i.e., the area of the rasterised focal zone)
            self.max_area = nonzero(distance_raster)[0].size * resolution**2
            
            # loop through all unique distances, calculate fab and load into output raster
            for d in unique(distance_raster):
                self.weighted[distance_raster == d] = self.fab(d, 
                    np_sum(where(logical_and(distance_raster > 0, distance_raster <= d), 1, 0)) * cell_area)

        # rescale so that the weights sum to 0-1 (not really necessary)
        self.weighted /= np_max(self.weighted)

        # set the centre cell to 1 for circle version
        if circle: 
            self.weighted[(int(self.weighted.shape[0]/2), int(self.weighted.shape[1]/2))] = 1

        # define the denominator (for working out proportion)
        self.denominator = np_sum(self.weighted)
        

    def fab(self, d, a):
        """
        Calculate the fab weight for any buffer (generalised / slower)
        """
        return 0 if d == 0 else (self.max_area * d) / (a * self.max_dist)
    

    def fab_circle(self, d, r):
        """
        Calculate the fab weight for a circular point buffer (specific / faster)
        """
        # TODO: should this be 1 if d=0??? That would save having to set to 1 later... 
        #   BUT would this then be normalised to something else...?
        return 0 if d == 0 else r / d


    def get_fab_correction(self, r=None):
        """
        Either weight a raster or return a raster containing the weights
        """
        try:
            return self.weighted if r is None else self.weighted * r
        except ValueError:
            print("ERROR: correction surface and provided surface are not the same size (correction function)")
            exit()


    def get_denominator(self):
        """
        Return the denominator for working out weighted proportions
        """
        return self.denominator


def euclidean_distance(arr, resolution, mask):
    """
    Calculate Euclidean distance surface
    """
    try:
        # get distance matrix in pixels (invert so it calculates distances from zero to non-zero rather than non-zero to zero)
        distance = distance_transform_edt(logical_not(arr), return_distances=True)

        # return distance array, converted to metres and masked
        return distance * resolution * mask

    except ValueError as e:
        print("ERROR: distance surface and mask surface are not the same size")
        exit()