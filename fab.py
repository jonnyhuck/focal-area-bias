import numpy as np
from scipy.ndimage import distance_transform_edt

class FAB:
    """
    Class for generic Focal Area Bias correction (vectorized).
    """

    def __init__(self, distance_raster, int_rounding=True, circle=True, resolution=None):
        """
        Constructor - pre-calculate all of the unique buffer areas required at runtime
        """
        # convert distances to integers if required (in-place if possible)
        if int_rounding:
            distance_raster = distance_raster.astype(int, copy=False)

        # keep shape for later
        self.shape = distance_raster.shape

        # max distance (note: zeros in masked area included, but max will be the largest non-zero)
        self.max_dist = np.max(distance_raster)

        if circle:
            # vectorized fab_circle values for all pixels
            with np.errstate(divide="ignore", invalid="ignore"):    # use errstate to ignore divide-by-zero warnings for d==0
                
                # r / d (r = self.max_dist); d==0 will produce inf, we will mask below
                raw = np.where(distance_raster == 0, 0.0, self.max_dist / distance_raster)

            # normalize so max becomes 1 (if raw is all zeros, avoid division by zero)
            max_raw = np.max(raw)
            if max_raw > 0:
                self.weighted = raw / max_raw
            else:
                self.weighted = raw  # all zeros case

            # now explicitly set centre cell to 1
            center = (int(self.shape[0] / 2), int(self.shape[1] / 2))
            self.weighted[center] = 1.0

        else:
            # polygon (generalised) branch
            if resolution is None:
                raise ValueError("resolution is a required argument for polygon FAB calculation")

            # get cell area
            cell_area = resolution ** 2
            self.max_area = np.count_nonzero(distance_raster) * cell_area

            # build histogram of counts per integer distance
            counts = np.bincount(distance_raster.ravel(), minlength=int(self.max_dist) + 1)

            # cumulative area up to distance d
            cumulative_area = np.cumsum(counts) * cell_area

            # fab values per integer distance index (0..max_dist)
            d_vals = np.arange(len(cumulative_area))
            fab_values = np.zeros_like(cumulative_area, dtype=float)
            mask = d_vals > 0

            # avoid division-by-zero when cumulative_area[mask] == 0
            fab_values[mask] = (self.max_area * d_vals[mask]) / (cumulative_area[mask] * self.max_dist)

            # map fab values back into raster (vectorized)
            self.weighted = fab_values[distance_raster]

            # normalize similarly
            max_w = np.max(self.weighted)
            if max_w > 0:
                self.weighted = self.weighted / max_w

        # denominator for working out proportions
        self.denominator = np.sum(self.weighted)

    def fab(self, d, a):
        """
        Calculate the fab weight for any buffer (generalised / slower) - kept for compatibility.
        """
        return 0 if d == 0 else (self.max_area * d) / (a * self.max_dist)

    def fab_circle(self, d, r):
        """
        Calculate the fab weight for a circular point buffer (specific / faster) - kept for compatibility.
        """
        return 0 if d == 0 else r / d

    def get_fab_correction(self, r=None):
        """
        Either return weight raster or apply weights to a raster
        """
        if r is None:
            return self.weighted
        if r.shape != self.weighted.shape:
            raise ValueError("correction surface and provided surface are not the same size")
        return self.weighted * r

    def get_denominator(self):
        """
        Return denominator for working out weighted proportions
        """
        return self.denominator


def euclidean_distance(arr, resolution, mask):
    """
    Calculate Euclidean distance surface.
        arr: binary raster where focal pixels are 1
        resolution: cell resolution in units of arr coordinates (e.g. metres)
        mask: raster of same shape used to clip (0/1)
    """
    if arr.shape != mask.shape:
        raise ValueError("distance surface and mask surface are not the same size")

    # distance_transform_edt expects True for background; compute distances from the focal 1s
    # We want distance from the point(s), so invert the boolean:
    distance_pixels = distance_transform_edt(~arr.astype(bool), return_distances=True)
    return distance_pixels * resolution * mask
