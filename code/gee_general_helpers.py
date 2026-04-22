import ee
import re

from typing import Literal, Any


def ensure_imagecollection(collection_id):
    try:
        asset = ee.data.getAsset(collection_id)
        
        if asset["type"] != "IMAGE_COLLECTION":
            raise ValueError(f"{collection_id} exists but is not an ImageCollection.")
        
        print(f"ImageCollection already exists: {collection_id}")
    
    except ee.EEException:
        print(f"Creating ImageCollection: {collection_id}")
        
        ee.data.createAsset(
            {"type": "IMAGE_COLLECTION"},
            collection_id
        )



def sanitize_name(x):
    x = str(x)
    x = re.sub(r"[^A-Za-z0-9_-]+", "_", x)
    x = re.sub(r"_+", "_", x).strip("_")
    return x


def filter_asset_ids(asset_ids, image_string=None):
    if image_string is None:
        return asset_ids
    return [a for a in asset_ids if image_string in a]


def list_child_assets(parent):
    """
    Return child asset IDs directly under a parent asset path.
    Works well for folders and ImageCollection asset paths.
    """
    response = ee.data.listAssets({"parent": parent})
    assets = response.get("assets", [])
    return [a["name"] for a in assets]


def get_assets_from_imagecollection(collection_path, image_string=None):
    """
    Prefer listing the collection's child assets directly.
    This gives actual asset IDs like:
    projects/.../assets/my_collection/2010
    """
    asset_ids = list_child_assets(collection_path)
    asset_ids = filter_asset_ids(asset_ids, image_string)
    return asset_ids


def get_assets_from_folder(folder_path, image_string=None):
    asset_ids = list_child_assets(folder_path)
    asset_ids = filter_asset_ids(asset_ids, image_string)
    return asset_ids


def get_export_name(image, asset_id, export_name_property):
    """
    Try requested property, then system:index, then basename of asset ID.
    """
    prop_val = image.get(export_name_property).getInfo()

    if prop_val is None:
        prop_val = image.get("system:index").getInfo()

    if prop_val is None:
        prop_val = asset_id.split("/")[-1]

    return sanitize_name(prop_val)





# Utility functions

# -------------------------------
# Reclassify an image and mask all values except a specific value (e.g., forest class)
def reclassify_image_binary(value):
    def _reclassify(image):
        binary = image.eq(value).selfMask()
        return binary.copyProperties(image, ['system:time_start'])
    return _reclassify

# -------------------------------
# Combine two collections using logical AND (per image pair), copying time from the first image
def combine_collections_with_and(collection1, collection2, and_name):
    list1 = collection1.toList(collection1.size())
    list2 = collection2.toList(collection2.size())
    combined = list1.zip(list2)

    def combine_pair(pair):
        image1 = ee.Image(ee.List(pair).get(0))
        image2 = ee.Image(ee.List(pair).get(1))
        and_result = image1.And(image2).rename(and_name)
        return and_result.copyProperties(image1, ['system:time_start'])

    return ee.ImageCollection(combined.map(combine_pair))

# -------------------------------
# Add a 'year' property to each image
def add_year_property(image):
    year = image.date().get('year')
    return image.set('year', year)



def moving_window(
    image: ee.Image,
    window_size: int | float,
    reducer: ee.Reducer,
    kernel_shape: Literal[
        'square', 'circle', 'rect', 'cross', 'gaussian', 'manhattan', 'chebyshev'
    ] = 'square',
    keep_original_bandnames: bool = False,
    **kwargs: Any
) -> ee.Image:
    """
    Apply a moving window operation on an image using a specified reducer, window size, and kernel shape.

    Parameters:
    - image: ee.Image
    - window_size: int or float, diameter in meters
    - reducer: ee.Reducer (e.g., ee.Reducer.mean())
    - kernel_shape: One of ['square', 'circle', 'rect', 'cross', 'gaussian', 'manhattan', 'chebyshev']
    - keep_original_bandnames: If True, retain the original band names after filtering
    - kwargs: Additional parameters:
        - 'xRadius', 'yRadius' for 'rect'
        - 'sigma' for 'gaussian'
        - 'normalize' (bool) for all kernels (default: False)

    Returns:
    - ee.Image: Result of the focal/moving window operation
    """
    normalize: bool = kwargs.get('normalize', False)

    if kernel_shape == 'square':
        kernel = ee.Kernel.square(radius=window_size / 2, units='meters', normalize=normalize)
    elif kernel_shape == 'circle':
        kernel = ee.Kernel.circle(radius=window_size / 2, units='meters', normalize=normalize)
    elif kernel_shape == 'rect':
        xRadius = kwargs.get('xRadius', window_size / 2)
        yRadius = kwargs.get('yRadius', window_size / 2)
        kernel = ee.Kernel.rect(xRadius=xRadius, yRadius=yRadius, units='meters', normalize=normalize)
    elif kernel_shape == 'cross':
        kernel = ee.Kernel.cross(radius=window_size / 2, units='meters', normalize=normalize)
    elif kernel_shape == 'gaussian':
        sigma = kwargs.get('sigma', window_size / 6)
        kernel = ee.Kernel.gaussian(radius=window_size / 2, sigma=sigma, units='meters', normalize=normalize)
    elif kernel_shape == 'manhattan':
        kernel = ee.Kernel.manhattan(radius=window_size / 2, units='meters', normalize=normalize)
    elif kernel_shape == 'chebyshev':
        kernel = ee.Kernel.chebyshev(radius=window_size / 2, units='meters', normalize=normalize)
    else:
        raise ValueError(f"Unsupported kernel shape: {kernel_shape}")

    result = image.reduceNeighborhood(reducer=reducer, kernel=kernel)

    if keep_original_bandnames:
        original_names = image.bandNames()
        result = result.rename(original_names)

    return result


def remove_to_bands_append(image):
    """
    Renames bands in the input Earth Engine Image by removing the first part
    (before the first underscore) from each band name.

    Args:
        image (ee.Image): The input image.

    Returns:
        ee.Image: The image with renamed bands.
    """
    old_names = image.bandNames()
    new_names = old_names.map(lambda name: ee.String(name).split('_').slice(1).join('_'))
    return image.rename(new_names)


def remove_to_bands_append2(image, og_ic):
    """
    Removes image ID prefixes (from toBands()) using IDs from the original ImageCollection.
    
    Args:
        image (ee.Image): The image with compound band names after toBands().
        og_ic (ee.ImageCollection): The original ImageCollection before stacking.
    
    Returns:
        ee.Image: Image with cleaned band names.
    """
    ids = ee.List(og_ic.aggregate_array('system:index'))

    def clean_name(name):
        name = ee.String(name)
        def iterate_fn(id_, acc):
            acc_str = ee.String(acc)
            id_str = ee.String(id_).cat('_')
            match = acc_str.slice(0, id_str.length()).equals(id_str)
            return ee.Algorithms.If(match, acc_str.slice(id_str.length()), acc_str)
        return ids.iterate(iterate_fn, name)

    old_names = image.bandNames()
    new_names = old_names.map(lambda b: clean_name(b))
    return image.rename(new_names)

def prepend_or_append_to_band_names(
    image: ee.Image,
    text: str,
    mode: Literal["APPEND", "PREPEND"]
) -> ee.Image:
    """
    Appends or prepends a string to all band names of the given Earth Engine image.

    Args:
        image (ee.Image): The input Earth Engine image.
        text (str): The string to append or prepend.
        mode (Literal["APPEND", "PREPEND"]): Specifies whether to append or prepend the string.

    Returns:
        ee.Image: Image with updated band names.
    """
    old_names = image.bandNames()

    if mode == "APPEND":
        new_names = old_names.map(lambda name: ee.String(name).cat(text))
    elif mode == "PREPEND":
        new_names = old_names.map(lambda name: ee.String(text).cat(name))
    else:
        # This should never be reached due to Literal type, but it's good practice
        raise ValueError("Mode must be either 'APPEND' or 'PREPEND'")

    return image.rename(new_names)


def insert_to_band_names(image: ee.Image, text: str, index: int) -> ee.Image:
    """
    Inserts a string into each band name of an Earth Engine image at the specified index.
    
    Args:
        image (ee.Image): The input image.
        text (str): The string to insert.
        index (int): The character index at which to insert the string.
                     Supports negative indexing from the end.
    
    Returns:
        ee.Image: Image with updated band names.
    """
    def insert_at(name):
        name_str = ee.String(name)
        name_len = name_str.length()
        
        # Compute final insertion index (accounting for negative values)
        insert_pos = ee.Algorithms.If(
            index >= 0,
            ee.Number(index),
            name_len.add(index)  # index is negative, so subtract from length
        )
        insert_pos = ee.Number(insert_pos).clamp(0, name_len)  # avoid out-of-bounds
        
        prefix = name_str.slice(0, insert_pos)
        suffix = name_str.slice(insert_pos)
        return prefix.cat(text).cat(suffix)

    old_names = image.bandNames()
    new_names = old_names.map(insert_at)
    return image.rename(new_names)

def toBands_with_projection(collection):
    collection = ee.ImageCollection(collection)
    image = collection.toBands()
    reference_img = ee.Image(collection.first())
    return image.setDefaultProjection(reference_img.projection())

def filter_bands_by_year(image: ee.Image, first_year: int, last_year: int) -> ee.Image:
    band_names = image.bandNames()

    # Create list of allowed years as strings
    allowed_years = ee.List.sequence(first_year, last_year).map(
        lambda y: ee.Number(y).format('%04d')
    )

    # Map over band names and keep only those that match an allowed year suffix
    def keep_if_valid(band):
        band_str = ee.String(band)
        year_suffix = band_str.slice(-4)
        return ee.Algorithms.If(
            allowed_years.contains(year_suffix),
            band_str,
            None
        )

    # Map and filter out None
    valid_bands = band_names.map(keep_if_valid).removeAll([None])

    return image.select(valid_bands)